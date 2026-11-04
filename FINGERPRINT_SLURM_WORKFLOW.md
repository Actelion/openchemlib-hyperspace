# Molecule fingerprints: local MCP jobs and Slurm arrays

This workflow generates graph-only 128D Deepspace7 fingerprints and their 16D
learned SkelSpheres projections for **enumerated molecule libraries**. It does
not import synthon spaces, perform pose-specific queries, or run exact PheSA.
The final cache is readable by `MoleculeSimilaritySearchCLI` without conversion.

## Build and deployment

Use Java 22 or newer. Build the CPU distribution for local testing:

```bash
mvn -pl openchemlib-hyperspace-3d,openchemlib-hyperspace-mcp -am clean package -DskipTests
```

For cluster GPU execution, build with `-Pcuda`. Always use a clean build when
switching profiles: do not put CPU and GPU ONNX Runtime jars in the same `lib/`.
Deploy the 3D module's JAR **and its entire `target/lib/` directory**, plus both
model bundles. Do not use the shaded hyperspace CLI JAR as the 3D distribution.
The MCP JAR remains a separate executable and does not load chemistry or ONNX.

The current runtime uses ONNX Runtime 1.21.0. CUDA native libraries, cuDNN, GPU
architecture support, driver compatibility and filesystem behavior must be
verified on the cluster. Exporting a bundle does not establish compatibility
and never installs system software. In particular, old CUDA 12.1/cuDNN 9.1
installations should not be assumed to support RTX 5080 GPUs.

## Local MCP preparation

Add optional paths to the existing MCP server configuration (relative paths
are relative to that file):

```json
{
  "workspace": "hyperspace-workspace",
  "cliJar": "openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar",
  "fingerprintJar": "openchemlib-hyperspace-3d/target/openchemlib-hyperspace-3d-3.1.0.jar",
  "fingerprintLibDirectory": "openchemlib-hyperspace-3d/target/lib",
  "modelBundle": "openchemlib-hyperspace-3d/model-bundles/deepspace7-v1",
  "compactBundle": "openchemlib-hyperspace-3d/model-bundles/deepspace7-skelspheres16",
  "threads": 4,
  "heap": "8G"
}
```

Old server configurations still work for synthon-space tools. Fingerprint
preparation reports missing configuration rather than changing their behavior.

1. `inspect_molecule_library` with `input` inspects at most 65,536 decoded
   characters and returns a small sample. This is not full input validation.
2. `prepare_molecule_fingerprints` accepts `input`, optional `outputDirectory`,
   `smilesColumn` (default `smiles`), `idColumn` (default `id`), `device`
   (`CPU` default or explicit `CUDA`), `threads`, `heap`, `encoderBatchSize`,
   `cudaDeviceId` and `sourceRowsPerShard` (default 1,000,000).
3. Poll `get_job_status`; use `get_job_log` with `log: "fingerprint"` for
   chemistry progress. `cancel_job` stops the computation and retains committed
   shards. Resubmit with the same `outputDirectory` to resume.
4. `list_molecule_libraries` lists successful local preparation jobs. These
   libraries are not synthon spaces and cannot be passed to `configure_gui`.

Only one local preparation job runs per MCP workspace. The worker survives an
MCP client disconnect. All inputs must be headered, tab-delimited UTF-8 tables,
optionally `.gz` or `.bz2`. CSV and headerless SMI are not accepted by this v1.
The model accepts 6-32 heavy atoms; unsuitable molecules are rejected,
not silently truncated. A library with no accepted molecules is not published
as a successful MCP job.

## Export a Slurm bundle

Copy and edit [the profile example](fingerprint-slurm-profile.example.json).
It contains `workflow` paths and settings plus `slurm` resource settings.
All cluster paths must be absolute and accessible from the allocated nodes;
the MCP machine does not need to mount them. Keep application distributions,
model bundles and any environment setup script immutable during a workflow.

Call `export_slurm_fingerprint_job` with:

```json
{
  "profilePath": "/local/path/cluster-profile.json",
  "outputDirectory": "/local/path/new-job-bundle"
}
```

This only writes files. It does not connect via SSH, transfer the library, or
submit anything. Copy the bundle to shared cluster storage using your usual
approved procedure. A trusted, user-managed `environmentSetup` shell script
can load modules or set library paths; scripts source it when executed on the
cluster. Do not put passwords or access tokens in the profile.

From the bundle on the submit node, execute these steps **in order**, waiting
for each stage to finish successfully:

```bash
bash submit_prepare.sh
bash submit_smoke.sh
bash submit_compute.sh
bash submit_finalize.sh
```

Each helper prints a Slurm job ID. Use your site's `squeue`/`sacct` commands and
the bundle's `logs/` directory to check completion. No automatic dependency
chain or cluster monitoring is performed by MCP. The compute helper refuses
submission without a successful smoke record matching the prepared workflow.
The smoke test encodes three molecules through the actual model and projection,
not merely `nvidia-smi`. It proves execution on the tested node, not every node
in a heterogeneous partition.

Default GPU resources are one GPU, eight CPUs, 24G Slurm memory, an 8G builder
heap, six CPU preparation workers, encoder batches of 512, and at most four
concurrent array tasks. The orchestration JVM has a separate 4G heap. Slurm
memory must cover both heaps, ONNX/native allocations and overhead. CPU-only
preparation/finalization use two CPUs, 8G memory and no GPU. Edit limits to fit
cluster policy and measured memory use; these defaults are not throughput
guarantees. `cudaDeviceId` is zero within Slurm's allocated GPU visibility;
the scripts do not overwrite `CUDA_VISIBLE_DEVICES`.

## Stages and restart behavior

**Prepare:** read/decompress the original library once; stream header-preserving
gzip chunks through scratch and publish each to shared `inputs/`. Record source
row intervals and SHA-256 hashes. Default chunk size is one million source rows
(maximum eight million). Duplicate IDs, blank rows and malformed records retain
their original positions; filtering happens in computation, not splitting.
An interrupted preparation rereads input and verifies previously published
chunks. Changed input/configuration/model data requires a new workflow directory.

**Compute:** each array task copies the application and model bundles onto a
private directory under `scratchRoot`, then downloads and processes one chunk
at a time. Fingerprints are written locally. Each completed partition uploads
to a temporary shared directory, is verified, and is atomically published under
`partitions/`. No tasks write a common builder checkpoint. If there are more
partitions than `maxArrayTasks` (default 1,000), task `i` processes partitions
`i`, `i + taskCount`, etc. This respects the configured array-size cap.

Re-submit compute to retry failures: verified completed partitions are skipped.
An incomplete partition is recomputed from its start; no partial scratch data
is required for recovery. Normal exit cleans task-owned scratch. A node crash
or SIGKILL may leave scratch files for your cluster's cleanup policy. Shared
partial uploads are never treated as completed partitions. Do not remove lock
files while a workflow is active.

**Finalize:** verify all expected partitions and publish `cache/manifest.json`.
All-rejected partitions are valid, but missing/corrupt partitions and an entirely
rejected library prevent publication. Source-row metadata is rewritten with
global offsets; IDs, SMILES and FP16 vectors are preserved. The final index
remains sharded. Vector/string files are hard-linked when possible, otherwise
copied with a warning. Finalization checkpoints each shard, so it can be rerun
after interruption. Source partitions are never edited or automatically deleted.

Published payload files are immutable, especially when hard-linked. New caches
record model-content identities so search can use a relocated identical model
bundle. Legacy caches without these identities retain their old path checks.

Shared storage must support atomic rename within a directory and advisory file
locks. Computation avoids repeated network reads of model and input data, but
preparation, upload verification and finalization necessarily perform sequential
shared-storage I/O. Verify available scratch space and shared quota first.
At 1.3 billion accepted molecules, vectors alone use about 374 GB (decimal),
plus strings, row metadata, chunks and temporary uploads. Copy-based finalization
needs another copy of vector/string data; it does not require new inference.

## Direct CLI and acceptance

The same workflow is usable without MCP or Slurm:

```bash
java -cp '/path/hyperspace-3d.jar:/path/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.MoleculeFingerprintWorkflowCLI \
  prepare --config /path/workflow-config.json
```

Replace `prepare` with `smoke` or `finalize`. For a worker, append
`--task 0 --tasks 1` after the config argument. Use `runtime.device: "CPU"`
for a local acceptance run; the exported GPU profile deliberately requires
`CUDA`. Direct worker execution is an advanced/manual path; the smoke gate is
enforced by the production submission helper.

Local verification covers splitting, publication, retry, consolidation, native
cache readability, MCP lifecycle and generated shell syntax. The remaining
cluster acceptance checks are scheduler policy, shared-filesystem semantics,
scratch access/capacity, CUDA/cuDNN/ONNX compatibility, and a representative
throughput/memory pilot before launching the full library.
