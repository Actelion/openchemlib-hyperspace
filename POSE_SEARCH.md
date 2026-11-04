# Supplied-conformation DeepSpace7 search

`pose_phesa_3d` exhaustively screens a supplied MOL/SDF query conformation against
an existing canonical 128D candidate fingerprint index. Candidate vectors describe
graph-only molecular ensembles. Outputs are learned directed pose→ensemble
predictions; hit SDFs contain 2D structures, with no aligned candidate conformations.

The selected upstream architecture is `SharedCandidatePosePredictor`. Its frozen
canonical candidate branch must match every parameter of the original predictor.
The deployment requires the protected six-layer 50/36 graph foundation (32 heavy
atoms, 16/14 latents), query distances in angstrom, a 128D pose query encoder, and
directed comparison heads. Total is calculated as `wPP * PP + (1-wPP) * shape`.
Weights in [0,1] are supported. The default is 0.5. Reflections preserve distance
matrices, so mirrored coordinates need not produce different predictions.

## Export and provenance

In `deepspace5`, use the ONNX-equipped environment:

```bash
python -m deepspace7.export.shared_candidate_pose_onnx \
  --foundation /models/foundation/selected.pt \
  --canonical /models/wpp050_seed17_d128/best.pt \
  --pose /models/shared_candidate/seed17_d128/best.pt \
  --latent-manifest /training/latents_fp16/manifest.json \
  --candidate-manifest /bundles/deepspace7-v1/manifest.json \
  --output /bundles/shared-candidate-pose
```

The export requires a new output directory. It verifies checkpoint type, the
canonical predictor SHA-256, the full immutable branch state, and the training
latent manifest/foundation identity. Evaluation mode disables dropout. ONNX opset
17 exports a dynamic batch with fixed atom dimensions. CPU parity is mandatory
at batches 1, 7, 256, and 1024, with absolute/relative tolerance 2e-5.

The external bundle contains `query.onnx`, `comparator.onnx`, `manifest.json`,
`feature-schema.json`, `verification.json`, and `checksums.sha256`. Its manifest
records tensor names/shapes, directed semantics, weights, feature schema, units,
foundation/canonical/pose checkpoint hashes, and the required canonical encoder
hash. Java checks the bundle against the supplied canonical graph bundle, then
the existing cache reader validates the library's checkpoint/encoder provenance
before any vector scan. Compatible indices are reused directly.

Seed 17 selected locally: pose checkpoint SHA-256
`5a0156760f2304f5002316c35a07f5912377ef945c31de4e6b36e80665e80045`,
canonical predictor `29b2d7770ae3740f7b58492e66d88c7b2743ffa9735f7a74543ad06c7182f2d2`,
foundation `713ec867036277d1197484e6124b876329d23a48cc7792d4fad57e6cb2c9d3e0`.
Models remain external to Git.

## CLI

Build CPU distribution with `mvn -pl openchemlib-hyperspace-3d clean package`.
Use `-Pcuda` for CUDA; clean the build when switching profiles so old native
runtime JARs cannot coexist. Keep the JAR and `target/lib` together. Native
ONNX Runtime is loaded from normal dependencies. CUDA selection is explicit.
CUDA disables TF32 (`use_tf32=0`) to preserve float32 deployment parity; see the
[ONNX Runtime CUDA options](https://onnxruntime.ai/docs/execution-providers/CUDA-ExecutionProvider.html#use_tf32).

Copy `openchemlib-hyperspace-3d/pose-search.example.json`, set paths, then:

```bash
java -cp 'openchemlib-hyperspace-3d/target/*:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.MoleculeSimilaritySearchCLI \
  --config pose-search.json
```

MOL/SDF record selection is zero based. Omit it only for single-record files.
Coordinates are preserved; no conformer generation or minimization occurs.
Explicit hydrogens are deleted with OCL's returned index mapping; every remaining
element and coordinate is verified against its source atom. Atom-map annotations
are cleared before graph featurization. Nonfinite/collapsed coordinates and
unsupported graphs fail with errors. Planar coordinates are accepted.

Java and DeepSpace7 supplier-cache readers provide bounded prefetch, batched ONNX
comparison, deterministic top-K ties, and deferred metadata resolution. TSV/SDF
contain total, shape, PP, selected objective, and score semantics. JSON contains
the source-file hash, record, source→heavy mapping, supplied coordinates, pose
and graph model identity, cache hash, and timings. `queryPreparationNanos` and
`queryEncodingNanos` are separate. Vector preparation is scan wall time less
comparator time; prefetch overlaps I/O and is not a standalone disk benchmark.

The first implementation provides exhaustive pose screening. Graph-based 16D
cascade selection remains available in the existing graph search modes.

## MCP

Configure `fingerprintJar` and `fingerprintLibDirectory` to the built 3D JAR and
its `lib` directory. The MCP package does not shade ONNX native libraries.
The server exposes:

- `inspect_3d_query(queryFile, record?)`: bounded OCL inspection through the 3D CLI.
- `search_pose_molecules(queryFile, index, modelBundle, poseBundle, record?,
  identifier?, objective?, ppWeight?, topK?, device?, cudaDeviceId?, batchSize?,
  prefetchDepth?, heap?)`: detached typed job submission.

`intraOpThreads` explicitly limits CPU inference threads (MCP default 1); CLI
value 0 retains the ONNX Runtime default. `heap` limits the detached search JVM.

Use `get_job_status`, `get_job_log` with `log="pose-search"`, `cancel_job`, and
the status artifact map. Searches survive client reconnects. Output paths are
allocated under the server workspace's job directory. `get_workflow_help` topic
`pose-search` exposes this guide in the packaged server.

## Verification and public fixture

Generate public hexane conformations, planar benzene, and a defined stereocenter:

```bash
python -m deepspace7.export.pose_golden --foundation ... --canonical ... \
  --pose ... --latent-manifest ... --output /deployment/golden
```

`golden.json.gz` contains RDKit tensors, masks, distances, PyTorch embeddings,
canonical candidate vectors and directed scores. MOL coordinates are rounded
before expected outputs are computed. The checked-in fixture contains synthetic
public examples only. With an external bundle at
`/deployment/shared-candidate-seed17-verified`, run:

```bash
mvn -pl openchemlib-hyperspace-3d test -Dpose.deployment=/deployment
```

The model-free tests check hydrogen mapping, invalid coordinates and invariance;
real-model tests verify exact RDKit/OCL tensor parity, ONNX/Java/PyTorch numerical
parity, conformational sensitivity, and ranking across batches/shards. Packaged
MCP tests additionally require `hyperspace.mcpJar`, `hyperspace3d.distributionJar`,
`hyperspace3d.modelBundle` and `hyperspace3d.compactBundle` system properties.

Explicit candidate-conformer databases, pose→pose inference, exact alignment,
and combinatorial optimization are subsequent workflows.

## Verified local deployment

CPU export verification over batches 1, 7, 256 and 1024 observed maximum absolute
errors of 6.49e-7 for embeddings and 3.28e-7 for scores. Java golden checks pass
at 2e-5 on CPU and an RTX 4070 CUDA provider with TF32 disabled, including atom
permutation, mixed comparison batches and directed roles. RDKit/OCL graph tensors
match exactly for the public molecular fixture; its two hexane conformations
differ by 0.00719 in the largest pose-embedding component.

The packaged MCP test launches the thin CPU distribution, reconnects its client,
inspects a real MOL, screens an existing two-record Java index, discovers results,
and checks failure/cancellation. A separate test screens a Python supplier cache
across shards and batch boundaries. The local two-record CLI example measured
36.7 ms query preparation, 73.5 ms cold encoding/session initialization, 9.2 ms
vector preparation/prefetch, 0.38 ms comparison and 5.6 ms metadata resolution.
These tiny cold-run measurements are a functionality check, not a throughput
benchmark. CUDA verification here covers golden inference, not a large-library
performance run.
