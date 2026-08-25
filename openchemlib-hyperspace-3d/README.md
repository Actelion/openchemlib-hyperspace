# Hyperspace3D

Hyperspace3D is the learned, graph-only 3D screening extension for
OpenChemLib Hyperspace. This initial module provides the deployable vertical
slice: OCL graph featurization, batched ONNX inference, a sampled-product
fingerprint index, streaming elite screening, and the batch contract needed
by a later local/beam optimizer.

## Data boundaries

The raw synthon space remains the authoritative source of reactions,
synthon sets, IDs, and structures. SkelSpheres downsampling/neighbor data
remains a separate layer used to choose representatives and local
substitutions. Complete sampled products and their 128D learned embeddings
live in sharded product-index files; they are never added to rawspace JSON.

## Selected model

The deployed model is the original compact Deepspace7 V1 similarity
architecture: the protected six-layer chemistry/geometry trunk (node latent
16, pair latent 14) followed by the graph-only, node-pair, normalized 128D
molecular encoder and the trained symmetric comparator (token hidden 64,
pair-comparison hidden 128). Inference uses no conformer, physical
coordinates, query-pose distance matrix, reconstruction heads, or loss code.

The primary model now uses the final 50D Deepspace7 V2 graph contract:

* atom tensor: float32 `[B,32,50]`
* ordered pair tensor: float32 `[B,32,32,36]`
* atom mask: bool `[B,32]`

All 50 atom channels and all 36 pair channels are checked against the
Python/RDKit golden fixture. Pharmacophore labels remain output supervision
for the protected trunk; they are not molecular graph inputs.

The bundled deployment uses the final `publication_v2_50d_ppaux_w025_2026`
foundation (epoch 29, pharmacophore auxiliary loss weight 0.25) and the canonical
balanced `wpp050_seed17_d128` similarity predictor (epoch 26, PheSA PP weight
0.5). These are separate training weights and both are recorded in the bundle.

## Model bundle

The canonical model used by this initial vertical slice is included at
`model-bundles/deepspace7-v1`:

```text
model-bundles/deepspace7-v1/
  manifest.json
  feature-schema.json
  encoder.onnx
  comparator.onnx
  checksums.sha256
  verification.json
```

The original training checkpoints and optimizer state are not included.
Future or substantially larger production bundles should remain external and
use the same versioned directory contract.

The loader verifies the manifest contract and SHA-256 hashes before opening
sessions. CPU and CUDA selection is explicit; a CUDA request fails if the
CUDA provider is unavailable. The default Maven dependency is CPU ONNX
Runtime. Build with `-Pcuda` for `onnxruntime_gpu`. Packaging is a normal
thin JAR plus `target/lib`, avoiding native-resource risks from shading.

The comparator preserves manifest target order:
`ffp_similarity`, `skelspheres_similarity`, `flexophore_similarity`,
`phesa_total`, `phesa_pharmacophore`, and `phesa_shape`.
`phesaPpWeight` is read from the bundle and is not treated as a universal
constant.

## Search flow

The query is featurized and encoded once. `FingerprintIndexScreener`
streams product embeddings in large batches, performs exact tuple
deduplication, and retains global and per-reaction top K results.
`phesa_total` can be selected directly; shape/pharmacophore composites are
explicitly separate objectives.

`BatchAssemblyScorer` establishes the later online search boundary. Its
initial implementation deduplicates tuples, assembles with the existing
`SynthonAssembler`, applies deterministic feature filters, batches the
encoder/comparator, routes results by stable candidate ID, and caches
embeddings per run. A production beam coordinator can add bounded CPU
workers, a shared bounded GPU queue, backpressure, and neighbor proposals
without replacing these interfaces. Exact PheSA, conformer generation,
force-field work, and ligand preparation remain downstream cluster jobs.

Candidate archives are streaming JSON Lines records with tuple identity,
total/shape/pharmacophore predictions, source, basin, round, bundle hash,
and query provenance.

## Build and verification

```bash
mvn -pl openchemlib-hyperspace-3d -am package -Dmaven.test.skip=true
mvn -pl openchemlib-hyperspace-3d test
mvn -pl openchemlib-hyperspace-3d test \
  -Dhyperspace3d.modelBundle=/path/to/deepspace7-v1
```

The final command enables Java ONNX/Python golden parity. Python export lives
in the sibling Deepspace7 repository:

```bash
python -m deepspace7.scripts.export_compact_v1_onnx \
  --foundation-checkpoint /path/to/selected.pt \
  --predictor-checkpoint /path/to/best.pt \
  --output-dir /path/to/deepspace7-v1 \
  --golden-fixture /path/to/golden.json.gz
```

The benchmark harness can be run from the packaged module:

```bash
java -cp 'target/classes:target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.benchmark.Hyperspace3DMicrobenchmark \
  /path/to/deepspace7-v1
```

These are first-order microbenchmarks, not JMH results. CPU preparation and
ONNX inference are reported separately.

## Building a sampled product fingerprint index

`ProductFingerprintIndexBuilder` creates the persistent data consumed by
`FingerprintIndexScreener`. It is query-independent: sampled products are
assembled from the downsampled representative synthon sets, featurized on CPU,
and encoded once into 128D vectors. The full rawspace remains authoritative for
structures and zero-based synthon ordinals.

Sampling is deterministic and without replacement inside each reaction. The
builder first attempts the configured minimum accepted coverage per reaction,
then fills the global record target using either the existing `EXPONENT` or
`BUCKETED_PRODUCT` reaction weighting. Rejected assemblies are counted by
stable reason and are not written.

Copy and edit `product-fingerprint-index-build.example.json`. All relative paths
are resolved against that JSON file. Build the normal JAR/lib-directory
distribution and run:

```bash
mvn -pl openchemlib-hyperspace-3d -am package -Dmaven.test.skip=true
java -cp 'openchemlib-hyperspace-3d/target/classes:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.ProductFingerprintIndexBuilderCLI \
  --config openchemlib-hyperspace-3d/product-fingerprint-index-build.example.json
```

Use `-Pcuda` while packaging when the JSON selects `CUDA`. CUDA requests never
fall back silently to CPU. Omitted encoder batch sizes default to 64 on CPU and
512 on CUDA.

An output directory contains `shard-NNNNN.h3di`, `build-state.json`, and, after
successful completion, `manifest.json`. Each shard is written to a temporary
file, SHA-256 checked, and committed with its checkpoint using atomic moves.
Interrupted runs resume at the last committed shard when `output.resume` is
true. Input, model, semantic configuration, and existing shard hashes must all
match; unrelated output is never overwritten.

The manifest records rawspace/downsampled hashes, model-bundle hash, sampling
algorithm and seed, reaction weighting, filters, shard metadata, rejection
counts, and separate CPU preparation, tensor packing, ONNX, writing, and
hashing timings. A typical float32 shard with one million records is roughly
0.5 GB plus tuple metadata.

## Learned 16D SkelSpheres search

The compact 2D path projects the canonical normalized 128D embedding through
the upstream seed-17 MLP (128 to 128 to 16) and L2 normalizes the result.
Similarity uses the trained sigmoid calibration from the bundle manifest.
The projection bundle is included at `model-bundles/deepspace7-skelspheres16`.

Compact indices use backward-compatible format version 2; version-1 128D
readers and files are unchanged. Each v2 shard separates sequential float16
vectors (`.vec`), fixed reaction/tuple references (`.rows`), and variable
synthon tuples (`.tuples`). Screening therefore reads only 32 vector bytes per
product plus fixed row metadata and resolves tuple payloads only for retained
global and per-reaction elites. All columns are independently checksummed and
the manifest binds both rawspaces, both model bundles, sampling, and filters.

Build directly from assembled graphs (there is no 128D-index conversion path):

```bash
java -cp 'openchemlib-hyperspace-3d/target/classes:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.CompactSkelSpheresIndexBuilderCLI \
  --config openchemlib-hyperspace-3d/compact-skelspheres-index-build.example.json
```

Run the standalone flat 2D screen:

```bash
java -cp 'openchemlib-hyperspace-3d/target/classes:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.CompactSkelSpheresSearchCLI \
  --config openchemlib-hyperspace-3d/compact-skelspheres-search.example.json
```

The query is encoded and projected once. Optional exact OCL SkelSpheres
reranking is disabled by default; when enabled, global and per-reaction learned
shortlists are unioned, assembled, and reranked on CPU. Hybrid 2D-to-3D search
and ANN remain explicit future extensions.

The Deepspace7 exporter is
`python -m deepspace7.scripts.export_skelspheres16_onnx`. It rejects
noncanonical dimensions, seed, provenance, or a failed quality gate; verifies
ONNX Runtime against PyTorch at batch sizes 1, 7, 256, and 4096; and records
float16 score error in `verification.json`.

## Building a flat supplier molecule fingerprint index

`MoleculeFingerprintIndexBuilder` is the high-throughput Java/ONNX equivalent
of the Python supplier-cache workflow. It accepts UTF-8 TSV, TSV.GZ, or TSV.BZ2
input with configurable SMILES and molecule-ID columns. OCL parsing and graph
featurization run in a bounded CPU pool while prepared batches queue ahead of
the ONNX lane. Each accepted graph is encoded once to the universal 128D
fingerprint and then projected to the 16D SkelSpheres metric fingerprint.

Copy and edit `molecule-fingerprint-index-build.example.json`, package with the
CUDA profile, and run:

```bash
mvn -pl openchemlib-hyperspace-3d -am -Pcuda package -Dmaven.test.skip=true
java -cp 'openchemlib-hyperspace-3d/target/classes:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.MoleculeFingerprintIndexBuilderCLI \
  --config openchemlib-hyperspace-3d/molecule-fingerprint-index-build.example.json
```

ONNX Runtime 1.21 GPU requires the CUDA 12 and cuDNN 9 runtime libraries.
CUDA provider initialization fails explicitly when either runtime is unavailable.

The output uses source-row-defined shard directories:

```text
shard-00000/
  vectors-128.fp16
  vectors-16.fp16
  rows.bin
  strings.bin
  shard.json
  .complete
```

Vector columns are little-endian, row-major float16 and can be memory mapped
independently. Fixed rows contain source-row and heavy-atom metadata; the
variable string table stores the original molecule ID and SMILES. The included
`MoleculeFingerprintIndexReader` exposes sequential batches, random rows, and
independent 128D/16D vector access.

A shard is first written to `.partial` and atomically renamed only after all
files and headers are closed. Restart removes abandoned partial directories and
continues from the last contiguous completed source interval. No index hashes or
active shard checksum scans are performed. A resumed BZip2 input must be
decompressed from the beginning to skip committed rows, but committed molecules
are not reparsed or re-encoded. The manifest reports rejection reasons and
separate preparation, tensor packing, ONNX, waiting, and writing timings.

## Screening a flat supplier molecule index

`MoleculeSkelSpheresSearchCLI` encodes one query into the learned 16D metric,
scans only the compact float16 vector columns, and resolves IDs and SMILES only
for the retained global shortlist. The shortlist is then scored with exact OCL
`DescriptorHandlerBinarySkelSpheres`. Exact ranks therefore apply within the
learned shortlist and are not global exact-search recall measurements.

Copy and edit `molecule-skelspheres-search.example.json`, then run:

```bash
java -cp 'openchemlib-hyperspace-3d/target/classes:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.MoleculeSkelSpheresSearchCLI \
  --config openchemlib-hyperspace-3d/molecule-skelspheres-search.example.json
```

The CLI writes a ranked TSV, an SDF with score fields, a concise Markdown
summary, and a JSON run manifest with scan and reranking timings. CPU is the
recommended query device because only one structure is encoded; the exhaustive
16D scan itself is CPU and does not require CUDA.
