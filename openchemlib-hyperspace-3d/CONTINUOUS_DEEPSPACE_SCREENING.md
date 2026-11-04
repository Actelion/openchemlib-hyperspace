# Continuous learned combinatorial screening

`ContinuousDeepSpaceScreeningCLI` samples product tuples, assembles and featurizes
them on CPU, encodes and compares batches through ONNX Runtime, and continuously
maintains global/per-reaction elites. Exact PheSA is an optional final pass.

```bash
mvn -pl openchemlib-hyperspace-3d -am clean package
java -cp 'openchemlib-hyperspace-3d/target/*:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.ContinuousDeepSpaceScreeningCLI \
  --config openchemlib-hyperspace-3d/continuous-deepspace-screening.example.json
```

For CUDA use a clean `-Pcuda` build and set `runtime.device` to `CUDA`; native
CUDA/cuDNN libraries must be available. The existing explicit-device runtime
rejects unavailable CUDA and disables TF32 for float32 parity.

## Inputs and sampling

Paths in JSON are relative to the configuration file. `inputs.rawFull` is the
authoritative rawspace. `inputs.rawDownsampled` optionally supplies representative
synthon sets; omit it to sample the full space directly. The existing
`ProductTupleSampler` verifies reduced synthon IDs, structures and full-space
ordinals. Weighted seeded affine traversal visits each eligible tuple once using
O(reaction-count) cursor state. This is deterministic scrambled enumeration,
not independent random draws or a claim of exhaustive chemical-space coverage.
Different tuples may assemble to the same chemical structure.

`reactionWeighting` uses the existing EXPONENT or BUCKETED_PRODUCT contracts;
weights are calculated from the chosen sampling-space synthon-set sizes. There
is no local optimizer or score-driven adaptation. Filters include atom count
(6–32), optional rotatable-bond maximum (0 disables it), and the model's standard
element/charge/fragment/stereochemistry compatibility checks.

The query is encoded once. For graph screening use `query.smiles`. For a supplied
pose replace smiles with `query.file` (MOL/SDF) and optional zero-based `record`,
and add `inputs.poseBundle`; multi-record files require explicit selection.
The [pose bundle/provenance contract](../POSE_SEARCH.md) applies. Candidates always
use the canonical graph-only encoder. `run.target` selects `phesa_total`,
`phesa_shape` or `phesa_pharmacophore`. Graph mode requires `run.ppWeight` to match
the graph manifest's trained total-head weight. Pose mode computes weighted total
from its directed shape/PP predictions. Scores are learned predictions.

## Execution, limits and recovery

One producer samples tuples and dispatches each preparation window to
`cpuWorkers`. A bounded `queueCapacity` supplies one encoder/comparator consumer.
`encoderBatchSize` limits graph batches; molecular rejection can make an actual
inference batch smaller. Inference sessions are reused. No growing embedding
cache or all-tuples seen-set is needed. Elite memory is bounded by global K plus
per-reaction K; graph-batch memory scales with batch size and queue capacity.
Equal-score ordering uses the reaction and ordered tuple identifiers.

Stop on the first reached `maxAccepted`, `maxProposals`, `maxRuntimeSeconds`, or
sampling-space exhaustion. Zero disables an individual numeric limit. Runtime
counts cumulative mining wall time across restarts and excludes the optional
finishing pass. Limits and stop requests are checked between inference batches;
an in-flight operation can exceed the wall-time limit. The accepted-count cap is
exact: remaining valid tuples in the final preparation window are counted as
`discardedAtAcceptedLimit`, not scored.

At `checkpointSeconds` intervals (0 means every consumed batch), atomically save
sampler position, metrics, global/per-reaction elites and any committed vector
shards. The cursor belongs to the consumed batch, so producer read-ahead is
replayed after a crash. Intermediate uncommitted work may be repeated; committed
tuples are not repeated. Ctrl-C requests a batch-boundary checkpoint and allows
up to 30 seconds for shutdown. Forced termination recovers from the most recent
durable checkpoint. A directory lock prevents concurrent writers.

Rerun the identical configuration with `output.resume=true`. Resume rejects
changed space/query/model/configuration identities and corrupted committed
fingerprint shards. Completed runs regenerate exports without further sampling;
interrupted runs continue. New budgets, resource settings, or queries currently
require a new output directory. If the final exact pass was interrupted it is
repeated, using the already saved learned shortlist. A matching completed exact
report is reused on subsequent invocations.

## Outputs

- `checkpoint.json`: authoritative versioned sampler/elite state.
- `run.json`: provenance, hashes, configuration, score-target order and metrics.
- `hits.tsv`, `reaction-hits.tsv`: periodically updated tuple IDs, full ordinals,
  objective and predicted total/shape/PP components.
- `summary.md`: periodic counters and timings.
- `hits.sdf`: final global elite structures, tuple references and learned objective.

CPU assembly/featurization times are summed worker times, not wall time. Tensor
packing, ONNX encoding and comparison are recorded separately. Query encoding
has a separate provenance timing. The supplied query conformation is preserved.

By default only the retained hits are saved. Set `persistFingerprints=true` to
write every scored product's canonical 128D vector using the existing float32
product-shard format. `fingerprints.json` is the standard versioned product-index
manifest, binding rawspace, sampling space, model, filters and reaction weights.
Only its committed shard paths should be read; `.partial` files and orphaned
unreferenced shards are not published data. This index is reusable for other
queries; the pose-query embedding is never stored as a candidate fingerprint.

## Optional exact PheSA post-filter

Set `exact.enabled=true`. From the saved global elites, reassemble and deduplicate
by OCL canonical structure identity, then label up to `exact.topK` unique products.
`exact.maxConformers` and `exact.threads` control the independent CPU workload.
Graph queries use generated query ensembles. Pose queries use only the supplied
heavy-atom conformation, adding hydrogen coordinates without moving heavy atoms;
planar queries are supported. This is rigid PheSA alignment over generated
candidate conformers, not force-field/beam optimization of the mining space.

`exact-hits.tsv` retains every selected observation, explicit failure status,
learned objective, exact total/shape/PP, and `passes_filter`. Rows are sorted by
the exact counterpart of `run.target`, with failures last. `exact.minScore`
filters on that exact target; `exact-hits.sdf` contains only successful passing
structures. `exact-report.json` records settings, query convention, counts and
timing. Hit SDFs contain 2D depictions, not aligned candidate poses. Exact scoring
does not feed back into sampling and does not hold up GPU batches during mining.

## Tests

```bash
mvn -pl openchemlib-hyperspace-3d -Dtest=ContinuousDeepSpaceMinerTest test
```

Tests cover deterministic bounded elites, filtering/exhaustion, exact proposal and
accepted limits, recovery with producer read-ahead after injected encoder failure,
committed-shard integrity/uniqueness, stop/resume, time limits, real CPU ONNX CLI
execution, exact finishing, and preservation of a planar fixed-pose query.
