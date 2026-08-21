# Initial implementation note

The existing Hyperspace design is reused rather than replaced:

* `RawSynthonSpace`, `RawSynthonSpaceIO`, and
  `RawSynthonSpaceAssembler` remain the chemistry-space authority.
* `DownsampledSynthonSpace` and `SynthonSetAccessor` supply representative
  synthons independently of product fingerprints.
* `SkelSpheresNeighborSampler` and the `NeighborSampler` interface are the
  local one-position substitution source.
* `SynthonAssembler.assembleSynthons_faster()` is used by the new batch
  scorer.
* `LocalOptimizationRequest` already defines beam width, neighbor pool and
  sample counts, per-position caps, rounds, patience, tolerance, and score
  threshold; those controls should be adapted into the later coordinator.

The current `AssemblyScorer`/`LocalBeamOptimizer` path is synchronous and
molecule-at-a-time. It was intentionally left unchanged. Hyperspace3D adds a
parallel batch boundary beside it, allowing multiple beams to share CPU
assembly workers and a bounded GPU queue in the next milestone.

The sampled-product index is deliberately a small versioned sequential
float32 format. Its manifest already binds rawspace, downsampled space, model
bundle, sampling, filters, dtype, count, and shards, leaving room for
float16/memory mapping without committing those optimizations prematurely.
