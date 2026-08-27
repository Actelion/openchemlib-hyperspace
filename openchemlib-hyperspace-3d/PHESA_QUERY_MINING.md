# Query-centric PheSA pair mining

Hyperspace3D can mine difficult learned-PheSA retrieval examples without
re-encoding a supplier library. The workflow keeps the canonical 50D graph
encoder and its normalized 128D molecular fingerprints frozen. Java performs
rank mining and exact OCL labeling; Python trains only symmetric comparison
heads.

## Artifact flow

1. `PheSAMiningUniverseBuilderCLI` extracts a deterministic bottom-hash sample
   from an existing `MoleculeFingerprintDataSource`. The intended first run is
   one million molecules with seed 17. Exclusion files contain one canonical
   SMILES per line and must include the locked benchmarks and the 6k query
   registry.
2. `PheSAQueryPairMiningCLI` accepts `query_uid`, `split`, and `smiles` TSV
   columns. It screens total, shape, and pharmacophore independently; samples
   six percentile ranges plus random and bounded flexibility, high-sp3,
   stereochemical, and size-tail reservoirs; deduplicates the union; and
   calculates exact OCL PheSA using eight conformers and `wPP=0.5`.

For a single frozen-baseline enrichment campaign, set `selectionPolicy` to
`fixed_baseline_1600_v1`. Each full one-million-candidate screen then selects
exactly 1,600 unique candidates:

- 400 from absolute ranks 1-400;
- 400 from the remainder of the top 1%;
- 400 from the 1-5% rank interval;
- 400 from the complete-library random background.

Each ranked group reserves 200 candidates from total rank, 100 from shape rank,
and 100 from pharmacophore rank. Earlier selections are excluded before later
channel quotas are filled, so the strata are disjoint.
3. Deepspace7 memory maps the resulting descriptors and fixed-width query
   shards. Mismatch labels are derived dynamically and are never baked into
   the cache.

Copy `phesa-mining-universe.example.json` and run:

```bash
java -cp 'openchemlib-hyperspace-3d/target/classes:openchemlib-hyperspace-3d/target/lib/*' \
  com.idorsia.research.chem.hyperspace3d.cli.PheSAMiningUniverseBuilderCLI \
  --config openchemlib-hyperspace-3d/phesa-mining-universe.example.json
```

Then copy `phesa-query-pair-mining.example.json` and run the query miner in the
same way. CUDA selection is explicit and never falls back to CPU.

## Cache format

The final directory contains `manifest.json`, `molecules.tsv.gz`, a raw
little-endian FP16 `[N,128]` `descriptors.f16`, `query-index.tsv`, and one
`queries/qNNNNNN.bin` file per query. Query files have a 64-byte header and
72-byte records containing candidate ID, three baseline predictions, three
full-screen ranks, three exact labels, selection mask, mining round, and exact
status. Every file and all model/source inputs are hash-bound in the manifest.

Exact failures remain explicit rows with stable status codes and NaN labels.
Completed outputs are immutable. A rerun must use a new output directory.

The first format deliberately favors auditability and streaming access over
Parquet dependencies. Exact labeling uses thread-confined worker handlers and checksum-committed
per-query work files. Interrupted runs validate and resume completed queries;
orphaned partial files are regenerated. The selection
mask already reserves a model-disagreement channel for a later multi-checkpoint
mining round.

## Initial experiment boundary

Use split A for mining-head training and split B for model selection. Split C
and the fixed external Enamine benchmark remain locked until a comparator has
been selected. After the 16/4 and 192/48 gates pass, the fixed-baseline
production campaign uses 800 diverse A queries and 200 diverse B queries. It
therefore produces exactly 1.6 million selected pair observations before
explicit exact-label failures are filtered during training. The same model
bundle, seed, selection policy, eight-conformer OCL settings, and `wPP=0.5`
remain fixed for all 1,000 queries; active-learning rescoring is intentionally
outside this dataset round.
