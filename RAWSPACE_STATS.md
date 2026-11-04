# Rawspace Statistics Export

`RawSynthonSpaceStatsCLI` creates a small report bundle for a `.rawspace` or `.rawspace.gz` file. It is intended for merged supplier spaces where exact full enumeration is impossible, but exact reaction sizes and sampled product-property summaries are useful.

For the rawspace JSON/gzip layout itself, see [`RAWSPACE_FORMAT.md`](RAWSPACE_FORMAT.md).

## Command

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceStatsCLI \
  --rawIn /mnt/cadd/VS/synthonspaces/merged/supplier_merged_2s.rawspace.gz \
  --outDir /mnt/cadd/VS/synthonspaces/stats/supplier_merged_2s_stats \
  --examplesPerReaction 10 \
  --productSamplesPerReaction 100 \
  --threads 12 \
  --seed 13
```

The same command is also available through the shaded launcher as `RAWSPACESTATS` with the same options.

Defaults:

- `--examplesPerReaction 10`
- `--productSamplesPerReaction 100`
- `--seed 13`
- `--threads <available processors>`
- `--maxReactions 0`, meaning all reactions

## Output Files

The output is a directory bundle:

- `summary.md`: human-readable overview and file descriptions.
- `space_summary.tsv`: one-row space totals plus rawspace metadata.
- `reaction_stats.tsv`: one row per reaction with source provenance, synthon-set sizes, exact combinatorial product count, synthon atom statistics, and sampled product atom/rotatable-bond statistics.
- `source_stats.tsv`: totals grouped by reaction metadata `source.spaceName` when present.
- `synthon_set_stats.tsv`: one row per reaction/set with synthon counts and atom statistics.
- `example_products.tsv`: deterministic random assembled examples, including `Structure [idcode]` and chosen fragment IDs.

## Notes

- Product counts are exact and computed from synthon-set sizes.
- Product atom and rotatable-bond distributions are sampled; full generated spaces are not enumerated.
- Example products are sampled deterministically from `--seed` and reaction id.
- The Synthon Space Explorer GUI exposes the same exporter under `Tools -> Export Rawspace Statistics...`.
