# SynthonSpaceMergeCLI Quick Workflow

`SynthonSpaceMergeCLI` merges multiple `RawSynthonSpace` files into:

1. one combined output, or
2. two outputs split by reaction synthon-set count (`2s` and `3s`).

It also writes provenance metadata per reaction so merged spaces remain traceable.

## Build

```bash
cd /home/liphath1/dev_3d/openchemlib-hyperspace
mvn -pl openchemlib-hyperspace-cli -am -DskipTests package
```

CLI jar:
`openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar`

## Input Options

- `--rawIn`: repeatable input paths, or comma-separated paths
- `--inputDirectory`: include all `*.rawspace.gz` files in one folder
- `--manifest`: text file with one input path per line
  - empty lines ignored
  - lines starting with `#` ignored
  - relative paths resolved relative to manifest file location

You can combine `--rawIn`, `--inputDirectory`, and `--manifest`.

## Mode A: Single Merged Output

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.SynthonSpaceMergeCLI \
  --rawIn /data/vendorA.rawspace.gz \
  --rawIn /data/vendorB.rawspace.gz \
  --synthonListOut /data/all_synthons.tsv \
  --out /data/merged_all.rawspace.gz \
  --spaceName merged_commercial \
  --version 1.0 \
  --duplicateStrategy rename
```

## Mode B: Split into 2s and 3s Outputs

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.SynthonSpaceMergeCLI \
  --inputDirectory /data/vendor_spaces \
  --splitBySetCount \
  --out2s /data/merged_2s.rawspace.gz \
  --out3s /data/merged_3s.rawspace.gz \
  --spaceName merged_commercial \
  --version 1.0 \
  --duplicateStrategy rename
```

In split mode:
- reactions with exactly 2 non-empty sets go to `--out2s`
- reactions with exactly 3 non-empty sets go to `--out3s`
- others are reported as skipped

## Duplicate Reaction IDs

`--duplicateStrategy`:
- `rename` (default): keep both reactions, rename conflicting IDs
- `fail`: stop on first duplicate reaction ID

## Provenance Metadata per Reaction

By default, metadata keys use prefix `source` (customizable with `--sourceMetadataPrefix`):

- `source.spaceName`
- `source.spacePath`
- `source.originalReactionId`
- `source.synthonSetCount`
- `source.mergedReactionId` (only when a rename happened)

This metadata is stored in reaction-level metadata (`Map<String,String>`).

## Optional Synthon List Export (TSV)

Use `--synthonListOut /path/file.tsv` to export all unique synthons found in loaded input spaces.

Columns:
- `Structure[idcode]`
- `NumAtoms`
- `NumSpaces`
- `NumSynthonSets`

## Manifest Example

`/data/rawspaces_manifest.txt`

```text
# Commercial raw spaces
/data/vendorA.rawspace.gz
/data/vendorB.rawspace.gz
relative/path/vendorC.rawspace
```
