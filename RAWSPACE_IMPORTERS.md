# RawSynthonSpace Importers

This page documents the import paths that create `RawSynthonSpace` JSON files. For the rawspace JSON layout itself, see [`RAWSPACE_FORMAT.md`](RAWSPACE_FORMAT.md).

## Import CLI

The main entry point is:

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceImportCLI \
  --format enamine \
  --input vendor.tsv \
  --spaceName my_space \
  --mode FragFp \
  --rawOut my_space.rawspace.gz
```

Supported `--format` values:

- `enamine`: one source file containing all reactions.
- `xtalpi`: one CSV source file using `SMILES`, `synton_id`, `synton_role`, and `reaction_id`.
- `csv`: directory of per-reaction CSV files; reaction identity is derived from filenames.

Both formats can write:

- `--rawOut`: required `RawSynthonSpace` JSON output (`.rawspace` or `.rawspace.gz`).
- `--synthonOut`: optional serialized `SynthonSpace`.
- `--similarityOut`: optional serialized `SynthonSimilaritySpace3`.

Descriptor modes accepted by the importer are:

- `FragFp`
- `PathFp`
- `mode_pfp` or `pfp`
- `mode_ffp` or `ffp`
- `mode_ppcore` or `ppc`

`--descriptorTag` adds entries to `metadata["descriptor.tags"]`. It is repeatable and also accepts comma-separated values. Multiple descriptor tags are only supported when producing raw output only.

`--metadata key=value` adds arbitrary rawspace metadata and is repeatable. Common supplier imports use keys such as `space.role=full`, `source.supplier`, `source.collection`, and `source.release`.

## Enamine TSV Import

Use:

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceImportCLI \
  --format enamine \
  --input idorsia_toy_space_a.txt \
  --spaceName toy \
  --mode FragFp \
  --threads 4 \
  --maxSets 3 \
  --rawOut toy.rawspace.gz
```

The parser skips the first line as a header, then expects tab-separated rows:

```text
SMILES<TAB>fragmentId<TAB>synthonSetIndex<TAB>reactionId
```

Example:

```text
[U]c([nH]nc1)c1-c(nn1)c[n]1[Np]    015119586-618282196    0    benzoimidazole_b-8
```

Notes:

- `--maxSets` filters out reactions with more than the configured number of non-empty synthon sets. The default is `3`.
- `--threads` controls parallel parsing batches.
- Import metadata includes `source.format=enamine-tsv`, `source.file`, `parser.mode`, `parser.maxSynthonSets`, `parser.threads`, `descriptor.shortName`, `descriptor.bits`, and `descriptor.tags`.

The Enamine importer can also read a table inside a zip archive:

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceImportCLI \
  --format enamine \
  --input "Molecule.One/1-step D2B space.zip" \
  --zipEntry synthons.txt \
  --reactionZipEntry reactions.txt \
  --spaceName molecule_one_1step_d2b \
  --sourceFormat molecule-one-zip \
  --metadata source.supplier=Molecule.One \
  --metadata "source.collection=1-step D2B" \
  --metadata space.role=full \
  --rawOut molecule_one_1step_d2b.rawspace.gz
```

`--zipEntry` points to the synthon table inside the archive. `--reactionZipEntry` is optional and imports per-reaction supplier metadata when the table contains columns such as `reaction_id`, `components`, `Reaction`, `Product`, `R1`, `R2`, `R3`, and `R4`.

## Unified CSV File Through Enamine Mode

If `--format enamine --input` points to a file whose name ends in `.csv`, the same parser switches to a single-file CSV mode with a header row.

Required columns:

| Column | Meaning |
| --- | --- |
| `SMILES` | Synthon structure as SMILES. |
| `synthon_id` | Fragment identifier. |
| `reaction_id` | Reaction identifier. |
| `synthon#` | Synthon set index. |

Optional column:

| Column | Meaning |
| --- | --- |
| `price` | Stored as fragment attribute `price`. |

This is still selected with `--format enamine`; the switch is based on the `.csv` file extension.

## Xtalpi CSV Import

Use:

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceImportCLI \
  --format xtalpi \
  --input VAST_synthon_2026_H1.csv \
  --spaceName xtalpi_vast_2026_h1 \
  --mode FragFp \
  --threads 4 \
  --rawOut xtalpi_vast_2026_h1.rawspace.gz
```

Default columns:

| Column | Meaning |
| --- | --- |
| `SMILES` | Synthon structure as SMILES. |
| `synton_id` | Fragment identifier. |
| `reaction_id` | Reaction identifier. |
| `synton_role` | Synthon set role, for example `synton_0`, `synton_1`, or `synton_2`. |

The role suffix is converted to the numeric synthon set index. Column names can be overridden with `--smilesColumn`, `--idColumn`, `--reactionColumn`, and `--synthonSetColumn`.

## Per-Reaction CSV Directory Import

Use:

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceImportCLI \
  --format csv \
  --directory /data/synple_or_vendor_csv \
  --spaceName vendor_csv \
  --mode FragFp \
  --smilesColumn SMILES \
  --idColumn bb1_parent_id \
  --priceColumn Price \
  --priceAttribute price \
  --synthonSetColumn SynthonSet \
  --rawOut vendor_csv.rawspace.gz
```

By default, only regular files ending in `.csv` are loaded. The Java API has an `includeAllFiles` option, but the current CLI does not expose it.

Required columns by default:

| CLI option | Default column | Meaning |
| --- | --- | --- |
| `--smilesColumn` | `SMILES` | Synthon structure as SMILES. |
| `--idColumn` | `bb1_parent_id` | Fragment identifier. |

Optional columns:

| CLI option | Default column | Meaning |
| --- | --- | --- |
| `--priceColumn` | `Price` | Price or other scalar metadata. If the column is absent, no price is stored. |
| `--priceAttribute` | `price` | Fragment attribute key used for price values. |
| `--synthonSetColumn` | unset | Synthon set index column. |
| `--synthonDefaultSet` | `0` | Synthon set index when no set column exists or the value is empty. |

The parser is case-insensitive for column-name lookup.

## Per-Reaction Filename Convention

For `--format csv`, the reaction ID comes from the filename:

1. Strip the final extension.
2. Strip a trailing `_synthon`, if present.
3. If the remaining name ends in `_<digits>`, use those digits as the synthon set index and strip the suffix from the reaction ID.
4. Otherwise use `--synthonSetColumn` if provided, then `--synthonDefaultSet`.

Examples:

| Filename | Reaction ID | Detected set index |
| --- | --- | --- |
| `benzoimidazole_b-8.csv` | `benzoimidazole_b-8` | none; use column/default |
| `benzoimidazole_b-8_0.csv` | `benzoimidazole_b-8` | `0` |
| `benzoimidazole_b-8_1_synthon.csv` | `benzoimidazole_b-8` | `1` |

Multiple files that resolve to the same reaction ID are merged. This supports layouts with one file per reaction/set, for example:

```text
amide_coupling_0.csv
amide_coupling_1.csv
amide_coupling_2.csv
```

It also supports one file per reaction if the synthon set index is present as a column:

```text
amide_coupling.csv
```

with rows containing `SynthonSet` values such as `0`, `1`, and `2`.

## Batch Importer

`RawSynthonSpaceBatchImporter` imports many sources from manifest files:

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceBatchImporter \
  --root /data/raw_sources \
  --enamine enamine_manifest.tsv \
  --synple synple_manifest.tsv \
  --out rawspaces \
  --mode FragFp \
  --threads 4 \
  --maxSets 3
```

Manifest format:

```text
spaceName<TAB>inputPath<TAB>optionalDescriptorTag
```

- `--enamine` manifest entries are imported as Enamine TSV/single-file inputs.
- `--synple` manifest entries are imported as per-reaction CSV directories when the input path is a directory. Consolidated Synple CSV files are imported through the single-file table parser.
- Manifest input paths are resolved relative to `--root`.
- Output files are written as `<spaceName>.json.gz` under `--out`.

## Supplier Import Commands

The checked supplier files under `supplier_data`, their supported import modes, output names, and exact commands are documented in [`SUPPLIER_SPACE_IMPORTS.md`](SUPPLIER_SPACE_IMPORTS.md).

## Validation and Skipped Reactions

After parsing, reactions are validated with `SynthonReactionValidator`. Reactions that fail validation are skipped with a message. Typical reasons include unsupported connector patterns or reactions with more synthon sets than `--maxSets` for Enamine imports.

Malformed individual CSV rows are reported and skipped; a completely empty per-reaction CSV file is skipped.
