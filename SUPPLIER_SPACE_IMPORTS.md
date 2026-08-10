# Supplier Rawspace Imports

This document records the supplier files checked under `supplier_data` and the commands used to import them as separate `RawSynthonSpace` files.

The helper script [`../import_supplier_spaces.sh`](../import_supplier_spaces.sh) runs all imports listed here.

```bash
HYPERSPACE_JAR=/path/to/openchemlib-hyperspace-cli.jar \
SUPPLIER_ROOT=/home/liphath1/dev_3d/supplier_data \
OUT_DIR=/mnt/cadd/VS/synthonspaces/imported \
THREADS=24 \
MODE=FragFp \
MAX_SETS=0 \
../import_supplier_spaces.sh
```

Defaults:

- `SUPPLIER_ROOT=/home/liphath1/dev_3d/supplier_data`
- `OUT_DIR=/mnt/cadd/VS/synthonspaces/imported`
- `REPO_DIR=/home/liphath1/dev_3d/openchemlib-hyperspace`
- `HYPERSPACE_JAR=$REPO_DIR/openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar`
- `MODE=FragFp`
- `THREADS=24`
- `MAX_SETS=0`, meaning reactions are not filtered by number of synthon sets before connector validation.
- Existing outputs are skipped. Set `OVERWRITE=1` to regenerate them.

## Checked Files

| Supplier | Source file | Import support | Output |
| --- | --- | --- | --- |
| Enamine | `Enamine/2024-02_REAL_synthons.smi` | Supported as Enamine-style TSV. | `enamine_real_2024_02.rawspace.gz` |
| ChemSpace | `ChemSpace/Freedom Space 5.0/2026-01_Freedom_5/2026-01_Freedom_5_synthons.txt` | Supported as Enamine-style TSV. | `chemspace_freedom_5_2026_01.rawspace.gz` |
| Xtalpi | `Xtalpi/VAST_synthon_2026_H1.csv` | Supported via `--format xtalpi`; `synton_role` values such as `synton_1` are mapped to set index `1`. | `xtalpi_vast_2026_h1.rawspace.gz` |
| Molecule.One | `Molecule.One/1-step D2B space.zip` with `synthons.txt` and `reactions.txt` | Supported as a zip table import. Reaction metadata from `reactions.txt` is stored per reaction. | `molecule_one_1step_d2b.rawspace.gz` |
| Molecule.One | `Molecule.One/1-step vials space.zip` with `synthons.txt` and `reactions.txt` | Supported as a zip table import. Reaction metadata from `reactions.txt` is stored per reaction. | `molecule_one_1step_vials.rawspace.gz` |
| Molecule.One | `Molecule.One/2+ step HTE space 2.zip` with `2+ step HTE space/synthons.txt` and `2+ step HTE space/reactions.txt` | Supported as a zip table import. Reaction metadata from `reactions.txt` is stored per reaction. | `molecule_one_2plus_step_hte.rawspace.gz` |
| Synple | `Synple/Synple_eMolecules_SMARTS_BBs_Jan26/Synthon Dataset/FileFor_RDKIt_SynthonSearch/SYNPLE_eXplore_Synthon_Library_Jan26.csv` | Supported as a consolidated CSV through `--format enamine`. | `synple_explore_jan26.rawspace.gz` |
| Toy/Test | `Test/idorsia_toy_space_a.txt` | Supported as Enamine-style TSV. | `idorsia_toy_space_a.rawspace.gz` |

All generated files are marked with `metadata["space.role"]="full"` and supplier metadata such as `source.supplier`, `source.collection`, and `source.release` where applicable.

## Individual Commands

The commands below use shell variables to keep paths readable:

```bash
JAR=/path/to/openchemlib-hyperspace-cli.jar
CLI=com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceImportCLI
SUPPLIER_ROOT=/home/liphath1/dev_3d/supplier_data
OUT_DIR=/mnt/cadd/VS/synthonspaces/imported
COMMON=(--mode FragFp --threads 24 --maxSets 0 --metadata space.role=full)
```

Enamine REAL:

```bash
java -cp "$JAR" "$CLI" \
  --format enamine \
  --input "$SUPPLIER_ROOT/Enamine/2024-02_REAL_synthons.smi" \
  --spaceName enamine_real_2024_02 \
  --sourceFormat enamine-real-tsv \
  --metadata source.supplier=Enamine \
  --metadata source.collection=REAL \
  --metadata source.release=2024-02 \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/enamine_real_2024_02.rawspace.gz"
```

ChemSpace Freedom 5.0:

```bash
java -cp "$JAR" "$CLI" \
  --format enamine \
  --input "$SUPPLIER_ROOT/ChemSpace/Freedom Space 5.0/2026-01_Freedom_5/2026-01_Freedom_5_synthons.txt" \
  --spaceName chemspace_freedom_5_2026_01 \
  --sourceFormat chemspace-freedom-tsv \
  --metadata source.supplier=ChemSpace \
  --metadata "source.collection=Freedom Space 5.0" \
  --metadata source.release=2026-01 \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/chemspace_freedom_5_2026_01.rawspace.gz"
```

Xtalpi VAST:

```bash
java -cp "$JAR" "$CLI" \
  --format xtalpi \
  --input "$SUPPLIER_ROOT/Xtalpi/VAST_synthon_2026_H1.csv" \
  --spaceName xtalpi_vast_2026_h1 \
  --metadata source.supplier=Xtalpi \
  --metadata source.collection=VAST \
  --metadata "source.release=2026 H1" \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/xtalpi_vast_2026_h1.rawspace.gz"
```

Molecule.One 1-step D2B:

```bash
java -cp "$JAR" "$CLI" \
  --format enamine \
  --input "$SUPPLIER_ROOT/Molecule.One/1-step D2B space.zip" \
  --zipEntry synthons.txt \
  --reactionZipEntry reactions.txt \
  --spaceName molecule_one_1step_d2b \
  --sourceFormat molecule-one-zip \
  --metadata source.supplier=Molecule.One \
  --metadata "source.collection=1-step D2B" \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/molecule_one_1step_d2b.rawspace.gz"
```

Molecule.One 1-step vials:

```bash
java -cp "$JAR" "$CLI" \
  --format enamine \
  --input "$SUPPLIER_ROOT/Molecule.One/1-step vials space.zip" \
  --zipEntry synthons.txt \
  --reactionZipEntry reactions.txt \
  --spaceName molecule_one_1step_vials \
  --sourceFormat molecule-one-zip \
  --metadata source.supplier=Molecule.One \
  --metadata "source.collection=1-step vials" \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/molecule_one_1step_vials.rawspace.gz"
```

Molecule.One 2+ step HTE:

```bash
java -cp "$JAR" "$CLI" \
  --format enamine \
  --input "$SUPPLIER_ROOT/Molecule.One/2+ step HTE space 2.zip" \
  --zipEntry "2+ step HTE space/synthons.txt" \
  --reactionZipEntry "2+ step HTE space/reactions.txt" \
  --spaceName molecule_one_2plus_step_hte \
  --sourceFormat molecule-one-zip \
  --metadata source.supplier=Molecule.One \
  --metadata "source.collection=2+ step HTE space" \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/molecule_one_2plus_step_hte.rawspace.gz"
```

Synple eXplore:

```bash
java -cp "$JAR" "$CLI" \
  --format enamine \
  --input "$SUPPLIER_ROOT/Synple/Synple_eMolecules_SMARTS_BBs_Jan26/Synthon Dataset/FileFor_RDKIt_SynthonSearch/SYNPLE_eXplore_Synthon_Library_Jan26.csv" \
  --spaceName synple_explore_jan26 \
  --sourceFormat synple-csv \
  --metadata source.supplier=Synple \
  --metadata source.collection=eXplore \
  --metadata source.release=Jan26 \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/synple_explore_jan26.rawspace.gz"
```

Toy/Test:

```bash
java -cp "$JAR" "$CLI" \
  --format enamine \
  --input "$SUPPLIER_ROOT/Test/idorsia_toy_space_a.txt" \
  --spaceName idorsia_toy_space_a \
  --sourceFormat toy-tsv \
  --metadata source.supplier=Idorsia \
  --metadata "source.collection=Toy/Test" \
  "${COMMON[@]}" \
  --rawOut "$OUT_DIR/idorsia_toy_space_a.rawspace.gz"
```

## Merge Imported Supplier Spaces

After importing the supplier spaces, merge the seven supplier rawspaces into split 2-set and 3-set files:

```bash
../merge_supplier_spaces.sh
```

Defaults:

- `IMPORT_DIR=/mnt/cadd/VS/synthonspaces/imported`
- `MERGED_DIR=/mnt/cadd/VS/synthonspaces/merged`
- `REPO_DIR=/home/liphath1/dev_3d/openchemlib-hyperspace`
- `HYPERSPACE_JAR=$REPO_DIR/openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar`
- Existing merged outputs are skipped when both split files already exist. Set `OVERWRITE=1` to regenerate them.

Outputs:

- `$MERGED_DIR/supplier_merged_2s.rawspace.gz`
- `$MERGED_DIR/supplier_merged_3s.rawspace.gz`
- `$MERGED_DIR/supplier_merged_synthons.tsv`

The merge script intentionally excludes `idorsia_toy_space_a.rawspace.gz`.

## Clean Merged Supplier Spaces

After merging, run LigPrep-like synthon cleaning on the split merged spaces:

```bash
../clean_2s_merged_space.sh
../clean_3s_merged_space.sh
```

Defaults:

- `MERGED_DIR=/mnt/cadd/VS/synthonspaces/merged`
- `CLEANED_DIR=/mnt/cadd/VS/synthonspaces/cleaned_merged`
- `NEON_DIR=/home/liphath1/dev_3d/neon`
- `NEON_JAR=$NEON_DIR/neon-drugforge-tools/target/neon-drugforge-tools.jar`
- `THREADS=12`
- `PH=7.4`, `MAX_IONS=6`, `MAX_STEREO=4`, `CONNECTIVITY_FALLBACK=keep_original`
- `MAX_SYNTHON_ATOMS=36`, passed to `--max-synthon-atoms` to drop oversized individual synthons before cleaning.

Outputs:

- `$CLEANED_DIR/supplier_merged_2s_cleaned.rawspace.gz`
- `$CLEANED_DIR/supplier_merged_3s_cleaned.rawspace.gz`

Set `OVERWRITE=1` to regenerate existing cleaned outputs.

## Notes

- The importers validate connector consistency per reaction. Invalid reactions are skipped with a diagnostic message.
- `--maxSets 0` keeps all reactions before validation. Use `--maxSets 3` if a workflow should intentionally retain only 2- and 3-component reactions.
- Molecule.One `reactions.txt` columns are imported as reaction metadata with keys such as `supplier.reaction.components`, `supplier.reaction.smarts`, `supplier.reaction.product`, and `supplier.reaction.R1`.
- The Synple consolidated CSV uses the same column shape as the single-file Enamine CSV path: `SMILES`, `synthon_id`, `synthon#`, and `reaction_id`.
