# RawSynthonSpace Format

`RawSynthonSpace` is the lightweight JSON format used by Hyperspace before or instead of materializing heavier serialized `SynthonSpace` data files. It stores synthon reactions, synthon sets, synthons, and optional metadata in plain JSON.

A `.rawspace` file is UTF-8 JSON. A `.rawspace.gz` file is the same JSON stream gzip-compressed; compression is selected only by the `.gz` filename suffix in `RawSynthonSpaceIO`.

The canonical Java model and reader/writer live in `com.idorsia.research.chem.hyperspace.rawspace`:

- `RawSynthonSpace`
- `RawSynthon`
- `RawSynthonSpaceIO`

## Concepts

- **Space**: a named collection of synthon reactions.
- **Reaction**: one combinatorial reaction, identified by `reactionId`.
- **Synthon set**: one reactant position in a reaction, keyed by an integer index such as `0`, `1`, or `2`.
- **Synthon** or **fragment**: one building-block candidate inside a synthon set.
- **Full rawspace**: a rawspace containing all imported synthons.
- **Downsampled rawspace**: a separate rawspace whose normal `fragmentSets` contain only retained representative synthons and whose metadata contains `space.role=downsampled`.

Old files with embedded `downsampledFragmentSets` are obsolete and should be regenerated as separate downsampled rawspace files.

## Minimal JSON

This is the smallest useful shape. Optional maps/lists may be absent when reading; the current writer emits them as empty structures.

```json
{
  "name": "toy_space",
  "version": "1.0",
  "metadata": {
    "space.role": "full"
  },
  "reactions": [
    {
      "reactionId": "amide_coupling",
      "fragmentSets": {
        "0": [
          {
            "reactionId": "amide_coupling",
            "fragIndex": 0,
            "fragmentId": "amine-001",
            "idcode": "gC...",
            "connectors": "AQ=="
          }
        ],
        "1": [
          {
            "reactionId": "amide_coupling",
            "fragIndex": 1,
            "fragmentId": "acid-001",
            "idcode": "fH...",
            "connectors": "AQ=="
          }
        ]
      }
    }
  ]
}
```

## Complete JSON Shape

The writer emits this full reaction shape:

```json
{
  "name": "supplier_merged_2s",
  "version": "1.0",
  "metadata": {
    "space.role": "full",
    "source.format": "merged-rawspace",
    "descriptor.tags": "FragFp"
  },
  "reactions": [
    {
      "reactionId": "amide_coupling",
      "fragmentSets": {
        "0": [
          {
            "reactionId": "amide_coupling",
            "fragIndex": 0,
            "fragmentId": "amine-001",
            "idcode": "gC...",
            "connectors": "AQ=="
          }
        ]
      },
      "exampleScaffolds": [],
      "partialAssemblies": {},
      "representativeCompounds": [],
      "descriptors": {},
      "reactionMetadata": {
        "source.spaceName": "enamine_real_2024_02",
        "source.originalReactionId": "amide_coupling"
      },
      "fragmentAttributes": {
        "amine-001": {
          "price": "12.50",
          "descriptor.FragFp": "..."
        }
      }
    }
  ]
}
```

## Top-Level Fields

| Field | Type | Required for writing | Reader behavior | Meaning |
| --- | --- | --- | --- | --- |
| `name` | string | Yes | Required by `RawSynthonSpace.builder` | Human-readable space name. |
| `version` | string | Yes | Required by current builder path | Producer/format version; builders default to `1.0`. |
| `metadata` | object string -> string | No | Missing becomes empty | Space-level metadata. |
| `reactions` | array | No | Missing becomes empty | Reaction records. |

All metadata values are strings. Use parseable scalar strings for numbers and booleans, for example `"24"`, `"0.75"`, or `"true"`.

## Reaction Fields

| Field | Type | Required for useful reaction | Reader behavior | Meaning |
| --- | --- | --- | --- | --- |
| `reactionId` | string | Yes | Used as the reaction key | Reaction identifier. Must match the containing fragments' reaction in normal writer output. |
| `fragmentSets` | object string -> array | Yes | Missing means no synthons | Map from stringified integer set index to fragment records. |
| `exampleScaffolds` | array of strings | No | Missing becomes empty | Optional assembled scaffold IDCodes for examples/inspection. |
| `partialAssemblies` | object string -> array of strings | No | Missing becomes empty | Optional partial assembly IDCodes keyed by missing set index. Keys are stringified integers. |
| `representativeCompounds` | array of strings | No | Missing becomes empty | Optional assembled representative compound IDCodes. |
| `descriptors` | object string -> string | No | Missing becomes empty | Reaction-level descriptors or provenance values. |
| `reactionMetadata` | object string -> string | No | Missing becomes empty | Reaction-level metadata such as merge provenance. |
| `fragmentAttributes` | object string -> object string -> string | No | Missing becomes empty | Fragment-level attributes keyed by `fragmentId`. |

`fragmentSets` and `partialAssemblies` are maps with integer keys in Java, but JSON object keys are strings. A valid producer should write keys like `"0"`, `"1"`, and `"2"`.

## Fragment Fields

| Field | Type | Required | Meaning |
| --- | --- | --- | --- |
| `reactionId` | string | Yes | Reaction this synthon belongs to. |
| `fragIndex` | integer | Yes | Synthon-set index, matching the enclosing `fragmentSets` key. |
| `fragmentId` | string | Yes | Stable fragment/building-block identifier inside the reaction. |
| `idcode` | string | Yes | OpenChemLib IDCode for the synthon structure. |
| `connectors` | string or null | No | Base64-encoded connector bitset. `null` or missing means no connector bits. |

`connectors` is encoded as `Base64.getEncoder().encodeToString(bitSet.toByteArray())` and decoded as `BitSet.valueOf(Base64.getDecoder().decode(connectors))`. External writers that do not know connector bits may write `null`, but files intended for building/searching should preserve connector information generated from the synthon structures.

## Producer Guidance

- Use OpenChemLib IDCodes, not SMILES, in `idcode` fields.
- Keep `fragmentId` stable and unique within a reaction; downstream attributes are keyed by this value.
- Store all metadata, descriptor values, and fragment attributes as strings.
- Keep `reactionId` and `fragIndex` in each fragment consistent with the enclosing reaction and synthon-set key.
- Use stringified integer keys for `fragmentSets` and `partialAssemblies`.
- Prefer explicit `metadata["space.role"]="full"` for full spaces, although unmarked full spaces remain accepted.
- For merged spaces, preserve reaction provenance in `reactionMetadata` using `source.spaceName`, `source.spacePath`, and `source.originalReactionId`.

## Full vs Downsampled Spaces

A full space is either unmarked or explicitly marked:

```json
"metadata": {
  "space.role": "full"
}
```

Unmarked full spaces are accepted for compatibility with imported rawspaces.

A downsampled space is a normal rawspace with fewer fragments in `fragmentSets`:

```json
"metadata": {
  "space.role": "downsampled",
  "downsampling.algorithm": "SkelSpheresKCentersRaw",
  "downsampling.maxCenters": "1000",
  "downsampling.minSimilarity": "0.75",
  "downsampling.seed": "13",
  "downsampling.enforceConnectorEquivalence": "true",
  "downsampling.sizeCapScale": "0.0",
  "downsampling.sizeCapOffset": "0.0",
  "downsampling.includeClusterMembers": "false"
}
```

The downsampled file does not carry the full input sets. Workflows that need both use two files: one full rawspace and one downsampled rawspace.

## Space-Level Metadata

Common keys:

| Key | Meaning |
| --- | --- |
| `space.role` | `full` or `downsampled`. Full spaces may be unmarked. |
| `source.format` | Import source format, e.g. `enamine-tsv`, `xtalpi-csv`, `csv-per-reaction`, or `molecule-one-zip`. |
| `source.file` | Original source file for single-file imports. |
| `source.directory` | Original source directory for per-reaction CSV imports. |
| `source.zipEntry` | Internal zip entry used for synthon tables. |
| `source.reactionZipEntry` | Internal zip entry used for reaction metadata tables. |
| `descriptor.shortName` | Primary descriptor used while importing/building, e.g. `FragFp`. |
| `descriptor.bits` | Descriptor bit count, when applicable. |
| `descriptor.tags` | Comma-separated descriptors available as fragment attributes or compatible with this rawspace. |
| `parser.mode` | Parser/build mode used by import. |
| `parser.maxSynthonSets` | Maximum synthon-set count accepted during import. |
| `parser.threads` | Import thread count. |
| `parser.smilesColumn` | CSV column used for synthon SMILES. |
| `parser.idColumn` | CSV column used for fragment IDs. |
| `parser.priceColumn` | CSV price column, if configured. |
| `parser.priceAttributeKey` | Fragment attribute key used for imported prices. |
| `parser.synthonSetColumn` | CSV column used for synthon set index, if configured. |
| `parser.defaultSynthonSet` | Default synthon set index when no set column is present. |
| `parser.includeAllFiles` | Whether CSV import included all files in the input directory. |

Downsampling keys:

| Key | Meaning |
| --- | --- |
| `downsampling.algorithm` | Downsampler implementation name. |
| `downsampling.maxCenters` | Requested maximum representatives per synthon set (`0` means unlimited unless size cap applies). |
| `downsampling.minSimilarity` | Similarity threshold used when assigning representatives. |
| `downsampling.seed` | Random seed used for shuffling/sampling. |
| `downsampling.enforceConnectorEquivalence` | Whether connector-equivalent synthons were clustered together only. |
| `downsampling.sizeCapScale` | Scale in `ceil(scale * sqrt(n) + offset)`. |
| `downsampling.sizeCapOffset` | Offset in `ceil(scale * sqrt(n) + offset)`. |
| `downsampling.includeClusterMembers` | Whether cluster-member assignments were requested. Usually `false` for rawspace output. |

Merge keys:

| Key | Meaning |
| --- | --- |
| `merge.mode` | `single`, `split_2s`, or `split_3s`. |
| `merge.sourceCount` | Number of rawspace inputs merged. |
| `merge.reactionCount` | Number of reactions written to this output. |
| `merge.source.N.name` | Name of source `N`. |
| `merge.source.N.path` | Path of source `N`. |

## Reaction Metadata

`reactionMetadata` is for metadata scoped to one reaction. Merge currently writes:

| Key | Meaning |
| --- | --- |
| `source.spaceName` | Source rawspace name for this reaction. |
| `source.spacePath` | Source rawspace path for this reaction. |
| `source.originalReactionId` | Reaction ID before merge/rename. |
| `source.synthonSetCount` | Number of non-empty synthon sets in the source reaction. |
| `source.mergedReactionId` | New reaction ID when duplicate handling renamed it. |

The prefix can be changed in `SynthonSpaceMergeCLI` with `--sourceMetadataPrefix`.

Supplier reaction metadata from zip/table imports may also be stored here, for example `supplier.reaction.components`, `supplier.reaction.smarts`, `supplier.reaction.product`, or `supplier.reaction.R1`.

## Fragment Attributes

`fragmentAttributes` is for metadata scoped to a single `fragmentId`.

Common keys:

| Key | Meaning |
| --- | --- |
| `price` | Default imported price attribute. |
| Custom price key | CSV import can set this with `--priceAttribute`. |
| `descriptor.<shortName>` | Precomputed descriptor value attached by `RawSynthonSpaceProcessorCLI`, e.g. `descriptor.FragFp` or `descriptor.SkelSpheres`. |

Only retained representatives keep their fragment attributes when creating a downsampled rawspace from a raw input.

## Compatibility Notes

- Readers ignore missing optional arrays/maps by treating them as empty.
- Full rawspaces may be unmarked for compatibility.
- Downsampled rawspaces should be marked with `space.role=downsampled`; seed finding and continuous screening expect this marker.
- `downsampledFragmentSets`, top-level `downsamplingAlgorithm`, and top-level `downsamplingRequest` are obsolete and should not be emitted.
- JSON is pretty-printed by the current writer, but readers do not require formatting.
- `.rawspace.gz` is gzip-compressed JSON, not a different binary format.
