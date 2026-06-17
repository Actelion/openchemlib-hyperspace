# RawSynthonSpace Format

`RawSynthonSpace` is the lightweight JSON representation used before or instead of materializing heavier `SynthonSpace` data files. It stores reactions, synthon sets, synthons, and optional metadata in a format that can be inspected as plain JSON (`.rawspace`) or compressed with gzip (`.rawspace.gz`).

The format is intentionally close to the Java model in `com.idorsia.research.chem.hyperspace.rawspace`.

## Concepts

- A **space** is a collection of synthon reactions.
- A **reaction** is identified by `reactionId` and contains one or more synthon sets.
- A **synthon set** is one position in a combinatorial reaction, keyed by an integer index such as `0`, `1`, or `2`.
- A **fragment** or **synthon** is one building block candidate in a synthon set.
- A **full rawspace** contains all imported synthons.
- A **downsampled rawspace** is a separate rawspace file whose normal `fragmentSets` contain only retained representatives. It is marked with `metadata["space.role"]="downsampled"`.

Old rawspace files with embedded `downsampledFragmentSets` are no longer a supported target format. Regenerate them as separate downsampled rawspace files.

## File Structure

Top-level JSON object:

```json
{
  "name": "toy",
  "version": "1.0",
  "metadata": {
    "source.format": "enamine-tsv",
    "descriptor.tags": "FragFp"
  },
  "reactions": [
    {
      "reactionId": "benzoimidazole_b-8",
      "fragmentSets": {
        "0": [
          {
            "reactionId": "benzoimidazole_b-8",
            "fragIndex": 0,
            "fragmentId": "BB-001",
            "idcode": "gC...",
            "connectors": "AQ=="
          }
        ]
      },
      "exampleScaffolds": [],
      "partialAssemblies": {},
      "representativeCompounds": [],
      "descriptors": {},
      "reactionMetadata": {},
      "fragmentAttributes": {
        "BB-001": {
          "price": "12.50",
          "descriptor.FragFp": "..."
        }
      }
    }
  ]
}
```

Fields:

- `name`: Human-readable space name.
- `version`: Rawspace format or producer version. Current builders default to `1.0`.
- `metadata`: Space-level `Map<String,String>`.
- `reactions`: List of reaction records.
- `fragmentSets`: Map from synthon set index, serialized as a string key, to fragment records.
- `idcode`: OpenChemLib IDCode for the synthon.
- `connectors`: Base64 encoding of the connector bitset. It can be `null` or absent-equivalent when the bitset is empty.
- `exampleScaffolds`, `partialAssemblies`, `representativeCompounds`: Optional reaction helper data used by import/build workflows.
- `descriptors`: Reaction-level descriptor/provenance map.
- `reactionMetadata`: Reaction-level `Map<String,String>`.
- `fragmentAttributes`: Map from `fragmentId` to arbitrary string attributes.

All metadata and attributes are stored as strings. Use parseable scalar strings for numbers and booleans.

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
  "downsampling.maxCenters": "8",
  "downsampling.minSimilarity": "0.75",
  "downsampling.seed": "13",
  "downsampling.enforceConnectorEquivalence": "true",
  "downsampling.sizeCapScale": "0.0",
  "downsampling.sizeCapOffset": "0.0"
}
```

The downsampled file does not carry the full input sets. Workflows that need both use two files: one full rawspace and one downsampled rawspace.

## Space-Level Metadata

Common keys:

| Key | Meaning |
| --- | --- |
| `space.role` | `full` or `downsampled`. Full spaces may be unmarked. |
| `source.format` | Import source format, e.g. `enamine-tsv` or `csv-per-reaction`. |
| `source.file` | Original source file for single-file imports. |
| `source.directory` | Original source directory for per-reaction CSV imports. |
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

## Fragment Attributes

`fragmentAttributes` is for metadata scoped to a single fragment ID.

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
