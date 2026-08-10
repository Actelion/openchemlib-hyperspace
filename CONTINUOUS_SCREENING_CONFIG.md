# Continuous Screening Config

This document is a quick-start reference for the `ContinuousScreeningCLI` JSON config. The workflow samples candidate products from a downsampled rawspace, optionally runs a small optimization on the downsampled space, then runs full local optimization against the full rawspace and writes hits to TSV.

## Command

```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.ContinuousScreeningCLI \
  continuous-screening-config.example.json
```

The config path can also be passed as `--config continuous-screening-config.example.json`.

## Minimal Working Config

```json
{
  "inputs": {
    "rawFull": "hyperspace.rawspace.gz",
    "rawDownsampled": "hyperspace_downsampled.rawspace.gz"
  },
  "query": {
    "phesaFile": "query.phesa"
  },
  "sampling": {
    "attemptsPerReaction": 20,
    "minAtoms": 1,
    "maxAtoms": 0,
    "maxRotatableBonds": 15,
    "minSimilarity": 0.65
  },
  "microOptimization": {
    "enabled": false,
    "request": {
      "beamSize": 5,
      "neighborPoolSize": 4,
      "sampledNeighbors": 4,
      "perPositionCap": 2,
      "maxRounds": 2,
      "patience": 1,
      "minPhesaSimilarity": 0.0,
      "minScoreThreshold": 0.0,
      "improvementTolerance": 1.0e-4,
      "reportAllCandidates": true,
      "logLevel": "NONE"
    }
  },
  "fullOptimization": {
    "request": {
      "beamSize": 10,
      "neighborPoolSize": 8,
      "sampledNeighbors": 4,
      "perPositionCap": 2,
      "maxRounds": 6,
      "patience": 3,
      "minPhesaSimilarity": 0.60,
      "minScoreThreshold": 0.0,
      "improvementTolerance": 1.0e-4,
      "reportAllCandidates": true,
      "logLevel": "SUMMARY"
    }
  },
  "orchestration": {
    "workerThreads": 8,
    "queueCapacity": 1000,
    "progressIntervalSeconds": 60,
    "reactionWeighting": {
      "mode": "BUCKETED_PRODUCT",
      "buckets": [
        { "maxProductExclusive": 2000, "weight": 0.01 },
        { "maxProductExclusive": 50000, "weight": 0.1 },
        { "weight": 1.0 }
      ]
    },
    "duplicateCacheSize": 200000
  },
  "run": {
    "maxRuntime": "24h",
    "randomSeed": null
  },
  "output": {
    "hitsTsv": "screening_hits.tsv",
    "minReportedSimilarity": 0.6
  }
}
```

A complete template is kept in [`continuous-screening-config.example.json`](continuous-screening-config.example.json). Split-query examples are [`continuous-screening-config.split.example.json`](continuous-screening-config.split.example.json) and [`continuous-screening-query.example.json`](continuous-screening-query.example.json).

## Config Blocks

### `inputs`

```json
"inputs": {
  "rawFull": "full.rawspace.gz",
  "rawDownsampled": "downsampled.rawspace.gz"
}
```

`rawFull` is the full merged rawspace used for final local optimization.

`rawDownsampled` is the reduced rawspace used for sampling, reaction weighting, and optional micro optimization. It should be marked in rawspace metadata with `space.role=downsampled`.

### `query` and `queryFile`

Specify exactly one query source. You can put it directly in the main config:

```json
"query": {
  "phesaFile": "query.phesa"
}
```

or put it in a separate file:

```json
"queryFile": "continuous-screening-query.example.json"
```

Supported query fields are:

- `smiles`: a SMILES string. The CLI generates one conformer and a PheSA descriptor.
- `idcode`: an OpenChemLib IDCode. The CLI generates one conformer and a PheSA descriptor.
- `sdfFile`: an SDF file. The selected record's 3D coordinates are used directly.
- `sdfRecordIndex`: zero-based SDF record index, default `0`.
- `phesaFile`: text file containing a string-serialized encoded PheSA descriptor.

For production screening, `phesaFile` is usually the most reproducible option.

### `sampling`

```json
"sampling": {
  "attemptsPerReaction": 20,
  "minAtoms": 1,
  "maxAtoms": 0,
  "maxRotatableBonds": 15,
  "minSimilarity": 0.65
}
```

`attemptsPerReaction` controls how many random sampled products are attempted per selected reaction.

`minAtoms`, `maxAtoms`, and `maxRotatableBonds` are hard filters on sampled assembled products. `maxAtoms: 0` means no upper atom limit.

`minSimilarity` is the first PheSA threshold. Candidates below this threshold do not enter downstream optimization.

### `microOptimization`

```json
"microOptimization": {
  "enabled": false,
  "request": {
    "beamSize": 5,
    "neighborPoolSize": 4,
    "sampledNeighbors": 4,
    "perPositionCap": 2,
    "maxRounds": 2,
    "patience": 1,
    "minPhesaSimilarity": 0.0,
    "minScoreThreshold": 0.0,
    "improvementTolerance": 1.0e-4,
    "reportAllCandidates": true,
    "logLevel": "NONE"
  }
}
```

This optional stage runs local optimization on the downsampled rawspace. It is intended as a cheap cleanup/refinement step before full-space optimization.

Defaults are intentionally small: low beam size, few neighbors, and only two rounds.

### `fullOptimization`

```json
"fullOptimization": {
  "request": {
    "beamSize": 10,
    "neighborPoolSize": 8,
    "sampledNeighbors": 4,
    "perPositionCap": 2,
    "maxRounds": 6,
    "patience": 3,
    "minPhesaSimilarity": 0.60,
    "minScoreThreshold": 0.0,
    "improvementTolerance": 1.0e-4,
    "reportAllCandidates": true,
    "logLevel": "SUMMARY"
  }
}
```

This is the expensive local optimization stage against the full rawspace.

Optimization request fields:

- `beamSize`: number of best assemblies retained during local search.
- `neighborPoolSize`: number of nearest/similar replacement synthons considered per position.
- `sampledNeighbors`: number of neighbors sampled/tried per position.
- `perPositionCap`: diversity cap for accepted beam entries per synthon position.
- `maxRounds`: maximum local-search rounds.
- `patience`: stop after this many rounds without meaningful improvement.
- `minPhesaSimilarity`: threshold for entries included in the optimizer result list.
- `minScoreThreshold`: scorer threshold during optimization. `0.0` is permissive; stricter values can reject seeds/candidates earlier.
- `improvementTolerance`: minimum score gain counted as an improvement.
- `reportAllCandidates`: if `true`, output/report all scored candidates above `minPhesaSimilarity`, not only the final beam.
- `randomSeed`: optional per-stage seed. If omitted/null, the top-level `run.randomSeed` is used.
- `logLevel`: `NONE`, `SUMMARY`, `IMPROVEMENTS`, or `VERBOSE`.

### `orchestration`

```json
"orchestration": {
  "workerThreads": 8,
  "queueCapacity": 1000,
  "progressIntervalSeconds": 60,
  "reactionWeighting": {
    "mode": "BUCKETED_PRODUCT",
    "buckets": [
      { "maxProductExclusive": 2000, "weight": 0.01 },
      { "maxProductExclusive": 50000, "weight": 0.1 },
      { "weight": 1.0 }
    ]
  },
  "duplicateCacheSize": 200000
}
```

`workerThreads` controls the number of parallel screening workers.

`queueCapacity` limits queued optimization work.

`progressIntervalSeconds` controls periodic `[Screening]` progress output. Use `0` to disable periodic progress output.

`duplicateCacheSize` controls how many assembled seed products are remembered for duplicate suppression.

`reactionWeighting` controls how often reactions are selected. The recommended mode is `BUCKETED_PRODUCT`, based on the product of downsampled synthon-set sizes.

For the example above:

```text
product < 2,000      -> weight 0.01
2,000 to < 50,000    -> weight 0.1
>= 50,000            -> weight 1.0
```

The older fields `reactionWeightExponent` and `reactionMinWeight` are still accepted when no explicit `reactionWeighting` is supplied.

### `run`

```json
"run": {
  "maxRuntime": "24h",
  "randomSeed": null
}
```

`maxRuntime` is the preferred run limit. It accepts compact values such as `4h`, `90m`, `3600s`, and ISO-8601 durations such as `PT4H`.

When the runtime is reached, the CLI stops queueing new jobs, drops queued-but-not-started jobs, waits for active jobs to finish, writes final output, and prints the final summary.

`iterations` is still accepted for bounded test runs, but wall-clock `maxRuntime` is usually better for production/cluster jobs.

`randomSeed: null` creates a time-derived seed. The effective seed is printed at startup. A fixed numeric seed can be used for more repeatable test runs, though exact reproducibility is still limited by parallel scheduling.

### `output`

```json
"output": {
  "hitsTsv": "screening_hits.tsv",
  "minReportedSimilarity": 0.6
}
```

`hitsTsv` is the output hit table.

`minReportedSimilarity` is the final global reporting threshold. Candidates below this value are not written, even if the optimizer evaluated them.

## Threshold Cheat Sheet

The most important score thresholds are separate:

- `sampling.minSimilarity`: early gate for sampled downsampled products.
- `microOptimization.request.minScoreThreshold`: scorer threshold during optional micro optimization.
- `fullOptimization.request.minScoreThreshold`: scorer threshold during full optimization.
- `fullOptimization.request.minPhesaSimilarity`: threshold for candidates included in full optimizer results.
- `output.minReportedSimilarity`: final TSV reporting threshold.

For broad exploration, keep early thresholds permissive and filter later. For faster stricter runs, raise `sampling.minSimilarity` and/or `output.minReportedSimilarity`.

## Logging

Default logging is quiet:

```json
"logLevel": "NONE"
```

Useful levels:

- `NONE`: no per-optimization output.
- `SUMMARY`: one concise `[FullOpt]` line after each full optimization in continuous screening.
- `IMPROVEMENTS`: summary plus local optimizer improvement lines.
- `VERBOSE`: summary plus verbose candidate-level local optimizer output.

Recommended continuous-screening setting while testing:

```json
"fullOptimization": {
  "request": {
    "logLevel": "SUMMARY"
  }
}
```

Example summary line:

```text
[FullOpt] job=123 rxn=m_abc sampled=0.6123 micro=0.6451 fullBest=0.7018 gain=0.0567 reported=3 candidates=84 rounds=4/5 stop=PATIENCE
```

If seed scoring fails and logging is enabled, the local optimizer reports it in structured form:

```text
[LocalOpt] rxn=m_abc seed=frag1;frag2 initial=0.4200 stop=SEED_SCORE_FAILED
```

## Output TSV

Continuous screening writes a TSV with columns similar to:

```text
rxnId  sourceSpace  fragIds  Structure [idcode]  atoms  rotatableBonds  phesaSimilarity  attemptIndex  seedFragIds
```

`sourceSpace` is filled from merged rawspace reaction metadata when available, typically `source.spaceName`; otherwise it is empty.

`Structure [idcode]` contains the assembled product IDCode.

`seedFragIds` records the sampled or post-micro seed tuple that started the full optimization.

## Practical Starting Points

For a first production-like run:

```json
"run": { "maxRuntime": "4h", "randomSeed": null },
"fullOptimization": { "request": { "logLevel": "SUMMARY" } }
```

For quieter cluster runs:

```json
"fullOptimization": { "request": { "logLevel": "NONE" } },
"orchestration": { "progressIntervalSeconds": 300 }
```

For stricter output without changing search behavior:

```json
"output": { "minReportedSimilarity": 0.7 }
```

For stricter full optimization scoring itself:

```json
"fullOptimization": {
  "request": {
    "minScoreThreshold": 0.6,
    "minPhesaSimilarity": 0.6
  }
}
```

Use the latter carefully: a high `minScoreThreshold` can reject the starting seed assembly and stop that optimization with `SEED_SCORE_FAILED`.
