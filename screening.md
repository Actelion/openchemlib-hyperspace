# Continuous Screening Run A

## Goal
- Run a continuous virtual screening pipeline in three stages:
  1. Sample candidate assemblies from a downsampled synthon space.
  2. Optionally run a fast micro-optimization on downsampled space.
  3. Run full local optimization on the full synthon space and write hits to TSV.

## Original command line (legacy)
```bash
--minAtoms 26 --maxAtoms 38 --rawFull /mnt/app/cadd/VS/temp/realspace2023.rawspace.gz --rawDownsampled /mnt/app/cadd/VS/temp/realspace2023_downsampled.rawspace.gz --querySmiles "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" --candidateThreshold 0.5 --candidateAttempts 20 --fullSampleNeighbors 20 --fullMinSimilarity 0.55 --microEnabled --fullThreads 11 --outputHits screening_hits_runA.tsv --iterations 10000
```

## JSON config equivalent
```json
{
  "inputs": {
    "rawFull": "/mnt/app/cadd/VS/temp/realspace2023.rawspace.gz",
    "rawDownsampled": "/mnt/app/cadd/VS/temp/realspace2023_downsampled.rawspace.gz"
  },
  "query": {
    "smiles": "xxxcxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
  },
  "sampling": {
    "attemptsPerReaction": 20,
    "minAtoms": 26,
    "maxAtoms": 38,
    "maxRotatableBonds": 15,
    "minSimilarity": 0.5
  },
  "microOptimization": {
    "enabled": true
  },
  "fullOptimization": {
    "request": {
      "sampledNeighbors": 20,
      "minPhesaSimilarity": 0.55
    }
  },
  "orchestration": {
    "workerThreads": 11
  },
  "run": {
    "iterations": 10000
  },
  "output": {
    "hitsTsv": "screening_hits_runA.tsv",
    "minReportedSimilarity": 0.55
  }
}
```

## Command with JSON config
```bash
java -cp openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  com.idorsia.research.chem.hyperspace.cli.ContinuousScreeningCLI \
  screening_runA.json
```

Query can also be externalized via `queryFile` (and `query` omitted), using a separate JSON with the same schema:
`{ "smiles": "...", "idcode": null, "sdfFile": null, "sdfRecordIndex": 0 }`.

The query object must contain exactly one source among `smiles`, `idcode`, and `sdfFile`.
When `sdfFile` is used, the selected SDF record conformation is used directly (single-conformation PheSA descriptor; no conformer generation), and implicit hydrogens are expanded to explicit 3D hydrogens first.

## Effective settings used by the CLI
- Candidate sampling:
  - `attemptsPerReaction=20`
  - `minAtoms=26`, `maxAtoms=38`
  - `maxRotatable=15` (default)
  - `candidateThreshold=0.50`
- Micro optimization (enabled):
  - `microBeam=5` (default)
  - `microTopL=4` (default)
  - `microSampleNeighbors=4` (default)
  - `microCap=2` (default)
  - `microRounds=2` (default)
  - `microPatience=1` (default)
  - No similarity gate in micro stage (scorer threshold fixed to `0.0`)
- Full optimization:
  - `fullBeam=10` (default)
  - `fullTopL=8` (default)
  - `fullSampleNeighbors=20` (provided)
  - `fullCap=2` (default)
  - `fullRounds=6` (default)
  - `fullPatience=3` (default)
  - `fullMinSimilarity=0.55`
- Global:
  - `workerThreads=11`
  - `iterations=10000`
  - `queueCapacity=1000` (default)
  - `progressSeconds=20` (default)
  - `dedupeMax=200000` (default)
  - `randomSeed=13` (default)
  - `output.minReportedSimilarity=0.55`
  - Output: `screening_hits_runA.tsv`

## Pipeline behavior summary
- Stage A (sampling on downsampled space):
  - Random assemblies are generated reaction-wise.
  - Hard filters: atom count, rotatable bonds, and `similarity >= 0.50`.
- Stage B (micro-opt on downsampled space, because `--microEnabled`):
  - Starts from sampled seed with no explicit similarity threshold gate.
  - Intention: sampled candidates can continue through micro search without additional similarity filtering.
- Stage C (full-opt on full space):
  - Runs only for candidates that survive Stage B.
  - Search uses scorer threshold `0.0` internally, but written results are filtered by both `fullMinSimilarity` and `output.minReportedSimilarity` (effective gate is the stricter one).

## What to monitor in this run
- `sampled` vs `submitted`:
  - A large drop now is less likely to come from micro similarity gating.
- `microOptComparisons`:
  - Near-zero values indicate micro stage is not effectively exploring.
- Final hit count and score distribution:
  - Compare with and without `--microEnabled` to quantify micro stage impact.
