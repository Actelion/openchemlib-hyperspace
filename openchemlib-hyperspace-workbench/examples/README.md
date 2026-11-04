# External assets and owner operations

These examples contain paths/placeholders, not datasets or credentials. They do
not load automatically. Replace paths with files available to the Java owner.
Read `GET /api/operations` or MCP `tools/list` for current validated arguments.

1. Build/start the owner following the module README. Configure your MCP client
   using `mcp-client.example.json` and call `wb_workspace_state`.
2. In Libraries, register paths from `assets.example.json`, or call
   `wb_asset_inspect` then `wb_asset_register`. New registrations get stable IDs.
   Link a common source snapshot only when it really represents the same coverage.
3. Explicitly load the canonical graph bundle and await its job. Load compact/pose
   bundles with `graphModelId` set to that loaded canonical asset. Load the cache.
   These steps validate hashes/contracts without converting the cache format.
4. Publish a concrete supported query, e.g. `CCCCCCO`. A graph-only query has no
   supplied pose; advanced substructure query features are rejected for embedding.
5. Adapt `learned-request.example.json` with returned IDs and the **current** revision.
   Call `wb_search_validate` before `wb_search_learned`. MCP arguments are the
   `arguments` object plus `expectedRevision` (the MCP bridge supplies actor).
   Poll the returned job; page retained hits with `wb_run_results`.

For HTTP, POST the full operation envelope to `/api/call` with `Content-Type:
application/json` and `X-Workbench-Token` from your private `connection.json`.
Every mutation needs a fresh expected revision. A stale response is a conflict;
inspect the new state before deciding what to publish. Do not put tokens in Git.

## Mode variants

- 2D: `mode=SKELSPHERES_2D`, `compactAssetId`, `shortlistTopK=2000`,
  `exactSkelSpheres=true`. Exact ranking is limited to that learned shortlist.
- Cascade: `mode=SKELSPHERES_THEN_PHESA`, `compactAssetId`, `shortlistTopK=2000`.
  Learned PheSA ranks only the shortlist.
- Pose: import a small MOL/SDF via `wb_query_import_pose` (explicit zero-based
  `record` for multi-record input), then `mode=POSE_PHESA_3D`, `poseAssetId`,
  `posePpWeight=0.25` if supported by that bundle, and its compatible canonical/cache
  asset IDs. Outputs are directed pose→ensemble predictions, not alignments.
- Bounded: call `wb_search_rank_candidates` with `request` (graph or pose mode,
  no index required), one `sourceRunId` or `collectionId`, and `maxCandidates`.
  Only the retained prefix is ranked; source tuples and rejection counts persist.

## Prepare a molecule cache instead

Configure runtime paths with `preparation-runtime.example.json`. Its canonical and
compact bundles must be compatible. Call `wb_preparation_inspect` with `molecules=true`,
then `wb_preparation_recipe` / `wb_preparation_submit` with:

```json
{
  "kind": "fingerprints",
  "options": {
    "input": "/data/public-molecules.tsv",
    "outputDirectory": "/data/public-cache",
    "smilesColumn": "smiles",
    "idColumn": "id",
    "device": "CPU",
    "threads": 4,
    "heap": "4G",
    "encoderBatchSize": 64,
    "sourceRowsPerShard": 100000
  }
}
```

A tiny public input is below; the three-carbon row demonstrates a stable model
size rejection. It does not prevent other rows from being indexed.

```tsv
id	smiles
hexanol	CCCCCCO
benzene	c1ccccc1
small	CCO
```

Successful outputs register without loading. Resume at the same output directory
using the job's **Resume / reattach** operation; do not restart a live worker.
Select CUDA explicitly only with a GPU distribution and compatible native libraries.
