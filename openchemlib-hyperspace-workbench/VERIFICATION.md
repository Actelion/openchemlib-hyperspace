# Workbench verification — 2026-11-04

The full requested local prototype is implemented and its acceptance scenarios
were exercised with real chemistry, search engines, model inference and detached
workers. This report distinguishes that evidence from GPU/production-scale claims.
The requirement-to-code audit is in `IMPLEMENTATION.md`.

## Final checks actually run

| Check | Result and coverage |
| --- | --- |
| Full isolated-source Maven reactor `package` | **171 tests reported across 61 suites: 166 passed, 5 existing desktop/display-gated skips, zero failures/errors.** All modules compiled/packaged. Final source includes workbench, MCP adaptations and core bridge fix. |
| Workbench suite | **32 passed, zero skipped/failed.** Query features/stale references, state/revision/undo/branching, saved-history pagination past 100 runs, explicit catalogue lifecycle/aliases/leases, HTTP authentication, cancellation/failure/persistence, real substructure/assembly, learned/cache/pose compatibility, native session cleanup, native/Python scores, preparation/recovery/exports. |
| Existing Hyperspace3D suite with available gates enabled | **63 passed, zero skipped/failed.** Includes graph golden features/ONNX, compact projection, native/product index builds, streamed Python cache, supplied-pose golden and existing continuous workflows. |
| Existing MCP suite with packaged integration gates enabled | **16 passed, zero skipped/failed.** Packaged standalone fingerprint/pose workflow and five-reaction toy preparation/reconnect/cancellation coexist with the new `--workbench` mode. |
| Java → OCL JS | **14 documents passed:** 3 original exact IDCode round trips plus 11 extended aromatic/ring/stereo/atom-exclusion/bond-alternative/bridge fixtures. Encoded features and decoded fields are checked. |
| TypeScript / Vite production build | Passed. OCL JS code/resources trigger the existing large-chunk advisory; no compile/build failure. |
| Workbench packaging | Passed as a normal application JAR with `target/lib`; packaged native CPU ONNX inference, detached worker classpaths and separate MCP attachment were exercised. |
| CUDA dependency selection | `-Pcuda dependency:tree -Dincludes=com.microsoft.onnxruntime` selects only `onnxruntime_gpu:1.21.0` for workbench. CPU-distribution negative test proves an explicit CUDA request fails with no CPU fallback and releases leases. No GPU execution claim. |
| Chromium end-to-end | **6 scenarios passed together**, including scientist toy flow, actual stdio MCP collaboration, every learned mode, real guided preparation, saved-run/collection edits and conflicting unpublished drawing. Final run took approximately 1.1 minutes while other CPU checks ran; this is not a throughput benchmark. |
| Browser visual review | Rendered light/dark structures, inspector, supplier/fingerprint preparation, learned supplied-pose results and agent collection inspected. Screenshots wait for theme transitions. |
| Change hygiene | `git diff --check`, example JSON parsing, documentation link checks. Generated frontend/build/cache/model outputs remain ignored; no new confidential data or production model added. |

## Browser acceptance scenarios

1. **Scientist toy workflow:** prepare/register/load the real public combinatorial
   index; publish `CCO`, search actual families, expand four tuples (two molecular
   identities), inspect grid/table/linked plot, save/export a collection and reload.
2. **Separate MCP agent:** discover this same owner, inspect the scientist-selected
   oxygen, publish typed N/O alternatives, run actual search/assembly, update shared
   hit selection, create a visible collection and navigate the scientist's browser.
3. **Learned search:** build a public native cache with real canonical/compact
   models, register/load assets via the UI, validate compatibility and run 2D with
   exact shortlist reranking, ensemble 3D and 16D→128D cascade. Open the actual
   substructure expansion from scenario 1, rank its products and inspect retained
   synthon provenance. Attach the available external pose model and supplied MOL
   conformation, select PP weight 0.25, and run directed pose→ensemble screening.
4. **Guided preparation:** inspect supplier input, preview resources and run detached
   import/index workers; use the prepared index for actual substructure search.
   Inspect a molecule table, preview CPU encoder recipe, build fingerprints and
   verify automatic output registration. Export a Slurm bundle without submission.
5. **Saved provenance:** rename the workspace, branch an original run query, inspect
   history, edit collection notes without changing membership/creation time, then
   explicitly replace members using shared selection. Reconnect and inspect notes.
6. **Concurrent drawing:** draw a real unpublished bond in OCL JS, have the MCP agent
   publish another graph, confirm the drawing survives, then publish it as an
   explicit branch of its original base. The agent revision also remains in history.

Playwright 1.64.1 used installed Chromium revision 1234 via
`PLAYWRIGHT_CHROMIUM_EXECUTABLE`. The final fresh owner workspace was
`/tmp/hyperspace-workbench-accepted-browser` on loopback port 8774, with
`java.awt.headless=true`. Screenshots/traces are generated under ignored
`frontend/test-results`. No Swing application was opened by workbench operations.

## Chemistry and numerical parity

- `GoldenFeatureParityTest`: Java 50/36 padded graph tensors match the existing
  approved Python/RDKit fixtures exactly (zero absolute tolerance).
- `JavaOnnxParityTest`: CPU Java ONNX embeddings/scores match the existing Python
  ONNX golden outputs at `1e-5` absolute tolerance; single-batch and atom-permutation
  checks pass. This reuses the approved export fixtures; no model retraining or
  fresh PyTorch export was performed for the workbench.
- Compact projection and real native index build tests pass with canonical and
  compact bundle properties explicitly supplied.
- Workbench graph screening matches existing `MoleculePheSAScreener` results at
  `1e-7`; changing comparator batch size preserves result order. Duplicate graphs
  retain distinct source occurrences.
- Workbench's native/Python-format test uses seven publicly specified molecules,
  actual encoder outputs and a generated supported supplier-cache wire fixture.
  Both readers produce equal retained identities/scores at `1e-7` for 2D, 3D and
  cascade. This tests the existing Python **format reader**; it does not regenerate
  production supplier fingerprints or claim arbitrary OCL/RDKit feature equivalence.
- Existing external `PoseDeploymentTest` checks feature/distance/embedding/comparator
  parity, coordinate/atom-order behavior and provenance rejection. Its model checks
  pass with the verified external deployment directory (embedding tolerance `2e-5`).
  Workbench pose tests additionally verify directed semantics and the configured
  shape/PP weighted total at `1e-6`.
- Native original-query tests catch zero-length bridge endpoint removal and two-bond
  index drift. Concrete acceptance/highlighting use native OCL, including ring and
  shortest-path behavior; path-expanded family discovery is labelled as candidates.

## Reproduce the stronger reactor checks

After building packages, optional upstream integration gates can be enabled:

```sh
mvn package \
  -Dhyperspace3d.modelBundle=/absolute/repo/openchemlib-hyperspace-3d/model-bundles/deepspace7-v1 \
  -Dhyperspace3d.compactBundle=/absolute/repo/openchemlib-hyperspace-3d/model-bundles/deepspace7-skelspheres16 \
  -Dhyperspace3d.distributionJar=/absolute/repo/openchemlib-hyperspace-3d/target/openchemlib-hyperspace-3d-3.1.0.jar \
  -Dhyperspace.cliJar=/absolute/repo/openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar \
  -Dhyperspace.mcpJar=/absolute/repo/openchemlib-hyperspace-mcp/target/openchemlib-hyperspace-mcp.jar \
  -Dhyperspace.toy=/absolute/repo/example/idorsia_toy_space_a.txt \
  -Dpose.deployment=/absolute/external/pose_deployment \
  -Dhyperspace.poseBundle=/absolute/external/pose_deployment/shared-candidate-seed17-verified
```

The existing packaged MCP regression expects the repository's five-reaction toy
file, including its diagnostics, rather than the workbench's smaller synthetic toy.
Without optional bundle/package properties, upstream tests retain their existing
assumption gates. Workbench's external pose test also skips when no external bundle
is configured; the local verification above used the available bundle and did not skip.

The full reactor was run in `/tmp/hyperspace-workbench-final-reactor`, copied from
current tracked source plus authored workbench/MCP files. This intentionally excluded
user-owned untracked datasets and the unrelated
`ContinuousScreeningOrchestratorTest.java`, which currently refers to an unavailable
`RawSynthonSpace.Builder.withDownsampledSpace` method. That file and its surrounding
user data were left untouched. Consequently the result proves the complete intended
source tree, not compilation of that unrelated untracked test in the original folder.

## Precision theme refresh (frontend only)

The later appearance update follows the supplied cyan/burnt-orange light/dark
reference: prism mark, sharper panels, cyan actions/match highlights and consistent
synthon-position accents. Dark SVG depictions use light ink and readable shades
of OCL's atom colors. The native canvas keeps its editor instance and uses a
presentation filter because its public API has no dark-palette setting.

Verification rerun for this change:

- TypeScript/Vite production build passes; the existing OCL chunk-size advisory
  remains. No Java, search/model contract or chemical identity change was made.
- All six real browser scenarios pass together (31.1 seconds), including an
  unpublished canvas drawing surviving a light/dark toggle and a competing agent
  edit. The existing 14 Java→OCL JS query-feature documents also pass.
- Light/dark structures, tables, library catalogue, saved collections and dialogs
  were visually inspected. Six actual public amide-coupling products were prepared
  in a separate review workspace to verify ring-containing molecule depictions.
- A focused browser check confirms that changing appearance leaves the published
  query and workspace revision unchanged. The explicit theme survives reload;
  a fresh browser follows the OS dark preference. No page runtime errors occurred.
- A 1024-pixel desktop viewport has no horizontal document overflow. Token contrast
  for body text, muted text, action text, primary buttons and orange chips is at
  least 4.5:1 in both themes. Chemical match colors and actual OCL red/blue atom
  paints are checked in the resulting dark SVGs.

Screenshots (`precision-light.png`, `precision-dark.png` and auxiliary views) are
ignored artifacts in `frontend/test-results`. The scientist's running workspace
was preserved; all preparation/interaction checks used a separate owner.

## Remaining practical limits

- Five existing desktop/display gates remain skipped; interactive Swing windows
  were not launched. Compilation/package and headless chemistry were verified.
- Actual CUDA provider initialization/inference and GPU throughput need a compatible
  GPU runtime. Only dependency selection and explicit-failure behavior were checked.
- The workbench was exercised with small public libraries and real models, not a
  multimillion-row production cache or huge supplier-space load. Resource/session
  ownership and bounded responses are tested; production capacity is not measured.
- Slurm exports are not submissions or cluster monitoring. Learned scores are not
  exact PheSA alignments, docking or prepared candidate conformer coordinates.
- No remaining implementation blocker for the requested local prototype was found.
  Production beam/exhaustive constrained ranking, candidate conformer databases,
  multi-user hosting and elaborate analytics remain outside this milestone.


## Later production-cache check — unresolved feature parity

The Enamine REAL January 2026 source directory contains 13.6M and 136M enumerated
molecule tables. A completed Python dual-vector cache was found for the larger
source: 136,430,449 source rows, 132,030,824 accepted molecules and 4,399,625
rejections in 137 shards. Its foundation/predictor/projection hashes match the
loaded canonical and compact bundles. Metadata/file-size validation passes and
registration/loading does not copy or scan the fingerprint payloads.

A seven-row real-source check exposed neutral imine-donor feature differences.
A targeted code correction plus an independently specified public 12-molecule
RDKit fixture resolves those cases exactly. On the first seven production rows,
maximum cache differences fall to approximately 0.000487 (128D) and 0.000955 (16D),
consistent with the original mixed-precision/FP16 generation. The updated 3D suite
reports 64 passed, zero failed/skipped with available integration gates enabled.

However, the broader 64-row check described in the 3D README identifies remaining
ring/aromaticity/conjugation discrepancies. The cache is **registered and enabled for exploratory screening** in the live
demo after discussing these limitations. Exact query-feature parity and the
impact on retrieval quality remain unverified; enabling it is not a claim that
all preprocessing differences have been resolved. No complete library screen, fingerprint regeneration or production-data
copy was started. No separate full 13.6M fingerprint cache was found.

The targeted featurizer correction is in the working tree and packaged 3D module;
it has not been deployed into the still-running workbench process. The existing
public-cache wire-format test remains valid for its scope; its Java-generated
vectors must not be taken as proof of arbitrary real RDKit/OCL feature parity.
The local detailed data report is in the demo workspace's
`data-notes/enamine-real-compatibility.json`; no supplier structures or identifiers
were added to the repository fixtures.
