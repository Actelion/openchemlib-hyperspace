# Workbench implementation and acceptance ledger

The requested local-first prototype is implemented. This audit maps its full scope
to code and verification; it does not expand the milestone into production search,
cluster submission or multi-user hosting. See `VERIFICATION.md` for actual checks
and environment limits, and `README.md` for reproducible workflows.

## Architecture and reuse

`WorkbenchApplication` owns the transport-independent operations; `WorkspaceStore`
owns persistent metadata and the exclusive owner lock. A thin JDK loopback HTTP
transport and the existing MCP server call the same operations. React/TypeScript
and native OCL JS provide the frontend. This fits the existing Maven reactor without
coupling chemistry to Hilla or changing the legacy server/desktop lifecycle.

Reused implementations include `SubstructureSearchHelper`, `SynthonAssembler`,
GUI 2's bounded expansion/highlighting concepts, `RawSynthonSpace` import/build
commands, native/Python fingerprint readers, learned screeners/exact 2D reranking,
ONNX model validation, supplied-pose runtime and the existing detached MCP workers/
Slurm exporter. Structurized was inspected separately: its immutable OCL inspection,
identity and occurrence-aware patterns informed the design; no PRISM stack was added.

The only core search algorithm change fixes two confirmed bridge expansion defects:
zero intermediate atoms must retain both endpoints, and two bridge bonds must be
expanded in descending original bond-index order. Final concrete-product matching
and highlights use the original native OCL query rather than assuming expanded
paths precisely capture every ring/shortest-path condition.

## Full acceptance audit

| Requirement | Implementation and evidence |
| --- | --- |
| [x] Local launch, persistent workspace and shared owner | Thin JAR/lib distribution, `WorkbenchServer`, file lock, token-protected loopback API and private connection file. Packaged browser runs, HTTP/reconnect/exclusive-lock tests. README includes Java/Node prerequisites, launch and MCP configuration. |
| [x] Logical libraries, source snapshots, separate assets and chemical coverage | `AssetService` records explicit provenance, capabilities, counts and coverage. All ten requested asset kinds are recognized. Raw/downsampled roles are preserved; standalone product manifests, native/Python caches and graph/compact/pose bundles are validated using existing contracts. Registration does not scan large payloads. `CatalogueFormatsTest`, preparation and learned compatibility tests. |
| [x] Explicit loading/unloading and active-job protection | Resource aliases share instances, leases protect queued/running work, native sessions reuse runtime options and close after final unload. Failed files become failed, restart resets actual loaded status, missing CUDA fails explicitly. `SearchLifecycleTest`, `LearnedWorkflowTest`, `CollaborationHistoryTest`. |
| [x] Register and prepare data, resource recipes, recovery and Slurm export | Add-data paths inspect bounded supplier/ZIP/table samples, validate runtime paths/models/CPU-GPU/threads/heap, call existing detached import/index/fingerprint workers and auto-register outputs. Real legacy similarity index preparation/loading is covered. Resume checks handles first, reuses committed shards/completed import, and registration is idempotent. Slurm remains export-only. `PreparationWorkflowTest` and browser scenario 4. Small query uploads are supported. |
| [x] Scientific interface and persistent exploration | Libraries, Screening, Saved destinations; native editor, families/grid/table/plot, collapsible inspector/jobs/activity, light/dark themes, keyboard shortcuts/focus-managed dialogs, loading/error/prerequisite states and paged results/history. Browser scenarios 1–6 plus visual review. Original run query/settings remain immutable. |
| [x] Java/JavaScript OCL query contract and real bridge semantics | Shared IDCode + separate coordinates/fragment state, typed atom/bond alternatives/exclusions, ring/aromaticity/stereo/bridge edits, readable inspector and example tests. Exact Java→JS IDCode/feature round trips cover 3 original + 11 extended fixtures. Native OCL bridge tests cover zero-length endpoints, two original bond indices and concrete product checks. |
| [x] Agent editing, stale refs, shared selection, undo and branching | `QueryService`, revision-scoped atom/bond refs and optimistic mutations; selected-feature editing, example matching, immutable publication, history/branch-from-run, graph/pose detachment and query/selection undo. Unpublished browser drawings survive competing agent publication and require explicit branch/discard. Java tests and real separate MCP/browser scenarios 2, 5 and 6. |
| [x] Real combinatorial search and bounded expansion | Existing headless search and assembly, bounded family DTOs with represented counts, deterministic globally capped expansion/tuple dedup, native original-query verification, optional aligned 2D depictions, reactions/synthons/source occurrences. `SearchLifecycleTest`, bridge tests, browser scenarios 1/2 and actual prepared supplier search. Counts explicitly distinguish families/combinations/products/identities. |
| [x] Real flat 2D, ensemble 3D, cascade and supplied-pose searches | `LearnedSearchService` adapts existing native/Python readers and screeners. Strict checkpoint/model/cache provenance is checked before scanning. Direct/composite objectives and shortlist/exhaustive scopes are explicit. Pose coordinates/record/hash persist; directed pose→ensemble predictions produce no computed alignments. Real CPU ONNX/native/Python tests and browser scenario 3 exercise all four modes with available external pose bundle. |
| [x] Bounded substructure-to-learned-ranking adapter | Frozen source prefix, unique-graph batched featurization/encoding, graph or pose comparator, occurrence/synthon provenance, expansion policy/settings, rejected/scored counts, objective and timing diagnostics. Real toy search→assembly→learned-ranking is tested in Java and browser scenario 3. This does not rank all represented combinations. |
| [x] Curation, occurrence identity and method scores | Grid, sortable table, selected inspector, retained-only property filters, linked plot, OCL FFP512 diversity, named/editable notes/memberships and canonical identity comparison. Canonical OCL IDCode includes charge/stereo with no salt/tautomer normalization; duplicate source occurrences remain separate. Original method scores remain in source provenance when reranked. Unit/integration and browser checks. |
| [x] SDF/TSV/Markdown/JSON provenance exports | `ResultService` exports run or collection occurrences, structures, synthons and raw scores plus versioned provenance. Export row/record/count/provenance assertions and real browser download. |
| [x] Shared MCP operations | Existing MCP `--workbench connection.json` discovers live `wb_` tools and proxies owner calls. No second loaded-index owner or direct state-file mutations. Catalogue/preparation/query/search/jobs/logs/cancellation/paged results/selection/collections/comparison/navigation/export/history all use the same validated service. Real stdio MCP/browser collaboration and packaged legacy MCP regressions verify both modes. |
| [x] Actor history and revision safety | Every meaningful scientist/agent mutation checks `expectedRevision` and records actor/activity. Background progress uses event sequence only; queued browser writes cross only their own proven mutations. Undo never reverses completed computation. Terminal job publication/cancellation share the metadata lock and preserve worker handles. Query/workspace/lifecycle tests and browser conflicts. |
| [x] Onboarding and reproducibility | In-app help covers registration, both preparation recipes, advanced/linker queries, learned/pose modes, bounded ranking, curation and agent attachment with scope/prerequisites. README, small external asset/runtime/search/MCP examples and public synthetic toy. Root README/MCP guide link the workbench. No new large/private datasets or production model files. |
| [x] Boundary tests and browser acceptance | Workbench Java suite, Java/OCL JS feature parity, reused upstream graph/ONNX/pose tests, full source-reactor package/tests, actual scientist/MCP/preparation/learned browser scenarios. Current numerical evidence and skipped desktop/GPU boundaries are reported in `VERIFICATION.md`. |

## Explicit limits and extension points

- Single-user, backend-local, Linux detached-worker prototype; no built-in chat or
  cloud/multi-user deployment. The browser has no filesystem picker for huge data.
- Legacy similarity and sampled-product assets are catalogued using their existing
  formats. Their legacy desktop/CLI and product-index workflows remain available;
  the web learned search form operates native/Python **flat-molecule** caches.
- Retention is bounded (20,000 products, 1,000 families, 200 rows/result page).
  Native original-query matching verifies retained bridge products; expanded path
  candidate-family discovery is explicitly labelled and can overgenerate families.
- CPU numerical and real workflow checks are completed. GPU dependencies/provider
  selection are implemented; actual GPU inference/throughput and production-scale
  caches are not validated here. External bundle provenance remains mandatory.
- Production beam coordination, exhaustive constrained 3D ranking over huge families,
  candidate-conformer databases, exact alignment/PheSA, docking/FEP, cluster job
  submission/monitoring and elaborate analytics are deferred as requested.
