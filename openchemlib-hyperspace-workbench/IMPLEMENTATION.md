# Workbench implementation and verification ledger

This ledger tracks the full requested prototype. An unchecked item remains in scope.
The existing desktop, CLI and external-model workflows must remain functional.

## Architecture

A new plain-Java application service owns one persistent workspace and its loaded
resources. Browser HTTP and MCP operations call that owner. Use a typed HTTP API
and React/TypeScript/OCL JS; this avoids coupling a new frontend to the existing
legacy Spring/WebSocket server or introducing Hilla-generated transport into the
chemistry domain. Local transport uses the JDK HTTP server, loopback binding and
validated DTOs. Native ONNX dependencies remain a thin JAR/lib distribution.

Inspected reuse boundaries: headless `SubstructureSearchHelper`,
`SynthonAssembler.expandCombinatorialHit`, GUI 2's bounded expansion/highlighting
policy, rawspace builders/importers, DeepSpace7 cache readers/model provenance,
`POSE_SEARCH.md`, and existing detached MCP preparation/pose workers. Structurized
was inspected at a shallow checkout: reuse its immutable molecular inspection,
OCL identity and occurrence-aware result patterns; do not import PRISM or treat
its SMARTS tools as a complete OCL query-feature editor.

## Full acceptance checklist

- [ ] Local launch command, persistent owner, browser connection, MCP shared owner.
- [ ] Catalogue: logical libraries/source snapshots, separate prepared assets,
  coverage/provenance/capabilities/counts, all requested asset and bundle types.
- [ ] Explicit load/unload/status and active-job protection; bounded registration
  without hashing huge fingerprint payloads.
- [ ] Add data: inspect/register prepared assets; inspect source and execute
  import/index/fingerprint recipes with explicit resources, resume, automatic
  output registration; retain Slurm export-only support; small query uploads.
- [ ] Three polished destinations, query/results/inspector/jobs layout, dark/light,
  keyboard and empty/error states, bounded pagination and persistent selections.
- [ ] Immutable query revisions, Java/OCL JS IDCode/query-feature round trips,
  atom/bond alternatives, ring/aromatic/stereo edits and true OCL bridge semantics.
- [ ] Agent inspection/edit/test/publish, revision-scoped references, stale-edit
  rejection, selected atoms/bonds, graph/pose detachment and undo/branching.
- [ ] Real combinatorial substructure search, family counts, bounded expansion,
  highlights/aligned 2D and synthon/source provenance.
- [ ] Real flat 2D, graph 3D, cascade and supplied-pose workflows with capability
  gating, strict provenance compatibility, correct score/ranking scope.
- [ ] Bounded substructure-to-learned-ranking adapter and expansion/rejection counts.
- [ ] Grid/table, inspector, retained-result property filters, linked plot, diversity,
  collections/notes, occurrence-preserving molecular identity and run comparison.
- [ ] SDF/TSV/Markdown/JSON exports; incompatible method scores stay separate.
- [ ] MCP workspace discovery/attach, catalogue/preparation, queries, jobs/logs/cancel,
  bounded results, shared selections/collections, browser navigation and exports.
- [ ] Actor activity, optimistic revisions and undo that never reverses computation.
- [ ] In-app guidance plus full workflow documentation; redistributable real toy
  setup and external-cache/model configuration examples; no large datasets in Git.
- [ ] Tests for all named boundaries and complete scientist/agent toy scenario.
- [ ] Real browser visual/interaction verification and concise verification report.

The original specification also defers production beam coordination, exhaustive
3D ranking over huge families, candidate-conformer databases, exact alignment,
docking/FEP, multi-user deployment and elaborate analytics. These are not substitutes
for any unchecked requirement above.
