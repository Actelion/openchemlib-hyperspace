# Hyperspace web workbench

A local Java owner, a React/TypeScript/OpenChemLib JS client and the existing MCP
server share one persistent screening workspace. The prototype runs real
combinatorial substructure searches, bounded assembly, learned 2D/ensemble 3D/
supplied-pose screening, detached preparation and scientific curation. Existing
CLI and desktop applications retain their entrypoints.

## Build and launch

Use JDK 22 or newer, Maven, and Node.js 22.12+; verified with JDK 23 and Node 24.
From the repository root:

```sh
mvn -pl openchemlib-hyperspace-workbench,openchemlib-hyperspace-cli -am package -Dmaven.test.skip=true
npm ci --prefix openchemlib-hyperspace-workbench/frontend
npm run build --prefix openchemlib-hyperspace-workbench/frontend
java -Djava.awt.headless=true -Xmx4g \
  -jar openchemlib-hyperspace-workbench/target/openchemlib-hyperspace-workbench-3.1.0.jar \
  --workspace /absolute/path/to/my-workspace --port 8765
```

Open the printed browser URL including its token fragment. The owner binds only
to `127.0.0.1`; all API requests require the token. The browser moves it to session
storage and removes it from the address bar. A mode-0600 `connection.json` holds
connection information for local agents. Only one owner can lock a workspace.
Foreign Origin/Host values are rejected. This is a single-user local application.

The owner serves `openchemlib-hyperspace-workbench/frontend/dist` relative to the
launch directory. For another working directory add
`--frontend /absolute/path/to/frontend/dist`. Keep the workbench JAR and its
`target/lib` directory together. ONNX native libraries remain in ordinary runtime
JARs rather than a new shaded workbench JAR. No model download occurs at startup.

CPU is the default distribution. To build a CUDA distribution, use
`mvn -pl openchemlib-hyperspace-workbench -am clean package -Pcuda -Dmaven.test.skip=true`.
The workbench profile selects `com.microsoft.onnxruntime:onnxruntime_gpu:1.21.0`
in place of the CPU runtime. Clean when switching profiles, especially if invoking
classes with a `lib/*` classpath. The inference provider is also selected explicitly
in the UI/request: requesting CUDA without a working provider fails. It never
silently uses CPU. Install the native CUDA/cuDNN dependencies required by that
ONNX Runtime version. GPU execution was not verified in this workbench milestone.

## Light and dark appearance

The Precision theme uses signal cyan for actions, selection and query-match
highlights, with burnt-orange accents for the second synthon position. Crisp
panels and a vector prism mark follow the supplied design inspiration. The dark
appearance uses deep blue-black surfaces and readable molecular ink, including
lighter shades for saturated atom colors.

Use the moon/sun button in the top bar. The choice persists in this browser;
without an override, the OS appearance preference is used. Changing appearance
preserves the current editor instance and any unpublished drawing. This is a
presentation change: IDCodes, atom/bond references and search scores are unchanged.

## Reproducible toy workflow

1. In **Libraries**, choose **Prepare toy library**. It creates real rawspace and
   FragFp/512 combinatorial index assets from two public synthetic reactions,
   each with two sets of two synthons. No licensed supplier structures are used.
2. Wait for the preparation job, then **Load resource** on **Toy · substructure**.
   Loading the rawspace separately is unnecessary for index search.
3. Choose **Use in screening**, paste `CCO` under **Paste or import a structure**,
   enable **Fragment query**, and publish.
4. Search. Select families and expand them with an explicit assembly-attempt cap,
   for example 4. Represented family counts may overlap; they are not deduplicated
   molecular-space totals. Expansion enumerates a deterministic bounded prefix.
5. Inspect highlights, reaction and synthon provenance. Switch between the structure
   grid, sortable table and linked molecular-weight/rotor plot. Property filters
   apply only to retained results; they cannot establish global filtered top-K.
6. Select occurrences, optionally choose **Diverse 20**, then save a collection with
   scientific notes. Edit notes or explicitly replace collection membership.
   Export SDF, TSV, Markdown and JSON provenance. Compare two saved runs by identity.
7. Reconnect or restart the owner. Original query revisions, run settings, results,
   collections and meaningful selections persist. Loaded resources require an
   explicit reload after owner restart. Saved histories and results are paged.

The same public data can be generated without the browser:

```sh
java -Djava.awt.headless=true \
  -cp 'openchemlib-hyperspace-workbench/target/openchemlib-hyperspace-workbench-3.1.0.jar:openchemlib-hyperspace-workbench/target/lib/*' \
  com.idorsia.research.chem.hyperspace.workbench.ToySetup /tmp/hyperspace-toy
```

## Catalogue and preparation

**Add data → Register prepared data** accepts a backend-local path, detects its
format and displays coverage/capabilities before registration. Supported kinds:
rawspace, downsampled rawspace, combinatorial substructure index, existing legacy
similarity index, sampled-product fingerprint index, native flat-molecule cache,
supported Python supplier cache, canonical graph bundle, compact projection and
supplied-pose bundle. A `manifest.json` path or its containing bundle/cache
directory is accepted. Legacy similarity assets can be prepared/loaded; their
existing desktop/CLI search engine has no additional web search form.

A logical library may have several source snapshots and prepared assets. Sharing
its name does not assert identical molecular coverage. Registration inspects
bounded metadata; it neither materializes all indexes nor scans/hashes enormous
fingerprint payloads. Explicit loading validates the actual resource and reports
ready/failed status, counts and model compatibility. Aliases share an instance;
queued/running jobs lease it and prevent premature unloading. Model checksums and
cache compatibility use the existing strict Hyperspace3D contracts.

**Add data → Prepare source data** supports two recipes:

- **Supplier synthons/rawspace:** inspect a bounded sample of Enamine, XtalPi,
  supported CSV/ZIP input or an existing rawspace; choose substructure and optionally
  legacy similarity indexes. Importer settings use the existing CLI schema.
- **Enumerated molecule table → fingerprints:** select source columns, canonical
  and compatible compact bundles, a stable output directory, CPU/CUDA, encoder batch,
  threads, child JVM heap and resumable source-shard size. Outputs contain accepted
  molecule rows and explicit rejection/model/source provenance.

Open **Runtime and model paths**, verify the builder CLI/fingerprint JAR/lib/model
paths, preview the recipe, then start. The owner calls existing detached MCP
workers, monitors their persistent handles and automatically registers successful
outputs without loading them. The browser can close during execution. **Jobs**
provides bounded logs, cancellation and **Resume / reattach**. Recovery checks the
recorded worker first; it does not launch a duplicate live computation. Fingerprint
resume reuses committed shards at the same output path; space recovery can reuse a
completed raw import. Unfinished owner jobs become interrupted on owner restart;
reattach/resume is explicit. Completion registration is idempotent.

The Slurm section reuses the existing profile/exporter: it writes scripts and
configuration only. It does not submit, contact or monitor a cluster.
See [examples](examples/README.md) for explicit runtime/asset/search requests and
[rawspace importer documentation](../RAWSPACE_IMPORTERS.md) and
[supplier preparation guide](../SUPPLIER_SPACE_IMPORTS.md);
existing CLI `--help` and MCP preparation schemas remain authoritative.

## Query studio and bridges

OCL IDCode is the shared Java/JavaScript representation; depiction coordinates are
stored separately. Import SMILES, SMARTS, IDCode or MOL, or draw using the native
OCL JS canvas. Publishing creates an immutable revision. An unpublished drawing is
preserved if an agent publishes a different query: explicitly branch the drawing
or discard it to load the shared revision. **Revision history & branches** and
**Branch run query** preserve the originals.

The decoded inspector and typed Java service support element alternatives/exclusion,
bond alternatives, ring/aromaticity/stereo conditions and bridges. Agents inspect
revision-scoped atom/bond references rather than inventing IDCodes. Editing or
canonicalization invalidates old references. Shared atom/bond selections can be
edited by either scientist or agent. **Test example molecules** uses native OCL
substructure semantics on up to 50 concrete SMILES.

Try **Ring scaffold**, select a ring atom and allow N/O (`7,8`). For a linker, use
**Linker endpoints**, select its bond and set intermediate-atom bounds. OCL stores
`min` directly (0–15) and `span = max - min` (0–15); `0..0` has an all-zero bridge
encoding and is rejected. Hyperspace supports at most two bridge bonds. The core
helper now preserves endpoints for zero-length paths and expands the higher-index
bridge first. Path-expanded searches produce candidate families; native original
bridge-query matching verifies every assembled retained product and its highlights.
The family scope warns that expanded-path discovery has this limitation.

**Supplied pose** accepts a backend MOL/SDF path (up to 16 MiB) or a small upload
(up to 1 MB). Multi-record SDF requires an explicit zero-based record. The selected
coordinates, source hash and record are retained in a workspace-owned snapshot.
Graph edits detach the pose; ordinary 2D depictions do not become query conformers.

## Learned search and bounded ranking

First register/load the canonical **50 atom / 36 pair / 32 atom-mask** graph bundle,
then its compatible cache. Compact and pose assets additionally require a loaded
canonical bundle. The selected bundles, hashes and settings are pinned per run.
Models/caches remain external; setup errors never silently select another model.

| Mode | Prerequisites | Meaning of ranking |
| --- | --- | --- |
| Learned 2D | Canonical bundle, compact 16D projection, compatible native/Python cache | Compact scan covers accepted cache rows. Optional exact SkelSpheres reranks only the learned shortlist. |
| Ensemble 3D | Canonical bundle and compatible 128D cache | Comparator scans accepted cache rows for the selected learned objective. |
| 2D → 3D | Canonical + compact bundles, compatible cache | 128D PheSA ranks only the selected compact shortlist. |
| Supplied pose → ensemble | Canonical + compatible external pose bundle/cache, one validated query conformation | Directed learned pose→graph-ensemble scores. No candidate conformers or computed 3D alignments are produced. |

**Validate model / cache** checks provenance before scanning. Choose direct total,
shape or pharmacophore for 3D. An explicit shape/PP composite is separately labelled
and is not claimed equivalent to the trained total head. Pose mode uses its supported
PP-weight contract. See [POSE_SEARCH.md](../POSE_SEARCH.md) for the external pose
bundle, supplied-coordinate and shared-candidate compatibility details.

For **bounded substructure → similarity**, expand selected families or open a
collection, choose ensemble/pose mode, then **Rank bounded candidates**. The adapter
freezes a bounded source prefix, encodes each unique accepted graph in batches and
preserves every source occurrence/synthon tuple. Diagnostics record expansion policy,
limits, rejected graphs/reasons and scored occurrences. This never claims exhaustive
ranking of every represented combinatorial combination. The toy's expanded products have six heavy atoms and can be ranked. First publish
a supported concrete query such as `CCCCCCO`: the initial `CCO` fragment is below
the learned model's six-heavy-atom minimum. The earlier expansion remains pinned
to its original fragment query.

## Connect an agent to the visible workspace

The build above also produces the existing MCP JAR. Start it with:

```sh
java -jar openchemlib-hyperspace-mcp/target/openchemlib-hyperspace-mcp.jar \
  --workbench /absolute/path/to/my-workspace/connection.json
```

Use [examples/mcp-client.example.json](examples/mcp-client.example.json) in a local
MCP client, replacing absolute paths. `wb_workspace_state` discovers the named
workspace, stable ID, revision and capabilities. `wb_workspace_history` pages older
runs/collections/jobs. `wb_` tools mirror the owner's operation catalogue: inspect,
register, prepare, edit/test queries, search/expand/rank, monitor/cancel/read logs,
page results (maximum 200), inspect a hit, share selections, curate/compare/export,
and navigate the browser. `query.history` pages immutable query revisions.

All mutations require `expectedRevision`; stale human/agent edits fail instead of
silently overwriting. Activity identifies the actor. Query/selection undo restores
context; it does not undo computation or delete original runs. Polling makes agent
changes visible promptly. Hover, plot camera and unpublished editor gestures remain
local. The separate MCP process calls the owner and never opens duplicate indexes
or independently changes workspace files. Existing `--config` MCP mode is preserved.

## Architecture, persistence and identity

`WorkbenchApplication` is the transport-independent service. A thin authenticated
JDK loopback HTTP server publishes typed operations used by React and MCP. This
avoids coupling chemistry to the legacy server or a new Spring/Hilla lifecycle.
`WorkspaceStore` uses an exclusive file lock and atomic copy-on-write metadata.
Domain mutations advance edit `revision`; all persisted events advance `sequence`.
Background progress does not invalidate unrelated edits. Browser writes are queued
and advance only across their own proven revisions, never an agent's changes.

`JobService` bounds workers/queue and retains resource leases until computation
returns after cancellation. ONNX sessions are shared by loaded model/runtime options
and closed after the last alias unloads. Results live in per-run files rather than
workspace metadata. Responses summarize recent history; paged operations retrieve
older records. Retained products are capped at 20,000, families at 1,000, response
pages at 200. Family snapshots retain only the per-position prefix sufficient for
permitted expansion while preserving original represented counts.

Molecular identity is canonical OCL IDCode including charge and defined stereo,
with no tautomer/salt/protonation normalization. Different supplier rows and synthon
tuples remain distinct occurrences. Method-specific raw scores and source ranking
scores stay separate. Diversity uses OCL FFP512 `DiversitySelector` over unique graphs
and returns occurrence references. Exports retain identity and source provenance.

## Verification

```sh
mvn -pl openchemlib-hyperspace-workbench test
npm run test:query-parity --prefix openchemlib-hyperspace-workbench/frontend
mvn -pl openchemlib-hyperspace-mcp test
npm run build --prefix openchemlib-hyperspace-workbench/frontend
```

Start an owner on a **fresh** temporary workspace before end-to-end tests:

```sh
cd openchemlib-hyperspace-workbench/frontend
npx playwright install chromium
WORKBENCH_CONNECTION=/tmp/my-fresh-workspace/connection.json npm run test:browser
```

`PLAYWRIGHT_CHROMIUM_EXECUTABLE` can select an installed Chromium binary. The real
learned browser fixture builds public molecule fingerprints with the canonical and
compact bundles. Pose coverage additionally requires the external deployment bundle;
set `WORKBENCH_POSE_BUNDLE` for the browser and `-Dhyperspace.poseBundle=...` for Java.
No confidential structures, large caches or new production models belong in Git.

See [VERIFICATION.md](VERIFICATION.md) for checks actually run and limitations,
[IMPLEMENTATION.md](IMPLEMENTATION.md) for the acceptance audit.
