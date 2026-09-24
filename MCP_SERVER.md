# Hyperspace MCP Server

A local, Linux-first MCP server for preparing supplier spaces and opening the Hyperspace desktop GUI. The agent can inspect inputs, import a rawspace, build substructure and similarity indexes, monitor independent jobs, and generate GUI configurations. No cloud service or API key is required by this server.

## Build and connect

Use a JDK 17 or newer and Linux `setsid` (util-linux). Build matching CLI and MCP artifacts:

```bash
mvn -pl openchemlib-hyperspace-cli,openchemlib-hyperspace-mcp -am -Dmaven.compiler.release=17 -DskipTests package
```

The executable JARs are in each module's `target` directory, named `openchemlib-hyperspace-cli.jar` and `openchemlib-hyperspace-mcp.jar`. Use matching builds: older CLI JARs do not support completion reports. Existing JVM serialization compatibility requirements still apply to `.data` files.

Create a server configuration based on `mcp-server.example.json`:

```json
{
  "workspace": "/data/hyperspace-workspace",
  "cliJar": "/opt/hyperspace/openchemlib-hyperspace-cli.jar",
  "javaExecutable": "/usr/bin/java",
  "threads": 4,
  "heap": "8G",
  "guiHeap": "8G"
}
```

Relative workspace/JAR paths are resolved against the server configuration directory. Input paths passed to tools should be absolute. Omit `javaExecutable` to use the server's Java installation. Four threads (capped by available processors), `8G` preparation heap, and `8G` GUI heap are defaults. These are not estimates for full supplier spaces: use explicit larger budgets when needed. A heap limit excludes native JVM memory, and simultaneous GUI and preparation processes need separate RAM.

Configure an MCP-capable client to run the following command over stdio:

```bash
java -jar /opt/hyperspace/openchemlib-hyperspace-mcp.jar --config /data/server.json
```

For clients accepting an `mcpServers` JSON configuration, see `mcp-client.example.json`. Other clients expose equivalent command/arguments settings. The MCP process runs on the machine holding the files. It does not upload the full spaces to the agent; tool results contain summaries, paths, and bounded samples/logs.

## First workflow

Ask the agent:

> Prepare `/absolute/path/to/idorsia_toy_space_a.txt` for substructure and similarity search, using two threads and 2G heap. Open the GUI when it is ready.

The agent should:

1. Call `get_environment` to check Java, CLI JAR, `setsid`, resource defaults, disk space, and DISPLAY.
2. Read `get_workflow_help` with `topic=workflow` and inspect the input with `inspect_input`.
3. Resolve the input format and any required column/ZIP-entry settings. Samples do not validate the whole space.
4. Call `prepare_space`, then retain the returned job ID. The call returns before chemistry starts.
5. Poll `get_job_status`; read `import` or `build` logs with `get_job_log` when useful.
6. Inspect counts and warnings on completion. Read `artifacts.guiConfig`; launch it only if requested.

Example preparation arguments:

```json
{
  "input": "/absolute/path/to/example/idorsia_toy_space_a.txt",
  "format": "enamine",
  "spaceName": "Toy space",
  "searchModes": ["substructure", "similarity"],
  "threads": 2,
  "heap": "2G",
  "guiHeap": "2G",
  "launchGui": true
}
```

`launchGui` defaults to false. The toy input is included in the repository's `example` directory. The default search profile is FragFp with 512 bits. Both search modes are enabled unless explicitly selected otherwise.

## Inputs and importer options

| `format` | Input | Options under `importOptions` |
| --- | --- | --- |
| `enamine` | Positional TSV, supported unified CSV, or ZIP table | `maxSets`, `zipEntry`, `reactionZipEntry`, `sourceFormat`, `metadata` |
| `xtalpi` | Xtalpi CSV | `maxSets`, `smilesColumn`, `idColumn`, `reactionColumn`, `synthonSetColumn`, `sourceFormat`, `metadata` |
| `csv` | Directory of per-reaction CSV files | `smilesColumn`, `idColumn`, `priceColumn`, `priceAttribute`, `synthonSetColumn`, `synthonDefaultSet`, `metadata` |
| `rawspace` | Existing `.rawspace` or `.rawspace.gz` | None; skips import |

These map to the existing importers; consult [RAWSPACE_IMPORTERS.md](RAWSPACE_IMPORTERS.md) for exact headers, filename conventions, and supplier examples. CSV directories are not arbitrary CSV collections. Enamine TSV remains positional; generic Enamine column overrides are not implemented by that CLI and are rejected rather than silently ignored.

For ZIP input, inspect the archive, then supply the exact `zipEntry`; optionally set `reactionZipEntry`. Archives are read directly, not extracted. Inspection lists at most 100 archive entries, samples at most five directory files and 20 records per table, and reads at most 1 MiB of decompressed sample content per file. Ambiguous input requires an explicit format.

`maxSets` defaults to 3 for Enamine/Xtalpi imports; 0 disables that count filter. Connector validation still applies. Filtering counts and available diagnostics are included in completion reports; logs describe excluded reactions. `metadata` is a string-to-string object. Newly imported spaces are marked `space.role=full`.

Import writes a rawspace first; the next JVM builds indexes from it. Existing rawspace input is never changed. Similarity-only preparation still builds and retains its underlying substructure index, but only similarity appears in the generated GUI configuration. Metadata descriptor tags are not substitutes for built indexes.

## Jobs and artifacts

Tools: `prepare_space`, `list_jobs`, `get_job_status`, `get_job_log`, `cancel_job`, `list_spaces`, `configure_gui`, `launch_gui`, `get_environment`, `inspect_input`, `get_workflow_help`.

Each preparation has a unique `workspace/jobs/<jobId>` directory with request/config snapshots, command argument arrays, CLI JAR SHA-256, status/timestamps, logs, reports, rawspace, indexes, and `gui.json`. Generated data belongs in that workspace, not in Git.

A detached worker owns the job. Closing the AI client or its MCP server does not intentionally stop the worker. Reconnect with the same workspace and use `list_jobs`. This does not guarantee survival of machine shutdown or an administrator/client that explicitly kills detached processes. There is no automatic restart after a worker or host crash.

One preparation runs per workspace; another submission returns `accepted=false` and the active job. There is no queue. Job states are `pending`, `running`, `succeeded`, `failed`, `cancelled`, or `interrupted`. Logs remain available after failure. A lost worker with a still-running recorded child is reported explicitly and can be cancelled before resubmitting.

`get_job_log` defaults to the worker log, offset 0, and 8192 bytes. Use `log=import` or `log=build` for chemistry progress; pass `nextOffset` on the next read. Maximum chunk size is 65536 bytes. Byte boundaries may split UTF-8 characters in display text; the saved log remains intact.

Success requires a successful process exit, completion report, nonempty expected outputs, and at least one valid reaction. Partial outputs remain in `staging` and are not advertised as searchable spaces. If import succeeded but building failed, `artifacts.rawspace` can be submitted as `format=rawspace` in a new job; no mid-stage resume is attempted.

Threads control import and similarity initialization. The underlying rawspace assembly stage retains its existing threading behavior; a requested thread count does not promise that every stage uses that many cores. No elapsed-time or RAM predictions are inferred from compressed file size.

## GUI handoff

`configure_gui` accepts `jobIds` for one or more successful preparations and creates a fresh config. Files use absolute index paths and explicit thread settings. `launch_gui` accepts its `configPath` and an optional `heap` override. A running GUI is independent from the MCP session.

The Linux desktop requires DISPLAY (including XWayland where applicable). Without it the tool returns the configuration, command argument array, and quoted shell command. A GUI launch failure does not invalidate prepared spaces. `launched` reports process startup, not successful index loading; use the GUI's process/provider states for readiness.

The GUI loads substructure and similarity providers, preserves both when saving configurations, and resolves relative index paths relative to the config file. Missing configs/indexes and load failures are reported in the GUI. Space loading happens in background threads; searches become available when ready.

## Troubleshooting

- Java/JAR missing or incompatible: check `get_environment`, use a supported Java runtime and matching CLI/MCP builds.
- Import format/columns: inspect a sample and consult importer documentation; supply the correct ZIP entry or directory column options.
- Zero retained reactions: inspect validation/filter diagnostics and input connector conventions.
- Memory failure: read the stage log and rerun with an appropriate heap. Reuse a completed rawspace to avoid repeating import.
- GUI startup failure: inspect `gui.log`, DISPLAY, paths, and GUI heap; preparation may still have succeeded.
- Unresponsive status: check heartbeat and worker identity; interrupted jobs are not silently declared successful.

## Scope

V1 is local, single-user, Linux-first. It provides preparation and GUI handoff, not direct agent-driven search, cluster scheduling, cleaning, downsampling, or 3D screening. The worker uses fixed CLI entry points and typed arguments; it is not a general shell tool. Use only trusted local serialized `.data` files in the GUI.

Documentation is also available as `hyperspace://docs/workflow`, `hyperspace://docs/rawspace`, and `hyperspace://docs/importers`. The `prepare_searchable_space` MCP prompt supplies the workflow recipe. Clients without resource/prompt UI can use `get_workflow_help`.

## Developer verification

Run the focused unit/regression tests:

```bash
mvn -pl openchemlib-hyperspace-cli,openchemlib-hyperspace-mcp -am -Dmaven.compiler.release=17 -Dtest=PreparationTest,WorkerTest,HyperspaceConfigurationTest,RawSynthonSpaceTest,RawSynthonSpaceImporterTest -Dsurefire.failIfNoSpecifiedTests=false test
```

After packaging both JARs, run `PackagedMcpTest` in the MCP module with absolute system properties `hyperspace.cliJar`, `hyperspace.mcpJar`, and `hyperspace.toy`. It tests SDK discovery, resources, input inspection, a toy preparation across client restart, cancellation, and headless GUI handoff. An optional `hyperspace.testManifest` output path records the completed job for subsequent verification.

`PreparedSpaceSmokeTest` in the GUI module accepts that `hyperspace.testManifest` and loads both indexes and searches for a known assembled product. Set `hyperspace.guiSmoke=true` to open and render the GUI. Optional `hyperspace.screenCapture=true` requests desktop screen-capture permission; normal tests do not capture the desktop. These artifact-dependent tests are skipped when their system properties are absent.

## Build threading

The MCP `threads` setting controls import, reaction-level search-index construction,
and similarity-index initialization. Reaction construction uses a fixed-size pool;
each running reaction creates its own temporary data. More workers can therefore
increase peak memory usage. A single large reaction can leave a single-core tail.
Final shared-index initialization and file serialization remain sequential.

For direct `RawSynthonSpaceBuildCLI` use, `--threads N` controls reaction workers
(default `1`, positive integers only). `--similarityThreads N` remains independent
(default `4`). Existing indexes do not need rebuilding. Restart the MCP server
after updating its jars; new jobs use the new settings.

## GUI selection

GUI launch remains opt-in. MCP now defaults to GUI 2, which supports **local
substructure search only**. Existing indexes and generated `gui.json` files are
reused without rebuilding. Preparation still builds the requested search modes.

For similarity search, or both search types in one window:
```json
{"configPath": "/absolute/path/to/gui.json", "guiVersion": "legacy"}
```
Pass this to `launch_gui`. Both `launch_gui` and `prepare_space` accept
`guiVersion: "gui2"` (default) or `"legacy"`. With `prepare_space`, launching
still requires `launchGui: true`.

Mixed configurations open their substructure providers in GUI 2 with a warning.
Similarity-only configurations require explicit selection of the legacy GUI;
there is no silent fallback. On headless hosts, the selected launch command is
returned instead.

Direct launch:
```bash
java -Xmx8G -jar openchemlib-hyperspace-cli.jar GUI2 /path/to/gui.json
java -Xmx8G -jar openchemlib-hyperspace-cli.jar GUI /path/to/gui.json
```

GUI 2 honors provider names and `MaxNumberOfThreads`, and resolves relative
index paths against the config directory. Standalone startup without arguments
still reads `spaces.conf` when present (one index path per line).
