# POST_DATA2 Architecture

## Public entry points

- `postdata_app` starts the desktop application.
- `postdata_startup` adds ordinary source folders to the MATLAB path.
- `postdata_run` is the only calculation entry point.

## Module boundaries

```text
app      -> core, postdata_run, plot, export
postdata_run -> core, io, analysis
analysis -> io
plot     -> core
export   -> core
io/core  -> MATLAB only
```

Reverse dependencies are not allowed. Analysis functions do not read GUI
controls or write exported files. The GUI always disables analyzer-owned
legacy figures and renders through `src/plot`; direct analyzer entry points
retain optional standalone figures only for script compatibility.

### `src/app`

Owns UI state and callbacks only. `PostDataApp` converts controls into request,
plot, and output structures, then delegates work to services.

### `src/core`

Owns versioned requests, validation, option catalogs, safe value parsing,
selection conversion, execution checkpoints, logging, configuration upgrades,
shared compatibility helpers, and UI/result-view localization. Executable
MATLAB sources remain ASCII-only; `pd_ui_text` explicitly decodes the UTF-8
catalog under `src/core/resources` so Windows 7 / MATLAB R2016b never relies on
the process code page. The option catalog is
the single source of truth for all editable calculation parameters.
Catalog metadata also defines allowed values, numeric constraints, required
values, vector length, and dependencies. The GUI and dispatcher both consume
the same metadata, so inactive parameters cannot leak into calculations.

### `src/io`

Owns input discovery, indexed chunk reading, Slurm parsing, and progress
reporting. Input preflight validates headers and required variables before
expensive calculation. Header parsing is centralized and rejects variable
names that collapse to the same MATLAB field. Indexed reads expose cooperative
cancellation, rebuild corrupt/stale sidecars, write caches atomically, and
return the in-memory index so analysis never depends on cache persistence. This
module does not contain analysis formulas.

### `src/analysis`

Owns numerical calculations. Each analyzer accepts explicit values and returns
a result structure. `network2d` has its own submodule because its grid,
geometry, topology, profile, and evolution concerns are substantially larger.
Mass-vx and mass-x both delegate density-column resolution, total-mass
reconstruction, sorting, and directional accumulation to
`pd_cumulative_areal_density`; analysis adapters only define their coordinate.
Shared descriptive statistics live in `analysis/common`; network evolution
selection and aggregation live in `network2d/pd_network_build_evolution`, while
the snapshot callback remains responsible for physical calculations.
Positive-valued histogram, PDF/CDF, empty-bin handling, automatic bin width,
and power-law/Gamma/lognormal curves are centralized in
`analysis/common/pd_distribution_build`. Analysis-specific adapters may rename
fields but must not reimplement distribution formulas.
Digital Betti/Euler topology, periodic 4-neighbor traversal, and connected
component winding detection are isolated in `pd_network_compute_topology`,
`pd_network_neighbor`, and `pd_network_label_components`.
Binary and PLIC cut-cell reconstruction is isolated in
`pd_network_build_geometry`. It owns phase fractions, centroids, interface
segments, perimeters, and explicit binary/error fallback behavior. The snapshot
analyzer consumes that contract and does not contain a second geometry pipeline.
Per-component measurements are isolated in
`pd_network_compute_component_stats`. It consumes labels and optional cut-cell
geometry, then returns the stable component table, size ranking, and explicit
open/periodic connectivity status used by the phase analyzer and plots.
Directional profile aggregation and directional connectivity classification are
isolated in `pd_network_build_directional_profiles` and
`pd_network_build_directional_connectivity`. Shared bin-edge and periodic-axis
semantics live in `pd_network_locate_bin` and `pd_network_is_periodic_axis`, so
geometry, topology, components, profiles, and plots use the same boundary rules.
Advanced paper-oriented result assembly lives in
`pd_network_build_advanced_stats`. Image Processing Toolbox-dependent skeleton,
branch/end-point, periodic tiling, and thickness calculations are isolated in
`pd_network_build_morphology_stats`, which exposes deterministic unavailable and
disabled states instead of leaking optional-toolbox failures into the analyzer.

### `src/plot`

Owns rendering into caller-provided axes. Plot conditions are parsed separately
from calculation parameters, so replotting does not repeat calculations.
`pd_result_plot_views` declares the views available from each result contract;
the application can render one selected view or arrange every available view
as a dashboard without adding graphics logic to analyzers.
`UpdateMode=replace|overlay` controls whether caller-owned axes are reset or
receive an additional labeled result layer. Network2d standalone snapshot and
evolution figures are implemented in `pd_render_network2d_snapshot` and
`pd_render_network2d_evolution`; the numerical analyzer contains no graphics
calls. Fixed plot options expose allowed values to the R2016b typed GUI editor.

### `src/export`

Owns result summaries, analysis-specific detail tables, MAT/CSV/PNG/FIG/PDF
output, reproducibility manifests, and collision-free file naming. A complete
export is first built in an isolated staging directory and committed only after
every requested format succeeds. Export settings are independent from analysis
and plotting settings.

## Extension rules

1. Add calculation options to `pd_option_catalog`.
2. Implement the calculation in `src/analysis`.
3. Add one dispatch case to `postdata_run` only when adding a new analysis type.
4. Add rendering in `pd_render_result` without changing calculation code.
5. Add output formats in `src/export` without changing analyzers or the GUI.
6. Add integration tests in `tests/test_postdata.m`.

All production code must remain compatible with MATLAB R2016b. Executable
`.m` files must remain ASCII-only; localized text belongs in the UTF-8 resource
catalog and is retrieved with `pd_ui_text`.
