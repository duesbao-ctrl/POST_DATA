# POST_DATA

MATLAB post-processing helpers for chunk, cluster, and velocity analysis workflows.

## Main entry point

Use `run_analysis(taskType, ...)` as the unified interface:

- `taskType = 'chunk'` for chunk field extraction and plotting
- `taskType = 'cluster'` for cluster statistics and post-processing
- `taskType = 'vx'` for cumulative velocity analysis
- `taskType = 'network2d'` for pore/matrix network analysis on 2D chunk grids

For coordinate-based chunk and `network2d` outputs, `CoordScale` can be used
to rescale coordinates before plotting, binning, and range filtering.

See `analysis_usage_examples.m` for end-to-end examples.

## Read progress

`read_chunk_step_fast.m` builds a `.stepidx.mat` cache on the first slow scan,
then reuses it on later reads. Use `ProgressMode` with `run_analysis(...)`,
`read_chunk_step_fast(...)`, or `read_bin_chunk(...)` to control progress
display:

- `auto` uses `waitbar` on desktop MATLAB and console output in batch mode
- `off` disables progress output
- `waitbar` forces GUI progress
- `console` forces command-window progress

## Included files

- `run_analysis.m` coordinates task selection and common options
- `read_chunk_step_fast.m` reads chunk data by index, timestep, or physical time
- `cluster_postprocess.m` computes cluster statistics
- `vx_chunk_cumulative.m` processes cumulative velocity outputs
- `analyze_chunk_network2d.m` reconstructs 2D pore/matrix networks from `Ncount`
- `selftest_generated/run_selftest.m` provides a lightweight self-test using bundled fixtures

## Cluster notes

`taskType='cluster'` reports equivalent-diameter statistics and distribution
plots. Use `DiameterRange = [min max]` to choose the diameter scale window used
by the statistics and distribution plots. Use `DiameterHistBinSize` to set the
physical bin width of the diameter histogram; leave it empty to use an automatic
data-driven bin width.
`DiameterEmptyBinMode` controls empty histogram bins: `zero` keeps zero counts,
`nan` marks empty-bin count/probability values as `NaN`, and `remove` drops
empty bins from the plotted/statistical distribution arrays while retaining the
raw bin fields. `DiameterPlotStyle` selects `bar` or `scatter` for the diameter
distribution plots.
`DiameterFitTypes` controls which fitted curves are overlaid on the distribution
plots: `powerlaw`, `gamma`, `lognormal`, `all`, or `none`. `HistScale` controls
the axis scaling and supports `linear`, `semilogx`, `semilogy`, and `loglog`.

## Network2d notes

`taskType='network2d'` classifies each valid 2D chunk cell by `Ncount < ThresholdN`:

- pore: `Ncount < ThresholdN`
- matrix: `Ncount >= ThresholdN`

The module reports 2D morphology statistics for both phases, including:

- porosity and matrix area fraction
- interface length and 2D specific interface measures
- connected-component counts, largest-component fractions, percolation/wrap flags
- hole count and Euler characteristic for open boundaries
- number-size distributions and mean-size-vs-position distributions
- directional profiles of porosity, specific interface, and perpendicular connectivity

By default, `network2d` now applies a cut-cell geometric correction for
area/interface/diameter statistics and phase-map plotting. The binary
classification above is still the source of topology, connectivity,
percolation, skeleton, and thickness masks. The cut-cell layer uses
`phi = ThresholdN - Ncount` and reconstructs the zero interface inside each
coarse chunk cell with piecewise-linear center-to-corner triangles. It does
not create particle subcells or assume that smaller display cells contain
particles; `CutCellPlotRefinement` only smooths the visual phase map.

Use `GeometryMode` to choose the geometric statistics/plotting mode:

- `GeometryMode = 'cutcell'` (default) uses the cut-cell correction
- `GeometryMode = 'original'` uses the original strict 0/1 grid cells
- `CutCellMethod = 'plic'` uses the piecewise-linear cut-cell reconstruction
- `CutCellFallback = 'binary'` falls back to the original 0/1 cell if local
  reconstruction is underdetermined; use `'error'` to fail instead
- `CutCellPlotRefinement = 4` controls display-only interpolation of the phase map

The fractional fields are returned in `out.cutCell.pore.fraction` and
`out.cutCell.matrix.fraction`. A value between 0 and 1 means that the original
coarse cell is cut by the reconstructed pore-matrix interface.

For paper-grade analysis, the module now also exposes a unified `out.stats`
structure that keeps the legacy outputs intact while organizing the main
topological, geometric, and evolution metrics for downstream scripts.

### `out.stats` overview

`out.stats.meta` records the analysis convention used by the snapshot:

- `boundary`
- `foregroundConnectivity = 4`
- `backgroundHoleConnectivity = 8`
- mean definition parameters `m` and `n`
- cut-cell settings and a note that topology/connectivity still use the binary mask
- toolbox availability notes for skeleton/thickness paths

`out.stats.topology` reports the 2D Betti numbers and Euler characteristic for
both phases:

- `stats.topology.pore.beta0`, `beta1`, `chi`
- `stats.topology.matrix.beta0`, `beta1`, `chi`

This version fixes the convention to foreground 4-connectivity and background
8-connectivity for the open-boundary hole interpretation. For periodic
boundaries the code switches to the wrapped 4-neighbor digital cell complex,
so `beta1` is still well defined and now counts all non-trivial 1D loops of
the wrapped structure, including both ordinary enclosed holes and periodic
wrapping loops.

`out.stats.connectivity` summarizes the publication-oriented connectivity
metrics for pore and matrix:

- `largestFraction`
- `largestArea`
- `componentCount`
- `percolatesX`, `percolatesY`
- `wrapsX`, `wrapsY` for periodic cases

`out.stats.geometry` collects the core interface measures:

- `phi`
- `interfaceLength`
- `specificInterface = interfaceLength / validArea`
- `poreArea`, `matrixArea`, `matrixFraction`, `validArea`

When cut-cell correction is enabled, these geometry fields use the fractional
cell areas and reconstructed interface segments. Disable it with
`GeometryMode='original'` to recover the strict binary grid-edge estimates.

`out.stats.size` stores per-component equivalent-diameter vectors and histogram
data for both phases:

- `diameter = 2 * sqrt(area/pi)`
- `hist.diameter` with `edges / centers / count / pdf`

Average diameter statistics use the generalized mean
`sum(d^m) / sum(d^n)` through `MeanPowerM` and `MeanPowerN` in
`NetworkOptions`. The default `m=1, n=0` reduces to the arithmetic mean.
The area-weighted mean diameter is reported separately as a fixed physical
quantity. `DiameterRange = [min max]` selects the component diameter window used
by size statistics, diameter histograms, and diameter distribution plots; global
topology/connectivity/porosity are still computed from the full network.
`DiameterHistBinSize` sets the physical histogram bin width; if it is empty, a
data-driven bin width is used. Empty diameter bins are controlled with
`DiameterEmptyBinMode = 'zero'`, `'nan'`, or `'remove'`; raw bin arrays remain
available as `rawCount`, `rawProbability`, and related fields when bins are
hidden or converted to `NaN`. Diameter distribution plots can be drawn as bars
or scatter points with `DiameterPlotStyle`, can overlay power-law, gamma, and
lognormal fits via `DiameterFitTypes`, and can be shown with
`HistScale = 'linear'`, `'semilogx'`, `'semilogy'`, or `'loglog'`.

`out.stats.network` and `out.stats.thickness` provide optional higher-level
descriptors when the Image Processing Toolbox is available:

- `stats.network.<phase>.skeletonLength`
- `stats.network.<phase>.branchPoints`
- `stats.network.<phase>.endPoints`
- `stats.network.<phase>.branchDensity`
- `stats.thickness.matrix.min`, `mean`, `p1`, `p5`

Skeleton and thickness analysis are enabled with:

- `EnableSkeletonGraph = true`

If the toolbox is missing, topology/connectivity/geometry/size still work and
the dependent fields return `NaN` together with an explanatory note.

`out.stats.fragmentation.matrix` is a dedicated alias for matrix breakup:

- `count`
- `largestFraction`
- `sizeDist`

This is meant to be used directly for fragmentation-stage analysis without
having to rebuild it from the generic topology fields.

### Physical interpretation

The most useful paper-level interpretation is:

- `beta0`: how many disconnected pieces a phase has
- `beta1`: how many loops / enclosed rings that phase contains
- `chi = beta0 - beta1`: a compact topology summary

Typical reading rules are:

- pore coalescence: pore `beta0` decreases and pore `largestFraction` increases
- network formation: pore `beta1` increases, often together with more branch points
- matrix breakup: matrix `beta0` increases sharply and matrix `largestFraction` drops
- complexification: `phi` increases while `specificInterface` also increases
- coarsening / smoothing: `phi` increases while `specificInterface` decreases

For skeleton metrics:

- many `endPoints` with small `beta1` usually indicate a tree-like structure
- large `branchPoints` together with large `beta1` usually indicate a mesh-like network

For matrix thickness:

- `min` is the most sensitive necking indicator
- `p1` and `p5` are more robust low-thickness indicators than a single-point minimum
- a drop in `min/p1/p5` is often an early warning for impending matrix rupture

### Why periodic `beta1` is different from ordinary hole count

Under open boundaries, a "hole" is easy to define: it is a bounded background
region enclosed by the phase. That is why `beta1` can be counted reliably from
the complemented image.

Under periodic boundaries, the left/right and/or top/bottom edges are glued
together. Geometrically the domain is no longer a flat rectangle with edges:

- periodic in one direction behaves like a cylinder
- periodic in both directions behaves like a torus

On such wrapped domains, two different situations must be distinguished:

- an ordinary enclosed hole inside the phase
- a non-contractible loop that wraps around the periodic domain

These are not the same object. A structure can cross a periodic seam and still
have no ordinary enclosed hole, or it can wrap around the torus without
surrounding a bounded cavity in the open-boundary sense.

That is why the periodic implementation does **not** reuse the open-boundary
"count bounded background islands" rule. Instead it computes Betti numbers on a
wrapped digital cell complex consistent with 4-neighbor foreground topology:

- 0-cells: occupied pixels
- 1-cells: occupied horizontal/vertical adjacencies
- 2-cells: occupied 2x2 plaquettes

With that convention:

- `beta0` is still the number of connected components
- `beta1` is the rank of all 1D loops in the wrapped domain
- `chi = beta0 - beta1 + beta2`

So under periodic boundaries, `beta1` is no longer just an "ordinary hole
count". It includes:

- ordinary enclosed holes
- loops that wrap around the periodic direction(s)

Examples:

- a fully occupied `periodic-x` strip behaves like a cylinder, so `beta1 = 1`
- a fully occupied `periodic-xy` field behaves like a torus, so `beta1 = 2`

This is the mathematically correct periodic interpretation and is the one used
by `out.stats.topology.*` and the legacy `holeCount/eulerCharacteristic` fields.

### Evolution scan

`network2d` can optionally scan a timestep/index/time interval and append the
result to `out.stats.evolution` while keeping the main output focused on the
selected snapshot.

Use these `NetworkOptions`:

- `EnableEvolution`
- `EvolutionSelectBy = 'Index' | 'TimeStep' | 'Time'`
- `EvolutionRange = [min max]`
- `EvolutionStride`

The evolution output includes:

- `geometry.phi`, `geometry.interfaceLength`, `geometry.specificInterface`
- `topology.pore/matrix.beta0`, `beta1`, `chi`
- `connectivity.pore/matrix.largestFraction`, `percolatesX`, `percolatesY`
- `thickness.matrix.min`, `mean`, `p1`, `p5`
- `fragmentation.matrix.count`, `largestFraction`
- transition deltas such as `deltaPhi`, `deltaSpecificInterface`,
  `deltaBeta0`, `deltaBeta1`, and `deltaLargestFraction`

The 2D "specific surface area" is defined here as pore-matrix interface length per
area, with bulk, pore-only, and matrix-only normalizations.

Directional profiles are built from slices along a chosen axis (`ProfileAxis`,
default `x`). The local porosity and specific-interface curves are computed from
cells/interfaces inside the requested profile range. The connectivity curve is
computed differently: the code first identifies globally connected pore
components along the perpendicular direction, then assigns each connected
component to a profile bin by its centroid. In each bin the module reports:

- local porosity
- local specific interface (bulk normalization)
- centroid-assigned pore connectivity across the perpendicular direction

`PositionAxis` and `ProfileAxis` answer different questions:

- `PositionAxis` bins whole connected components by their centroids and reports
  mean component size versus position.
- `ProfileAxis` slices the domain spatially and reports local field-like
  quantities such as porosity, interface density, and perpendicular
  connectivity versus position.

For periodic boundaries, note that topological loops and directional
connectivity are not the same thing. A component may merge across a periodic
seam without being directionally connected; `wrapsX` / `wrapsY` still require a
true non-zero winding around the periodic domain.

Coordinate ranges in `network2d` are intentionally split by purpose:

- `PositionRangeX/Y` filters only the mean-size-vs-position statistics
- `ProfileRangeX/Y` filters only the directional profile curves
- `PlotRangeX/Y` crops only the 2D phase/label/connectivity figures

This means changing `PositionRange*` or `ProfileRange*` will not zoom the 2D
phase map. Use `PlotRangeX/Y` when you want the displayed x/y limits of the
2D network figures to change.

## Quick start

1. Open MATLAB in this repository root.
2. Review `analysis_usage_examples.m`.
3. Run the examples or call `run_analysis(...)` directly with your case directory.
