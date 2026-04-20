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

For periodic boundaries, hole count and Euler characteristic are returned as `NaN`
because the current version does not define topological holes on the wrapped domain.
Also note that a component may merge across a periodic seam without being treated
as directionally connected; `wrapsX` / `wrapsY` require a true non-zero winding
around the periodic domain.

## Quick start

1. Open MATLAB in this repository root.
2. Review `analysis_usage_examples.m`.
3. Run the examples or call `run_analysis(...)` directly with your case directory.
