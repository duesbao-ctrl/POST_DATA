# POST_DATA

MATLAB post-processing helpers for chunk, cluster, and velocity analysis workflows.

## Main entry point

Use `run_analysis(taskType, ...)` as the unified interface:

- `taskType = 'chunk'` for chunk field extraction and plotting
- `taskType = 'cluster'` for cluster statistics and post-processing
- `taskType = 'vx'` for cumulative velocity analysis

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
- `selftest_generated/run_selftest.m` provides a lightweight self-test using bundled fixtures

## Quick start

1. Open MATLAB in this repository root.
2. Review `analysis_usage_examples.m`.
3. Run the examples or call `run_analysis(...)` directly with your case directory.
