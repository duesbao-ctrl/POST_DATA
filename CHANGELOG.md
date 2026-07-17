# Changelog

## 0.6.0 - 2026-07-16

- Added current SPID chunk parsing for extensible metadata headers and
  per-frame `# Time` values while retaining legacy three-line input support.
- Added SPID `spatial`, `field`, and `cluster` compatibility, including 3D
  spatial field inspection and old/new cluster coordinate aliases.
- Added metadata-driven conversions for SPID length, velocity, density, and
  areal-density units; legacy files retain their manual conversion defaults.
- Added direct physical-time selection without Slurm when SPID embeds time.
- Added regression fixtures for SPID spatial 1D/3D, field mass-v, and cluster
  output.

## 0.5.0 - 2026-07-16

- Replaced mass-x areal-density-column input with the SPH `Ncount` physical
  model based on initial density, particle spacing, raw length unit, and y
  normalization width.
- Added bin1d full-width and bin2d full-width/multi-slice mass-x for 2D/3D
  SPH, including z-width normalization for 3D, fractional y-bin overlap,
  particle-count plots, and 10 um coordinate defaults.
- Kept SPH raw fields and mass-x independent from the MD `dV` derived-field path.
- Added a Plot Data copy option that can omit column headers while keeping CSV
  export headers unchanged.

## 0.4.0 - 2026-07-16

- Added a paged Plot Data tab for every available result view.
- Added selected-row and full-table clipboard copy plus current-view CSV export.
- Added `SavePlotDataCSV` to export the exact numeric table behind every view.
- Centralized derived histogram, CDF, mean-profile, mass-distribution, and
  network-grid values in `pd_result_plot_data`; rendering and export now share
  the same numeric contract.
- Moved generated verification dashboards under ignored `outputs/verification`
  and removed committed/distributed verification output.

## 0.3.1 - 2026-07-15

- Removed non-ASCII literals from executable MATLAB source files.
- Added an explicitly decoded UTF-8 Chinese UI resource with an English fallback.
- Added regression checks for Chinese Unicode values and ASCII-only M-files on MATLAB R2016b.

## 0.3.0 - 2026-07-14

- Unified programmatic and GUI calculation-option validation.
- Added ambiguous chunk-header and first-row preflight checks.
- Added cooperative cancellation to indexed chunk reads.
- Rebuilt corrupt index caches safely and changed cache writes to atomic replacement.
- Removed network evolution's dependency on a writable sidecar cache.
- Added transactional multi-format export and an isolated `outputs/` default.
- Added GUI preflight, environment diagnostics, detailed error logging, and log rotation.
- Retained MATLAB R2016b compatibility and expanded regression coverage.

## 0.2.0

- Introduced the modular request, analysis, rendering, export, and GUI architecture.
- Added chunk, cluster, mass-v, mass-x, and network2d workflows with large fixtures.
