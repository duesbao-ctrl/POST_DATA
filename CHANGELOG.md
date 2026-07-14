# Changelog

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
