# RECURRING ISSUES

## Environment pitfalls

- On this Windows workspace, `rg.exe` may fail with `Access is denied` even when the repo is readable. When that happens, use PowerShell-native fallbacks such as `Select-String`, `Get-ChildItem`, and `Get-Content` instead of retrying `rg`.

## Validation traps

- In `network2d`, `PositionRangeX/Y` and `ProfileRangeX/Y` are analysis-only filters. They do not automatically crop 2D phase/label/connectivity figures unless dedicated plot-range handling is wired in.
