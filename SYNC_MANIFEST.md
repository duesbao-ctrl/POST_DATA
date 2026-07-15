# Non-GUI synchronization

`POST_DATA2` is the canonical superset. Run this command after a verified core
change to update the sibling script-only project:

```powershell
powershell -ExecutionPolicy Bypass -File tools\sync_non_gui_to_post_data.ps1
```

The whitelist is `postdata_run.m`, `postdata_startup.m`, `run_analysis.m`,
`CHANGELOG.md`, and
`src/{analysis,core,export,io,plot}` plus examples, fixtures, and tests. The
script deliberately excludes `postdata_app.m` and `src/app`, removes the old
root-level duplicate implementations, and cleans MATLAB index/backup files.
The recursive `src/core` copy includes `pd_ui_text.m` and the UTF-8
`resources/ui_zh_CN.tsv` catalog required by localized result-view labels.
