param(
    [string]$TargetRoot = ''
)

$ErrorActionPreference = 'Stop'
$sourceRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
if ([string]::IsNullOrWhiteSpace($TargetRoot)) {
    $TargetRoot = Join-Path (Split-Path $sourceRoot -Parent) 'POST_DATA'
}
$targetRoot = [System.IO.Path]::GetFullPath($TargetRoot)

if ((Split-Path $sourceRoot -Leaf) -ne 'POST_DATA2') {
    throw "Refusing to sync from unexpected source: $sourceRoot"
}
if ((Split-Path $targetRoot -Leaf) -ne 'POST_DATA') {
    throw "Refusing to sync to unexpected target: $targetRoot"
}
if ((Split-Path $sourceRoot -Parent) -ne (Split-Path $targetRoot -Parent)) {
    throw 'POST_DATA2 and POST_DATA must be sibling directories.'
}
if (-not (Test-Path -LiteralPath (Join-Path $targetRoot '.git'))) {
    throw "Target is not the expected Git project: $targetRoot"
}

function Assert-TargetChild([string]$PathValue) {
    $full = [System.IO.Path]::GetFullPath($PathValue)
    $prefix = $targetRoot.TrimEnd('\') + '\'
    if (-not $full.StartsWith($prefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing filesystem operation outside POST_DATA: $full"
    }
    return $full
}

$moduleDirectories = @('analysis', 'core', 'export', 'io', 'plot')
$replaceDirectories = @('src', 'examples', 'fixtures', 'tests')
foreach ($name in $replaceDirectories) {
    $destination = Assert-TargetChild (Join-Path $targetRoot $name)
    if (Test-Path -LiteralPath $destination) {
        Remove-Item -LiteralPath $destination -Recurse -Force
    }
}

$targetSrc = Assert-TargetChild (Join-Path $targetRoot 'src')
New-Item -ItemType Directory -Path $targetSrc | Out-Null
foreach ($name in $moduleDirectories) {
    Copy-Item -LiteralPath (Join-Path (Join-Path $sourceRoot 'src') $name) -Destination $targetSrc -Recurse
}
foreach ($name in @('examples', 'fixtures', 'tests')) {
    Copy-Item -LiteralPath (Join-Path $sourceRoot $name) -Destination (Join-Path $targetRoot $name) -Recurse
}
foreach ($name in @('.gitignore', 'postdata_run.m', 'postdata_startup.m',
        'run_analysis.m', 'RECURRING_ISSUES.md', 'CHANGELOG.md')) {
    Copy-Item -LiteralPath (Join-Path $sourceRoot $name) -Destination (Join-Path $targetRoot $name) -Force
}

$legacyFiles = @(
    'analysis_usage_examples.m', 'analyze_chunk_field.m',
    'analyze_chunk_network2d.m', 'cluster_postprocess.m',
    'compute_temp_stress_chunk.m', 'demo_network2d_periodic_topology.m',
    'fix_rdp_3389.ps1', 'get_slurm_txt_fullpath.m', 'make_file_progress.m',
    'matlab_probe.txt', 'plot_cloud2d.m', 'plot_line1d.m',
    'read_bin_chunk.m', 'read_chunk_step_fast.m', 'read_slurm_stepcpu.m',
    'vx_chunk_cumulative.m', 'postdata_app.m'
)
foreach ($name in $legacyFiles) {
    $pathValue = Assert-TargetChild (Join-Path $targetRoot $name)
    if (Test-Path -LiteralPath $pathValue) {
        Remove-Item -LiteralPath $pathValue -Force
    }
}
foreach ($name in @('demo_case_generated', 'selftest_generated')) {
    $pathValue = Assert-TargetChild (Join-Path $targetRoot $name)
    if (Test-Path -LiteralPath $pathValue) {
        Remove-Item -LiteralPath $pathValue -Recurse -Force
    }
}

Get-ChildItem -LiteralPath $targetRoot -Recurse -File |
    Where-Object { $_.Name -like '*.stepidx.mat' -or
        $_.Extension -in @('.asv', '.tmp', '.bak') } |
    ForEach-Object {
        $safePath = Assert-TargetChild $_.FullName
        Remove-Item -LiteralPath $safePath -Force
    }

Write-Output "Synced non-GUI modules from $sourceRoot"
Write-Output "Target: $targetRoot"
