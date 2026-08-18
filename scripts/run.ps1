# Run GAMES102 (sets the Qt DLL path, reuses the build recorded by build.ps1)
# Usage: .\scripts\run.ps1
#
# NOTE: keep this file ASCII-only so it parses on any Windows locale / editor.
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot

$meta = Join-Path $root "build\last_build.json"
if (-not (Test-Path $meta)) {
    throw "Not built yet. Run .\scripts\build.ps1 first."
}
$info = Get-Content $meta | ConvertFrom-Json
if (-not (Test-Path $info.exe)) {
    throw "Program missing: $($info.exe)`nPlease re-run .\scripts\build.ps1"
}

$env:PATH = "$($info.qtDir)\bin;" + $env:PATH
Write-Host "Starting: $($info.exe)" -ForegroundColor Cyan
& $info.exe
