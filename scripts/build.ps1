# GAMES102 one-click build script
# Usage:
#   .\scripts\build.ps1                       # Debug + C:\Qt\5.15.2\msvc2019_64
#   .\scripts\build.ps1 -Config Release
#   .\scripts\build.ps1 -QtDir D:\Qt\5.15.2\msvc2019_64
#   .\scripts\build.ps1 -ConfigureOnly        # configure only
#
# NOTE: keep this file ASCII-only so it parses on any Windows locale / editor.
param(
    [string]$QtDir = "C:\Qt\5.15.2\msvc2019_64",
    [ValidateSet("Debug", "Release")][string]$Config = "Debug",
    [switch]$ConfigureOnly
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot

# ---- Validate Qt ----
if (-not (Test-Path "$QtDir\lib\cmake\Qt5Widgets")) {
    throw "Qt 5.15 (msvc2019_64) not found at: $QtDir`nPass -QtDir, e.g. .\scripts\build.ps1 -QtDir C:\Qt\5.15.2\msvc2019_64"
}

# ---- Detect Visual Studio and pick the CMake generator ----
$vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $vswhere)) {
    throw "vswhere.exe not found. Install Visual Studio 2019/2022 with the 'Desktop development with C++' workload."
}
$vsProps = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -format json | ConvertFrom-Json
if (-not $vsProps) {
    throw "No Visual Studio with C++ toolchain detected. Install VS2019/VS2022 'Desktop development with C++'."
}
$vsVersion = $vsProps.installationVersion
$vsPath = $vsProps.installationPath
$major = [int]($vsVersion.Split('.')[0])
if ($major -ge 17)      { $preset = "msvc2022"; $buildDir = "build\msvc2022" }
elseif ($major -eq 16)  { $preset = "msvc2019"; $buildDir = "build\msvc2019" }
else                    { throw "Unsupported Visual Studio version: $vsVersion (need 2019 or 2022)" }

# CMakePresets.json reads the Qt path from this env var
$env:QT_PREFIX = $QtDir

Write-Host "[1/3] CMake configure (preset=$preset, Qt=$QtDir)" -ForegroundColor Cyan
Push-Location $root
try {
    cmake --preset $preset
    if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }
} finally { Pop-Location }

if ($ConfigureOnly) {
    Write-Host "Configure done. Build directory: $buildDir" -ForegroundColor Green
    exit 0
}

Write-Host "[2/3] Build ($Config)" -ForegroundColor Cyan
cmake --build (Join-Path $root $buildDir) --config $Config
if ($LASTEXITCODE -ne 0) { throw "Build failed" }

$exe = Join-Path $root "$buildDir\$Config\GAMES102.exe"
if (-not (Test-Path $exe)) { throw "Output not found: $exe" }

# ---- Deploy Qt/VC runtime DLLs next to the exe with windeployqt ----
Write-Host "[3/3] windeployqt" -ForegroundColor Cyan
$wdq = Join-Path $QtDir "bin\windeployqt.exe"
if (Test-Path $wdq) {
    # Let windeployqt find the VC runtime (msvcp140.dll etc.)
    if ($vsPath) { $env:VCINSTALLDIR = Join-Path $vsPath "VC\" }
    & $wdq --no-translations --compiler-runtime $exe | Out-Null
}

# Remember this build so run.ps1 can use it
New-Item -ItemType Directory -Force -Path (Join-Path $root "build") | Out-Null
@{ exe = $exe; qtDir = $QtDir; config = $Config } | ConvertTo-Json |
    Set-Content (Join-Path $root "build\last_build.json") -Encoding UTF8

Write-Host ""
Write-Host "Build OK: $exe" -ForegroundColor Green
Write-Host "Run it with:  .\scripts\run.ps1"
