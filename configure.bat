@echo off
setlocal
rem ============================================================================
rem  GAMES102 - one-click CMake configure
rem
rem  Generator : Visual Studio 17 2022  (x64)
rem  Toolset   : v142                   (MSVC2019, matches Qt 5.15.2 msvc2019_64)
rem
rem  Usage:
rem      configure.bat                      -> uses C:\Qt\5.15.2\msvc2019_64
rem      configure.bat D:\Qt\5.15.2\msvc2019_64
rem
rem  Output: build\msvc2022-v142\GAMES102.sln
rem  NOTE: keep this file ASCII-only so it works on any Windows locale.
rem ============================================================================

rem Make this script location-independent (works when double-clicked)
cd /d "%~dp0"

set "QT_DIR=%~1"
if "%QT_DIR%"=="" set "QT_DIR=C:\Qt\5.15.2\msvc2019_64"

if not exist "%QT_DIR%\lib\cmake\Qt5Widgets" (
    echo [ERROR] Qt 5.15 msvc2019_64 not found at: "%QT_DIR%"
    echo.
    echo Usage: configure.bat [path-to-qt-msvc2019_64]
    echo Example: configure.bat D:\Qt\5.15.2\msvc2019_64
    exit /b 1
)

echo [1/2] CMake configure ^(VS2022 generator, v142 toolset, Qt=%QT_DIR%^)
set "QT_PREFIX=%QT_DIR%"
cmake --preset msvc2022-v142
if errorlevel 1 (
    echo.
    echo [ERROR] CMake configure failed.
    echo Hint: if it reports the v142 toolset is missing, open "Visual Studio
    echo        Installer", click Modify on VS2022, and install the component
    echo        "MSVC v142 - VS 2019 C++ x64/x86 build tools".
    exit /b 1
)

echo.
echo [OK] Configured: build\msvc2022-v142\GAMES102.sln
echo.
echo Next steps:
echo     cmake --build build\msvc2022-v142 --config Debug
echo     set PATH=%QT_DIR%\bin;%%PATH%%
echo     build\msvc2022-v142\Debug\GAMES102.exe
echo.
echo Or open build\msvc2022-v142\GAMES102.sln in Visual Studio and press F5.
endlocal
