@echo off
rem Double-click, or run from PowerShell/cmd, to start the AstraFlow GUI through WSLg.
rem Works from the repository (\wsl.localhost\Ubuntu-24.04\...) or from a copy or shortcut anywhere.
rem Arguments go to astraflow_gui, e.g.  AstraFlow.cmd --config examples/rocket_nozzle/cooled_wall.json
setlocal
set "P=%~dp0"
set "P=%P:\wsl.localhost\Ubuntu-24.04=%"
set "P=%P:\=/%"
if "%P:~0,1%" neq "/" set "P=/home/arnav/Arnavs Folder/dev/projects/portfolio/AstraFlow/"
wsl.exe -d Ubuntu-24.04 -e bash -lc "cd '%P%' && exec ./build/astraflow_gui %*"
