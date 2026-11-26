# Regenerate everything in Results_Main:  data/stats CSVs (Python), then figures (MATLAB).
#   powershell -ExecutionPolicy Bypass -File "Export_scripts\run_export.ps1"
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$matlab = 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe'

python -W ignore (Join-Path $here 'export_results_main_data.py')
if ($LASTEXITCODE -ne 0) { throw 'export_results_main_data.py failed' }

& $matlab -batch "cd('$here'); make_results_main_figures"
if ($LASTEXITCODE -ne 0) { throw 'make_results_main_figures.m failed' }
