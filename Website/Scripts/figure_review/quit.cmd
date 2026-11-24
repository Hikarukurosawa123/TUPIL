@echo off
py -B "%~dp0background.py" quit
set "result=%errorlevel%"
if /i not "%~1"=="--no-pause" pause
exit /b %result%
