@echo off
py -B "%~dp0background.py" start
set "result=%errorlevel%"
if /i not "%~1"=="--no-pause" pause
exit /b %result%
