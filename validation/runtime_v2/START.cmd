@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
set "COV_CONFIG=%~dp0config.json"
if exist "%~dp0..\..\..\..\active.json" set "COV_CONFIG=%~dp0..\..\..\..\config.json"
if not exist "%COV_CONFIG%" (
  copy /y "config.example.json" "%COV_CONFIG%" >nul
  echo Created config.json. Check its three paths before choosing Import or Run.
)
if exist "F:\Dev\Python312\python.exe" (
  "F:\Dev\Python312\python.exe" -X utf8 "%~dp0resume.py" --config "%COV_CONFIG%" menu
) else (
  where py >nul 2>nul
  if errorlevel 1 (
    python -X utf8 "%~dp0resume.py" --config "%COV_CONFIG%" menu
  ) else (
    py -3.12 -X utf8 "%~dp0resume.py" --config "%COV_CONFIG%" menu
  )
)
if errorlevel 1 pause
endlocal
