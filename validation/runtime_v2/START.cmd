@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
if not exist "config.json" (
  copy /y "config.example.json" "config.json" >nul
  echo Created config.json. Check its three paths before choosing Import or Run.
)
if exist "F:\Dev\Python312\python.exe" (
  "F:\Dev\Python312\python.exe" -X utf8 "%~dp0resume.py" menu
) else (
  where py >nul 2>nul
  if errorlevel 1 (
    python -X utf8 "%~dp0resume.py" menu
  ) else (
    py -3.12 -X utf8 "%~dp0resume.py" menu
  )
)
if errorlevel 1 pause
endlocal
