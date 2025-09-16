@echo off
setlocal

REM Usage:
REM   git-quick-push.bat
REM   git-quick-push.bat Your commit message here

set "MSG=%*"
if "%MSG%"=="" set "MSG=latest updates"

echo Staging changes...
git add -A

echo Committing ("%MSG%")...
git commit -m "%MSG%" >NUL 2>&1
if errorlevel 1 (
  echo Nothing to commit or commit skipped; continuing.
)

echo Pushing...
git push -u origin main

pause

endlocal

