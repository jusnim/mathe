@echo off
rem Startet die Rechenwerkzeuge im Browser. Beim ersten Start wird alles eingerichtet.
cd /d "%~dp0"
if not exist .venv (
  echo Erster Start: richte Python-Umgebung ein ...
  py -3.12 -m venv .venv 2>nul || python -m venv .venv
  .venv\Scripts\pip install -q -r requirements.txt
)
.venv\Scripts\python -m skripte.web %*
