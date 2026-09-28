#!/usr/bin/env bash
# Startet die Rechenwerkzeuge im Browser. Beim ersten Start wird alles eingerichtet.
set -e
cd "$(dirname "$0")"
PY=python3.12
command -v "$PY" >/dev/null 2>&1 || PY=python3
if [ ! -d .venv ]; then
  echo "Erster Start: richte Python-Umgebung ein …"
  "$PY" -m venv .venv
  .venv/bin/pip install -q -r requirements.txt
fi
exec .venv/bin/python -m skripte.web "$@"
