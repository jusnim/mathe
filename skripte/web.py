"""Web-Oberfläche für alle Rechenwerkzeuge.

Zweck
-----
Die Werkzeuge rechnen Aufgaben aus Zahlentheorie, Kryptographie,
Codierungstheorie und Kombinatorik. Sie zeigen jeden Rechenschritt.
Diese Seite ist die Bedienoberfläche: Zahlen eingeben, rechnen,
den Lösungsweg lesen. Abschnitte lassen sich auf- und zuklappen.

Die Seite läuft auf dem eigenen Rechner. Sie braucht kein Internet.

Aufruf
------
    python -m skripte.web                   # öffnet http://localhost:8000
    python -m skripte.web --port 8080       # anderer Port
    python -m skripte.web --host 0.0.0.0    # im lokalen Netz freigeben (für das Team)
    python -m skripte.web --kein-browser

Wie es funktioniert
-------------------
- Der Server liest die Optionen jedes Werkzeugs aus dessen argparse-Parser.
  Daraus baut die Seite die Formulare. Neue Optionen erscheinen also von selbst.
- Beim Klick auf „Rechnen“ startet der Server das Werkzeug als eigenen Prozess.
  Es sind nur die Werkzeuge aus der Liste erlaubt. Es gibt keine Shell.
"""

from __future__ import annotations

import argparse
import importlib
import json
import shlex
import subprocess
import sys
import threading
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

WURZEL = Path(__file__).resolve().parent.parent
STATISCH = Path(__file__).resolve().parent / "web_static"

# Gruppen und Reihenfolge in der Seitenleiste.
GRUPPEN: list[tuple[str, list[tuple[str, str]]]] = [
    ("Zahlentheorie", [
        ("euklid", "Euklidischer Algorithmus"),
        ("primfaktor_phi", "Primfaktoren und φ(n)"),
        ("schnell_potenzieren", "Modulares Potenzieren"),
        ("crt", "Chinesischer Restsatz"),
        ("einheitengruppe", "Einheitengruppe ℤₘ*"),
        ("teilbarkeit", "Teilbarkeit"),
    ]),
    ("Kryptographie und Prüfziffern", [
        ("rsa", "RSA"),
        ("isbn10", "ISBN-10"),
    ]),
    ("Codierungstheorie", [
        ("linearer_code", "Lineare Codes"),
        ("mod11_code", "Modulo-11-Code"),
        ("zyklischer_code", "Zyklische Codes"),
        ("schranken", "Schranken für Codes"),
    ]),
    ("Kombinatorik", [
        ("kombinatorik", "Zählen und Anordnen"),
        ("inklusion_exklusion", "Inklusion–Exklusion"),
        ("turnierplan", "Turnierplan"),
    ]),
]
SKRIPTE: list[tuple[str, str]] = [s for _, liste in GRUPPEN for s in liste]
GRUPPE_VON = {name: gruppe for gruppe, liste in GRUPPEN for name, _ in liste}
ERLAUBT = {name for name, _ in SKRIPTE}


# ---------------------------------------------------------------------------
# Parser der Werkzeuge auslesen
# ---------------------------------------------------------------------------

class _ParserGefunden(Exception):
    """Trägt den Parser nach außen. Wird nur intern benutzt."""

    def __init__(self, parser: argparse.ArgumentParser):
        self.parser = parser


def parser_von(modulname: str) -> argparse.ArgumentParser:
    """Holt den argparse-Parser eines Skripts.

    Die Skripte bauen ihren Parser in main(). Deshalb ersetzen wir parse_args
    kurz durch eine Funktion, die den Parser zurückgibt statt zu rechnen.
    """
    modul = importlib.import_module(f"skripte.{modulname}")
    original = argparse.ArgumentParser.parse_args

    def abfangen(self, *args, **kwargs):
        raise _ParserGefunden(self)

    argparse.ArgumentParser.parse_args = abfangen
    try:
        modul.main([])
    except _ParserGefunden as gefunden:
        return gefunden.parser
    finally:
        argparse.ArgumentParser.parse_args = original
    raise RuntimeError(f"Kein Parser in skripte.{modulname} gefunden.")


def _metavar(aktion: argparse.Action) -> list[str] | str | None:
    mv = aktion.metavar
    if isinstance(mv, tuple):
        return list(mv)
    return mv


def schema(parser: argparse.ArgumentParser) -> dict:
    """Beschreibt alle Optionen eines Parsers als JSON-fähiges Wörterbuch."""
    felder = []
    unterbefehle = None
    for a in parser._actions:
        if isinstance(a, argparse._HelpAction):
            continue
        if isinstance(a, argparse._SubParsersAction):
            unterbefehle = {
                "dest": a.dest,
                "befehle": [
                    {"name": name, "hilfe": sp.description or "", **schema(sp)}
                    for name, sp in a.choices.items()
                ],
            }
            continue
        if isinstance(a, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            art = "schalter"
        elif isinstance(a, argparse._AppendAction):
            art = "mehrfach"
        elif a.choices is not None:
            art = "auswahl"
        else:
            art = "wert"
        nargs = a.nargs
        felder.append({
            "art": art,
            "optionen": list(a.option_strings),
            "dest": a.dest,
            "hilfe": (a.help or "").replace("%(default)s", str(a.default)),
            "standard": a.default if isinstance(a.default, (str, int, float, bool, type(None))) else None,
            "auswahl": list(a.choices) if a.choices is not None else None,
            "nargs": nargs if isinstance(nargs, (int, str, type(None))) else None,
            "pflicht": bool(a.required) or not a.option_strings,
            "metavar": _metavar(a),
        })
    return {
        "beschreibung": parser.description or "",
        "epilog": parser.epilog or "",
        "felder": felder,
        "unterbefehle": unterbefehle,
    }


def alle_schemata() -> list[dict]:
    ergebnis = []
    for name, titel in SKRIPTE:
        try:
            ergebnis.append({"name": name, "titel": titel, "gruppe": GRUPPE_VON[name], **schema(parser_von(name))})
        except Exception as fehler:  # ein kaputtes Skript soll die Seite nicht blockieren
            ergebnis.append({"name": name, "titel": titel, "gruppe": GRUPPE_VON[name], "fehler": str(fehler),
                             "beschreibung": "", "epilog": "", "felder": [], "unterbefehle": None})
    return ergebnis


# ---------------------------------------------------------------------------
# Werkzeug ausführen
# ---------------------------------------------------------------------------

def ausfuehren(name: str, argv: list[str]) -> dict:
    """Startet ein Werkzeug als eigenen Prozess und gibt die Ausgabe zurück."""
    if name not in ERLAUBT:
        return {"code": 2, "stdout": "", "stderr": f"Unbekanntes Werkzeug: {name}", "befehl": ""}
    befehl = [sys.executable, "-m", f"skripte.{name}", *argv]
    try:
        lauf = subprocess.run(befehl, cwd=WURZEL, capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired:
        return {"code": -1, "stdout": "", "stderr": "Abbruch: Die Rechnung lief länger als 60 Sekunden.",
                "befehl": ""}
    anzeige = "python -m skripte." + name + ("" if not argv else " " + shlex.join(argv))
    return {"code": lauf.returncode, "stdout": lauf.stdout, "stderr": lauf.stderr, "befehl": anzeige}


# ---------------------------------------------------------------------------
# HTTP-Server
# ---------------------------------------------------------------------------

class Handler(BaseHTTPRequestHandler):
    schemata_json: bytes = b"[]"

    def log_message(self, *args):  # keine Zeile pro Anfrage im Terminal
        pass

    def _senden(self, status: int, inhalt: bytes, typ: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", typ)
        self.send_header("Content-Length", str(len(inhalt)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(inhalt)

    def _json(self, daten, status: int = 200) -> None:
        self._senden(status, json.dumps(daten, ensure_ascii=False).encode(), "application/json; charset=utf-8")

    def do_GET(self):
        pfad = self.path.split("?", 1)[0]
        if pfad in ("/", "/index.html"):
            self._senden(200, (STATISCH / "index.html").read_bytes(), "text/html; charset=utf-8")
        elif pfad == "/api/skripte":
            self._senden(200, self.schemata_json, "application/json; charset=utf-8")
        else:
            self._senden(HTTPStatus.NOT_FOUND, b"Nicht gefunden", "text/plain; charset=utf-8")

    def do_POST(self):
        if self.path != "/api/ausfuehren":
            self._senden(HTTPStatus.NOT_FOUND, b"Nicht gefunden", "text/plain; charset=utf-8")
            return
        laenge = int(self.headers.get("Content-Length", "0"))
        try:
            daten = json.loads(self.rfile.read(laenge) or b"{}")
            name = str(daten["skript"])
            argv = [str(x) for x in daten.get("argv", [])]
        except (ValueError, KeyError) as fehler:
            self._json({"code": 2, "stdout": "", "stderr": f"Ungültige Anfrage: {fehler}", "befehl": ""}, 400)
            return
        self._json(ausfuehren(name, argv))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m skripte.web",
                                 description="Startet die Web-Oberfläche für alle Rechenwerkzeuge.")
    ap.add_argument("--port", type=int, default=8000, help="Port (Standard: 8000)")
    ap.add_argument("--host", default="127.0.0.1",
                    help="Adresse (Standard: 127.0.0.1, nur dieser Rechner). "
                         "0.0.0.0 gibt die Seite im lokalen Netz frei.")
    ap.add_argument("--kein-browser", action="store_true", help="Browser nicht automatisch öffnen")
    args = ap.parse_args(argv)

    Handler.schemata_json = json.dumps(alle_schemata(), ensure_ascii=False).encode()
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    adresse = f"http://localhost:{args.port}"
    print(f"Rechenwerkzeuge laufen: {adresse}")
    if args.host not in ("127.0.0.1", "localhost"):
        print(f"Im Netz erreichbar unter http://<Name-dieses-Rechners>:{args.port}")
    print("Beenden mit Strg+C.")
    if not args.kein_browser:
        threading.Timer(0.5, webbrowser.open, args=(adresse,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nBeendet.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
