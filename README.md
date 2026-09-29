# Rechenwerkzeuge

Rechenwerkzeuge für Zahlentheorie, Kryptographie, Codierungstheorie und Kombinatorik.
Jedes Werkzeug zeigt den ganzen Lösungsweg Schritt für Schritt.
Die Bedienung läuft über eine Seite im Browser.

## Starten

Du brauchst nur Python 3.12 (3.11 geht auch).

**Linux / macOS**

```bash
./start.sh
```

**Windows**

Doppelklick auf `start.bat`.

Beim ersten Start richtet das Skript eine eigene Python-Umgebung ein (Ordner `.venv`).
Danach öffnet sich der Browser unter http://localhost:8000.
Die Seite braucht kein Internet. Beenden mit Strg + C im Terminalfenster.

Weitere Optionen:

```bash
./start.sh --port 8080           # anderer Port
./start.sh --host 0.0.0.0        # im lokalen Netz für das Team freigeben
./start.sh --kein-browser        # Browser nicht automatisch öffnen
```

Bei `--host 0.0.0.0` erreichen andere die Seite unter `http://<Name-des-Rechners>:8000`.
Gib die Seite nur in vertrauenswürdigen Netzen frei.

## Bedienung

1. Links ein Werkzeug wählen. Die Suche hilft beim Finden.
2. Zahlen in die Felder eingeben. Pflichtfelder haben einen roten Stern.
3. „Rechnen“ klicken oder Strg + Enter drücken.
4. Der Lösungsweg erscheint darunter.

Weitere Funktionen:

- **Beispiele:** Ein Klick füllt die Felder mit einem Beispiel und rechnet sofort.
- **Varianten und Optionen:** Einige Aufgaben haben je nach Konvention verschiedene Lösungen,
  zum Beispiel die Form der Kontrollmatrix. Dieser Bereich ist aufklappbar.
  „Alle Varianten“ zeigt alle Versionen nacheinander.
- **Klappbare Abschnitte:** Jeder Abschnitt des Lösungswegs lässt sich auf- und zuklappen.
  Zugeklappt zeigt die Kopfzeile das Ergebnis des Abschnitts.
- **Ergebnisse auf einen Blick:** oben im Lösungsweg.
- **Farben:** grün = Ergebnis, gelb = Warnung (z. B. „57 ist keine Primzahl“), grau = Annahme.
- **A− / A+** ändert die Schriftgröße, **◐** wechselt zwischen hell und dunkel.
- **Text kopieren** und **Drucken** für den ganzen Lösungsweg.

Die Seite merkt sich die letzten Eingaben in diesem Browser.

## Werkzeuge

| Bereich | Werkzeug | Was es rechnet |
|---|---|---|
| Zahlentheorie | Euklidischer Algorithmus | ggT, kgV, Bezout, ax + by = c, Inverse mod m |
| | Primfaktoren und φ(n) | Primfaktorzerlegung, Eulersche φ-Funktion |
| | Modulares Potenzieren | aᵏ mod m, Ordnung, primitive Elemente |
| | Chinesischer Restsatz | simultane Kongruenzen |
| | Einheitengruppe ℤₘ* | Einheiten, Inverse, Gruppentafel, Ordnungen |
| | Teilbarkeit | Regeln für 9, 11, 7; Teileranzahl, Teilersumme, vollkommene Zahlen |
| Kryptographie und Prüfziffern | RSA | Schlüssel, Verschlüsseln, Entschlüsseln, Probe |
| | ISBN-10 | Prüfziffer, Prüfung, Zahlendreher |
| Codierungstheorie | Lineare Codes | Erzeuger- und Kontrollmatrix, Hamming-Codes, Syndromdecodierung |
| | Modulo-11-Code | Codieren, Syndrom, Fehlerkorrektur |
| | Zyklische Codes | Kreisteilungsklassen, Generator- und Kontrollpolynom |
| | Schranken für Codes | Hamming- und Singleton-Schranke |
| Kombinatorik | Zählen und Anordnen | Binomialkoeffizient, Gitterwege, Derangements, Partitionen |
| | Inklusion–Exklusion | Siebformel |
| | Turnierplan | Rundenturnier |

## Für Entwickler

- Die Werkzeuge liegen in `skripte/`. Jedes Werkzeug ist ein Python-Modul mit argparse.
- Die Seite (`skripte/web.py`, `skripte/web_static/index.html`) liest die Optionen aus dem Parser.
  Eine neue Option erscheint deshalb ohne Änderung an der Seite.
- Ein neues Werkzeug trägst du in `GRUPPEN` in `skripte/web.py` ein.
- Die Werkzeuge laufen auch ohne Seite, mit Tab-Vervollständigung:
  `python3 skripte/<name>.py --help` (oder `python -m skripte.<name> --help`).
  Vorher die Umgebung aktivieren: `source .venv/bin/activate` (Windows: `.venv\Scripts\activate`).
- Tests: `python -m pytest -q` (in der Umgebung `.venv`).
- Abhängigkeiten: `requirements.txt` (sympy, pytest).

Die Ordner `latex/` und die Datei `PRIORISIERUNG.md` sind Lernmaterial.
Die Werkzeuge brauchen sie nicht.
