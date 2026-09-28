# mathe – Diskrete Mathematik (HTW Berlin)

Dieses Repository hilft bei der Klausur.
Es hat drei Teile.

## 1. Wissen und Übungen als LaTeX (`latex/`)

- `latex/wissen/` enthält Definitionen, Sätze, Beweise und Verfahren.
- `latex/uebungen/` enthält Aufgaben mit Lösungen.
- `latex/main.tex` bindet alle Dateien ein.
- Die Texte sind 1:1 aus den Originalen übertragen.
- Fehler im Original sind mit `% UNSICHER` markiert.

## 2. Priorisierung (`PRIORISIERUNG.md`)

Die Datei sortiert Themen und Aufgaben.
Oben steht, was am wahrscheinlichsten in der Klausur vorkommt.

## 3. Python-Skripte (`skripte/`)

Jedes Skript gibt den Lösungsweg Schritt für Schritt aus.
Man schreibt ihn in der Klausur auf Papier ab.

### Einrichten (einmal, vor der Klausur)

```bash
python3.12 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m pytest -q              # alle Tests müssen bestehen
```

### Benutzen im Browser (empfohlen)

```bash
python -m skripte.web
```

Der Browser öffnet sich unter http://localhost:8000.
Die Seite braucht kein Internet.

- Links steht die Liste der Skripte. Die Reihenfolge ist die Priorität für die Klausur.
- Für jedes Skript gibt es ein Formular. Ausfüllen und „Rechnen“ klicken (oder Strg + Enter).
- Unter „Beispiele aus den Quellen“ startet ein Klick eine Aufgabe aus Übung oder Klausur.
- Die Ausgabe ist gegliedert: Abschnitte, grüne Ergebnis-Kästen, gelbe Warnungen.
- Oben steht „Alle Ergebnisse auf einen Blick“.
- Mit A− und A+ ändert man die Schriftgröße. ◐ wechselt zwischen hell und dunkel.
- Der Reiter „Priorisierung“ zeigt `PRIORISIERUNG.md` als Seite.

Beenden mit Strg + C im Terminal.

### Benutzen im Terminal

Jedes Skript hat eine Hilfe mit Beispiel: `python -m skripte.<name> --help`.
Führen Konventionen zu verschiedenen Ergebnissen, dann gibt es `--variante`.
`--alle-varianten` zeigt alle Versionen nacheinander.

Die Reihenfolge entspricht der Priorität für die Klausur:

| Nr. | Skript | Thema | Beispiel |
|---:|---|---|---|
| 1 | `euklid` | Euklidischer Algorithmus, ggT, Bezout, ax+by=c, Inverse | `python -m skripte.euklid --a 2406 --b 654 --c 24` |
| 2 | `primfaktor_phi` | Primfaktoren, φ(n) | `python -m skripte.primfaktor_phi --n 2451` |
| 3 | `schnell_potenzieren` | aᵏ mod m, Ordnung, primitive Elemente | `python -m skripte.schnell_potenzieren --a 12 --k 100 --m 34` |
| 4 | `rsa` | RSA: verschlüsseln, d berechnen, Probe | `python -m skripte.rsa --p 43 --q 57 --e 221 --w 1511` |
| 5 | `mod11_code` | [10,8,3]₁₁-Code: codieren, Syndrom, decodieren | `python -m skripte.mod11_code --a 11500005 --r 1150000511` |
| 6 | `inklusion_exklusion` | Siebformel, Derangements | `python -m skripte.inklusion_exklusion --n 720 --teiler 3 4 5` |
| 7 | `crt` | Chinesischer Restsatz | `python -m skripte.crt -k 5 13 -k 4 15` |
| 8 | `isbn10` | ISBN-10-Prüfziffer, Zahlendreher | `python -m skripte.isbn10 --nummer 3-05-501517` |
| 9 | `linearer_code` | G, H, Hamming-Codes, Syndromdecodierung | `python -m skripte.linearer_code --q 3 --typ hamming --m 2 --r 2200` |
| 10 | `schranken` | Hamming- und Singleton-Schranke | `python -m skripte.schranken --n 10 --k 8 --d 3 --q 11` |
| 11 | `turnierplan` | Rundenturnier | `python -m skripte.turnierplan --mannschaften 8 --gegner` |
| 12 | `einheitengruppe` | ℤₘ*, Inverse, Gruppentafel | `python -m skripte.einheitengruppe --m 17 --modus alles` |
| 13 | `zyklischer_code` | Kreisteilungsklassen, g(X), h(X), Syndrom | `python -m skripte.zyklischer_code --p 2 --n 7 --g "X^3+X+1" --y "X^6+X+1"` |
| 14 | `kombinatorik` | Binomialkoeffizient, Gitterwege, Partitionen | `python -m skripte.kombinatorik gitter 4 4` |
| 15 | `teilbarkeit` | Regeln für 9, 11, 7; d(n), σ(n) | `python -m skripte.teilbarkeit regeln 299792` |

### Wichtig

- Das Skript prüft die Eingaben. Ein Beispiel: Ist p keine Primzahl, erscheint eine Warnung.
- Lies immer die Zeile „Annahme: …“. Sie sagt, welche Konvention gilt.
- Vergleiche die Konvention mit der Aufgabenstellung. Wähle sonst eine andere `--variante`.
