"""Schneidet Teile aus den übertragenen LaTeX-Dateien aus, damit das Lernskript sie 1:1 einbindet.

Aufruf (im Ordner latex/):  python3 skript/auszug.py
Ergebnis: Dateien in latex/skript/auszug/
"""
import re
from pathlib import Path

WURZEL = Path(__file__).resolve().parent.parent
ZIEL = WURZEL / "skript" / "auszug"
ZIEL.mkdir(exist_ok=True)


def lies(pfad):
    return (WURZEL / pfad).read_text(encoding="utf-8")


def schreibe(name, text):
    (ZIEL / f"{name}.tex").write_text(text.strip() + "\n", encoding="utf-8")


def wissen(name, pfad, titel_anfang):
    """Unterabschnitt, dessen Titel mit titel_anfang beginnt, bis zum nächsten Unterabschnitt."""
    s = lies(pfad)
    m = re.search(r"^\\subsection\{" + re.escape(titel_anfang) + r".*$", s, re.M)
    assert m, (pfad, titel_anfang)
    ende = re.search(r"^\\(sub)?section\{", s[m.end():], re.M)
    rumpf = s[m.end(): m.end() + ende.start()] if ende else s[m.end():]
    titel = m.group(0)[len("\\subsection{"):-1]
    schreibe(name, f"\\subsubsection*{{{titel}}}\n{rumpf}")


def aufgabe(name, pfad, titel, neuer_titel):
    """Aufgabe mit [titel] und die direkt folgende Lösung."""
    s = lies(pfad)
    a = s.index(f"\\begin{{aufgabe}}[{titel}]")
    a_ende = s.index("\\end{aufgabe}", a) + len("\\end{aufgabe}")
    text = s[a:a_ende].replace(f"[{titel}]", f"[{neuer_titel}]", 1)
    schreibe(name + "_a", text)
    rest = s[a_ende:]
    l = rest.find("\\begin{loesung}")
    naechste = rest.find("\\begin{aufgabe}")
    if l != -1 and (naechste == -1 or l < naechste):
        l_ende = rest.index("\\end{loesung}", l) + len("\\end{loesung}")
        schreibe(name + "_l", f"\\textbf{{Lösung zu {neuer_titel}}}\n\n" + rest[l:l_ende])


def klausurloesung(name, nr, neuer_titel):
    s = lies("klausur1_loesungen.tex")
    m = re.search(r"^\\section\*\{Aufgabe " + str(nr) + r":.*$", s, re.M)
    ende = re.search(r"^(\\section\*\{|\\end\{document\})", s[m.end():], re.M)
    schreibe(name + "_l", f"\\textbf{{Lösung zu {neuer_titel}}}\n" + s[m.end(): m.end() + ende.start()])


def rezept(name, titel_anfang):
    s = lies("skript/klausurskript.tex")
    m = re.search(r"^\\section\{" + re.escape(titel_anfang) + r".*$", s, re.M)
    rest = s[m.end():]
    ende = min(x for x in (rest.find("\\subsection*{Übungen}"), rest.find("\n\\section{"),
                           rest.find("% ====")) if x != -1)
    schreibe(name, rest[:ende])


# Wissen (1:1 aus der Vorlesung)
wissen("w_vl3_31", "wissen/03_potenzen_zm.tex", "3.1")
wissen("w_vl3_32", "wissen/03_potenzen_zm.tex", "3.2")
wissen("w_vl3_33", "wissen/03_potenzen_zm.tex", "3.3")
wissen("w_vl4_41", "wissen/04_lineare_codes.tex", "4.1")
wissen("w_vl4_42", "wissen/04_lineare_codes.tex", "4.2")
wissen("w_vl6_63", "wissen/06_decodierung.tex", "6.3")
wissen("w_vl1_ea", "wissen/01_ring_ganze_zahlen.tex", "Division mit Rest")
wissen("w_vl1_crt", "wissen/01_ring_ganze_zahlen.tex", "Der Chinesische")
wissen("w_vl1_phi", "wissen/01_ring_ganze_zahlen.tex", "Eulersche")
wissen("w_vl1_isbn", "wissen/01_ring_ganze_zahlen.tex", "Der ISBN")
wissen("w_vl2_25", "wissen/02_zaehlprinzipien.tex", "2.5")

# Aufgaben mit Lösungen des Dozenten
KV = "uebungen/91_klausurvorbereitung.tex"
for nr in range(1, 9):
    aufgabe(f"kv{nr}", KV, f"Aufgabe {nr}", f"Klausurvorbereitung {nr}")
for nr in (1, 2, 3, 4, 7, 8, 9):
    aufgabe(f"u4_{nr}", "uebungen/04_uebung.tex", f"Aufgabe {nr}", f"Übung 4, Aufgabe {nr}")
for nr in (1, 2, 4, 5):
    aufgabe(f"u2_{nr}", "uebungen/02_uebung.tex", f"Aufgabe {nr}", f"Übung 2, Aufgabe {nr}")
for nr in (1, 2):
    aufgabe(f"u3_{nr}", "uebungen/03_uebung.tex", f"Aufgabe {nr}", f"Übung 3, Aufgabe {nr}")
aufgabe("u6_9", "uebungen/06_uebung.tex", "6\\_9", "Übung 6, Aufgabe 9")

# Klausur vom 29.07.2026 und Gedächtnisprotokoll
for nr in range(1, 6):
    aufgabe(f"k{nr}", "uebungen/92_klausur_pzr1.tex", f"Aufgabe {nr}", f"Klausur 29.07.2026, Aufgabe {nr}")
    klausurloesung(f"k{nr}", nr, f"Klausur 29.07.2026, Aufgabe {nr}")
for nr in (1, 2, 4):
    aufgabe(f"gp{nr}", "uebungen/90_gedaechtnisprotokoll.tex", f"Aufgabe {nr}", f"Gedächtnisprotokoll, Aufgabe {nr}")

# Rezepte und Musterrechnungen
rezept("r_rsa", "RSA")
rezept("r_sieb", "Siebformel")
rezept("r_euklid", "Euklid")
rezept("r_crt", "Chinesischer")
rezept("r_mod11", "Modulo-11")
rezept("r_formeln", "Formelsammlung")
print("fertig:", len(list(ZIEL.glob("*.tex"))), "Dateien")
