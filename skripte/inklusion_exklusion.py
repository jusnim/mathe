"""Inklusion und Exklusion (Siebformel) mit Lösungsweg.

Zweck
-----
Das Skript rechnet Aufgaben zur Siebformel aus der Vorlesung (VL2, Abschnitt 2.5).
Es gibt jeden Schritt so aus, wie er in den Musterlösungen steht.
Die Siebformel lautet:

    |Ω \\ A₁ ∪ … ∪ A_r| = |Ω| − α₁ + α₂ − … + (−1)^r·α_r

Dabei ist α_k die Summe der Anzahlen aller Schnitte von k Mengen.

Es gibt drei Modi:

1. Teiler: Wie viele Zahlen 1 … N sind durch keine der Zahlen t₁, t₂, … teilbar?
       python -m skripte.inklusion_exklusion --n 720 --teiler 3 4 5
2. Mengen: Man kennt |Ω|, |Aᵢ|, |Aᵢ ∩ Aⱼ|, … direkt (zum Beispiel Sportverein).
       python -m skripte.inklusion_exklusion --omega 55 --stufe 35 27 12 --stufe 13 7 5 --stufe 2
   Jede Option --stufe enthält die Anzahlen für eine Stufe k.
   Die Schnitte stehen in der Reihenfolge A₁∩A₂, A₁∩A₃, A₂∩A₃ (usw.).
3. Derangement: Wie viele fixpunktfreie Permutationen gibt es in S_n?
       python -m skripte.inklusion_exklusion --derangement 4

Varianten (nur Modus Teiler)
----------------------------
--variante kgv      (Standard) |A_a ∩ A_b| = ⌊N / kgV(a, b)⌋. Immer richtig.
--variante produkt  |A_a ∩ A_b| = ⌊N / (a·b)⌋. Nur richtig bei teilerfremden Teilern.
--alle-varianten    gibt beide Varianten nacheinander aus.

Schreibweise der Indizes: --indizes teiler (A₂, A₃, A₅ wie Ü4 A1)
oder --indizes nummer (A₁, A₂, A₃ wie VL2 und KV A6).
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import combinations, permutations

VARIANTEN = ("kgv", "produkt")
STANDARD_VARIANTE = "kgv"

_TIEF = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def tief(zahl: int | str) -> str:
    """Schreibt eine Zahl als tiefgestellten Index, zum Beispiel 12 → ₁₂."""
    return str(zahl).translate(_TIEF)


# ---------------------------------------------------------------------------
# Kleine Zahlentheorie (lokal, damit das Skript allein läuft)
# ---------------------------------------------------------------------------

def ggt(a: int, b: int) -> int:
    """Größter gemeinsamer Teiler mit dem Euklidischen Algorithmus."""
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def kgv(a: int, b: int) -> int:
    """Kleinstes gemeinsames Vielfaches: kgV(a, b) = a·b / ggT(a, b)."""
    return a * b // ggt(a, b)


def kgv_liste(zahlen) -> int:
    ergebnis = 1
    for z in zahlen:
        ergebnis = kgv(ergebnis, z)
    return ergebnis


def faktorisiere(n: int) -> dict[int, int]:
    """Primfaktorzerlegung durch Probedivision. Gibt {Primzahl: Exponent} zurück."""
    faktoren: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            faktoren[p] = faktoren.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        faktoren[n] = faktoren.get(n, 0) + 1
    return faktoren


def zerlegung_text(n: int) -> str:
    hoch = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    teile = []
    for p, e in faktorisiere(n).items():
        teile.append(str(p) if e == 1 else f"{p}{str(e).translate(hoch)}")
    return "·".join(teile) if teile else "1"


# ---------------------------------------------------------------------------
# Modus Teiler
# ---------------------------------------------------------------------------

@dataclass
class Schnitt:
    """Ein Schnitt A_a ∩ A_b ∩ … im Modus Teiler."""

    teiler: tuple[int, ...]
    modul: int          # kgV (oder Produkt bei Variante "produkt")
    anzahl: int         # ⌊N / modul⌋
    kgv: int            # immer das richtige kgV
    richtig: int        # ⌊N / kgV⌋


@dataclass
class TeilerErgebnis:
    n: int
    teiler: list[int]
    variante: str
    stufen: list[list[Schnitt]]
    alphas: list[int]
    ergebnis: int
    richtiges_ergebnis: int
    produktformel: Fraction
    produktformel_exakt: bool
    paarweise_teilerfremd: bool
    probe: int | None
    hinweise: list[str] = field(default_factory=list)


def pruefe_teiler(n: int, teiler: list[int]) -> tuple[list[int], list[str]]:
    """Prüft die Eingaben. Gibt bereinigte Teiler und Hinweise zurück."""
    hinweise: list[str] = []
    if n < 1:
        raise ValueError(f"N = {n} ist nicht positiv. Ω = {{1, …, N}} braucht N ≥ 1.")
    if not teiler:
        raise ValueError("Es fehlen die Teiler. Beispiel: --teiler 2 3 5")
    for t in teiler:
        if t < 1:
            raise ValueError(f"Der Teiler {t} ist nicht positiv.")
    bereinigt: list[int] = []
    for t in teiler:
        if t in bereinigt:
            hinweise.append(f"Hinweis: Der Teiler {t} steht doppelt in der Eingabe. Ich zähle ihn nur einmal.")
        else:
            bereinigt.append(t)
    if 1 in bereinigt:
        hinweise.append("Hinweis: Jede Zahl ist durch 1 teilbar. Dann ist A₁ = Ω und das Ergebnis ist 0.")
    for a, b in combinations(bereinigt, 2):
        klein, gross = min(a, b), max(a, b)
        if klein > 1 and gross % klein == 0:
            hinweise.append(
                f"Hinweis: {gross} ist ein Vielfaches von {klein}. Also gilt A{tief(gross)} ⊆ A{tief(klein)}. "
                f"Nur die kgV-Variante rechnet hier richtig."
            )
        elif ggt(a, b) > 1:
            hinweise.append(
                f"Hinweis: {a} und {b} sind nicht teilerfremd (ggT = {ggt(a, b)}). "
                f"Der Schnitt braucht kgV({a}, {b}) = {kgv(a, b)}, nicht {a}·{b} = {a * b}."
            )
    for t in bereinigt:
        if t > 1 and sum(faktorisiere(t).values()) > 1:
            hinweise.append(f"Hinweis: {t} = {zerlegung_text(t)} ist keine Primzahl. "
                            f"Das ist erlaubt (K A2 hat den Teiler 4).")
    return bereinigt, hinweise


def berechne_teiler(n: int, teiler: list[int], variante: str = STANDARD_VARIANTE,
                    probe_grenze: int = 10_000_000) -> TeilerErgebnis:
    """Rechnet die Siebformel für Zahlen 1 … n, die durch keinen Teiler teilbar sind."""
    if variante not in VARIANTEN:
        raise ValueError(f"Unbekannte Variante {variante!r}. Erlaubt: {', '.join(VARIANTEN)}")
    teiler, hinweise = pruefe_teiler(n, list(teiler))
    r = len(teiler)
    stufen: list[list[Schnitt]] = []
    alphas: list[int] = []
    richtige_alphas: list[int] = []
    for k in range(1, r + 1):
        stufe = []
        for auswahl in combinations(teiler, k):
            richtiges_kgv = kgv_liste(auswahl)
            modul = richtiges_kgv if variante == "kgv" else math.prod(auswahl)
            stufe.append(Schnitt(auswahl, modul, n // modul, richtiges_kgv, n // richtiges_kgv))
        stufen.append(stufe)
        alphas.append(sum(s.anzahl for s in stufe))
        richtige_alphas.append(sum(s.richtig for s in stufe))
    ergebnis = n + sum((-1) ** k * a for k, a in enumerate(alphas, start=1))
    richtig = n + sum((-1) ** k * a for k, a in enumerate(richtige_alphas, start=1))

    produktformel = Fraction(n)
    for t in teiler:
        produktformel *= 1 - Fraction(1, t)
    teilerfremd = all(ggt(a, b) == 1 for a, b in combinations(teiler, 2))
    exakt = teilerfremd and n % math.prod(teiler) == 0

    probe = None
    if n <= probe_grenze:
        probe = sum(1 for k in range(1, n + 1) if all(k % t for t in teiler))
    return TeilerErgebnis(n, teiler, variante, stufen, alphas, ergebnis, richtig,
                          produktformel, exakt, teilerfremd, probe, hinweise)


def _bruch_text(b: Fraction) -> str:
    return str(b.numerator) if b.denominator == 1 else f"{b.numerator}/{b.denominator}"


def _dezimal(b: Fraction, stellen: int = 4) -> str:
    return f"{float(b):.{stellen}f}".replace(".", ",")


def _vorzeichen_summe(start: int, alphas: list[int]) -> str:
    text = str(start)
    for k, a in enumerate(alphas, start=1):
        text += f" {'−' if k % 2 else '+'} {a}"
    return text


def _alpha_formel(r: int) -> str:
    teile = ["|Ω|"]
    for k in range(1, r + 1):
        teile.append(f"{'−' if k % 2 else '+'} α{tief(k)}")
    return " ".join(teile)


def text_teiler(e: TeilerErgebnis, indizes: str = "teiler") -> str:
    """Erzeugt den Lösungsweg im Modus Teiler."""
    nummer = {t: i for i, t in enumerate(e.teiler, start=1)}

    def a(t: int) -> str:
        return f"A{tief(t if indizes == 'teiler' else nummer[t])}"

    def schnitt_name(auswahl) -> str:
        return " ∩ ".join(a(t) for t in auswahl)

    r = len(e.teiler)
    vereinigung = " ∪ ".join(a(t) for t in e.teiler)
    z: list[str] = []
    if e.variante == "kgv":
        z.append("Annahme: Variante 'kgv' (Standard). |A_a ∩ A_b| = ⌊N / kgV(a, b)⌋.")
    else:
        z.append("Annahme: Variante 'produkt'. |A_a ∩ A_b| = ⌊N / (a·b)⌋.")
        z.append("Diese Variante stimmt nur, wenn die Teiler paarweise teilerfremd sind.")
    z.append("Schreibweise: ⌊x⌋ heißt abrunden auf die nächste ganze Zahl.")
    z.append("")
    for h in e.hinweise:
        z.append(h)
    if e.hinweise:
        z.append("")
    z.append(f"Ω = {{1 ≤ k ≤ {e.n}}},   |Ω| = {e.n}")
    if indizes == "teiler":
        z.append(f"Sei A_i = {{k ∈ Ω | i ist Teiler von k}} für i = {', '.join(map(str, e.teiler))}.")
    else:
        for t in e.teiler:
            z.append(f"Sei {a(t)} = {{k ∈ Ω | {t} ist Teiler von k}}.")
    z.append("(Achtung: Ü4 A1 schreibt „i ist Teiler von 1000“. Das ist ein Fehler der Quelle.")
    z.append(" Richtig ist „i teilt k“, also die Zahlen in Ω, die durch i teilbar sind.)")
    z.append("")
    z.append("Die gesuchte Anzahl ist dann")
    z.append(f"|Ω \\ {vereinigung}| = {_alpha_formel(r)}")
    z.append("mit")
    for k, stufe in enumerate(e.stufen, start=1):
        z.append(f"α{tief(k)} = " + " + ".join(f"|{schnitt_name(s.teiler)}|" for s in stufe))
    z.append("")
    for k, stufe in enumerate(e.stufen, start=1):
        teile = []
        for s in stufe:
            if e.n % s.modul == 0:
                teile.append(f"|{schnitt_name(s.teiler)}| = {e.n}/{s.modul} = {s.anzahl}")
            else:
                teile.append(f"|{schnitt_name(s.teiler)}| = ⌊{e.n}/{s.modul}⌋ = {s.anzahl}")
        for t in teile:
            z.append("    " + t)
        if k >= 2:
            erklaerung = []
            for s in stufe:
                args = ", ".join(map(str, s.teiler))
                if e.variante == "kgv":
                    prod = math.prod(s.teiler)
                    if s.kgv == prod:
                        erklaerung.append(f"kgV({args}) = {s.kgv}")
                    else:
                        erklaerung.append(f"kgV({args}) = {s.kgv} (nicht {prod})")
                else:
                    erklaerung.append(f"{'·'.join(map(str, s.teiler))} = {s.modul}")
            z.append("    (" + "; ".join(erklaerung) + ")")
        if len(stufe) == 1:
            z.append(f"⇒ α{tief(k)} = {e.alphas[k - 1]}")
        else:
            z.append(f"⇒ α{tief(k)} = " + " + ".join(str(s.anzahl) for s in stufe) + f" = {e.alphas[k - 1]}")
    z.append("")
    z.append(f"⇒ |Ω \\ {vereinigung}| = {_alpha_formel(r)}")
    z.append(f"    = {_vorzeichen_summe(e.n, e.alphas)} = {e.ergebnis}")
    z.append("")

    if e.variante == "produkt" and e.ergebnis != e.richtiges_ergebnis:
        z.append(f"ACHTUNG: Diese Variante ist hier FALSCH. Mit kgV kommt {e.richtiges_ergebnis} heraus.")
        for stufe in e.stufen:
            for s in stufe:
                if s.modul != s.kgv:
                    z.append(f"    {schnitt_name(s.teiler)}: Zahlen durch {s.kgv} teilbar, nicht nur durch {s.modul}. "
                             f"Richtig: ⌊{e.n}/{s.kgv}⌋ = {s.richtig} statt {s.anzahl}.")
        z.append("")
    elif e.variante == "produkt":
        z.append("Hier sind die Teiler paarweise teilerfremd. Also ist a·b = kgV(a, b) und das Ergebnis stimmt.")
        z.append("")

    # Kontrolle mit Produktformel
    faktoren = "".join(f"(1 − 1/{t})" for t in e.teiler)
    z.append("Kontrolle mit der Produktformel:")
    z.append(f"    {e.n} · {faktoren} = {_bruch_text(e.produktformel)}"
             + ("" if e.produktformel.denominator == 1 else f" ≈ {_dezimal(e.produktformel)}"))
    if e.produktformel_exakt:
        z.append(f"    Die Formel ist hier exakt: Die Teiler sind paarweise teilerfremd und {e.n} ist durch "
                 f"{'·'.join(map(str, e.teiler))} = {math.prod(e.teiler)} teilbar.")
    elif not e.paarweise_teilerfremd:
        z.append("    Die Formel passt hier nicht: Die Teiler sind nicht paarweise teilerfremd.")
    else:
        z.append(f"    Die Formel gilt hier nur näherungsweise: {e.n} ist nicht durch "
                 f"{'·'.join(map(str, e.teiler))} = {math.prod(e.teiler)} teilbar.")
        z.append("    (So steht es auch in Ü4 A1: 266,66 statt 266.)")
    z.append("")
    if e.probe is not None:
        status = "stimmt" if e.probe == e.richtiges_ergebnis else "WEICHT AB"
        z.append(f"Probe durch Abzählen aller k = 1 … {e.n}: {e.probe} Zahlen sind durch keinen Teiler teilbar ({status}).")
    else:
        z.append("Probe durch Abzählen: übersprungen, weil N sehr groß ist.")
    z.append("")
    z.append(f"Ergebnis: {e.ergebnis} Zahlen zwischen 1 und {e.n} sind durch keine der Zahlen "
             f"{', '.join(map(str, e.teiler))} teilbar." + ("" if e.ergebnis == e.richtiges_ergebnis
                                                          else f" (Variante 'produkt' – richtig ist {e.richtiges_ergebnis}.)"))
    return "\n".join(z)


def liste_teiler(n: int, teiler: list[int]) -> list[int]:
    """Alle k in 1 … n, die durch keinen der Teiler teilbar sind."""
    return [k for k in range(1, n + 1) if all(k % t for t in teiler)]


# ---------------------------------------------------------------------------
# Modus Mengen
# ---------------------------------------------------------------------------

@dataclass
class MengenErgebnis:
    omega: int
    r: int
    stufen: list[list[int]]
    schnitte: list[list[tuple[int, ...]]]
    alphas: list[int]
    ergebnis: int
    namen: list[str]
    hinweise: list[str] = field(default_factory=list)


def berechne_mengen(omega: int, stufen: list[list[int]], namen: list[str] | None = None) -> MengenErgebnis:
    """Siebformel, wenn |Ω|, |Aᵢ|, |Aᵢ ∩ Aⱼ|, … direkt gegeben sind."""
    if omega < 0:
        raise ValueError(f"|Ω| = {omega} ist negativ.")
    if not stufen or not stufen[0]:
        raise ValueError("Es fehlen die Anzahlen |Aᵢ|. Beispiel: --stufe 35 27 12")
    r = len(stufen[0])
    hinweise: list[str] = []
    if len(stufen) > r:
        raise ValueError(f"Es gibt nur {r} Mengen, aber {len(stufen)} Stufen.")
    schnitte = [list(combinations(range(1, r + 1), k)) for k in range(1, r + 1)]
    stufen = [list(s) for s in stufen]
    for k, s in enumerate(stufen, start=1):
        soll = math.comb(r, k)
        if len(s) != soll:
            raise ValueError(
                f"Stufe {k} braucht C({r}, {k}) = {soll} Zahlen, es sind aber {len(s)}. "
                f"Reihenfolge: " + ", ".join("∩".join(f"A{tief(i)}" for i in c) for c in schnitte[k - 1])
            )
        for wert in s:
            if wert < 0:
                raise ValueError(f"In Stufe {k} steht die negative Zahl {wert}.")
    for k in range(len(stufen) + 1, r + 1):
        hinweise.append(f"Hinweis: Stufe {k} fehlt. Ich nehme an, alle Schnitte von {k} Mengen sind leer (α{tief(k)} = 0).")
        stufen.append([0] * math.comb(r, k))
    # Plausibilität: ein Schnitt ist nie größer als eine seiner Teilmengen.
    wert_von = {c: stufen[len(c) - 1][i] for k in range(r) for i, c in enumerate(schnitte[k])}
    for c, w in wert_von.items():
        if w > omega:
            hinweise.append(f"Warnung: |{'∩'.join(f'A{tief(i)}' for i in c)}| = {w} ist größer als |Ω| = {omega}.")
        for teil in combinations(c, len(c) - 1):
            if teil and w > wert_von[teil]:
                hinweise.append(
                    f"Warnung: |{'∩'.join(f'A{tief(i)}' for i in c)}| = {w} ist größer als "
                    f"|{'∩'.join(f'A{tief(i)}' for i in teil)}| = {wert_von[teil]}. Das ist unmöglich."
                )
    alphas = [sum(s) for s in stufen]
    ergebnis = omega + sum((-1) ** k * a for k, a in enumerate(alphas, start=1))
    if ergebnis < 0:
        hinweise.append(f"Warnung: Das Ergebnis {ergebnis} ist negativ. Die Eingaben passen nicht zusammen.")
    namen = list(namen or [])
    if namen and len(namen) != r:
        hinweise.append(f"Hinweis: {len(namen)} Namen für {r} Mengen. Ich ignoriere die Namen.")
        namen = []
    return MengenErgebnis(omega, r, stufen, schnitte, alphas, ergebnis, namen, hinweise)


def text_mengen(e: MengenErgebnis) -> str:
    """Erzeugt den Lösungsweg im Modus Mengen (wie VL2 Musikschule und KV A6)."""
    z: list[str] = []
    z.append("Annahme: Die Anzahlen der Mengen und Schnitte sind direkt gegeben.")
    z.append("")
    for h in e.hinweise:
        z.append(h)
    if e.hinweise:
        z.append("")
    z.append(f"|Ω| = {e.omega},   r = {e.r}")
    if e.namen:
        for i, name in enumerate(e.namen, start=1):
            z.append(f"A{tief(i)} = {{a ∈ Ω | a gehört zu „{name}“}}")
    else:
        z.append("A_i = {a ∈ Ω | a besitzt Eigenschaft i}")
    vereinigung = " ∪ ".join(f"A{tief(i)}" for i in range(1, e.r + 1))

    def name(c) -> str:
        return " ∩ ".join(f"A{tief(i)}" for i in c)

    z.append(f"Gesucht: |Ω \\ {vereinigung}|")
    z.append("")
    z.append(f"|Ω \\ {vereinigung}| = {_alpha_formel(e.r)}")
    z.append("mit")
    for k in range(e.r):
        z.append(f"α{tief(k + 1)} = " + " + ".join(f"|{name(c)}|" for c in e.schnitte[k]))
    z.append("")
    for k in range(e.r):
        z.append("    " + ",  ".join(f"|{name(c)}| = {w}" for c, w in zip(e.schnitte[k], e.stufen[k])))
        if len(e.stufen[k]) == 1:
            z.append(f"⇒ α{tief(k + 1)} = {e.alphas[k]}")
        else:
            z.append(f"⇒ α{tief(k + 1)} = " + " + ".join(map(str, e.stufen[k])) + f" = {e.alphas[k]}")
    z.append("")
    # Zeile wie in VL2: 73 − (20 + 25 + 52) + (7 + 12 + 17) − 1
    teile = [str(e.omega)]
    for k in range(e.r):
        klammer = " + ".join(map(str, e.stufen[k]))
        if len(e.stufen[k]) > 1:
            klammer = f"({klammer})"
        teile.append(f"{'−' if (k + 1) % 2 else '+'} {klammer}")
    z.append(f"|Ω \\ {vereinigung}| = {_alpha_formel(e.r)}")
    z.append("    = " + " ".join(teile))
    z.append(f"    = {_vorzeichen_summe(e.omega, e.alphas)}")
    z.append(f"    = {e.ergebnis}")
    z.append("")
    # Probe: Anzahl der Elemente in genau j Mengen muss zusammen |Ω| ergeben.
    genau = genau_in_mengen(e)
    z.append("Probe: Anzahl der Elemente in genau j Mengen (aus den α berechnet):")
    z.append("    " + ",  ".join(f"genau {j}: {w}" for j, w in enumerate(genau)))
    summe = sum(genau)
    z.append(f"    Summe = {' + '.join(map(str, genau))} = {summe}"
             + (" = |Ω| (stimmt)" if summe == e.omega else f" ≠ |Ω| (FEHLER)"))
    if any(w < 0 for w in genau):
        z.append("    Warnung: Eine Anzahl ist negativ. Die Eingaben passen nicht zusammen.")
    z.append("")
    z.append(f"Ergebnis: {e.ergebnis} Elemente von Ω liegen in keiner der Mengen A₁ … A{tief(e.r)}.")
    return "\n".join(z)


def genau_in_mengen(e: MengenErgebnis) -> list[int]:
    """Anzahl der Elemente in genau j Mengen: E_j = Σ_k (−1)^(k−j)·C(k, j)·α_k (α₀ = |Ω|)."""
    alle = [e.omega] + e.alphas
    return [sum((-1) ** (k - j) * math.comb(k, j) * alle[k] for k in range(j, e.r + 1))
            for j in range(e.r + 1)]


# ---------------------------------------------------------------------------
# Modus Derangement
# ---------------------------------------------------------------------------

@dataclass
class DerangementErgebnis:
    n: int
    alphas: list[int]
    d_n: int
    wahrscheinlichkeit: Fraction
    reihe: list[Fraction]
    liste: list[tuple[int, ...]] | None


def berechne_derangement(n: int, liste_bis: int = 6) -> DerangementErgebnis:
    """Anzahl der fixpunktfreien Permutationen d_n (VL2, Derangement-Problem)."""
    if n < 1:
        raise ValueError(f"n = {n} ist zu klein. Man braucht n ≥ 1.")
    fak = math.factorial(n)
    alphas = [math.comb(n, k) * math.factorial(n - k) for k in range(1, n + 1)]
    d_n = fak + sum((-1) ** k * a for k, a in enumerate(alphas, start=1))
    reihe = [Fraction((-1) ** k, math.factorial(k)) for k in range(n + 1)]
    liste = None
    if n <= liste_bis:
        liste = [p for p in permutations(range(1, n + 1)) if all(p[i] != i + 1 for i in range(n))]
    return DerangementErgebnis(n, alphas, d_n, Fraction(d_n, fak), reihe, liste)


def text_derangement(e: DerangementErgebnis) -> str:
    n = e.n
    z: list[str] = []
    z.append("Annahme: Ω = S_n, die Menge aller Permutationen von n Objekten.")
    z.append("Ein Fixpunkt ist ein i mit π(i) = i. Ein Derangement ist eine Permutation ohne Fixpunkt.")
    z.append("")
    z.append(f"Ω = S{tief(n)},   |Ω| = {n}! = {math.factorial(n)},   r = {n}")
    z.append("A_i = {π ∈ Ω | π(i) = i}")
    z.append("|Schnitt von k Mengen A_i| = (n − k)!")
    z.append("α_k = C(n, k)·(n − k)! = n!/k!")
    z.append("")
    for k, a in enumerate(e.alphas, start=1):
        z.append(f"    α{tief(k)} = C({n}, {k})·{n - k}! = {math.comb(n, k)}·{math.factorial(n - k)} = {a}")
    z.append("")
    z.append(f"d{tief(n)} = {n}! − α₁ + α₂ − … + (−1)^{n}·α{tief(n)}")
    z.append(f"   = {_vorzeichen_summe(math.factorial(n), e.alphas)} = {e.d_n}")
    z.append("")
    reihe_text = " ".join(
        ("1" if k == 0 else f"{'−' if k % 2 else '+'} 1/{k}!") for k in range(n + 1)
    )
    summe = sum(e.reihe)
    z.append("Mit der Formel aus VL2:")
    z.append(f"d{tief(n)} = {n}!·({reihe_text}) = {n}!·{_bruch_text(summe)} = {e.d_n}")
    z.append("")
    z.append(f"P({n}) = d{tief(n)}/{n}! = {e.d_n}/{math.factorial(n)} = {_bruch_text(e.wahrscheinlichkeit)}"
             f" ≈ {_dezimal(e.wahrscheinlichkeit)}      (Grenzwert 1/e ≈ 0,3679)")
    if e.liste is not None:
        z.append("")
        z.append(f"Die {len(e.liste)} fixpunktfreien Permutationen (Bild von 1, 2, …, {n}):")
        for i in range(0, len(e.liste), 6):
            z.append("    " + ",   ".join(" ".join(map(str, p)) for p in e.liste[i:i + 6]))
        status = "stimmt" if len(e.liste) == e.d_n else "WEICHT AB"
        z.append(f"Probe durch Abzählen: {len(e.liste)} Permutationen ({status}).")
    z.append("")
    z.append(f"Ergebnis: d{tief(n)} = {e.d_n}")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

BEISPIELE = """Beispiele aus den Quellen:
  Klausur A2 (Zahlen ≤ 720, weder durch 3, 4 noch 5 teilbar):
    python -m skripte.inklusion_exklusion --n 720 --teiler 3 4 5
  Ü4 A1 (1 … 1000, weder durch 2, 3 noch 5 teilbar; Ergebnis 266):
    python -m skripte.inklusion_exklusion --n 1000 --teiler 2 3 5
  KV A6 (Sportverein, 55 Athleten; Ergebnis 4):
    python -m skripte.inklusion_exklusion --omega 55 --stufe 35 27 12 --stufe 13 7 5 --stufe 2 \\
        --namen Fußball Leichtathletik Judo
  Ü4 A2 (fixpunktfreie Permutationen in S₄; Ergebnis 9):
    python -m skripte.inklusion_exklusion --derangement 4
"""


def baue_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m skripte.inklusion_exklusion",
        description="Siebformel (Inklusion und Exklusion) mit vollständigem Lösungsweg.",
        epilog=BEISPIELE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    g = p.add_argument_group("Modus Teiler")
    g.add_argument("--n", type=int, help="Obere Grenze N. Gezählt werden die Zahlen 1 … N.")
    g.add_argument("--teiler", type=int, nargs="+", help="Die Teiler, zum Beispiel 3 4 5.")
    g.add_argument("--variante", choices=VARIANTEN, default=STANDARD_VARIANTE,
                   help="Wie man Schnitte zählt: kgv (Standard, immer richtig) oder produkt (nur bei teilerfremden Teilern richtig).")
    g.add_argument("--alle-varianten", action="store_true", help="Gibt alle Varianten nacheinander aus.")
    g.add_argument("--indizes", choices=("teiler", "nummer"), default="teiler",
                   help="Namen der Mengen: teiler (A₂, A₃, A₅ wie Ü4 A1, Standard) oder nummer (A₁, A₂, A₃ wie VL2).")
    g.add_argument("--liste", action="store_true", help="Gibt zusätzlich alle gezählten Zahlen aus.")
    g2 = p.add_argument_group("Modus Mengen")
    g2.add_argument("--omega", type=int, help="Größe der Grundmenge |Ω|, zum Beispiel 55.")
    g2.add_argument("--stufe", type=int, nargs="+", action="append",
                    help="Anzahlen einer Stufe. Erst |A₁| |A₂| |A₃|, dann |A₁∩A₂| |A₁∩A₃| |A₂∩A₃|, dann |A₁∩A₂∩A₃|. "
                         "Die Option mehrmals angeben.")
    g2.add_argument("--namen", nargs="+", help="Optionale Namen der Mengen, zum Beispiel Flöte Geige Klavier.")
    g3 = p.add_argument_group("Modus Derangement")
    g3.add_argument("--derangement", type=int, metavar="N",
                    help="Anzahl der fixpunktfreien Permutationen in S_N berechnen.")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = baue_parser()
    args = parser.parse_args(argv)
    modi = [m for m, aktiv in (("teiler", args.teiler is not None or args.n is not None),
                               ("mengen", args.stufe is not None or args.omega is not None),
                               ("derangement", args.derangement is not None)) if aktiv]
    if len(modi) != 1:
        parser.error("Bitte genau einen Modus wählen: --n mit --teiler, --omega mit --stufe, oder --derangement.")
    try:
        if modi[0] == "teiler":
            if args.n is None or args.teiler is None:
                parser.error("Im Modus Teiler braucht man --n und --teiler.")
            varianten = VARIANTEN if args.alle_varianten else (args.variante,)
            for i, v in enumerate(varianten):
                if len(varianten) > 1:
                    if i:
                        print()
                    print(f"===== Variante: {v} =====")
                erg = berechne_teiler(args.n, args.teiler, v)
                print(text_teiler(erg, args.indizes))
            if args.liste:
                zahlen = liste_teiler(args.n, erg.teiler)
                print()
                print(f"Die {len(zahlen)} Zahlen: " + ", ".join(map(str, zahlen)))
        elif modi[0] == "mengen":
            if args.omega is None or args.stufe is None:
                parser.error("Im Modus Mengen braucht man --omega und mindestens eine --stufe.")
            print(text_mengen(berechne_mengen(args.omega, args.stufe, args.namen)))
        else:
            print(text_derangement(berechne_derangement(args.derangement)))
    except ValueError as fehler:
        print(f"Fehler in der Eingabe: {fehler}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
