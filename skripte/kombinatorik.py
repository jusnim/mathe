"""Kombinatorik: Zählformeln mit Lösungsweg.

Zweck
-----
Das Skript löst Zählaufgaben: Auswahlen, Permutationen, Gitterwege,
Derangements, Partitionen und Summen mit Binomialkoeffizienten.
Es zeigt jeden Rechenschritt. Am Ende steht "Ergebnis: ...".
Danach folgt meist eine Probe.

Unterbefehle
------------
  binom        Binomialkoeffizient C(n, k) ("n über k")
  geordnet     geordnete Auswahl von k aus n (ohne oder mit Wiederholung)
  permutation  Anzahl der Permutationen n!
  kombination  ungeordnete Auswahl von k aus n (ohne oder mit Wiederholung)
  gitter       kürzeste Gitterwege im m×n-Gitter
  derangement  fixpunktfreie Permutationen dₙ
  partition    Partitionen p(n) einer Zahl, auch mit Einschränkung
  ferrers      Ferrers-Diagramm und konjugierte Partition
  lehrsatz     binomischer Lehrsatz (a + b)ⁿ
  summe        Summen mit Binomialkoeffizienten
  pascal       Pascalsches Dreieck

Aufruf (Beispiele)
------------------
  python -m skripte.kombinatorik gitter 4 4              # 4×4-Gitter: 70
  python -m skripte.kombinatorik gitter 3 5              # 3×5-Gitter
  python -m skripte.kombinatorik derangement 4 --liste   # d₄ = 9
  python -m skripte.kombinatorik partition 5             # p(5) = 7
  python -m skripte.kombinatorik summe alternierend 6    # Ergebnis 1
  python -m skripte.kombinatorik summe vandermonde 3 4 --r 2
  python -m skripte.kombinatorik kombination 5 7 --wiederholung --auswahl 1 1 1 2 4 4 5

Varianten
---------
Nur der Unterbefehl "gitter" hat Varianten. Sie betreffen die Frage:
Was bedeutet "m×n-Gitter"?
  kaestchen (Standard): m×n Kästchen, also (m+1)×(n+1) Gitterpunkte.
                        Ein 4×4-Gitter hat dann Wege aus 8 Strecken.
  punkte:               m×n Gitterpunkte, also (m−1)×(n−1) Kästchen.
Mit --alle-varianten gibt das Skript beide Rechnungen aus.
Alle anderen Formeln haben keine abweichenden Konventionen.
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
from dataclasses import dataclass, field
from fractions import Fraction

# ---------------------------------------------------------------------------
# Kleine Hilfen für die Ausgabe
# ---------------------------------------------------------------------------


def komma(x: Fraction | float, stellen: int = 4) -> str:
    """Schreibt eine Zahl mit deutschem Dezimalkomma."""
    text = f"{float(x):.{stellen}f}".rstrip("0").rstrip(".")
    return text.replace(".", ",")


def produkt_text(faktoren: list[int]) -> str:
    """Schreibt Faktoren als 8·7·6·5. Eine leere Liste ist 1."""
    return "·".join(str(f) for f in faktoren) if faktoren else "1"


def C(n: int, k: int) -> str:
    """Schreibweise für den Binomialkoeffizienten."""
    return f"C({n},{k})"


def binom(n: int, k: int) -> int:
    """Binomialkoeffizient. Für k < 0 oder k > n ist er 0."""
    if k < 0 or n < 0 or k > n:
        return 0
    return math.comb(n, k)


def pruefe_nichtnegativ(**werte: int) -> list[str]:
    """Gibt Hinweise zurück, wenn eine Zahl negativ ist."""
    return [f"Hinweis: {name} = {wert} ist negativ. Erlaubt sind nur Zahlen ≥ 0."
            for name, wert in werte.items() if wert < 0]


# ---------------------------------------------------------------------------
# Binomialkoeffizient
# ---------------------------------------------------------------------------


@dataclass
class BinomDaten:
    n: int
    k: int
    k_rechnung: int          # benutztes k (evtl. n − k wegen Symmetrie)
    zaehler: list[int]       # n·(n−1)·…·(n−k+1)
    nenner: list[int]        # 1·2·…·k
    gekuerzt: list[int]      # Zähler nach dem Kürzen
    wert: int
    hinweise: list[str] = field(default_factory=list)


def kuerze(zaehler: list[int], nenner: list[int]) -> tuple[list[int], list[int]]:
    """Kürzt den Bruch Faktor für Faktor, wie man es auf Papier macht.

    Jeder Nennerfaktor wird (vom größten her) gegen Zählerfaktoren gekürzt.
    Faktoren 1 fallen am Ende weg.
    """
    z = list(zaehler)
    rest_nenner = []
    for d in sorted(nenner, reverse=True):
        for i, f in enumerate(z):
            if d == 1:
                break
            g = math.gcd(d, f)
            if g > 1:
                z[i] //= g
                d //= g
        if d != 1:
            rest_nenner.append(d)
    return [f for f in z if f != 1], rest_nenner


def berechne_binom(n: int, k: int, symmetrie: bool = True) -> BinomDaten:
    """Berechnet C(n, k) mit Zähler, Nenner und gekürzter Form."""
    hinweise = pruefe_nichtnegativ(n=n, k=k)
    if k > n >= 0:
        hinweise.append(f"Hinweis: k = {k} > n = {n}. Man kann nicht mehr Elemente "
                        f"auswählen, als da sind. Also C({n},{k}) = 0.")
    if hinweise:
        return BinomDaten(n, k, k, [], [], [], 0, hinweise)
    kr = n - k if (symmetrie and n - k < k) else k
    zaehler = list(range(n, n - kr, -1))
    nenner = list(range(1, kr + 1))
    gekuerzt, rest = kuerze(zaehler, nenner)
    if rest:  # sollte nicht vorkommen, da C(n,k) ganzzahlig ist
        gekuerzt = [math.prod(zaehler) // math.prod(nenner)]
    return BinomDaten(n, k, kr, zaehler, nenner, gekuerzt, math.comb(n, k))


def binom_zeilen(d: BinomDaten) -> list[str]:
    """Rechenzeile für C(n, k), z. B. C(8,4) = 8·7·6·5/(1·2·3·4) = 7·2·5 = 70."""
    if d.hinweise:
        return list(d.hinweise) + [f"{C(d.n, d.k)} = 0"]
    zeilen = []
    if d.k_rechnung != d.k:
        zeilen.append(f"Symmetrie: {C(d.n, d.k)} = {C(d.n, d.k_rechnung)} "
                      f"(weniger Faktoren)")
    links = C(d.n, d.k) if d.k_rechnung == d.k else C(d.n, d.k_rechnung)
    if d.k_rechnung == 0:
        zeilen.append(f"{links} = 1   (Festlegung C(n,0) = 1)")
        return zeilen
    teile = [links,
             f"{produkt_text(d.zaehler)} / ({produkt_text(d.nenner)})"]
    gek = produkt_text(d.gekuerzt)
    if len(d.gekuerzt) > 1:
        teile.append(gek)
    teile.append(str(d.wert))
    zeilen.append(" = ".join(teile))
    return zeilen


def text_binom(d: BinomDaten) -> str:
    z = ["Binomialkoeffizient C(n,k) = \"n über k\"",
         "= Anzahl der k-elementigen Teilmengen einer n-elementigen Menge.",
         "Formel: C(n,k) = n(n−1)···(n−k+1) / (1·2···k) = n! / (k!(n−k)!)",
         ""]
    z += binom_zeilen(d)
    z += ["", f"Ergebnis: {C(d.n, d.k)} = {d.wert}"]
    if not d.hinweise:
        z.append(f"Probe mit Fakultäten: {d.n}! / ({d.k}!·{d.n - d.k}!) = "
                 f"{math.factorial(d.n)} / ({math.factorial(d.k)}·"
                 f"{math.factorial(d.n - d.k)}) = {d.wert}")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Geordnete Auswahl, Permutationen
# ---------------------------------------------------------------------------


@dataclass
class GeordnetDaten:
    n: int
    k: int
    wiederholung: bool
    faktoren: list[int]
    wert: int
    hinweise: list[str] = field(default_factory=list)


def berechne_geordnet(n: int, k: int, wiederholung: bool = False) -> GeordnetDaten:
    """Anzahl geordneter Auswahlen von k aus n Elementen."""
    hinweise = pruefe_nichtnegativ(n=n, k=k)
    if hinweise:
        return GeordnetDaten(n, k, wiederholung, [], 0, hinweise)
    if wiederholung:
        return GeordnetDaten(n, k, True, [n] * k, n ** k)
    if k > n:
        hinweise.append(f"Hinweis: Ohne Wiederholung kann man nicht k = {k} > n = {n} "
                        "Elemente auswählen. Die Anzahl ist 0.")
        return GeordnetDaten(n, k, False, [], 0, hinweise)
    faktoren = list(range(n, n - k, -1))
    return GeordnetDaten(n, k, False, faktoren, math.prod(faktoren))


def text_geordnet(d: GeordnetDaten, ueberschrift: str | None = None) -> str:
    art = "mit Wiederholung" if d.wiederholung else "ohne Wiederholung"
    z = [ueberschrift or f"Geordnete Auswahl von k = {d.k} aus n = {d.n} Elementen, {art}.",
         "Geordnet heißt: Die Reihenfolge zählt."]
    if d.hinweise:
        return "\n".join(z + d.hinweise + ["", "Ergebnis: 0"])
    z.append("")
    z.append("Wahl des 1. Elementes, 2. Elementes, …, k-ten Elementes:")
    z.append("  " + ", ".join(str(f) for f in d.faktoren) + " Möglichkeiten"
             if d.faktoren else "  keine Wahl (k = 0): 1 Möglichkeit")
    if d.wiederholung:
        z.append(f"Formel: n^k = {d.n}^{d.k} = {d.wert}")
    else:
        z.append("Formel (fallende Faktorielle): n(n−1)···(n−k+1) = n!/(n−k)!")
        z.append(f"  = {produkt_text(d.faktoren)} = {d.wert}")
        z.append(f"Probe: {d.n}!/{d.n - d.k}! = {math.factorial(d.n)}/"
                 f"{math.factorial(d.n - d.k)} = "
                 f"{math.factorial(d.n) // math.factorial(d.n - d.k)}")
    z += ["", f"Ergebnis: {d.wert}"]
    return "\n".join(z)


def berechne_permutation(n: int) -> GeordnetDaten:
    """Permutationen = geordnete Auswahl ohne Wiederholung mit k = n."""
    return berechne_geordnet(n, n, False)


def text_permutation(d: GeordnetDaten) -> str:
    kopf = (f"Permutationen von n = {d.n} Elementen.\n"
            "Eine Permutation ist eine geordnete Auswahl ohne Wiederholung mit k = n.")
    z = [kopf]
    if d.hinweise:
        return "\n".join(z + d.hinweise + ["", "Ergebnis: 0"])
    z.append("")
    if d.n == 0:
        z.append("0! = 1 (Festlegung)")
    else:
        z.append(f"{d.n}! = {produkt_text(list(range(1, d.n + 1)))} = {d.wert}")
    z += ["", f"Ergebnis: {d.n}! = {d.wert}"]
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kombinationen (ungeordnet), auch mit Wiederholung
# ---------------------------------------------------------------------------


@dataclass
class KombinationDaten:
    n: int
    k: int
    wiederholung: bool
    binom: BinomDaten
    wert: int
    bijektion: dict | None = None
    hinweise: list[str] = field(default_factory=list)


def bijektion_wiederholung(n: int, auswahl: list[int]) -> dict:
    """Bijektion A → A' (mit n−1 Nullen) → A'' (Plätze der Nullen).

    Vor m₁ stehen m₁−1 Nullen, zwischen mᵢ₋₁ und mᵢ stehen mᵢ−mᵢ₋₁ Nullen,
    hinter m_k stehen n−m_k Nullen.
    """
    a = sorted(auswahl)
    wort: list[int] = []
    vorher = 1
    for m in a:
        wort += [0] * (m - vorher)
        wort.append(m)
        vorher = m
    wort += [0] * (n - vorher)
    plaetze = [i + 1 for i, x in enumerate(wort) if x == 0]
    return {"A": a, "A1": wort, "A2": plaetze}


def wort_text(wort: list[int]) -> str:
    """Gruppiert gleiche Nachbarn, z. B. 111 0 2 00 44 0 5."""
    gruppen = ["".join(str(x) for x in g) for _, g in itertools.groupby(wort)]
    return " ".join(gruppen)


def berechne_kombination(n: int, k: int, wiederholung: bool = False,
                         auswahl: list[int] | None = None) -> KombinationDaten:
    hinweise = pruefe_nichtnegativ(n=n, k=k)
    if wiederholung:
        if n == 0 and k > 0:
            hinweise.append("Hinweis: Aus einer leeren Menge kann man nichts auswählen.")
        b = berechne_binom(n + k - 1, k) if n + k - 1 >= 0 else berechne_binom(0, 0)
        wert = 0 if hinweise else binom(n + k - 1, k)
        bij = None
        if auswahl is not None:
            falsch = [m for m in auswahl if not 1 <= m <= n]
            if falsch:
                hinweise.append(f"Hinweis: Die Auswahl enthält {falsch}. "
                                f"Erlaubt sind nur Zahlen von 1 bis n = {n}.")
            elif len(auswahl) != k:
                hinweise.append(f"Hinweis: Die Auswahl hat {len(auswahl)} Elemente, "
                                f"aber k = {k}.")
            else:
                bij = bijektion_wiederholung(n, auswahl)
        return KombinationDaten(n, k, True, b, wert, bij, hinweise)
    b = berechne_binom(n, k)
    return KombinationDaten(n, k, False, b, b.wert, None, hinweise + b.hinweise)


def text_kombination(d: KombinationDaten) -> str:
    if not d.wiederholung:
        z = [f"Kombinationen von k = {d.k} aus n = {d.n} Elementen ohne Wiederholung.",
             "Kombination heißt: ungeordnete Auswahl. Die Reihenfolge zählt nicht.",
             "Anzahl = C(n,k).", ""]
        z += binom_zeilen(d.binom)
        z += ["", f"Ergebnis: {C(d.n, d.k)} = {d.wert}"]
        return "\n".join(z)
    n, k = d.n, d.k
    z = [f"Kombinationen von k = {k} aus n = {n} Elementen mit Wiederholung.",
         "Idee (Bijektion): Man ordnet die Auswahl m₁ ≤ m₂ ≤ … ≤ m_k",
         "und fügt n−1 Nullen ein. Die Plätze der Nullen bilden eine",
         f"(n−1)-elementige Teilmenge von N_(n+k−1) = N_{n + k - 1}."]
    for h in d.hinweise:
        z.append(h)
    if d.bijektion:
        b = d.bijektion
        z += ["", "Bijektion für die gegebene Auswahl:",
              f"  A   = {wort_text(b['A'])}",
              f"  A'  = {wort_text(b['A1'])}   ({n - 1} Nullen)",
              "  A'' = {" + ", ".join(str(p) for p in b["A2"]) + "}   (Plätze der Nullen)"]
    if d.hinweise and d.wert == 0:
        return "\n".join(z + ["", "Ergebnis: 0"])
    z += ["", f"Anzahl = C(n+k−1, k) = C(n+k−1, n−1) = {C(n + k - 1, k)} "
              f"= {C(n + k - 1, n - 1)}", ""]
    z += binom_zeilen(d.binom)
    z += ["", f"Ergebnis: {C(n + k - 1, k)} = {d.wert}"]
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Gitterwege
# ---------------------------------------------------------------------------

GITTER_VARIANTEN = {
    "kaestchen": "m×n bedeutet m×n Kästchen, also (m+1)×(n+1) Gitterpunkte.",
    "punkte": "m×n bedeutet m×n Gitterpunkte, also (m−1)×(n−1) Kästchen.",
}


@dataclass
class GitterDaten:
    m: int                   # Zeilen (Schritte nach oben)
    n: int                   # Spalten (Schritte nach rechts)
    variante: str
    oben: int
    rechts: int
    binom: BinomDaten
    wert: int
    wege: list[str] | None = None
    beispiel: str | None = None
    hinweise: list[str] = field(default_factory=list)


def berechne_gitter(m: int, n: int, variante: str = "kaestchen",
                    liste: bool = False, beispiel: str | None = None) -> GitterDaten:
    """Anzahl der kürzesten Wege von links unten nach rechts oben.

    m ist die Höhe (Schritte nach oben), n die Breite (Schritte nach rechts).
    """
    if variante not in GITTER_VARIANTEN:
        raise ValueError(f"Unbekannte Variante: {variante}")
    hinweise = pruefe_nichtnegativ(m=m, n=n)
    if variante == "punkte":
        oben, rechts = m - 1, n - 1
        if not hinweise and (m < 1 or n < 1):
            hinweise.append("Hinweis: Bei der Variante 'punkte' braucht das Gitter "
                            "mindestens 1×1 Punkte.")
    else:
        oben, rechts = m, n
    if hinweise:
        return GitterDaten(m, n, variante, oben, rechts, berechne_binom(0, 0), 0,
                           None, None, hinweise)
    b = berechne_binom(oben + rechts, rechts, symmetrie=False)
    wege = None
    if liste:
        wege = []
        for plaetze in itertools.combinations(range(oben + rechts), rechts):
            w = ["o"] * (oben + rechts)
            for p in plaetze:
                w[p] = "r"
            wege.append(" ".join(w))
    if beispiel is not None:
        w = beispiel.replace(" ", "").lower()
        if w.count("r") != rechts or w.count("o") != oben or len(w) != oben + rechts:
            hinweise.append(f"Hinweis: Der Beispielweg '{beispiel}' passt nicht. "
                            f"Er braucht genau {rechts}× r und {oben}× o.")
            beispiel = None
        else:
            beispiel = " ".join(w)
    return GitterDaten(m, n, variante, oben, rechts, b, b.wert, wege, beispiel, hinweise)


def text_gitter(d: GitterDaten) -> str:
    z = [f"Gitterwege im {d.m}×{d.n}-Gitter von links unten nach rechts oben.",
         "Erlaubt sind nur Schritte nach rechts (r) und nach oben (o).",
         f"Annahme (Variante '{d.variante}'): {GITTER_VARIANTEN[d.variante]}"]
    if d.hinweise and d.wert == 0:
        return "\n".join(z + d.hinweise + ["", "Ergebnis: 0"])
    z += d.hinweise
    s = d.oben + d.rechts
    z += ["",
          f"Diese Wege setzen sich aus insgesamt {s} Strecken zusammen,",
          f"wobei {d.oben} Strecken nach oben (o) und {d.rechts} Strecken "
          f"nach rechts (r) verlaufen."]
    if d.beispiel:
        z.append(f"Ein Weg kann geschrieben werden als:  {d.beispiel}")
    z += [f"Man hat also eine Bijektion zwischen den Wegen und der Menge der",
          f"Kombinationen von {d.rechts} Strecken aus {s}. "
          f"(Man wählt die {d.rechts} Plätze, auf denen r steht.)",
          "Das sind keine Kombinationen mit Wiederholung: Jeder Platz wird "
          "höchstens einmal gewählt.",
          "",
          f"Anzahl = C(m+n, n) = C(m+n, m) = {C(s, d.rechts)}"
          + (f" = {C(s, d.oben)}" if d.oben != d.rechts else "")]
    z += ["⇒ " + zeile for zeile in binom_zeilen(d.binom)]
    if d.wege is not None:
        z += ["", f"Alle {len(d.wege)} Wege:"]
        z += [f"  {i:>3}. {w}" for i, w in enumerate(d.wege, 1)]
    z += ["", f"Ergebnis: {d.wert} Wege"]
    z.append(f"Probe: {C(s, d.oben)} = {binom(s, d.oben)} (Symmetrie C(m+n,m) = C(m+n,n))")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Derangements (Beispiel: vertauschte Hüte)
# ---------------------------------------------------------------------------


@dataclass
class DerangementDaten:
    n: int
    alphas: list[int]            # α_k = n!/k! für k = 0 … n (α₀ = n!)
    summe_brueche: list[Fraction]  # (−1)^k / k!
    zaehler: int                 # d_n (Zähler über n!)
    wert: int
    wahrscheinlichkeit: Fraction
    liste: list[tuple[int, ...]] | None = None
    hinweise: list[str] = field(default_factory=list)


def derangement_zahl(n: int) -> int:
    """dₙ = Σ (−1)^k n!/k!."""
    return sum((-1) ** k * (math.factorial(n) // math.factorial(k)) for k in range(n + 1))


def alle_derangements(n: int) -> list[tuple[int, ...]]:
    """Alle fixpunktfreien Permutationen von 1…n in Kurzschreibweise π(1)…π(n)."""
    return [p for p in itertools.permutations(range(1, n + 1))
            if all(p[i] != i + 1 for i in range(n))]


def berechne_derangement(n: int, liste: bool = False) -> DerangementDaten:
    hinweise = pruefe_nichtnegativ(n=n)
    if hinweise:
        return DerangementDaten(n, [], [], 0, 0, Fraction(0), None, hinweise)
    fak = math.factorial(n)
    alphas = [fak // math.factorial(k) for k in range(n + 1)]
    brueche = [Fraction((-1) ** k, math.factorial(k)) for k in range(n + 1)]
    wert = derangement_zahl(n)
    lst = None
    if liste:
        if n > 8:
            hinweise.append(f"Hinweis: Für n = {n} gibt es {wert} Derangements. "
                            "Die Liste wird nur bis n = 8 ausgegeben.")
        else:
            lst = alle_derangements(n)
    return DerangementDaten(n, alphas, brueche, wert, wert, Fraction(wert, fak),
                            lst, hinweise)


def _bruch_glied(k: int) -> str:
    return "1" if k == 0 else f"1/{k}!"


def text_derangement(d: DerangementDaten) -> str:
    n = d.n
    z = [f"Derangements (fixpunktfreie Permutationen) von n = {n} Elementen.",
         "Ein Fixpunkt von π ist ein i mit π(i) = i.",
         "Ein Derangement ist eine Permutation ohne Fixpunkt."]
    if d.hinweise and not d.alphas:
        return "\n".join(z + d.hinweise + ["", "Ergebnis: 0"])
    fak = math.factorial(n)
    # Zeile 1: dₙ = n!(1 − 1/1! + 1/2! − …)
    glieder = ""
    for k in range(n + 1):
        zeichen = "" if k == 0 else (" − " if k % 2 else " + ")
        glieder += zeichen + _bruch_glied(k)
    zeile = f"d_{n} = {n}!({glieder})"
    # Zeile 2: 1 − 1/1! fällt weg, Rest mit Brüchen 1/2 − 1/6 + 1/24 …
    if n >= 2:
        rest = ""
        for k in range(2, n + 1):
            zeichen = "" if k == 2 else (" − " if k % 2 else " + ")
            rest += zeichen + f"1/{math.factorial(k)}"
        zeile += f" = {n}!({rest})"
    zeile += f" = {n}!·{d.zaehler}/{fak} = {d.wert}"
    z += ["", "Formel (Inklusion-Exklusion): dₙ = n!(1 − 1/1! + 1/2! − … + (−1)ⁿ/n!)", zeile]
    p = d.wahrscheinlichkeit
    ptext = f"P({n}) = {d.wert}/{fak}"
    if p.denominator != fak and p != 0:
        ptext += f" = {p.numerator}/{p.denominator}"
    ptext += f" = {komma(p)}"
    z += ["", "Wahrscheinlichkeit, dass keiner seinen eigenen Hut bekommt:",
          f"{ptext}      Grenzwert: P = 1/e = 0,3678…"]
    for h in d.hinweise:
        z.append(h)
    if d.liste is not None:
        z += ["", f"Die {len(d.liste)} Derangements (Kurzschreibweise π(1) π(2) … π({n})):"]
        texte = [" ".join(str(x) for x in p) for p in d.liste]
        for i in range(0, len(texte), 6):
            z.append("  " + ",   ".join(texte[i:i + 6]))
    z += ["", f"Ergebnis: d_{n} = {d.wert}"]
    # Probe mit α_k = n!/k!
    teile = []
    for k, a in enumerate(d.alphas):
        teile.append(str(a) if k == 0 else (f" − {a}" if k % 2 else f" + {a}"))
    z.append(f"Probe mit α_k = C(n,k)·(n−k)! = n!/k!:  d_{n} = {''.join(teile)} = "
             f"{sum((-1) ** k * a for k, a in enumerate(d.alphas))}")
    if d.liste is not None:
        z.append(f"Probe durch Abzählen der Liste: {len(d.liste)}")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Partitionen
# ---------------------------------------------------------------------------

EINSCHRAENKUNGEN = {
    "keine": "",
    "ungerade": "Teile sind ungerade",
    "verschieden": "Teile sind pw. verschieden",
    "max-teil": "Größe der Teile ist ≤ k",
    "max-anzahl": "Anzahl der Teile ist ≤ k",
    "anzahl": "genau k Teile",
}


def partitionen(n: int, max_teil: int | None = None) -> list[tuple[int, ...]]:
    """Alle Partitionen von n, absteigend sortiert (5, 4+1, 3+2, …)."""
    if max_teil is None:
        max_teil = n
    if n == 0:
        return [()]
    ergebnis = []
    for erster in range(min(n, max_teil), 0, -1):
        for rest in partitionen(n - erster, erster):
            ergebnis.append((erster,) + rest)
    return ergebnis


def erfuellt(p: tuple[int, ...], art: str, k: int | None) -> bool:
    if art == "keine":
        return True
    if art == "ungerade":
        return all(t % 2 == 1 for t in p)
    if art == "verschieden":
        return len(set(p)) == len(p)
    if art == "max-teil":
        return all(t <= k for t in p)
    if art == "max-anzahl":
        return len(p) <= k
    if art == "anzahl":
        return len(p) == k
    raise ValueError(art)


def partitionszahl(n: int) -> int:
    """p(n) ohne Liste (dynamische Programmierung), für große n."""
    tab = [1] + [0] * n
    for teil in range(1, n + 1):
        for s in range(teil, n + 1):
            tab[s] += tab[s - teil]
    return tab[n]


def partition_text(p: tuple[int, ...]) -> str:
    return "+".join(str(t) for t in p) if p else "(leer)"


@dataclass
class PartitionDaten:
    n: int
    art: str
    k: int | None
    wert: int
    liste: list[tuple[int, ...]] | None
    probe: tuple[str, int] | None = None
    hinweise: list[str] = field(default_factory=list)


def berechne_partition(n: int, art: str = "keine", k: int | None = None,
                       liste: bool = True) -> PartitionDaten:
    hinweise = pruefe_nichtnegativ(n=n)
    if art in ("max-teil", "max-anzahl", "anzahl") and k is None:
        hinweise.append(f"Hinweis: Die Einschränkung '{art}' braucht --k.")
    if hinweise:
        return PartitionDaten(n, art, k, 0, None, None, hinweise)
    if n > 40:
        # Liste wäre zu lang, nur p(n) ohne Einschränkung ist schnell
        if art != "keine":
            hinweise.append("Hinweis: Mit Einschränkung rechnet das Skript nur bis n = 40.")
            return PartitionDaten(n, art, k, 0, None, None, hinweise)
        return PartitionDaten(n, art, k, partitionszahl(n), None, None,
                              ["Hinweis: n > 40. Die Liste wird nicht ausgegeben."])
    alle = [p for p in partitionen(n) if erfuellt(p, art, k)]
    probe = None
    gegen = {"ungerade": "verschieden", "verschieden": "ungerade",
             "max-teil": "max-anzahl", "max-anzahl": "max-teil"}
    if art in gegen:
        andere = gegen[art]
        probe = (andere, sum(1 for p in partitionen(n) if erfuellt(p, andere, k)))
    return PartitionDaten(n, art, k, len(alle), alle if liste else None, probe)


def einschraenkung_text(art: str, k: int | None) -> str:
    t = EINSCHRAENKUNGEN[art]
    return t.replace("k", str(k)) if k is not None and "k" in t else t


def text_partition(d: PartitionDaten) -> str:
    bed = einschraenkung_text(d.art, d.k)
    name = f"p({d.n})" if d.art == "keine" else f"p({d.n} | {bed})"
    z = [f"Partitionen der Zahl n = {d.n}.",
         "Eine Partition von n ist eine Summe n = n₁ + n₂ + … + n_k mit Teilen nᵢ > 0.",
         "Die Reihenfolge der Teile zählt nicht. Man schreibt die Teile absteigend."]
    if bed:
        z.append(f"Einschränkung: {bed}")
    z += d.hinweise
    if d.hinweise and d.wert == 0:
        return "\n".join(z + ["", f"Ergebnis: {name} = 0 (nicht berechnet)"])
    if d.liste is not None:
        z += ["", f"Alle Partitionen ({len(d.liste)} Stück):"]
        z += ["  " + partition_text(p) for p in d.liste]
    z += ["", f"Ergebnis: {name} = {d.wert}"]
    if d.probe:
        andere, wert = d.probe
        satz = ("Euler-Identität" if andere in ("ungerade", "verschieden")
                else "Satz (Ferrers-Diagramm spiegeln)")
        z.append(f"Probe ({satz}): p({d.n} | {einschraenkung_text(andere, d.k)}) = {wert}"
                 + ("  ✓" if wert == d.wert else "  ✗ (Abweichung!)"))
    elif d.art == "keine" and d.liste is not None:
        z.append(f"Probe (Rekursion ohne Liste): p({d.n}) = {partitionszahl(d.n)}")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Ferrers-Diagramm
# ---------------------------------------------------------------------------


@dataclass
class FerrersDaten:
    teile: list[int]
    n: int
    konjugiert: list[int]
    hinweise: list[str] = field(default_factory=list)


def berechne_ferrers(teile: list[int]) -> FerrersDaten:
    hinweise = []
    if any(t <= 0 for t in teile):
        hinweise.append("Hinweis: Alle Teile müssen > 0 sein. Nicht positive Teile "
                        "werden weggelassen.")
        teile = [t for t in teile if t > 0]
    sortiert = sorted(teile, reverse=True)
    if sortiert != list(teile):
        hinweise.append("Hinweis: Die Teile werden absteigend sortiert.")
    konj = [sum(1 for t in sortiert if t >= j) for j in range(1, (sortiert[0] if sortiert else 0) + 1)]
    return FerrersDaten(sortiert, sum(sortiert), konj, hinweise)


def text_ferrers(d: FerrersDaten) -> str:
    z = ["Ferrers-Diagramm einer Partition.",
         "Jede Zeile ist ein Teil. Die Zeilenlänge ist die Größe des Teils."]
    z += d.hinweise
    z += ["", f"{d.n} = {' + '.join(map(str, d.teile))}", "",
          "Größe der Teile →"]
    z += ["  " + " ".join("•" * t) for t in d.teile]
    z.append("↓ Anzahl der Teile")
    z += ["",
          "Das Diagramm von oben nach unten (spaltenweise) gelesen ergibt:",
          f"{d.n} = {' + '.join(map(str, d.konjugiert))}",
          "",
          f"Ergebnis: konjugierte Partition {d.n} = {' + '.join(map(str, d.konjugiert))}",
          f"Probe: Summe = {sum(d.konjugiert)}; Anzahl der Teile {len(d.konjugiert)} "
          f"= größter Teil vorher ({d.teile[0] if d.teile else 0})"]
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Binomischer Lehrsatz und Summen
# ---------------------------------------------------------------------------


@dataclass
class LehrsatzDaten:
    n: int
    koeffizienten: list[int]
    a: int | None
    b: int | None
    summanden: list[int] | None
    wert: int | None
    hinweise: list[str] = field(default_factory=list)


def berechne_lehrsatz(n: int, a: int | None = None, b: int | None = None) -> LehrsatzDaten:
    hinweise = pruefe_nichtnegativ(n=n)
    if hinweise:
        return LehrsatzDaten(n, [], a, b, None, None, hinweise)
    koeff = [math.comb(n, k) for k in range(n + 1)]
    if a is None or b is None:
        return LehrsatzDaten(n, koeff, a, b, None, None)
    summanden = [koeff[k] * a ** (n - k) * b ** k for k in range(n + 1)]
    return LehrsatzDaten(n, koeff, a, b, summanden, sum(summanden))


def _potenz(basis: str, e: int) -> str:
    if e == 0:
        return ""
    return basis if e == 1 else f"{basis}^{e}"


def text_lehrsatz(d: LehrsatzDaten) -> str:
    n = d.n
    z = ["Binomischer Lehrsatz:",
         "(a+b)ⁿ = C(n,0)aⁿ + C(n,1)aⁿ⁻¹b + … + C(n,k)aⁿ⁻ᵏbᵏ + … + C(n,n)bⁿ",
         "Grund: Aus jedem der n Faktoren (a+b) wählt man a oder b.",
         "Der Term aⁿ⁻ᵏbᵏ entsteht, wenn man k-mal b wählt. Das geht auf C(n,k) Arten."]
    if d.hinweise:
        return "\n".join(z + d.hinweise)
    glieder = []
    for k, c in enumerate(d.koeffizienten):
        teil = (_potenz("a", n - k) + _potenz("b", k)) or "1"
        glieder.append(f"{c}{teil}" if c != 1 or teil == "1" else teil)
    z += ["", f"(a+b)^{n} = " + " + ".join(glieder)]
    z.append("Koeffizienten: " + ", ".join(f"C({n},{k}) = {c}"
                                           for k, c in enumerate(d.koeffizienten)))
    if d.summanden is not None:
        a, b = d.a, d.b
        z += ["", f"Einsetzen a = {a}, b = {b}:"]
        teile = [f"{c}·({a})^{n - k}·({b})^{k}" for k, c in enumerate(d.koeffizienten)]
        z.append(f"({a}+({b}))^{n} = " + " + ".join(teile))
        z.append("= " + " + ".join(f"({s})" if s < 0 else str(s) for s in d.summanden))
        z += ["", f"Ergebnis: ({a}+({b}))^{n} = {d.wert}",
              f"Probe: ({a + b})^{n} = {(a + b) ** n}"]
    else:
        z += ["", f"Ergebnis: (a+b)^{n} = " + " + ".join(glieder),
              f"Probe: Summe der Koeffizienten = {sum(d.koeffizienten)} = 2^{n} = {2 ** n}"]
    return "\n".join(z)


SUMMEN_ARTEN = {
    "alle": "C(n,0) + C(n,1) + … + C(n,n) = 2ⁿ",
    "wechselnd": "C(n,0) − C(n,1) + … + (−1)ⁿC(n,n) = 0",
    "alternierend": "C(k,1) − C(k,2) + … + (−1)^(k−1)C(k,k) = 1",
    "quadrate": "C(n,0)² + C(n,1)² + … + C(n,n)² = C(2n,n)",
    "vandermonde": "C(m,0)C(n,r) + C(m,1)C(n,r−1) + … + C(m,r)C(n,0) = C(m+n,r)  (Vandermonde)",
}


@dataclass
class SummeDaten:
    art: str
    parameter: dict
    summanden: list[int]
    wert: int
    formelwert: int
    hinweise: list[str] = field(default_factory=list)


def berechne_summe(art: str, n: int, m: int | None = None, r: int | None = None) -> SummeDaten:
    """Rechnet eine Summe aus Binomialkoeffizienten Glied für Glied aus.

    Bei 'vandermonde' sind m, n, r die Zahlen aus C(m+n, r).
    Bei 'alternierend' ist n das k aus C(k,1) − C(k,2) + … .
    """
    hinweise = pruefe_nichtnegativ(n=n, m=m or 0, r=r or 0)
    if hinweise:
        return SummeDaten(art, {"n": n, "k": n, "m": m, "r": r}, [], 0, 0, hinweise)
    if art == "alle":
        s = [math.comb(n, k) for k in range(n + 1)]
        formel = 2 ** n
        par = {"n": n}
    elif art == "wechselnd":
        s = [(-1) ** k * math.comb(n, k) for k in range(n + 1)]
        formel = 0 if n > 0 else 1
        par = {"n": n}
        if n == 0:
            hinweise.append("Hinweis: Die Formel gilt nur für n > 0. Für n = 0 ist die Summe 1.")
    elif art == "alternierend":
        s = [(-1) ** (j - 1) * math.comb(n, j) for j in range(1, n + 1)]
        formel = 1 if n > 0 else 0
        par = {"k": n}
        if n == 0:
            hinweise.append("Hinweis: Die Formel setzt k ≥ 1 voraus. Für k = 0 ist die Summe leer (0).")
    elif art == "quadrate":
        s = [math.comb(n, k) ** 2 for k in range(n + 1)]
        formel = math.comb(2 * n, n)
        par = {"n": n}
    elif art == "vandermonde":
        if m is None or r is None:
            raise ValueError("vandermonde braucht m, n und r")
        if not (m >= r and n >= r):
            hinweise.append(f"Hinweis: Die Formel setzt meist m ≥ r und n ≥ r voraus. "
                            "Die Formel gilt trotzdem, wenn man C(a,b) = 0 für b > a setzt.")
        s = [binom(m, k) * binom(n, r - k) for k in range(r + 1)]
        formel = binom(m + n, r)
        par = {"m": m, "n": n, "r": r}
    else:
        raise ValueError(f"Unbekannte Art: {art}")
    return SummeDaten(art, par, s, sum(s), formel, hinweise)


def text_summe(d: SummeDaten) -> str:
    z = [f"Summe mit Binomialkoeffizienten: {SUMMEN_ARTEN[d.art]}"]
    z += d.hinweise
    if not d.summanden and d.hinweise:
        return "\n".join(z + ["", "Ergebnis: nicht berechnet"])
    p = d.parameter
    z.append("")
    if d.art in ("alle", "wechselnd"):
        n = p["n"]
        a, b = (1, 1) if d.art == "alle" else (1, -1)
        z.append(f"Setzen im binomischen Lehrsatz a = {a}, b = {b}:")
        z.append(f"({a} + ({b}))^{n} = Σ C({n},k)·{a}^({n}−k)·({b})^k")
        glieder = " ".join(("+ " if x >= 0 else "− ") + str(abs(x)) for x in d.summanden)
        z.append(f"= {glieder.lstrip('+ ')}")
        z.append(f"= {d.wert}")
        z += ["", f"Ergebnis: {d.wert}",
              f"Probe: ({a + b})^{n} = {(a + b) ** n}"]
    elif d.art == "alternierend":
        k = p["k"]
        z.append("Setzen im binomischen Lehrsatz (a+b)^k = … a = 1, b = −1:")
        z.append(f"⇒ 1 − C({k},1) + C({k},2) − + … + (−1)^{k}·C({k},{k}) = (1 − 1)^{k} = 0")
        z.append(f"⇒ C({k},1) − C({k},2) + … + (−1)^{k - 1}·C({k},{k}) = 1")
        glieder = " ".join(("+ " if x >= 0 else "− ") + str(abs(x)) for x in d.summanden)
        z += ["", "Nachrechnen Glied für Glied:",
              f"{glieder.lstrip('+ ')} = {d.wert}",
              "", f"Ergebnis: {d.wert}",
              f"Probe: Formelwert 1 {'✓' if d.wert == d.formelwert else '✗'}"]
    elif d.art == "quadrate":
        n = p["n"]
        z.append(f"Koeffizient von x^{n} in (1+x)^{n}·(1+x)^{n} = (1+x)^{2 * n}:")
        z.append("Σ C(n,k)·C(n,n−k) = Σ C(n,k)²")
        z.append(" + ".join(f"{math.comb(n, k)}²" for k in range(n + 1))
                 + " = " + " + ".join(map(str, d.summanden)) + f" = {d.wert}")
        z += ["", f"Ergebnis: {d.wert}",
              f"Probe: C({2 * n},{n}) = {d.formelwert} "
              f"{'✓' if d.wert == d.formelwert else '✗'}"]
    else:
        m, n, r = p["m"], p["n"], p["r"]
        z.append(f"Koeffizient von x^{r} in (1+x)^{m}·(1+x)^{n} = (1+x)^{m + n}:")
        z.append("Die Potenz x^r entsteht, wenn C(m,k)x^k mit C(n,r−k)x^(r−k) "
                 "multipliziert wird, 0 ≤ k ≤ r.")
        z.append(" + ".join(f"C({m},{k})·C({n},{r - k})" for k in range(r + 1)))
        z.append("= " + " + ".join(f"{binom(m, k)}·{binom(n, r - k)}" for k in range(r + 1)))
        z.append("= " + " + ".join(map(str, d.summanden)) + f" = {d.wert}")
        z += ["", f"Ergebnis: {d.wert}",
              f"Probe: C({m + n},{r}) = {d.formelwert} "
              f"{'✓' if d.wert == d.formelwert else '✗'}"]
    return "\n".join(z)


@dataclass
class PascalDaten:
    zeilen: list[list[int]]


def berechne_pascal(n: int) -> PascalDaten:
    zeilen = [[1]]
    for _ in range(max(n, 0)):
        vorher = zeilen[-1]
        zeilen.append([1] + [vorher[i] + vorher[i + 1] for i in range(len(vorher) - 1)] + [1])
    return PascalDaten(zeilen)


def text_pascal(d: PascalDaten) -> str:
    z = ["Pascalsches Dreieck. Zeile n enthält C(n,0), C(n,1), …, C(n,n).",
         "Rekursion: C(n+1,k+1) = C(n,k) + C(n,k+1).", ""]
    breite = len("  ".join(map(str, d.zeilen[-1])))
    for n, zeile in enumerate(d.zeilen):
        z.append(f"n={n:<3}" + "  ".join(map(str, zeile)).center(breite))
    letzte = d.zeilen[-1]
    n = len(d.zeilen) - 1
    z += ["", f"Ergebnis: Zeile {n}: " + ", ".join(map(str, letzte)),
          f"Probe: Summe der Zeile = {sum(letzte)} = 2^{n}"]
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

BEISPIELE = """\
Beispiele:
  python -m skripte.kombinatorik gitter 4 4 --beispielweg rorooorr   (4×4-Gitter mit Beispielweg)
  python -m skripte.kombinatorik gitter 4 4 --alle-varianten   (beide Gitter-Varianten)
  python -m skripte.kombinatorik derangement 4 --liste   (Derangements von 4 Elementen)
  python -m skripte.kombinatorik partition 5   (alle Partitionen von 5)
  python -m skripte.kombinatorik partition 5 --einschraenkung ungerade   (nur ungerade Teile)
  python -m skripte.kombinatorik ferrers 8 6 6 3 2 1 1   (Ferrers-Diagramm)
  python -m skripte.kombinatorik summe alternierend 6   (alternierende Summe)
  python -m skripte.kombinatorik summe vandermonde 3 4 --r 2   (Vandermonde-Identität)
  python -m skripte.kombinatorik kombination 5 7 --wiederholung --auswahl 1 1 1 2 4 4 5   (Auswahl mit Wiederholung)
"""


def baue_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m skripte.kombinatorik",
        description="Zählformeln mit Lösungsweg (Kombinatorik). "
                    "Wählen Sie einen Unterbefehl. Hilfe zu einem Unterbefehl: "
                    "python -m skripte.kombinatorik UNTERBEFEHL --help",
        epilog=BEISPIELE, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="befehl", required=True, metavar="UNTERBEFEHL")

    def neu(name, hilfe, beispiel):
        return sub.add_parser(name, help=hilfe, description=hilfe,
                              epilog=f"Beispiel:\n  {beispiel}",
                              formatter_class=argparse.RawDescriptionHelpFormatter)

    p = neu("binom", "Binomialkoeffizient C(n,k) berechnen.",
            "python -m skripte.kombinatorik binom 8 4   (C(8,4) = 70)")
    p.add_argument("n", type=int, help="Größe der Menge n")
    p.add_argument("k", type=int, help="Anzahl der ausgewählten Elemente k")

    p = neu("geordnet", "Geordnete Auswahl von k aus n Elementen zählen.",
            "python -m skripte.kombinatorik geordnet 5 3   (3 aus 5 geordnet)")
    p.add_argument("n", type=int, help="Größe der Menge n")
    p.add_argument("k", type=int, help="Anzahl der ausgewählten Elemente k")
    p.add_argument("--wiederholung", action="store_true",
                   help="Elemente dürfen mehrfach gewählt werden (dann n^k)")

    p = neu("permutation", "Anzahl der Permutationen n! berechnen.",
            "python -m skripte.kombinatorik permutation 4   (4! = 24)")
    p.add_argument("n", type=int, help="Anzahl der Elemente n")

    p = neu("kombination", "Ungeordnete Auswahl von k aus n Elementen zählen.",
            "python -m skripte.kombinatorik kombination 5 7 --wiederholung "
            "--auswahl 1 1 1 2 4 4 5   (Auswahl mit Wiederholung)")
    p.add_argument("n", type=int, help="Größe der Menge n")
    p.add_argument("k", type=int, help="Anzahl der ausgewählten Elemente k")
    p.add_argument("--wiederholung", action="store_true",
                   help="Elemente dürfen mehrfach gewählt werden (dann C(n+k−1,k))")
    p.add_argument("--auswahl", type=int, nargs="+", metavar="M",
                   help="eine Auswahl m₁ … m_k aus 1 … n; zeigt die Bijektion A → A' → A''")

    p = neu("gitter", "Kürzeste Gitterwege von links unten nach rechts oben zählen.",
            "python -m skripte.kombinatorik gitter 4 4   (4×4-Gitter: 70 Wege)")
    p.add_argument("m", type=int, help="Höhe m des Gitters (Schritte nach oben)")
    p.add_argument("n", type=int, help="Breite n des Gitters (Schritte nach rechts)")
    p.add_argument("--variante", choices=list(GITTER_VARIANTEN), default="kaestchen",
                   help="Bedeutung von m×n: 'kaestchen' (Standard, m×n Kästchen) "
                        "oder 'punkte' (m×n Gitterpunkte)")
    p.add_argument("--alle-varianten", action="store_true",
                   help="beide Varianten nacheinander ausgeben")
    p.add_argument("--liste", action="store_true",
                   help="alle Wege als Folge von r und o ausgeben (nur bis 500 Wege)")
    p.add_argument("--beispielweg", metavar="WORT",
                   help="einen Weg als Wort aus r und o zeigen, z. B. rorooorr")

    p = neu("derangement", "Fixpunktfreie Permutationen dₙ zählen (vertauschte Hüte).",
            "python -m skripte.kombinatorik derangement 4 --liste   (d₄ = 9)")
    p.add_argument("n", type=int, help="Anzahl der Elemente n")
    p.add_argument("--liste", action="store_true",
                   help="alle Derangements auflisten (bis n = 8)")

    p = neu("partition", "Partitionen p(n) einer Zahl zählen und auflisten.",
            "python -m skripte.kombinatorik partition 5   (p(5) = 7)")
    p.add_argument("n", type=int, help="die Zahl n")
    p.add_argument("--einschraenkung", choices=list(EINSCHRAENKUNGEN), default="keine",
                   help="Einschränkung: ungerade Teile, verschiedene Teile, "
                        "größter Teil ≤ k, höchstens k Teile, genau k Teile")
    p.add_argument("--k", type=int, help="Grenze k für max-teil, max-anzahl, anzahl")
    p.add_argument("--ohne-liste", action="store_true",
                   help="nur die Anzahl ausgeben, ohne Liste")

    p = neu("ferrers", "Ferrers-Diagramm zeichnen und die konjugierte Partition ablesen.",
            "python -m skripte.kombinatorik ferrers 8 6 6 3 2 1 1   (konjugierte Partition)")
    p.add_argument("teile", type=int, nargs="+", help="die Teile der Partition")

    p = neu("lehrsatz", "Binomischer Lehrsatz: (a+b)ⁿ ausmultiplizieren.",
            "python -m skripte.kombinatorik lehrsatz 4 --a 2 --b -1   (a = 2, b = −1)")
    p.add_argument("n", type=int, help="Exponent n")
    p.add_argument("--a", type=int, help="Zahl für a (optional)")
    p.add_argument("--b", type=int, help="Zahl für b (optional)")

    p = neu("summe", "Summen mit Binomialkoeffizienten Glied für Glied ausrechnen.",
            "python -m skripte.kombinatorik summe alternierend 6   (Ergebnis 1)\n"
            "  python -m skripte.kombinatorik summe vandermonde 3 4 --r 2   (Vandermonde-Identität)")
    p.add_argument("art", choices=list(SUMMEN_ARTEN),
                   help="welche Summe: " + "; ".join(f"{k}: {v}" for k, v in SUMMEN_ARTEN.items()))
    p.add_argument("zahlen", type=int, nargs="+",
                   help="n (bei 'alternierend' das k); bei 'vandermonde' m und n")
    p.add_argument("--r", type=int, help="r bei 'vandermonde'")

    p = neu("pascal", "Pascalsches Dreieck bis Zeile n ausgeben.",
            "python -m skripte.kombinatorik pascal 6   (Zeilen 0 bis 6)")
    p.add_argument("n", type=int, help="letzte Zeile n")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = baue_parser().parse_args(argv)
    b = args.befehl
    if b == "binom":
        print(text_binom(berechne_binom(args.n, args.k)))
    elif b == "geordnet":
        print(text_geordnet(berechne_geordnet(args.n, args.k, args.wiederholung)))
    elif b == "permutation":
        print(text_permutation(berechne_permutation(args.n)))
    elif b == "kombination":
        print(text_kombination(berechne_kombination(args.n, args.k, args.wiederholung,
                                                    args.auswahl)))
    elif b == "gitter":
        varianten = list(GITTER_VARIANTEN) if args.alle_varianten else [args.variante]
        for i, v in enumerate(varianten):
            if len(varianten) > 1:
                print(("\n" if i else "") + "=" * 60 + f"\nVariante: {v}\n" + "=" * 60)
            liste = args.liste
            if liste:
                mm, nn = (args.m, args.n) if v == "kaestchen" else (args.m - 1, args.n - 1)
                if min(mm, nn) >= 0 and binom(mm + nn, nn) > 500:
                    print("Hinweis: Mehr als 500 Wege. Die Liste wird nicht ausgegeben.")
                    liste = False
            print(text_gitter(berechne_gitter(args.m, args.n, v, liste, args.beispielweg)))
    elif b == "derangement":
        print(text_derangement(berechne_derangement(args.n, args.liste)))
    elif b == "partition":
        print(text_partition(berechne_partition(args.n, args.einschraenkung, args.k,
                                                not args.ohne_liste)))
    elif b == "ferrers":
        print(text_ferrers(berechne_ferrers(args.teile)))
    elif b == "lehrsatz":
        if (args.a is None) != (args.b is None):
            print("Hinweis: Geben Sie a und b zusammen an. Das Skript zeigt nur die Formel.")
            args.a = args.b = None
        print(text_lehrsatz(berechne_lehrsatz(args.n, args.a, args.b)))
    elif b == "summe":
        if args.art == "vandermonde":
            if len(args.zahlen) != 2 or args.r is None:
                print("Hinweis: 'vandermonde' braucht zwei Zahlen m n und --r. "
                      "Beispiel: summe vandermonde 3 4 --r 2")
                return 2
            m, n = args.zahlen
            print(text_summe(berechne_summe("vandermonde", n, m, args.r)))
        else:
            if len(args.zahlen) != 1:
                print(f"Hinweis: '{args.art}' braucht genau eine Zahl. "
                      f"Das Skript nimmt die erste: {args.zahlen[0]}.")
            print(text_summe(berechne_summe(args.art, args.zahlen[0])))
    elif b == "pascal":
        print(text_pascal(berechne_pascal(args.n)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
