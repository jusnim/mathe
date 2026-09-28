"""Teilbarkeitsregeln (9, 11, 7) und Teilerfunktionen d(n), σ(n).

Zweck
-----
Das Skript zeigt den Rechenweg für diese Fälle:

* Regel für 9: Die Quersumme ist die Summe aller Ziffern.
  x ist durch 9 teilbar, wenn die Quersumme durch 9 teilbar ist.
* Regel für 11: Die alternierende Quersumme wechselt das Vorzeichen
  von Ziffer zu Ziffer: x₀ − x₁ + x₂ − …
* Regel für 7 (auch 11 und 13): Man teilt x von rechts in Dreierblöcke.
  Wegen 1000 ≡ −1 mod 1001 und 1001 = 7 · 11 · 13 gilt
  x ≡ y₀ − y₁ + y₂ − …
* Teileranzahl d(n) und Teilersumme σ(n) mit Formel.
* Vollkommene Zahlen: σ(n) = 2n.
* Mersenne-Zahlen 2ᵏ − 1 und die vollkommene Zahl 2ᵏ⁻¹(2ᵏ − 1).

Aufruf (Unterbefehle)
---------------------
    python -m skripte.teilbarkeit regeln 299792          # 9, 11 und 7
    python -m skripte.teilbarkeit neun 299792
    python -m skripte.teilbarkeit elf 299792
    python -m skripte.teilbarkeit sieben "10^19 + 1"
    python -m skripte.teilbarkeit teiler 28              # d(n), σ(n)
    python -m skripte.teilbarkeit vollkommen 6 28 496
    python -m skripte.teilbarkeit mersenne 2 3 5 11      # oder: mersenne --bis 13

Zahlen darf man als Ausdruck eingeben, zum Beispiel "10^19 + 1" oder "2^11 - 1".
Leerzeichen und Unterstriche in Zahlen werden ignoriert ("299 792").

Varianten (nur bei elf, sieben, regeln)
---------------------------------------
Die Vorzeichen der alternierenden Summe kann man von zwei Seiten her setzen.

* ``einer`` (Standard): Die Einerstelle (bzw. der rechte
  Dreierblock) bekommt „+“. Die Summe ist dann ≡ x.
* ``links``: Die erste Ziffer links (bzw. der linke Block) bekommt „+“.
  Die Summe ist dann ≡ (−1)ⁿ · x. Die Teilbarkeit ist gleich,
  der Rest kann aber ein anderer sein.

``--alle-varianten`` gibt beide Versionen nacheinander aus.
"""

from __future__ import annotations

import argparse
import ast
import sys
from dataclasses import dataclass, field

from sympy import factorint, isprime

VARIANTEN = {
    "einer": "Vorzeichen beginnt rechts (Einerstelle bzw. rechter Block) mit +. "
             "Die Summe ist dann ≡ x. Das ist der Standard.",
    "links": "Vorzeichen beginnt links (erste Ziffer bzw. linker Block) mit +. "
             "Die Summe ist dann ≡ (−1)ⁿ · x. Teilbar ja/nein bleibt gleich.",
}
STANDARD_VARIANTE = "einer"

_HOCH = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
_TIEF = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def hoch(k: int) -> str:
    """Schreibt k als hochgestellte Zahl (Exponent)."""
    return str(k).translate(_HOCH)


def tief(k: int) -> str:
    """Schreibt k als tiefgestellte Zahl (Index)."""
    return str(k).translate(_TIEF)


def zahl_mit_leerzeichen(x: int) -> str:
    """Gruppiert die Ziffern in Dreierblöcke: 299792 -> '299 792'."""
    return f"{x:,}".replace(",", " ")


def z(k: int) -> str:
    """Schreibt eine ganze Zahl mit echtem Minuszeichen: -2 -> '−2'."""
    return f"−{-k}" if k < 0 else str(k)


def vorzeichen_summe(vorzeichen: list[int], werte: list[int]) -> str:
    """Schreibt Terme mit Vorzeichen: ([-1, 1, -1], [2, 9, 9]) -> '−2 + 9 − 9'.

    Vorzeichen und Wert sind getrennt, damit auch „− 0“ richtig erscheint.
    """
    teile = []
    for i, (v, w) in enumerate(zip(vorzeichen, werte)):
        if i == 0:
            teile.append(f"−{w}" if v < 0 else f"{w}")
        else:
            teile.append(f"− {w}" if v < 0 else f"+ {w}")
    return " ".join(teile)


# ---------------------------------------------------------------------------
# Eingabe
# ---------------------------------------------------------------------------

_MAX_EXPONENT = 100_000


def zahl_lesen(text: str) -> int:
    """Liest eine ganze Zahl oder einen Ausdruck wie '10^19 + 1'.

    Erlaubt sind ganze Zahlen, +, −, ·/*, ^/** und Klammern.
    """
    bereinigt = text.replace(" ", "").replace("_", "").replace(" ", "")
    bereinigt = bereinigt.replace("^", "**").replace("·", "*").replace("−", "-")
    try:
        baum = ast.parse(bereinigt, mode="eval")
    except SyntaxError as fehler:
        raise ValueError(f"„{text}“ ist keine gültige Zahl.") from fehler

    def auswerten(knoten: ast.AST) -> int:
        if isinstance(knoten, ast.Expression):
            return auswerten(knoten.body)
        if isinstance(knoten, ast.Constant) and isinstance(knoten.value, int):
            return knoten.value
        if isinstance(knoten, ast.UnaryOp) and isinstance(knoten.op, (ast.USub, ast.UAdd)):
            wert = auswerten(knoten.operand)
            return -wert if isinstance(knoten.op, ast.USub) else wert
        if isinstance(knoten, ast.BinOp):
            links, rechts = auswerten(knoten.left), auswerten(knoten.right)
            if isinstance(knoten.op, ast.Add):
                return links + rechts
            if isinstance(knoten.op, ast.Sub):
                return links - rechts
            if isinstance(knoten.op, ast.Mult):
                return links * rechts
            if isinstance(knoten.op, ast.Pow):
                if rechts < 0 or rechts > _MAX_EXPONENT:
                    raise ValueError(f"Exponent {rechts} ist nicht erlaubt (0 … {_MAX_EXPONENT}).")
                return links ** rechts
        raise ValueError(f"„{text}“ enthält etwas, das keine ganze Zahl ist.")

    return auswerten(baum)


# ---------------------------------------------------------------------------
# Rechnen: Teilbarkeitsregeln
# ---------------------------------------------------------------------------

def ziffern(x: int) -> list[int]:
    """Ziffern von links nach rechts: 299792 -> [2, 9, 9, 7, 9, 2]."""
    return [int(z) for z in str(abs(x))]


def bloecke(x: int, breite: int = 3) -> list[int]:
    """Blöcke von links nach rechts, von rechts abgeteilt: 299792 -> [299, 792]."""
    x = abs(x)
    ergebnis = []
    basis = 10 ** breite
    while True:
        ergebnis.append(x % basis)
        x //= basis
        if x == 0:
            break
    return list(reversed(ergebnis))


@dataclass
class Regel9:
    x: int
    ziffern: list[int]
    quersumme: int
    rest: int                 # Quersumme mod 9 = x mod 9
    ohne_neunen: list[int]    # Ziffern ohne 9 (Kurzweg „Neunen streichen“)
    teilbar: bool


def regel_9(x: int) -> Regel9:
    """Regel für 9: x ≡ Quersumme (mod 9)."""
    z = ziffern(x)
    qs = sum(z)
    rest = qs % 9
    assert rest == abs(x) % 9
    return Regel9(abs(x), z, qs, rest, [d for d in z if d != 9], rest == 0)


@dataclass
class AlternierendeSumme:
    """Ergebnis für die 11-Regel (Ziffern) oder 7/11/13-Regel (Dreierblöcke)."""
    x: int
    art: str                  # "ziffern" oder "bloecke"
    variante: str
    teile: list[int]          # Ziffern bzw. Blöcke, von links nach rechts
    vorzeichen: list[int]     # +1 oder −1 je Teil, von links nach rechts
    terme: list[int]          # Teil mal Vorzeichen, von links nach rechts
    summe: int
    faktor: int               # summe ≡ faktor · x (faktor = ±1)
    module: list[int]         # geprüfte Teiler, z. B. [11] oder [7, 11, 13]
    reste_summe: dict[int, int] = field(default_factory=dict)   # summe mod p
    reste_x: dict[int, int] = field(default_factory=dict)       # x mod p
    teilbar: dict[int, bool] = field(default_factory=dict)


def _alternierend(x: int, teile: list[int], art: str, module: list[int],
                  variante: str) -> AlternierendeSumme:
    if variante not in VARIANTEN:
        raise ValueError(f"Unbekannte Variante „{variante}“.")
    n = len(teile) - 1  # höchster Index
    if variante == "einer":
        # Teil mit Index i (von rechts gezählt) bekommt (−1)^i.
        vorzeichen = [(-1) ** (n - pos) for pos in range(n + 1)]
        faktor = 1
    else:
        vorzeichen = [(-1) ** pos for pos in range(n + 1)]
        faktor = (-1) ** n
    terme = [v * t for v, t in zip(vorzeichen, teile)]
    summe = sum(terme)
    reste_summe = {p: summe % p for p in module}
    reste_x = {p: abs(x) % p for p in module}
    for p in module:  # Kontrolle: summe ≡ faktor · x
        assert (summe - faktor * abs(x)) % p == 0
    return AlternierendeSumme(abs(x), art, variante, teile, vorzeichen, terme, summe, faktor,
                              module, reste_summe, reste_x,
                              {p: reste_x[p] == 0 for p in module})


def regel_11(x: int, variante: str = STANDARD_VARIANTE) -> AlternierendeSumme:
    """Regel für 11: alternierende Quersumme der Ziffern."""
    return _alternierend(x, ziffern(x), "ziffern", [11], variante)


def regel_1001(x: int, variante: str = STANDARD_VARIANTE,
               module: tuple[int, ...] = (7, 11, 13)) -> AlternierendeSumme:
    """Regel für 7, 11, 13: alternierende Summe der Dreierblöcke (1000 ≡ −1 mod 1001)."""
    return _alternierend(x, bloecke(x), "bloecke", list(module), variante)


# ---------------------------------------------------------------------------
# Rechnen: Teilerfunktionen
# ---------------------------------------------------------------------------

@dataclass
class TeilerDaten:
    n: int
    faktoren: dict[int, int]      # p -> e
    teiler: list[int]
    d: int
    sigma: int
    d_formel_faktoren: list[int]  # (1 + eᵢ)
    sigma_formel_faktoren: list[int]  # (p^(1+e) − 1)/(p − 1)
    art: str                      # "vollkommen", "defizient", "abundant"


def teiler_daten(n: int) -> TeilerDaten:
    """Berechnet Teiler, d(n) und σ(n) direkt und mit der Formel."""
    if n < 1:
        raise ValueError("n muss mindestens 1 sein.")
    faktoren = dict(sorted(factorint(n).items()))
    teiler = [1]
    for p, e in faktoren.items():
        teiler = [t * p ** k for t in teiler for k in range(e + 1)]
    teiler.sort()
    d_f = [1 + e for e in faktoren.values()]
    s_f = [(p ** (1 + e) - 1) // (p - 1) for p, e in faktoren.items()]
    d = _produkt(d_f)
    sigma = _produkt(s_f)
    assert d == len(teiler) and sigma == sum(teiler)
    if sigma == 2 * n:
        art = "vollkommen"
    elif sigma < 2 * n:
        art = "defizient"
    else:
        art = "abundant"
    return TeilerDaten(n, faktoren, teiler, d, sigma, d_f, s_f, art)


def _produkt(werte: list[int]) -> int:
    ergebnis = 1
    for w in werte:
        ergebnis *= w
    return ergebnis


@dataclass
class MersenneDaten:
    k: int
    m: int                        # 2^k − 1
    ist_prim: bool
    faktoren: dict[int, int]
    k_zerlegung: tuple[int, int] | None  # (a, b) mit k = a·b, a, b > 1
    vollkommen: int | None        # 2^(k−1)(2^k − 1), falls m prim


def mersenne_daten(k: int) -> MersenneDaten:
    """Prüft 2ᵏ − 1 auf Primzahl und bildet die vollkommene Zahl."""
    if k < 1:
        raise ValueError("k muss mindestens 1 sein.")
    m = 2 ** k - 1
    prim = bool(isprime(m))
    faktoren = {} if m == 1 else ({m: 1} if prim else dict(sorted(factorint(m).items())))
    zerlegung = None
    for a in range(2, k):
        if k % a == 0:
            zerlegung = (a, k // a)
            break
    return MersenneDaten(k, m, prim, faktoren, zerlegung,
                         2 ** (k - 1) * m if prim else None)


# ---------------------------------------------------------------------------
# Ausgabe (Lösungsweg als Text)
# ---------------------------------------------------------------------------

def zerlegung_text(faktoren: dict[int, int]) -> str:
    if not faktoren:
        return "1"
    return " · ".join(f"{p}{hoch(e)}" if e > 1 else f"{p}" for p, e in faktoren.items())


def _rest_text(summe: int, p: int) -> str:
    """Zeigt summe = p · q + r mit 0 ≤ r < p."""
    q, r = divmod(summe, p)
    if 0 <= summe < p:
        return f"{summe} < {p}"
    return f"{z(summe)} = {p} · {q} + {r}" if q >= 0 else f"{z(summe)} = {p} · ({z(q)}) + {r}"


def text_9(r: Regel9) -> str:
    z = r.ziffern
    zeilen = [
        "Regel für 9:",
        "Wegen 10ʳ ≡ 1 (mod 9) gilt x ≡ Quersumme (mod 9).",
        "Die Quersumme ist die Summe aller Ziffern.",
        "",
        f"x = {zahl_mit_leerzeichen(r.x)}",
        f"x ≡ x{tief(len(z) - 1)} + … + x₁ + x₀ (mod 9)" if len(z) > 2
        else "x ≡ Summe der Ziffern (mod 9)",
        f"  = {' + '.join(map(str, z))}",
        f"  = {r.quersumme} ≡ {r.rest} (mod 9)",
    ]
    if 9 in z and r.ohne_neunen:
        zeilen.append(f"  (Kurzweg: Neunen streichen, denn 9 ≡ 0: "
                      f"{' + '.join(map(str, r.ohne_neunen))} = {sum(r.ohne_neunen)} "
                      f"≡ {r.rest} (mod 9))")
    zeilen.append("")
    if r.teilbar:
        zeilen.append(f"Ergebnis: x = {zahl_mit_leerzeichen(r.x)} ist durch 9 teilbar.")
    else:
        zeilen.append(f"Ergebnis: x = {zahl_mit_leerzeichen(r.x)} ist nicht durch 9 teilbar "
                      f"(Rest {r.rest}).")
    zeilen.append(f"Probe: {r.x} mod 9 = {r.x % 9}.")
    return "\n".join(zeilen)


def _annahme_text(variante: str) -> str:
    return f"Annahme (Variante „{variante}“): {VARIANTEN[variante]}"


def text_11(r: AlternierendeSumme) -> str:
    n = len(r.teile) - 1
    zeilen = [
        "Regel für 11:",
        "Wegen 10 ≡ −1 (mod 11) gilt 10ʳ ≡ (−1)ʳ (mod 11).",
        "Die alternierende Quersumme wechselt das Vorzeichen von Ziffer zu Ziffer.",
        _annahme_text(r.variante),
        "",
        f"x = {zahl_mit_leerzeichen(r.x)}",
    ]
    if r.variante == "einer":
        zeilen.append(f"x ≡ (−1){hoch(n)} x{tief(n)} + … − x₁ + x₀ (mod 11)")
    else:
        zeilen.append(f"(−1){hoch(n)} · x ≡ x{tief(n)} − x{tief(max(n - 1, 0))} "
                      f"+ … + (−1){hoch(n)} x₀ (mod 11)")
    rest = r.reste_summe[11]
    zeilen.append(f"  = {vorzeichen_summe(r.vorzeichen, r.teile)}")
    zeilen.append(f"  = {z(r.summe)}" + (f" ≡ {rest} (mod 11)" if r.summe != rest else " (mod 11)"))
    if r.faktor == -1:
        zeilen.append(f"Achtung: Diese Summe ist ≡ −x. Also x ≡ −({z(r.summe)}) ≡ {r.reste_x[11]} (mod 11).")
    zeilen.append("")
    if r.teilbar[11]:
        zeilen.append(f"Ergebnis: x = {zahl_mit_leerzeichen(r.x)} ist durch 11 teilbar.")
    else:
        zeilen.append(f"Ergebnis: x = {zahl_mit_leerzeichen(r.x)} ist nicht durch 11 teilbar "
                      f"(x ≡ {r.reste_x[11]} mod 11).")
    zeilen.append(f"Probe: {r.x} mod 11 = {r.x % 11}.")
    return "\n".join(zeilen)


def text_1001(r: AlternierendeSumme) -> str:
    m = len(r.teile) - 1
    blocktext = " ".join(f"{b:03d}" if i > 0 else str(b) for i, b in enumerate(r.teile))
    zeilen = [
        "Regel für 7 mit 1001:",
        "Es gilt 1001 = 7 · 11 · 13, also 1000 ≡ −1 (mod 7), (mod 11) und (mod 13).",
        "Man teilt x von rechts in Dreierblöcke y₀, y₁, …, y_m.",
        "x = yₘ · 1000ᵐ + … + y₁ · 1000 + y₀",
        "  ≡ (−1)ᵐ yₘ + … − y₁ + y₀ (mod 7)",
        _annahme_text(r.variante),
        "",
        f"x = {blocktext}   (Dreierblöcke, m = {m})",
        f"x ≡ {vorzeichen_summe(r.vorzeichen, r.teile)}",
    ]
    if len(r.terme) > 1:
        zeilen.append(f"  = {z(r.summe)}")
    for p in r.module:
        zeilen.append(f"  {p}: {_rest_text(r.summe, p)}"
                      + (f"   → Summe ≡ {r.reste_summe[p]} (mod {p})"))
    if r.faktor == -1:
        zeilen.append("Achtung: Diese Summe ist ≡ −x. Der Rest von x ist das Negative.")
    zeilen.append("")
    teile = []
    for p in r.module:
        if r.teilbar[p]:
            teile.append(f"durch {p} teilbar")
        else:
            teile.append(f"nicht durch {p} teilbar (x ≡ {r.reste_x[p]} mod {p})")
    haupt = teile[0]
    zeilen.append(f"Ergebnis: x = {zahl_mit_leerzeichen(r.x)} ist {haupt}.")
    if len(teile) > 1:
        zeilen.append("Nebenbei (gleiche Summe): x ist " + "; ".join(teile[1:]) + ".")
    zeilen.append("Probe: " + ", ".join(f"{r.x} mod {p} = {r.x % p}" for p in r.module) + ".")
    return "\n".join(zeilen)


def text_teiler(t: TeilerDaten) -> str:
    zeilen = [
        "Teileranzahl d(n) = Σ_{d|n} 1 und Teilersumme σ(n) = Σ_{d|n} d.",
        "",
        f"n = {t.n} = {zerlegung_text(t.faktoren)}" if t.n > 1 else "n = 1 (keine Primfaktoren)",
        f"Teiler von {t.n}: {', '.join(map(str, t.teiler))}",
        "",
        f"d({t.n}) = Anzahl der Teiler = {t.d}",
    ]
    if t.faktoren:
        zeilen.append("Formel: d(n) = (1 + e₁) · … · (1 + e_r)")
        zeilen.append(f"d({t.n}) = " + " · ".join(f"(1 + {e})" for e in t.faktoren.values())
                      + f" = {' · '.join(map(str, t.d_formel_faktoren))} = {t.d}")
    summe = " + ".join(map(str, t.teiler))
    zeilen += ["", f"σ({t.n}) = {summe} = {t.sigma}" if len(t.teiler) > 1 else f"σ({t.n}) = {t.sigma}"]
    if t.faktoren:
        zeilen.append("Formel: σ(n) = (p₁^(1+e₁) − 1)/(p₁ − 1) · … · (p_r^(1+e_r) − 1)/(p_r − 1)")
        brueche = " · ".join(f"({p}{hoch(1 + e)} − 1)/({p} − 1)" for p, e in t.faktoren.items())
        zeilen.append(f"σ({t.n}) = {brueche}")
        zeilen.append(f"       = {' · '.join(map(str, t.sigma_formel_faktoren))} = {t.sigma}")
    zeilen.append("")
    zeilen.append(f"Ergebnis: d({t.n}) = {t.d}, σ({t.n}) = {t.sigma}.")
    zeilen.append(f"Probe: Liste und Formel geben dieselben Werte ({len(t.teiler)} Teiler, "
                  f"Summe {sum(t.teiler)}).")
    return "\n".join(zeilen)


def text_vollkommen(t: TeilerDaten) -> str:
    zeilen = [
        "Eine Zahl n heißt vollkommen, wenn σ(n) = 2n gilt.",
        f"n = {t.n} = {zerlegung_text(t.faktoren)}",
        f"σ({t.n}) = {' + '.join(map(str, t.teiler))} = {t.sigma}",
        f"2 · {t.n} = {2 * t.n}",
    ]
    if t.art == "vollkommen":
        zeilen.append(f"σ({t.n}) = {t.sigma} = 2 · {t.n}")
        zeilen.append(f"Ergebnis: {t.n} ist vollkommen.")
    else:
        vergleich = "<" if t.art == "defizient" else ">"
        zeilen.append(f"σ({t.n}) = {t.sigma} {vergleich} {2 * t.n} = 2 · {t.n}")
        zeilen.append(f"Ergebnis: {t.n} ist nicht vollkommen "
                      f"({'zu kleine' if t.art == 'defizient' else 'zu große'} Teilersumme).")
    return "\n".join(zeilen)


def text_mersenne(md: MersenneDaten) -> str:
    k = md.k
    zeilen = [f"k = {k}: 2{hoch(k)} − 1 = {md.m}"]
    if md.m == 1:
        zeilen.append("Ergebnis: 1 ist keine Primzahl.")
        return "\n".join(zeilen)
    if md.k_zerlegung:
        a, b = md.k_zerlegung
        q = 2 ** a
        zeilen.append(f"  k = {a} · {b} ist keine Primzahl. Setze q = 2{hoch(a)} = {q}:")
        zeilen.append(f"  q{hoch(b)} − 1 = (q − 1)(1 + q + … + q{hoch(b - 1)}), "
                      f"also teilt {q - 1} die Zahl {md.m}.")
    if md.ist_prim:
        n = md.vollkommen
        zeilen.append(f"  {md.m} ist eine Primzahl, also eine Mersenne-Primzahl.")
        zeilen.append(f"  n = 2{hoch(k - 1)} · (2{hoch(k)} − 1) = {2 ** (k - 1)} · {md.m} = {n}")
        zeilen.append(f"  σ(n) = σ(2{hoch(k - 1)}) · σ({md.m}) = (2{hoch(k)} − 1) · 2{hoch(k)} "
                      f"= {md.m} · {2 ** k} = {md.m * 2 ** k} = 2 · {n}")
        zeilen.append(f"Ergebnis: 2{hoch(k)} − 1 = {md.m} ist prim; {n} ist vollkommen.")
    else:
        zeilen.append(f"  {md.m} = {zerlegung_text(md.faktoren)}")
        zeilen.append(f"Ergebnis: 2{hoch(k)} − 1 = {md.m} = {zerlegung_text(md.faktoren)} "
                      f"ist keine Primzahl.")
        zeilen.append(f"Probe: {zerlegung_text(md.faktoren)} = {md.m}.")
    return "\n".join(zeilen)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

BEISPIELE = """Beispiele:
  python -m skripte.teilbarkeit regeln 299792  (Regeln für 9, 11 und 7)
  python -m skripte.teilbarkeit sieben "10^19 + 1"  (Regel für 7 mit Ausdruck)
  python -m skripte.teilbarkeit teiler 28  (Teileranzahl und Teilersumme)
  python -m skripte.teilbarkeit vollkommen 6 28 496  (drei vollkommene Zahlen)
  python -m skripte.teilbarkeit mersenne 2 3 5 11  (vier Mersenne-Zahlen)
"""


def _varianten_liste(args: argparse.Namespace) -> list[str]:
    return list(VARIANTEN) if args.alle_varianten else [args.variante]


def _zahl_pruefen(x: int, hinweise: list[str]) -> int:
    if x < 0:
        hinweise.append(f"Hinweis: x = {x} ist negativ. Ich rechne mit |x| = {-x}.")
        return -x
    return x


def _mit_varianten(args, funktion) -> str:
    teile = []
    varianten = _varianten_liste(args)
    for v in varianten:
        if len(varianten) > 1:
            teile.append(f"=== Variante „{v}“ ===")
        teile.append(funktion(v))
        teile.append("")
    return "\n".join(teile).rstrip()


def ausfuehren(args: argparse.Namespace) -> str:
    hinweise: list[str] = []
    ausgabe: list[str] = []
    befehl = args.befehl

    if befehl in ("neun", "elf", "sieben", "regeln"):
        x = _zahl_pruefen(zahl_lesen(args.x), hinweise)
        if befehl in ("neun", "regeln"):
            ausgabe.append(text_9(regel_9(x)))
        if befehl in ("elf", "regeln"):
            ausgabe.append(_mit_varianten(args, lambda v: text_11(regel_11(x, v))))
        if befehl in ("sieben", "regeln"):
            ausgabe.append(_mit_varianten(args, lambda v: text_1001(regel_1001(x, v))))
    elif befehl == "teiler":
        for text in args.n:
            n = zahl_lesen(text)
            if n < 1:
                hinweise.append(f"Hinweis: n = {n} ist nicht positiv. d(n) und σ(n) gibt es nur für n ≥ 1.")
                continue
            ausgabe.append(text_teiler(teiler_daten(n)))
    elif befehl == "vollkommen":
        for text in args.n:
            n = zahl_lesen(text)
            if n < 1:
                hinweise.append(f"Hinweis: n = {n} ist nicht positiv. Übersprungen.")
                continue
            ausgabe.append(text_vollkommen(teiler_daten(n)))
    elif befehl == "mersenne":
        ks = [zahl_lesen(t) for t in args.k]
        if args.bis:
            ks += list(range(1, args.bis + 1))
        if not ks:
            hinweise.append("Hinweis: Bitte k angeben oder --bis K nutzen.")
        for k in ks:
            if k < 1:
                hinweise.append(f"Hinweis: k = {k} ist zu klein. k muss mindestens 1 sein.")
                continue
            if k > 256:
                hinweise.append(f"Hinweis: k = {k} ist sehr groß. Die Zerlegung kann lange dauern.")
            ausgabe.append(text_mersenne(mersenne_daten(k)))
        prim = [k for k in ks if k >= 1 and isprime(2 ** k - 1)]
        if len(ks) > 1:
            ausgabe.append("Zusammenfassung: 2ᵏ − 1 ist prim für k = "
                           + (", ".join(map(str, prim)) if prim else "keines") + ".")
    return "\n\n".join(hinweise + ausgabe)


def parser_bauen() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m skripte.teilbarkeit",
        description="Teilbarkeitsregeln für 9, 11 und 7 und Teilerfunktionen d(n), σ(n). "
                    "Das Skript zeigt jeden Rechenschritt.",
        epilog=BEISPIELE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    unter = parser.add_subparsers(dest="befehl", required=True, metavar="BEFEHL",
                                  help="Art der Rechnung. Hilfe zu einem Befehl: BEFEHL --help")

    def varianten_optionen(p: argparse.ArgumentParser) -> None:
        p.add_argument("--variante", choices=list(VARIANTEN), default=STANDARD_VARIANTE,
                       help="Wo beginnt das Vorzeichen „+“? einer = rechts bei der "
                            "Einerstelle (Standard). links = bei der ersten Ziffer links.")
        p.add_argument("--alle-varianten", action="store_true",
                       help="Gibt die Rechnung für alle Varianten nacheinander aus.")

    hilfe_x = "Die Zahl x. Auch als Ausdruck möglich, z. B. \"10^19 + 1\"."
    p = unter.add_parser("neun", help="Regel für 9 mit der Quersumme.",
                         epilog="Beispiele:\n  python -m skripte.teilbarkeit neun 299792  (Quersumme)",
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("x", help=hilfe_x)
    p = unter.add_parser("elf", help="Regel für 11 mit der alternierenden Quersumme.",
                         epilog="Beispiele:\n  python -m skripte.teilbarkeit elf 299792  (alternierende Quersumme)",
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("x", help=hilfe_x)
    varianten_optionen(p)
    p = unter.add_parser("sieben", help="Regel für 7 (und 11, 13) mit Dreierblöcken und 1001.",
                         epilog="Beispiele:\n  python -m skripte.teilbarkeit sieben \"10^19 + 1\"  (Dreierblöcke)",
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("x", help=hilfe_x)
    varianten_optionen(p)
    p = unter.add_parser("regeln", help="Alle drei Regeln (9, 11, 7) für eine Zahl.",
                         epilog="Beispiele:\n  python -m skripte.teilbarkeit regeln 299792  (alle drei Regeln)",
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("x", help=hilfe_x)
    varianten_optionen(p)
    p = unter.add_parser("teiler", help="Teileranzahl d(n) und Teilersumme σ(n) mit Formel.",
                         epilog="Beispiele:\n  python -m skripte.teilbarkeit teiler 28  (d(28) und σ(28))",
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("n", nargs="+", help="Eine oder mehrere Zahlen n ≥ 1.")
    p = unter.add_parser("vollkommen", help="Prüft σ(n) = 2n (vollkommene Zahl).",
                         epilog="Beispiele:\n  python -m skripte.teilbarkeit vollkommen 6 28 496  (drei Zahlen prüfen)",
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("n", nargs="+", help="Eine oder mehrere Zahlen n ≥ 1.")
    p = unter.add_parser("mersenne", help="Prüft, ob 2ᵏ − 1 prim ist, und bildet 2ᵏ⁻¹(2ᵏ − 1).",
                         epilog="Beispiele:\n  python -m skripte.teilbarkeit mersenne 2 3 5 11  (vier Exponenten)",
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("k", nargs="*", help="Ein oder mehrere Exponenten k ≥ 1.")
    p.add_argument("--bis", type=int, metavar="K",
                   help="Prüft alle k von 1 bis K.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = parser_bauen()
    args = parser.parse_args(argv)
    try:
        print(ausfuehren(args))
    except ValueError as fehler:
        print(f"Fehler: {fehler}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
