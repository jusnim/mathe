"""Chinesischer Restsatz: simultane Kongruenzen mit Lösungsweg.

Zweck
-----
Das Skript löst Systeme simultaner Kongruenzen.
Simultan heißt: Eine Zahl z soll alle Kongruenzen gleichzeitig erfüllen.
Beispiel: z ≡ 8 (mod 17) und z ≡ 5 (mod 15).
Das Skript zeigt jeden Rechenschritt wie in den Musterlösungen.
Sie können den Lösungsweg direkt auf Papier abschreiben.

Lösungsweg (VL1, Beweis „Bestimmung von f⁻¹“)
---------------------------------------------
Aus z ≡ a (mod m) und z ≡ b (mod n) folgt der Ansatz
    z = m·x + a = n·y + b   ⇒   m·x − n·y = b − a   (Lemma von Bezout).
Der Euklidische Algorithmus (EA) liefert x und y.
Die Gleichung ist genau dann lösbar, wenn ggT(m, n) die Zahl b − a teilt.
Alle Lösungen sind z = z₀ + k·kgV(m, n) mit k ∈ ℤ.
Bei teilerfremden Moduln ist kgV(m, n) = m·n.
Bei mehr als zwei Kongruenzen fasst das Skript schrittweise zusammen:
erst die ersten beiden, dann das Ergebnis mit der dritten, und so weiter (wie VL1).

Aufruf (Beispiele aus den Quellen)
----------------------------------
    python -m skripte.crt -k 5 13 -k 4 15                      (KV A8: −86 ≡ 109 mod 195)
    python -m skripte.crt -k 8 17 -k 5 15 --vorgabe 8 9        (K A4 mit Ergebnis aus K A3)
    python -m skripte.crt -k 3 5 -k 4 7                        (Ü3 A1 a: 18 mod 35)
    python -m skripte.crt -k -1 12 -k 2 5                      (Ü3 A1 b: −73 ≡ 47 mod 60)
    python -m skripte.crt -k 3 11 -k 6 8 -k 1 15               (VL1: −74 ≡ 1246 mod 1320)
Dabei bedeutet „-k A M“ die Kongruenz z ≡ A (mod M).

Varianten
---------
--variante (welche Kongruenz ist „m, a“):
- erste (Standard): m, a kommen aus der ersten Kongruenz (VL1 Z. 595, KV A8, Ü3 A1).
- zweite: m, a kommen aus der zweiten Kongruenz.
--weg (wie der EA die Bezout-Lösung liefert):
- rueckwaerts (Standard): EA rückwärts einsetzen (KV A8, Ü3 A1 Variante 1).
- matrix: Matrix Q (Ü3 A1 Variante 2).
--vorgabe X Y: Eine Lösung von m·x − n·y = ggT(m, n) ist schon bekannt
  (K A4: „mit Hilfe des Ergebnisses aus Aufgabe 3“).
--alle-varianten gibt alle Kombinationen nacheinander aus.
Alle Varianten sind richtig. Sie geben oft ein anderes z₀.
Alle z₀ sind aber kongruent modulo kgV. Die kleinste positive Lösung ist gleich.

Nutzung als Modul
-----------------
    from skripte.crt import loese_system
    erg = loese_system([(5, 13), (4, 15)])
    erg.z0, erg.modul, erg.kleinste_positive    # −86, 195, 109
    erg.text                                    # Lösungsweg
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field

from skripte.euklid import (
    EAErgebnis, MatrixErgebnis, RueckwaertsErgebnis,
    ea_matrix, ea_rueckwaerts, euklid, repraesentant, text_ea,
)

__all__ = [
    "PaarErgebnis", "CRTErgebnis", "loese_paar", "loese_system",
    "text_paar", "text_system", "kleinste_positive", "main",
    "VARIANTEN", "WEGE",
]

VARIANTEN = ("erste", "zweite")
VARIANTEN_TEXT = {
    "erste": "m, a aus der ersten Kongruenz (wie VL1 Z. 595, KV A8, Ü3 A1)",
    "zweite": "m, a aus der zweiten Kongruenz (Reihenfolge getauscht)",
}
WEGE = ("rueckwaerts", "matrix")
WEGE_TEXT = {
    "rueckwaerts": "EA rückwärts einsetzen",
    "matrix": "Matrix Q",
    "vorgabe": "vorgegebene Lösung von m·x − n·y = ggT(m, n)",
}


# ---------------------------------------------------------------------------
# Hilfsfunktionen für die Schreibweise
# ---------------------------------------------------------------------------

def _k(n: int) -> str:
    """Negative Zahl in Klammern: −5 → '(−5)'."""
    return str(n) if n >= 0 else f"({n})"


def _plus(z: int) -> str:
    """„+ 5“ oder „− 5“."""
    return f"+ {z}" if z >= 0 else f"− {-z}"


def _mx_ny(m: int, x: int, n: int, y: int) -> str:
    """Schreibt m·x − n·y mit eingesetzten Zahlen."""
    return f"{m}·{_k(x)} − {n}·{_k(y)}"


def kleinste_positive(z: int, modul: int) -> int:
    """Kleinste positive Zahl, die ≡ z (mod modul) ist (Wert in 1 … modul)."""
    return (z - 1) % modul + 1


# ---------------------------------------------------------------------------
# Zwei Kongruenzen zusammenfassen
# ---------------------------------------------------------------------------

@dataclass
class PaarErgebnis:
    """Lösung von z ≡ a (mod m), z ≡ b (mod n) über m·x − n·y = b − a.

    start: Zahl, bei der die Bezout-Rechnung ankommt (ggT oder |b − a|).
    x_s, y_s: m·x_s − n·y_s = start.
    faktor: (b − a) / start.  x = x_s·faktor, y = y_s·faktor.
    z0 = m·x + a = n·y + b.  Alle Lösungen: z0 + k·modul mit modul = kgV(m, n).
    """
    a: int
    m: int
    b: int
    n: int
    d: int                      # b − a
    ggt: int
    loesbar: bool
    modul: int                  # kgV(m, n)
    weg: str                    # "rueckwaerts", "matrix" oder "vorgabe"
    ea: EAErgebnis
    rueckwaerts: RueckwaertsErgebnis | None = None
    matrix: MatrixErgebnis | None = None
    start: int | None = None
    x_s: int | None = None
    y_s: int | None = None
    faktor: int | None = None
    x: int | None = None
    y: int | None = None
    z0: int | None = None
    hinweise: list[str] = field(default_factory=list)
    text: str = ""


def _auf_m_n(ea: EAErgebnis, kg: int, kk: int) -> tuple[int, int]:
    """Koeffizienten für (größere, kleinere Zahl) → Koeffizienten für (m, n)."""
    return (kk, kg) if ea.getauscht else (kg, kk)


def loese_paar(a: int, m: int, b: int, n: int, weg: str = "rueckwaerts",
               vorgabe: tuple[int, int] | None = None) -> PaarErgebnis:
    """Fasst z ≡ a (mod m) und z ≡ b (mod n) zu einer Kongruenz zusammen.

    weg = "rueckwaerts" oder "matrix".
    vorgabe = (x, y) mit m·x − n·y = ggT(m, n). Dann entfällt die Bezout-Rechnung.
    Die Moduln m, n müssen ≥ 1 sein. Sie dürfen einen gemeinsamen Teiler haben.
    """
    if m < 1 or n < 1:
        raise ValueError(f"Die Moduln müssen ≥ 1 sein (hier m = {m}, n = {n}).")
    if weg not in WEGE:
        raise ValueError(f"Unbekannter Weg {weg!r}. Erlaubt: {', '.join(WEGE)}.")
    ea = euklid(m, n)
    g = ea.ggt
    d = b - a
    erg = PaarErgebnis(a, m, b, n, d, g, d % g == 0, m * n // g, weg, ea)
    if not erg.loesbar:
        erg.text = text_paar(erg)
        return erg

    if vorgabe is not None:
        vx, vy = vorgabe
        if m * vx - n * vy == g:
            erg.weg = "vorgabe"
            erg.start, erg.x_s, erg.y_s = g, vx, vy
        else:
            erg.hinweise.append(
                f"Hinweis: Die Vorgabe (x, y) = ({vx}, {vy}) erfüllt m·x − n·y = ggT nicht: "
                f"{_mx_ny(m, vx, n, vy)} = {m * vx - n * vy} ≠ {g}. "
                f"Das Skript rechnet deshalb mit dem EA.")
    if erg.start is None:
        if weg == "matrix":
            mx = ea_matrix(m, n)
            erg.matrix = mx
            cm, cn = _auf_m_n(ea, mx.koef_gross, mx.koef_klein)
            erg.start = g
        else:
            # Ist |b − a| selbst ein Rest im EA, startet das Einsetzen dort (wie VL1).
            rw = ea_rueckwaerts(m, n, ziel=abs(d) if d != 0 else None)
            erg.rueckwaerts = rw
            cm, cn = _auf_m_n(ea, rw.koef_gross, rw.koef_klein)
            erg.start = rw.ziel
        # start = cm·m + cn·n = m·cm − n·(−cn)
        erg.x_s, erg.y_s = cm, -cn
    erg.faktor = d // erg.start
    erg.x, erg.y = erg.x_s * erg.faktor, erg.y_s * erg.faktor
    erg.z0 = m * erg.x + a
    assert erg.z0 == n * erg.y + b
    erg.text = text_paar(erg)
    return erg


def text_paar(erg: PaarErgebnis, nr: int | None = None) -> str:
    """Lösungsweg für zwei Kongruenzen."""
    a, m, b, n, d, g = erg.a, erg.m, erg.b, erg.n, erg.d, erg.ggt
    z: list[str] = []
    kopf = f"Schritt {nr}: " if nr is not None else ""
    z.append(f"{kopf}z ≡ {a} (mod {m}),  z ≡ {b} (mod {n})")
    z.append(f"  m = {m}, a = {a};   n = {n}, b = {b}")
    z.append("  Ansatz:  z = m·x + a = n·y + b")
    z.append("  ⇒  m·x − n·y = b − a        (Lemma von Bezout)")
    z.append(f"  ⇒  {m}·x − {n}·y = {b} − {_k(a)} = {d}")
    z.append("")
    z.append(text_ea(erg.ea))
    z.append("")
    z.append("Lösbarkeit:")
    if g == 1:
        z.append(f"  ggT({m}, {n}) = 1. Die Moduln sind teilerfremd.")
        z.append(f"  1 teilt b − a = {d}. Das System ist lösbar (Chinesischer Restsatz).")
    else:
        z.append(f"  ggT({m}, {n}) = {g} ≠ 1. Die Moduln sind nicht teilerfremd.")
        z.append("  Der Chinesische Restsatz gilt hier nicht direkt.")
        z.append(f"  Das System ist genau dann lösbar, wenn {g} die Zahl b − a = {d} teilt.")
        if not erg.loesbar:
            z.append(f"  {d} = {d // g}·{g} + {d % g}. Also teilt {g} die Zahl {d} nicht.")
            z.append(f"  Grund: Jede Lösung erfüllt z ≡ a und z ≡ b (mod {g}).")
            z.append(f"  Aber {a} mod {g} = {a % g} und {b} mod {g} = {b % g} sind verschieden.")
            z.append(f"  ⇒ Das System z ≡ {a} (mod {m}), z ≡ {b} (mod {n}) hat keine Lösung.")
            return "\n".join(z)
        z.append(f"  {d} = {d // g}·{g}. Also teilt {g} die Zahl {d}. Das System ist lösbar.")
    z.append("")
    z.extend(erg.hinweise)

    if d == 0:
        z.append("  b − a = 0. Also passt x = 0, y = 0.")
    elif erg.weg == "vorgabe":
        z.append(f"Bezout-Lösung aus der Vorgabe (vorherige Aufgabe):")
        z.append(f"  {_mx_ny(m, erg.x_s, n, erg.y_s)} = {m * erg.x_s} − {_k(n * erg.y_s)} = {g}")
    elif erg.weg == "matrix":
        mx = erg.matrix
        z.append(mx.text)
        z.append(f"  ⇒ {_mx_ny(m, erg.x_s, n, erg.y_s)} = {g}")
    else:
        rw = erg.rueckwaerts
        if rw.ziel != g:
            z.append(f"  |b − a| = {rw.ziel} ist selbst ein Rest im EA (r_{rw.k}).")
            z.append("  Das Einsetzen startet deshalb dort.")
        z.append(rw.text)
        z.append(f"  ⇒ {erg.start} = {_mx_ny(m, erg.x_s, n, erg.y_s)}")
    if d != 0:
        if erg.faktor != 1:
            z.append(f"  Mit b − a = {d} = {_k(erg.faktor)}·{erg.start} multiplizieren:   | ·{_k(erg.faktor)}")
        z.append(f"  ⇒ b − a = {d} = {_mx_ny(m, erg.x, n, erg.y)}")
        z.append(f"  Also x = {erg.x}, y = {erg.y}.")
    z.append(f"  ⇒ z = m·x + a = {m}·{_k(erg.x)} {_plus(a)} = {erg.z0}")
    z.append(f"      = n·y + b = {n}·{_k(erg.y)} {_plus(b)} = {erg.z0}")
    if g == 1:
        z.append(f"  Alle Lösungen: z ≡ {erg.z0} (mod m·n = {m}·{n} = {erg.modul}),")
    else:
        z.append(f"  kgV({m}, {n}) = {m}·{n} / {g} = {erg.modul}")
        z.append(f"  Alle Lösungen: z ≡ {erg.z0} (mod kgV = {erg.modul}),")
    z.append(f"  also z = {erg.z0} + k·{erg.modul} mit k ∈ ℤ.")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# System mit beliebig vielen Kongruenzen
# ---------------------------------------------------------------------------

@dataclass
class CRTErgebnis:
    """Lösung des Systems z ≡ a_i (mod m_i).

    Bei Lösbarkeit: alle Lösungen z = z0 + k·modul.
    """
    kongruenzen: list[tuple[int, int]]      # Paare (a_i, m_i)
    variante: str
    weg: str
    rep: str
    schritte: list[PaarErgebnis]
    loesbar: bool
    z0: int | None = None
    modul: int | None = None
    kleinste_positive: int | None = None
    repraesentant: int | None = None
    text: str = ""


def loese_system(kongruenzen: list[tuple[int, int]], variante: str = "erste",
                 weg: str = "rueckwaerts", vorgabe: tuple[int, int] | None = None,
                 rep: str = "standard") -> CRTErgebnis:
    """Löst z ≡ a_i (mod m_i) für alle Paare (a_i, m_i) in „kongruenzen“.

    variante = "erste": m, a stammen aus der ersten (bzw. bisher zusammengefassten) Kongruenz.
    variante = "zweite": m, a stammen aus der neu hinzukommenden Kongruenz.
    weg = "rueckwaerts" oder "matrix".
    vorgabe: Bezout-Lösung (x, y) mit m·x − n·y = ggT(m, n). Sie gilt nur im ersten Schritt.
    rep = "standard" (0 … M−1) oder "symmetrisch" (−M/2 … M/2) für den Repräsentanten.
    """
    if variante not in VARIANTEN:
        raise ValueError(f"Unbekannte Variante {variante!r}. Erlaubt: {', '.join(VARIANTEN)}.")
    if rep not in ("standard", "symmetrisch"):
        raise ValueError(f"Unbekannter Repräsentant {rep!r}. Erlaubt: standard, symmetrisch.")
    kongruenzen = [(int(a), int(m)) for a, m in kongruenzen]
    if len(kongruenzen) < 2:
        raise ValueError("Bitte mindestens zwei Kongruenzen angeben.")
    for a, m in kongruenzen:
        if m < 1:
            raise ValueError(f"Der Modul in z ≡ {a} (mod {m}) muss ≥ 1 sein.")
    erg = CRTErgebnis(kongruenzen, variante, weg, rep, [], True)
    akt_a, akt_m = kongruenzen[0]
    for i, (b, n) in enumerate(kongruenzen[1:]):
        vg = vorgabe if i == 0 else None
        if variante == "erste":
            p = loese_paar(akt_a, akt_m, b, n, weg, vg)
        else:
            p = loese_paar(b, n, akt_a, akt_m, weg, vg)
        erg.schritte.append(p)
        if not p.loesbar:
            erg.loesbar = False
            break
        akt_a, akt_m = p.z0, p.modul
    if erg.loesbar:
        erg.z0, erg.modul = akt_a, akt_m
        erg.kleinste_positive = kleinste_positive(akt_a, akt_m)
        erg.repraesentant = repraesentant(akt_a, akt_m, rep)
    erg.text = text_system(erg)
    return erg


def _probe_zeile(z: int, a: int, m: int) -> str:
    q, r = divmod(z, m)
    ok = (z - a) % m == 0
    zeichen = "✓" if ok else "FEHLER"
    if 0 <= a < m:
        return f"  {z} = {q}·{m} + {r}  ≡ {r} (mod {m})  {zeichen}"
    diff = z - a
    return f"  {z} − {_k(a)} = {diff} = {diff // m}·{m}  ⇒  {z} ≡ {a} (mod {m})  {zeichen}"


def text_system(erg: CRTErgebnis) -> str:
    """Lösungsweg für das ganze System."""
    z: list[str] = ["Gegeben: System simultaner Kongruenzen"]
    for a, m in erg.kongruenzen:
        z.append(f"  z ≡ {a} (mod {m})")
    z.append("")
    mehrere = len(erg.kongruenzen) > 2
    if mehrere:
        z.append("Mehr als zwei Kongruenzen: Wir fassen schrittweise zusammen.")
        z.append("Erst die ersten beiden. Dann das Ergebnis mit der nächsten Kongruenz.")
        z.append("")
    for i, p in enumerate(erg.schritte, start=1):
        z.append(text_paar(p, i if mehrere else None))
        z.append("")
    if not erg.loesbar:
        p = erg.schritte[-1]
        z.append(f"Ergebnis: Das System hat keine Lösung, weil ggT({p.m}, {p.n}) = {p.ggt} "
                 f"die Zahl b − a = {p.d} nicht teilt.")
        return "\n".join(z)
    z0, M, kp = erg.z0, erg.modul, erg.kleinste_positive
    z.append("Alle Lösungen:")
    z.append(f"  z = {z0} + k·{M}  (k ∈ ℤ),  also  z ≡ {z0} (mod {M})")
    if kp == z0:
        z.append(f"Kleinste positive Lösung: z = {kp}")
    else:
        k = (kp - z0) // M
        schritt = f"{z0} {_plus(k * M)}" if abs(k) == 1 else f"{z0} {'+' if k > 0 else '−'} {abs(k)}·{M}"
        z.append(f"Kleinste positive Lösung: z = {schritt} = {kp}")
    if erg.rep == "symmetrisch":
        z.append(f"Symmetrischer Repräsentant (−{M}/2 … {M}/2): z ≡ {erg.repraesentant} (mod {M})")
    z.append("")
    z.append(f"Probe mit z = {kp}:")
    for a, m in erg.kongruenzen:
        z.append(_probe_zeile(kp, a, m))
    z.append("")
    rep_txt = f" ≡ {erg.repraesentant}" if erg.rep == "symmetrisch" and erg.repraesentant != kp else ""
    z.append(f"Ergebnis: z ≡ {z0}{'' if kp == z0 else f' ≡ {kp}'}{rep_txt} (mod {M}); "
             f"alle Lösungen z = {kp} + k·{M} (k ∈ ℤ); kleinste positive Lösung z = {kp}.")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

BEISPIELE = """Beispiele aus den Quellen („-k A M“ heißt z ≡ A (mod M)):
  KV A8:   python -m skripte.crt -k 5 13 -k 4 15               (−86 ≡ 109 mod 195)
  K A4:    python -m skripte.crt -k 8 17 -k 5 15 --vorgabe 8 9 (−400 ≡ 110 mod 255)
           (Vorgabe aus K A3: 17·8 − 15·9 = 1)
  Ü3 A1 a: python -m skripte.crt -k 3 5 -k 4 7 --weg matrix    (18 mod 35)
  Ü3 A1 b: python -m skripte.crt -k -1 12 -k 2 5               (−73 ≡ 47 mod 60)
  VL1:     python -m skripte.crt -k 3 11 -k 6 8 -k 1 15        (−74 ≡ 1246 mod 1320)
  Nicht teilerfremd: python -m skripte.crt -k 2 4 -k 4 6      (10 mod 12)
"""


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m skripte.crt",
        description="Chinesischer Restsatz: löst simultane Kongruenzen z ≡ a (mod m) mit Lösungsweg. "
                    "Auch mehr als zwei Kongruenzen und nicht teilerfremde Moduln.",
        epilog=BEISPIELE, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-k", "--kongruenz", nargs=2, type=int, action="append", metavar=("A", "M"),
                   required=True,
                   help="Eine Kongruenz z ≡ A (mod M). Mehrmals angeben, "
                        "zum Beispiel -k 5 13 -k 4 15.")
    p.add_argument("--variante", choices=VARIANTEN, default="erste",
                   help="Welche Kongruenz ist „m, a“? erste (Standard, wie VL1 und KV A8) "
                        "oder zweite.")
    p.add_argument("--weg", choices=WEGE, default="rueckwaerts",
                   help="Wie der EA die Bezout-Lösung liefert: rueckwaerts (Standard, EA rückwärts) "
                        "oder matrix (Matrix Q).")
    p.add_argument("--vorgabe", nargs=2, type=int, metavar=("X", "Y"),
                   help="Bekannte Lösung von m·x − n·y = ggT(m, n) für den ersten Schritt, "
                        "zum Beispiel aus einer vorherigen Aufgabe (K A4: --vorgabe 8 9).")
    p.add_argument("--alle-varianten", action="store_true",
                   help="Gibt alle Varianten (Reihenfolge und Weg) nacheinander aus.")
    p.add_argument("--rep", choices=("standard", "symmetrisch"), default="standard",
                   help="Repräsentant des Ergebnisses: standard (0 … M−1) "
                        "oder symmetrisch (−M/2 … M/2). Die kleinste positive Lösung steht immer da.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    kongr = [tuple(k) for k in args.kongruenz]
    vorgabe = tuple(args.vorgabe) if args.vorgabe else None
    if args.alle_varianten:
        laeufe = [(v, w, None) for v in VARIANTEN for w in WEGE]
        if vorgabe is not None:
            laeufe.insert(0, ("erste", "rueckwaerts", vorgabe))
    else:
        laeufe = [(args.variante, args.weg, vorgabe)]
    aus: list[str] = []
    ergebnisse = []
    try:
        for i, (v, w, vg) in enumerate(laeufe):
            weg_txt = WEGE_TEXT["vorgabe"] if vg is not None else WEGE_TEXT[w]
            if len(laeufe) > 1:
                aus.append(f"\n=== Variante {i + 1}: {v}, {weg_txt} ===")
            aus.append(f"Annahme: {VARIANTEN_TEXT[v]}; Bezout-Lösung: {weg_txt}"
                       + ("" if vg is None else f" (x, y) = {vg}") + ".")
            e = loese_system(kongr, v, w, vg, args.rep)
            ergebnisse.append(e)
            aus.append(e.text)
    except ValueError as f:
        print(f"Hinweis: {f}")
        return 2
    geloest = [e for e in ergebnisse if e.loesbar]
    if len(ergebnisse) > 1 and geloest:
        z0s = sorted({e.z0 for e in geloest})
        aus.append("")
        aus.append(f"Vergleich: Die Varianten geben z₀ = {', '.join(map(str, z0s))}. "
                   f"Alle sind ≡ {geloest[0].kleinste_positive} (mod {geloest[0].modul}). "
                   f"Alle Varianten sind richtig.")
    print("\n".join(aus).lstrip("\n"))
    return 0 if all(e.loesbar for e in ergebnisse) else 1


if __name__ == "__main__":
    sys.exit(main())
