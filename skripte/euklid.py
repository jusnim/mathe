"""Euklidischer Algorithmus (EA), ggT, kgV, Bezout, Inverse, diophantische Gleichung.

Zweck
-----
Das Skript rechnet Aufgaben zum Euklidischen Algorithmus.
Es zeigt jeden Rechenschritt.
Sie können den Lösungsweg direkt übernehmen.

Begriffe
--------
- ggT: größter gemeinsamer Teiler.
- kgV: kleinstes gemeinsames Vielfaches.
- EA: wiederholte Division mit Rest. Der letzte Rest ungleich 0 ist der ggT.
- Bezout: Es gibt ganze Zahlen x, y mit a·x + b·y = ggT(a, b).
- Diophantische Gleichung: Gleichung, deren Lösungen ganze Zahlen sein müssen.
- Matrix Q: Produkt der Matrizen Q_i = (q_i 1; 1 0) aus dem EA.

Aufruf (Beispiele)
------------------
    python -m skripte.euklid --a 2406 --b 654 --c 24 --form -     (Gleichung a·x − b·y = c)
    python -m skripte.euklid --a 1001 --b 840 --c 98 --form +     (Gleichung a·x + b·y = c)
    python -m skripte.euklid --a 17 --b 15 --c 1 --form -         (Gleichung mit c = 1)
    python -m skripte.euklid --a 1965 --b 225                     (nur ggT und Bezout)
    python -m skripte.euklid --a 3 --mod 10 --rep symmetrisch     (Inverse, symmetrisch)
    python -m skripte.euklid --a 12 --b 18 --kgv                  (kgV)

Varianten (--variante)
----------------------
- beide (Standard): EA rückwärts und Matrix Q.
- rueckwaerts: nur „EA rückwärts“ (Einsetzen von unten nach oben).
- matrix: nur Matrix Q = Q_0·…·Q_n.
Beide Wege sind richtig. Sie können aber verschiedene spezielle Lösungen geben.
Beispiel 2406·x − 654·y = 24: EA rückwärts gibt (3, 11), Matrix Q gibt (112, 412).
--alle-varianten gibt beide Wege nacheinander aus.

Weitere Optionen
----------------
- --form + oder -: Gleichung a·x + b·y = c oder a·x − b·y = c. Standard: −.
- --rep standard oder symmetrisch: Repräsentant der Inversen (0 … m−1 oder −m/2 … m/2).

Nutzung als Modul (für andere Skripte)
--------------------------------------
    from skripte.euklid import ggt, kgv, bezout, inverse, loese_diophantisch
    erg = loese_diophantisch(2406, 654, 24, form="-")
    erg.x, erg.y          # spezielle Lösung
    erg.text              # Lösungsweg als Text

Jede Rechenfunktion gibt ein Ergebnis-Objekt zurück.
Das Objekt hat die Werte als Felder und den Lösungsweg im Feld „text“.
Die Textfunktionen (text_…) erzeugen den Lösungsweg aus den Daten.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field

__all__ = [
    "EAZeile", "EAErgebnis", "ErwEAZeile", "RueckwaertsErgebnis", "MatrixErgebnis",
    "BezoutErgebnis", "InverseErgebnis", "DiophantErgebnis", "KgVErgebnis",
    "KeinInversesFehler",
    "euklid", "erweiterter_ea", "ggt", "ggt_mit_weg", "kgv", "kgv_mit_weg",
    "ea_rueckwaerts", "ea_matrix", "bezout", "bezout_mit_weg",
    "inverse", "inverse_mit_weg", "loese_diophantisch", "repraesentant",
    "text_ea", "text_ea_tabelle", "text_rueckwaerts", "text_matrix", "text_bezout",
    "text_inverse", "text_diophantisch", "text_kgv",
    "VARIANTEN",
]

VARIANTEN = ("beide", "rueckwaerts", "matrix")
VARIANTEN_TEXT = {
    "beide": "EA rückwärts und Matrix Q",
    "rueckwaerts": "EA rückwärts (Einsetzen von unten nach oben)",
    "matrix": "Matrix Q = Q_0·…·Q_n",
}


# ---------------------------------------------------------------------------
# Hilfsfunktionen für die Schreibweise
# ---------------------------------------------------------------------------

def _k(n: int) -> str:
    """Zahl in Klammern, wenn sie negativ ist: 5 → '5', −5 → '(−5)'."""
    return str(n) if n >= 0 else f"({n})"


def _gl(a: int, x: int, b: int, y: int, form: str) -> str:
    """Schreibt a·x ± b·y mit eingesetzten Zahlen."""
    return f"{a}·{_k(x)} {'+' if form == '+' else '−'} {b}·{_k(y)}"


def _term(koef: int, wert: str, erster: bool) -> str:
    """Ein Summand koef·wert mit Vorzeichen. Koeffizient 1 wird weggelassen."""
    betrag = abs(koef)
    kern = wert if betrag == 1 else f"{betrag}·{wert}"
    if erster:
        return kern if koef >= 0 else f"−{kern}"
    return f" + {kern}" if koef >= 0 else f" − {kern}"


def _summe(terme: list[tuple[int, str]]) -> str:
    return "".join(_term(k, w, i == 0) for i, (k, w) in enumerate(terme)) or "0"


def _lin(p: int, kp: int, q: int, kq: int) -> str:
    """Schreibt p·kp ± q·|kq|, zum Beispiel 2406·28 − 654·103."""
    return f"{p}·{_k(kp)} {'+' if kq >= 0 else '−'} {q}·{abs(kq)}"


def repraesentant(x: int, m: int, rep: str = "standard") -> int:
    """Repräsentant von x mod m.

    rep = "standard": Wert in 0 … m−1.
    rep = "symmetrisch": Wert in (−m/2, m/2], zum Beispiel −5 … 5 für m = 11.
    """
    r = x % m
    if rep == "symmetrisch" and r > m // 2:
        r -= m
    return r


# ---------------------------------------------------------------------------
# Euklidischer Algorithmus
# ---------------------------------------------------------------------------

@dataclass
class EAZeile:
    """Eine Zeile r_(i−1) = q_i · r_i + r_(i+1) des EA."""
    i: int
    dividend: int   # r_(i−1)
    q: int          # q_i
    divisor: int    # r_i
    rest: int       # r_(i+1)


@dataclass
class EAErgebnis:
    """Ergebnis des EA für zwei Zahlen.

    Der EA rechnet mit den Beträgen und mit der größeren Zahl zuerst.
    r_(−1) = groß, r_0 = klein.
    """
    a: int                      # Eingabe 1 (wie übergeben)
    b: int                      # Eingabe 2 (wie übergeben)
    gross: int                  # r_(−1)
    klein: int                  # r_0
    getauscht: bool             # True, wenn |a| < |b| und der EA mit (b, a) läuft
    zeilen: list[EAZeile]
    ggt: int
    text: str = ""

    @property
    def quotienten(self) -> list[int]:
        return [z.q for z in self.zeilen]

    @property
    def reste(self) -> list[int]:
        """r_(−1), r_0, r_1, …, r_n, 0. Index j in der Liste entspricht r_(j−1)."""
        if not self.zeilen:
            return [self.gross, self.klein]
        return [self.gross, self.klein] + [z.rest for z in self.zeilen]

    def r(self, i: int) -> int:
        """Rest r_i (i ≥ −1)."""
        return self.reste[i + 1]

    @property
    def n(self) -> int:
        """Index n des letzten Rests ungleich 0 (r_n = ggT)."""
        return len(self.zeilen) - 1


def euklid(a: int, b: int) -> EAErgebnis:
    """EA für a und b. Gibt alle Divisionszeilen und den ggT zurück.

    Mindestens eine der Zahlen muss ungleich 0 sein.
    Negative Zahlen sind erlaubt. Der EA rechnet dann mit den Beträgen.
    """
    if a == 0 and b == 0:
        raise ValueError("ggT(0, 0) ist nicht definiert. Mindestens eine Zahl muss ≠ 0 sein.")
    A, B = abs(a), abs(b)
    getauscht = A < B
    gross, klein = (B, A) if getauscht else (A, B)
    zeilen: list[EAZeile] = []
    x, y, i = gross, klein, 0
    while y != 0:
        q, r = divmod(x, y)
        zeilen.append(EAZeile(i, x, q, y, r))
        x, y, i = y, r, i + 1
    erg = EAErgebnis(a, b, gross, klein, getauscht, zeilen, x)
    erg.text = text_ea(erg)
    return erg


def text_ea(ea: EAErgebnis, titel: bool = True) -> str:
    """Lösungsweg des EA: Divisionszeilen und ggT."""
    z: list[str] = []
    if titel:
        z.append(f"Euklidischer Algorithmus (EA) für {ea.gross}, {ea.klein}")
        if ea.getauscht:
            z.append(f"(Die größere Zahl steht vorne. Deshalb rechnet der EA mit {ea.gross}, {ea.klein}.)")
        if ea.a < 0 or ea.b < 0:
            z.append("(Der EA rechnet mit den Beträgen. Der ggT ist immer positiv.)")
    if ea.klein == 0:
        z.append(f"Eine Zahl ist 0. Es gilt ggT({ea.gross}, 0) = {ea.gross}.")
        return "\n".join(z)
    breite = max(len(str(ea.gross)), 1)
    for zl in ea.zeilen:
        zeile = f"  {zl.dividend:>{breite}} = {zl.q}·{zl.divisor} + {zl.rest}"
        z.append(zeile)
    z.append(f"Der letzte Rest ≠ 0 ist r_{ea.n} = {ea.ggt}.")
    z.append(f"⇒ ggT({ea.gross}, {ea.klein}) = {ea.ggt}")
    return "\n".join(z)


@dataclass
class ErwEAZeile:
    """Zeile der erweiterten EA-Tabelle: r_i = x_i·gross + y_i·klein."""
    i: int
    q: int | None
    r: int
    x: int
    y: int


def erweiterter_ea(a: int, b: int) -> tuple[EAErgebnis, list[ErwEAZeile]]:
    """Erweiterter EA. Gibt den EA und die Tabelle i | q_i | r_i | x_i | y_i zurück.

    Für jede Zeile gilt r_i = x_i·r_(−1) + y_i·r_0.
    Dabei ist r_(−1) die größere und r_0 die kleinere Zahl (Beträge).
    """
    ea = euklid(a, b)
    reste = ea.reste
    qs = ea.quotienten
    x0, y0, x1, y1 = 1, 0, 0, 1
    tab = [ErwEAZeile(-1, None, reste[0], x0, y0)]
    for j in range(1, len(reste)):
        i = j - 1
        q = qs[i] if i < len(qs) else None
        tab.append(ErwEAZeile(i, q, reste[j], x1, y1))
        if q is not None:
            x0, y0, x1, y1 = x1, y1, x0 - q * x1, y0 - q * y1
    return ea, tab


def text_ea_tabelle(ea: EAErgebnis, tabelle: list[ErwEAZeile] | None = None,
                    erweitert: bool = False) -> str:
    """EA als Tabelle: Spalten i | q_i | r_i (optional x_i, y_i)."""
    if tabelle is None:
        _, tabelle = erweiterter_ea(ea.a, ea.b)
    kopf = ["i", "q_i", "r_i"] + (["x_i", "y_i"] if erweitert else [])
    zeilen = []
    for t in tabelle:
        zl = [str(t.i), "" if t.q is None else str(t.q), str(t.r)]
        if erweitert:
            zl += ([str(t.x), str(t.y)] if t.r != 0 or t.i == -1 else ["", ""])
        zeilen.append(zl)
    br = [max(len(kopf[s]), *(len(z[s]) for z in zeilen)) for s in range(len(kopf))]

    def fmt(z: list[str]) -> str:
        return " " + str(z[0]).rjust(br[0]) + " | " + " ".join(c.rjust(w) for c, w in zip(z[1:], br[1:]))

    out = [fmt(kopf), "-" * len(fmt(kopf))] + [fmt(z) for z in zeilen]
    if erweitert:
        out.append(f"(Jede Zeile erfüllt r_i = x_i·{ea.gross} + y_i·{ea.klein}.)")
    return "\n".join(out)


def ggt(a: int, b: int) -> int:
    """ggT(a, b) als Zahl (ohne Lösungsweg)."""
    return euklid(a, b).ggt


def ggt_mit_weg(a: int, b: int) -> tuple[int, str]:
    """ggT(a, b) und der Lösungsweg (EA) als Text."""
    ea = euklid(a, b)
    return ea.ggt, ea.text


@dataclass
class KgVErgebnis:
    a: int
    b: int
    ggt: int
    kgv: int
    ea: EAErgebnis
    text: str = ""


def kgv_mit_weg(a: int, b: int) -> KgVErgebnis:
    """kgV(a, b) über kgV·ggT = a·b. Gibt Wert und Lösungsweg zurück."""
    ea = euklid(a, b)
    wert = 0 if a == 0 or b == 0 else abs(a * b) // ea.ggt
    erg = KgVErgebnis(a, b, ea.ggt, wert, ea)
    erg.text = text_kgv(erg)
    return erg


def kgv(a: int, b: int) -> int:
    """kgV(a, b) als Zahl."""
    return kgv_mit_weg(a, b).kgv


def text_kgv(erg: KgVErgebnis) -> str:
    A, B = abs(erg.a), abs(erg.b)
    z = [text_ea(erg.ea), "", "kgV mit der Formel kgV(a, b)·ggT(a, b) = a·b:"]
    if erg.kgv == 0:
        z.append("Eine Zahl ist 0. Dann ist kgV = 0.")
    else:
        z.append(f"  kgV({A}, {B}) = {A}·{B} / {erg.ggt} = {A * B} / {erg.ggt} = {erg.kgv}")
        z.append(f"Probe: {erg.kgv} / {A} = {erg.kgv // A}, {erg.kgv} / {B} = {erg.kgv // B}.")
    z.append(f"Ergebnis: kgV({A}, {B}) = {erg.kgv}")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Bezout: EA rückwärts und Matrix Q
# ---------------------------------------------------------------------------

@dataclass
class RueckwaertsErgebnis:
    """EA rückwärts: ziel = koef_gross·gross + koef_klein·klein."""
    ea: EAErgebnis
    ziel: int          # der Rest, bei dem das Einsetzen startet
    k: int             # Index dieses Rests (ziel = r_k)
    koef_gross: int
    koef_klein: int
    zeilen: list[str]  # Rechenzeilen ohne "ziel = "
    text: str = ""


def ea_rueckwaerts(a: int, b: int, ziel: int | None = None) -> RueckwaertsErgebnis:
    """EA rückwärts: stellt einen Rest r_k als Kombination von a und b dar.

    Ohne „ziel“ startet das Einsetzen beim ggT r_n.
    Ist „ziel“ selbst ein Rest r_k (k ≥ 1) im EA, startet es dort.
    Rückgabe: Koeffizienten für die größere und die kleinere Zahl (Beträge).
    """
    ea = euklid(a, b)
    k = ea.n
    if ziel is not None:
        for i in range(1, ea.n + 1):
            if ea.r(i) == ziel:
                k = i
                break
    zeilen: list[str] = []
    if ea.klein == 0:
        erg = RueckwaertsErgebnis(ea, ea.gross, -1, 1, 0, [f"1·{ea.gross} + 0·0"])
        erg.text = text_rueckwaerts(erg)
        return erg
    if k == 0:
        # b teilt a: ggT = b = 0·a + 1·b
        erg = RueckwaertsErgebnis(ea, ea.klein, 0, 0, 1, [f"0·{ea.gross} + 1·{ea.klein}"])
        erg.text = text_rueckwaerts(erg)
        return erg
    # Terme als Liste [koeffizient, index]; der Wert eines Index j ist r_j.
    j = k - 1
    terme = [[1, j - 1], [-ea.zeilen[j].q, j]]
    zeilen.append(_summe([(c, str(ea.r(i))) for c, i in terme]))
    while j >= 1:
        # Ersetze r_j = r_(j−2) − q_(j−1)·r_(j−1)
        pos = 0 if terme[0][1] == j else 1
        c = terme[pos][0]
        q = ea.zeilen[j - 1].q
        innen = _summe([(1, str(ea.r(j - 2))), (-q, str(ea.r(j - 1)))])
        teile = []
        for p, (cc, ii) in enumerate(terme):
            if p == pos:
                betrag = "" if abs(c) == 1 else f"{abs(c)}·"
                vz = ("" if c >= 0 else "−") if p == 0 else (" + " if c >= 0 else " − ")
                teile.append(f"{vz}{betrag}({innen})")
            else:
                teile.append(_term(cc, str(ea.r(ii)), p == 0))
        zeilen.append("".join(teile))
        # Vereinfachen: Der ersetzte Platz bekommt r_(j−2), der andere r_(j−1).
        andere = 1 - pos
        terme[andere][0] += -q * c
        terme[pos] = [c, j - 2]
        zeilen.append(_summe([(cc, str(ea.r(ii))) for cc, ii in terme]))
        j -= 1
    kg = sum(cc for cc, ii in terme if ii == -1)
    kk = sum(cc for cc, ii in terme if ii == 0)
    erg = RueckwaertsErgebnis(ea, ea.r(k), k, kg, kk, zeilen)
    erg.text = text_rueckwaerts(erg)
    return erg


def text_rueckwaerts(erg: RueckwaertsErgebnis) -> str:
    """EA rückwärts in drei Schritten: Zeilen nummerieren, nach dem Rest umstellen, einsetzen.

    Jeder Einsetz-Schritt zeigt: welcher Rest ersetzt wird und aus welcher Zeile,
    dann das Ausmultiplizieren, dann das Zusammenfassen.
    """
    ea, k, ziel = erg.ea, erg.k, erg.ziel
    z = [f"EA rückwärts (von r_{k} = {ziel} nach oben einsetzen):"]
    if k < 1:  # Sonderfälle: eine Zahl ist 0 oder b teilt a
        z += [f"  {ziel} = {zl}" for zl in erg.zeilen]
        return "\n".join(z)

    br = len(str(ea.gross))
    z.append("")
    z.append("1) EA-Zeilen nummerieren")
    for i, zl in enumerate(ea.zeilen, start=1):
        z.append(f"  ({i})  {zl.dividend:>{br}} = {zl.q}·{zl.divisor} + {zl.rest}")

    # Zeile (i) liefert den Rest r_i = r_(i−2) − q_(i−1)·r_(i−1).
    z.append("")
    z.append("2) Zeilen nach dem Rest umstellen (von unten nach oben, ohne Null-Zeile)")
    for i in range(k, 0, -1):
        zl = ea.zeilen[i - 1]
        z.append(f"  ({i}')  {zl.rest} = {zl.dividend} − {zl.q}·{zl.divisor}")

    z.append("")
    z.append("3) Rückwärts einsetzen")
    rechts: list[tuple[str, str]] = []   # (Rechnung, Kommentar)
    j = k - 1
    terme = [[1, j - 1], [-ea.zeilen[j].q, j]]
    rechts.append((_summe([(c, str(ea.r(i))) for c, i in terme]), f"Start: ({k}')"))
    while j >= 1:
        pos = 0 if terme[0][1] == j else 1
        c = terme[pos][0]
        q = ea.zeilen[j - 1].q
        rest = ea.r(j)
        # a) Rest r_j durch die umgestellte Zeile (j') ersetzen, in Klammern
        innen = f"{ea.r(j - 2)} − {q}·{ea.r(j - 1)}"
        teile = []
        for p, (cc, ii) in enumerate(terme):
            if p == pos:
                betrag = "" if abs(c) == 1 else f"{abs(c)}·"
                vz = ("" if c >= 0 else "−") if p == 0 else (" + " if c >= 0 else " − ")
                teile.append(f"{vz}{betrag}({innen})")
            else:
                teile.append(_term(cc, str(ea.r(ii)), p == 0))
        rechts.append(("".join(teile), f"{rest} ersetzen durch ({j}')"))
        # b) ausmultiplizieren
        aus: list[tuple[int, str]] = []
        for p, (cc, ii) in enumerate(terme):
            if p == pos:
                aus += [(c, str(ea.r(j - 2))), (-c * q, str(ea.r(j - 1)))]
            else:
                aus.append((cc, str(ea.r(ii))))
        rechts.append((_summe(aus), "ausmultiplizieren"))
        # c) gleiche Zahlen zusammenfassen
        andere = 1 - pos
        alt = terme[andere][0]
        terme[andere][0] += -q * c
        terme[pos] = [c, j - 2]
        kommentar = (f"{ea.r(j - 1)}er zusammenfassen: {_k(alt)} {'−' if q * c >= 0 else '+'} "
                     f"{abs(q * c)} = {terme[andere][0]}")
        rechts.append((_summe([(cc, str(ea.r(ii))) for cc, ii in terme]), kommentar))
        j -= 1

    breite = max(len(r) for r, _ in rechts)
    for n, (r, kom) in enumerate(rechts):
        links = f"{ziel} =" if n == 0 else " " * len(str(ziel)) + " ="
        z.append(f"  {links} {r:<{breite}}   | {kom}")
    return "\n".join(z)


def _matmul(A: tuple, B: tuple) -> tuple:
    (a, b), (c, d) = A
    (e, f), (g, h) = B
    return ((a * e + b * g, a * f + b * h), (c * e + d * g, c * f + d * h))


def _mat_text(ms: list[tuple], trenner: str = "·") -> list[str]:
    """Mehrere 2×2-Matrizen nebeneinander in zwei Zeilen."""
    oben, unten = [], []
    for M in ms:
        w = max(len(str(v)) for zeile in M for v in zeile)
        oben.append(f"⎛{M[0][0]:>{w}} {M[0][1]:>{w}}⎞")
        unten.append(f"⎝{M[1][0]:>{w}} {M[1][1]:>{w}}⎠")
    return [trenner.join(oben), (" " * len(trenner)).join(unten)]


@dataclass
class MatrixErgebnis:
    """Matrix Q = Q_0·…·Q_n = (s u; t v) mit Det Q = (−1)^(n+1).

    Es gilt gross·v − klein·u = Det Q · ggT.
    koef_gross, koef_klein: ggT = koef_gross·gross + koef_klein·klein.
    """
    ea: EAErgebnis
    faktoren: list[tuple]
    zwischen: list[tuple]   # Q_0, Q_0·Q_1, …, Q
    Q: tuple
    det: int
    koef_gross: int
    koef_klein: int
    text: str = ""


def ea_matrix(a: int, b: int) -> MatrixErgebnis:
    """Bezout-Koeffizienten mit der Matrix Q (Beweis des Lemmas von Bezout)."""
    ea = euklid(a, b)
    if ea.klein == 0:
        raise ValueError("Mit einer 0 gibt es keinen EA-Schritt und keine Matrix Q.")
    faktoren = [((q, 1), (1, 0)) for q in ea.quotienten]
    zwischen = [faktoren[0]]
    for M in faktoren[1:]:
        zwischen.append(_matmul(zwischen[-1], M))
    Q = zwischen[-1]
    (s, u), (t, v) = Q
    det = s * v - t * u
    # gross·v − klein·u = det·g  ⇒  g = (det·v)·gross + (−det·u)·klein
    erg = MatrixErgebnis(ea, faktoren, zwischen, Q, det, det * v, -det * u)
    erg.text = text_matrix(erg)
    return erg


def text_matrix(erg: MatrixErgebnis, form: str = "-") -> str:
    ea = erg.ea
    (s, u), (t, v) = erg.Q
    n = ea.n
    z = [f"Matrix Q (Q_i = (q_i 1; 1 0), Quotienten q_0 … q_{n} = {', '.join(map(str, ea.quotienten))}):"]
    kopf = "  Q = "
    zeilen = _mat_text(erg.faktoren)
    z.append(kopf + zeilen[0])
    z.append(" " * len(kopf) + zeilen[1])
    for i, M in enumerate(erg.zwischen[1:], start=1):
        name = "Q_0·…·Q_" + str(i) if i > 1 else "Q_0·Q_1"
        zl = _mat_text([M])
        pre = f"  {name} = "
        z.append(pre + zl[0])
        z.append(" " * len(pre) + zl[1])
    z.append(f"  Q = (s u; t v) mit s = {s}, u = {u}, t = {t}, v = {v}")
    z.append(f"  Kontrolle: (s, t)·ggT = ({s}·{ea.ggt}, {t}·{ea.ggt}) = ({s * ea.ggt}, {t * ea.ggt})"
             f" ✓")
    z.append(f"  Det Q = s·v − t·u = {s}·{v} − {t}·{u} = {erg.det}"
             f"   (= (−1)^(n+1) mit n = {n})")
    return "\n".join(z)


@dataclass
class BezoutErgebnis:
    """ggT(a, b) = x·a + y·b für die Eingaben a, b (mit Vorzeichen)."""
    a: int
    b: int
    ggt: int
    x: int
    y: int
    weg: str                                   # "rueckwaerts" oder "matrix"
    ea: EAErgebnis
    rueckwaerts: RueckwaertsErgebnis | None = None
    matrix: MatrixErgebnis | None = None
    text: str = ""


def _auf_eingabe(ea: EAErgebnis, kg: int, kk: int) -> tuple[int, int]:
    """Koeffizienten für (gross, klein) → Koeffizienten für (a, b) mit Vorzeichen."""
    ka, kb = (kk, kg) if ea.getauscht else (kg, kk)
    if ea.a < 0:
        ka = -ka
    if ea.b < 0:
        kb = -kb
    return ka, kb


def bezout_mit_weg(a: int, b: int, weg: str = "rueckwaerts") -> BezoutErgebnis:
    """Bezout-Koeffizienten x, y mit a·x + b·y = ggT(a, b), mit Lösungsweg.

    weg = "rueckwaerts" (EA rückwärts) oder "matrix" (Matrix Q).
    """
    if weg not in ("rueckwaerts", "matrix"):
        raise ValueError(f"Unbekannter Weg: {weg!r}. Erlaubt: rueckwaerts, matrix.")
    ea = euklid(a, b)
    rw = mx = None
    if weg == "matrix" and ea.klein != 0:
        mx = ea_matrix(a, b)
        kg, kk = mx.koef_gross, mx.koef_klein
    else:
        rw = ea_rueckwaerts(a, b)
        kg, kk = rw.koef_gross, rw.koef_klein
        weg = "rueckwaerts"
    x, y = _auf_eingabe(ea, kg, kk)
    erg = BezoutErgebnis(a, b, ea.ggt, x, y, weg, ea, rw, mx)
    erg.text = text_bezout(erg)
    return erg


def bezout(a: int, b: int, weg: str = "rueckwaerts") -> tuple[int, int, int]:
    """(ggT, x, y) mit a·x + b·y = ggT(a, b)."""
    e = bezout_mit_weg(a, b, weg)
    return e.ggt, e.x, e.y


def text_bezout(erg: BezoutErgebnis) -> str:
    ea = erg.ea
    z = [text_ea(ea), ""]
    if erg.rueckwaerts is not None:
        z.append(erg.rueckwaerts.text)
    if erg.matrix is not None:
        mx = erg.matrix
        (s, u), (t, v) = mx.Q
        z.append(mx.text)
        z.append(f"  Det Q = {mx.det}  | ·({mx.det * ea.ggt})")
        z.append(f"  ⇒ {_lin(ea.gross, mx.koef_gross, ea.klein, mx.koef_klein)} = {ea.ggt}")
    z.append(f"Probe: {erg.a}·{_k(erg.x)} + {erg.b}·{_k(erg.y)} = "
             f"{erg.a * erg.x} + {_k(erg.b * erg.y)} = {erg.a * erg.x + erg.b * erg.y}")
    z.append(f"Ergebnis: ggT({erg.a}, {erg.b}) = {erg.ggt} = {erg.a}·{_k(erg.x)} + {erg.b}·{_k(erg.y)}"
             f"  (Bezout-Koeffizienten x = {erg.x}, y = {erg.y})")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Inverse modulo m
# ---------------------------------------------------------------------------

class KeinInversesFehler(ValueError):
    """a hat kein Inverses mod m, weil ggT(a, m) ≠ 1."""

    def __init__(self, a: int, m: int, g: int, text: str = ""):
        self.a, self.m, self.ggt, self.text = a, m, g, text
        super().__init__(f"{a} hat kein Inverses modulo {m}: ggT({a}, {m}) = {g} ≠ 1.")


@dataclass
class InverseErgebnis:
    a: int
    m: int
    inverse: int          # Repräsentant nach rep
    roh: int              # Wert direkt aus dem EA (z. B. −3)
    rep: str
    weg: str
    bezout: BezoutErgebnis
    text: str = ""


def inverse_mit_weg(a: int, m: int, rep: str = "standard", weg: str = "rueckwaerts") -> InverseErgebnis:
    """a⁻¹ mod m über Bezout: a·x + m·y = 1 ⇒ a·x ≡ 1 (mod m).

    Wirft KeinInversesFehler, wenn ggT(a, m) ≠ 1. Der Fehler enthält den EA als Text.
    rep: "standard" (0 … m−1) oder "symmetrisch" (−m/2 … m/2).
    """
    if m < 2:
        raise ValueError(f"Der Modul m muss mindestens 2 sein (m = {m}).")
    a_red = a % m
    if a_red == 0:
        raise KeinInversesFehler(a, m, m, f"{a} ≡ 0 (mod {m}). 0 hat kein Inverses.")
    bz = bezout_mit_weg(m, a_red, weg)
    if bz.ggt != 1:
        text = (text_ea(bz.ea) + "\n" + f"Hinweis: ggT({a_red}, {m}) = {bz.ggt} ≠ 1. "
                f"Deshalb hat {a} kein Inverses in ℤ_{m}.")
        raise KeinInversesFehler(a, m, bz.ggt, text)
    roh = bz.y
    erg = InverseErgebnis(a, m, repraesentant(roh, m, rep), roh, rep, bz.weg, bz)
    erg.text = text_inverse(erg)
    return erg


def inverse(a: int, m: int, rep: str = "standard") -> int:
    """a⁻¹ mod m als Zahl. Wirft KeinInversesFehler, wenn ggT(a, m) ≠ 1."""
    return inverse_mit_weg(a, m, rep).inverse


def text_inverse(erg: InverseErgebnis) -> str:
    a, m = erg.a, erg.m
    a_red = a % m
    bz = erg.bezout
    z = [f"Gesucht: {a}⁻¹ in ℤ_{m}.",
         f"Idee (Bezout): Für ggT({a_red}, {m}) = 1 hat {a_red}·x + {m}·y = 1 eine Lösung.",
         f"Dann gilt {a_red}·x ≡ 1 (mod {m}), also x = {a_red}⁻¹ in ℤ_{m}."]
    if a_red != a:
        z.append(f"Zuerst reduzieren: {a} ≡ {a_red} (mod {m}).")
    z.append("")
    z.append(f"EA für {m}, {a_red}")
    z.append(text_ea(bz.ea, titel=False))
    if bz.rueckwaerts is not None:
        z.append(bz.rueckwaerts.text)
    if bz.matrix is not None:
        mx = bz.matrix
        z.append(mx.text)
        z.append(f"  Det Q = {mx.det}  | ·({mx.det})")
        z.append(f"  ⇒ {_lin(m, mx.koef_gross, a_red, mx.koef_klein)} = 1")
    z.append(f"⇒ 1 = {a_red}·{_k(erg.roh)} + {m}·{_k(bz.x)} ≡ {a_red}·{_k(erg.roh)} (mod {m})")
    rep_text = "0 … m−1" if erg.rep == "standard" else "−m/2 … m/2 (symmetrisch)"
    if erg.roh != erg.inverse:
        z.append(f"⇒ {a_red}⁻¹ = {erg.roh} ≡ {erg.inverse} (mod {m})   (Repräsentant aus {rep_text})")
    else:
        z.append(f"⇒ {a_red}⁻¹ = {erg.inverse} in ℤ_{m}   (Repräsentant aus {rep_text})")
    p = a_red * erg.inverse
    z.append(f"Probe: {a_red}·{_k(erg.inverse)} = {p} = {p // m}·{m} + {p % m} ≡ 1 (mod {m})")
    z.append(f"Ergebnis: {a}⁻¹ = {erg.inverse} in ℤ_{m}")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Lineare diophantische Gleichung a·x ± b·y = c
# ---------------------------------------------------------------------------

@dataclass
class DiophantErgebnis:
    """Lösung von a·x + b·y = c (form "+") oder a·x − b·y = c (form "-").

    Allgemeine Lösung: x = x + dx·t, y = y + dy·t (t ganz).
    """
    a: int
    b: int
    c: int
    form: str
    ggt: int
    loesbar: bool
    x: int | None = None
    y: int | None = None
    dx: int | None = None
    dy: int | None = None
    x_klein: int | None = None     # spezielle Lösung mit kleinstem x > 0
    y_klein: int | None = None
    weg: str = "rueckwaerts"
    ea: EAErgebnis | None = None
    rueckwaerts: RueckwaertsErgebnis | None = None
    matrix: MatrixErgebnis | None = None
    faktor: int | None = None       # c / (Startwert des Einsetzens)
    x_start: int | None = None      # Lösung für die rechte Seite „Startwert“ (vor dem Hochmultiplizieren)
    y_start: int | None = None
    text: str = ""

    def loesung(self, t: int) -> tuple[int, int]:
        """Lösung für den Parameter t."""
        return self.x + self.dx * t, self.y + self.dy * t


def loese_diophantisch(a: int, b: int, c: int, form: str = "-", weg: str = "rueckwaerts") -> DiophantErgebnis:
    """Löst a·x + b·y = c (form "+") oder a·x − b·y = c (form "-") in ganzen Zahlen.

    weg: "rueckwaerts" (EA rückwärts) oder "matrix" (Matrix Q).
    Ist die Gleichung nicht lösbar (ggT ∤ c), ist loesbar = False.
    """
    form = _form(form)
    if weg not in ("rueckwaerts", "matrix"):
        raise ValueError(f"Unbekannter Weg: {weg!r}. Erlaubt: rueckwaerts, matrix.")
    if a == 0 or b == 0:
        raise ValueError("a und b müssen ≠ 0 sein.")
    ea = euklid(a, b)
    g = ea.ggt
    erg = DiophantErgebnis(a, b, c, form, g, c % g == 0, weg=weg, ea=ea)
    if not erg.loesbar:
        erg.text = text_diophantisch(erg)
        return erg
    B = b if form == "+" else -b          # Gleichung als a·x + B·y = c
    if weg == "rueckwaerts":
        rw = ea_rueckwaerts(a, B, ziel=c)
        erg.rueckwaerts = rw
        kg, kk, start = rw.koef_gross, rw.koef_klein, rw.ziel
    else:
        mx = ea_matrix(a, B)
        erg.matrix = mx
        kg, kk, start = mx.koef_gross, mx.koef_klein, g
    ka, kB = _auf_eingabe(euklid(a, B), kg, kk)   # a·ka + B·kB = start
    erg.faktor = c // start
    # Mit B = ±b ist a·x + B·y genau die Gleichung der Aufgabe.
    erg.x_start, erg.y_start = ka, kB
    erg.x, erg.y = ka * erg.faktor, kB * erg.faktor
    # Homogene Lösung: a·(−B/g) + B·(a/g) = 0.
    # Form "−": (b/g, a/g). Form "+": (−b/g, a/g).
    erg.dx, erg.dy = -B // g, a // g
    if erg.dx < 0 and form == "-":
        erg.dx, erg.dy = -erg.dx, -erg.dy
    h = abs(erg.dx)
    xk = (erg.x - 1) % h + 1          # kleinstes x > 0
    tk = (xk - erg.x) // erg.dx
    erg.x_klein, erg.y_klein = erg.loesung(tk)
    erg.text = text_diophantisch(erg)
    return erg


def _form(form: str) -> str:
    f = {"+": "+", "plus": "+", "-": "-", "−": "-", "minus": "-"}.get(str(form).strip().lower())
    if f is None:
        raise ValueError(f"Unbekannte Form {form!r}. Erlaubt: + oder -.")
    return f


def text_diophantisch(erg: DiophantErgebnis) -> str:
    a, b, c, g = erg.a, erg.b, erg.c, erg.ggt
    op = "+" if erg.form == "+" else "−"
    z = [f"Gleichung: {a}·x {op} {b}·y = {c}   (Form a·x {op} b·y = c)", "",
         "a) ggT mit dem EA", text_ea(erg.ea), ""]
    z.append("b) Lösbarkeit")
    if not erg.loesbar:
        z.append(f"  Die Gleichung ist genau dann lösbar, wenn ggT({abs(a)}, {abs(b)}) = {g} die Zahl {c} teilt.")
        z.append(f"  {c} = {c // g}·{g} + {c % g}. Also teilt {g} die Zahl {c} nicht.")
        z.append(f"Ergebnis: Die Gleichung {a}·x {op} {b}·y = {c} hat keine ganzzahlige Lösung.")
        return "\n".join(z)
    z.append(f"  {c} = {c // g}·{g} ist ein Vielfaches von ggT = {g}. Die Gleichung ist daher lösbar.")
    z.append("")
    z.append(f"c) Eine spezielle Lösung. Weg: {VARIANTEN_TEXT[erg.weg]}")
    ea2 = euklid(a, b if erg.form == "+" else -b)
    G, K = ea2.gross, ea2.klein
    if erg.rueckwaerts is not None:
        rw = erg.rueckwaerts
        if rw.ziel == c and rw.k != ea2.n:
            z.append(f"  {c} ist selbst ein Rest im EA (r_{rw.k}). Wir setzen direkt ab {c} ein.")
        z.append(rw.text)
        start = rw.ziel
        z.append(f"  ⇒ {start} = {_lin(G, rw.koef_gross, K, rw.koef_klein)}")
    else:
        mx = erg.matrix
        (s, u), (t, v) = mx.Q
        z.append(mx.text)
        z.append(f"  Det Q = {mx.det}   | ·({mx.det * g})")
        z.append(f"  ⇒ {_lin(G, mx.koef_gross, K, mx.koef_klein)} = {g}")
        start = g
    xs, ys = erg.x_start, erg.y_start
    z.append(f"  In der Form der Aufgabe: {_gl(a, xs, b, ys, erg.form)} = {start}")
    if erg.faktor != 1:
        z.append(f"  Hochmultiplizieren mit {c}/{start} = {erg.faktor}:   | ·{_k(erg.faktor)}")
        z.append(f"  {_gl(a, erg.x, b, erg.y, erg.form)} = {c}")
    z.append(f"  Spezielle Lösung: x₀ = {erg.x}, y₀ = {erg.y}")
    z.append("")
    # Homogene Gleichung
    z.append("d) Alle Lösungen")
    z.append(f"  Für alle weiteren Lösungen muss gelten: {a}·x {op} {b}·y = 0.")
    hx, hy = erg.dx, erg.dy
    z.append(f"  Es ist {_gl(a, hx * g, b, hy * g, erg.form)} = 0   | : {g}")
    z.append(f"  also {_gl(a, hx, b, hy, erg.form)} = 0   (homogene Lösung ({hx}, {hy}))")
    z.append(f"  Allgemeine Lösung (t ∈ ℤ):  x = {erg.x} {'+' if hx >= 0 else '−'} {abs(hx)}·t,"
             f"  y = {erg.y} {'+' if hy >= 0 else '−'} {abs(hy)}·t")
    z.append(f"  Vektorform: (x, y) = ({erg.x}, {erg.y}) + ({hx}, {hy})·t")
    z.append(f"  Lösung mit kleinstem x > 0: x = {erg.x_klein}, y = {erg.y_klein}")
    z.append("")
    wert = a * erg.x + (b * erg.y if erg.form == "+" else -b * erg.y)
    z.append(f"Probe: {_gl(a, erg.x, b, erg.y, erg.form)} = {a * erg.x} {op} {_k(b * erg.y)} = {wert}"
             + ("  ✓" if wert == c else "  FEHLER"))
    z.append(f"Ergebnis: ggT = {g}; spezielle Lösung x = {erg.x}, y = {erg.y}; "
             f"alle Lösungen x = {erg.x} {'+' if hx >= 0 else '−'} {abs(hx)}t, "
             f"y = {erg.y} {'+' if hy >= 0 else '−'} {abs(hy)}t (t ∈ ℤ)")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

BEISPIELE = """Beispiele:
  python -m skripte.euklid --a 2406 --b 654 --c 24 --form -   (Gleichung a·x − b·y = c, beide Wege)
  python -m skripte.euklid --a 1001 --b 840 --c 98 --form + --variante matrix   (Gleichung a·x + b·y = c mit Matrix Q)
  python -m skripte.euklid --a 17 --b 15 --c 1 --form -   (Gleichung mit c = 1)
  python -m skripte.euklid --a 1965 --b 225 --variante matrix   (ggT und Bezout mit Matrix Q)
  python -m skripte.euklid --a 6 --mod 11   (Inverse modulo 11)
  python -m skripte.euklid --a 13 --mod 100 --variante matrix   (Inverse modulo 100 mit Matrix Q)
  python -m skripte.euklid --a 12 --b 18 --kgv   (ggT und kgV)
"""


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m skripte.euklid",
        description="Euklidischer Algorithmus mit Lösungsweg: ggT, kgV, Bezout, "
                    "Inverse modulo m, diophantische Gleichung a·x ± b·y = c.",
        epilog=BEISPIELE, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", type=int, required=True, help="Erste Zahl a.")
    p.add_argument("--b", type=int, help="Zweite Zahl b.")
    p.add_argument("--c", type=int, help="Rechte Seite c der Gleichung a·x ± b·y = c.")
    p.add_argument("--form", default="-",
                   help="Form der Gleichung: + für a·x + b·y = c, - für a·x − b·y = c. "
                        "Standard: -.")
    p.add_argument("--mod", type=int, help="Modul m. Das Skript berechnet dann a⁻¹ in ℤ_m.")
    p.add_argument("--variante", choices=VARIANTEN, default="beide",
                   help="Lösungsweg: beide (Standard), rueckwaerts (EA rückwärts) oder matrix (Matrix Q).")
    p.add_argument("--alle-varianten", action="store_true",
                   help="Gibt alle Lösungswege nacheinander aus.")
    p.add_argument("--rep", choices=("standard", "symmetrisch"), default="standard",
                   help="Repräsentant der Inversen: standard (0 … m−1) oder symmetrisch (−m/2 … m/2).")
    p.add_argument("--kgv", action="store_true", help="Berechnet zusätzlich kgV(a, b).")
    p.add_argument("--erweitert", action="store_true",
                   help="Zeigt die EA-Tabelle mit den Spalten x_i, y_i (erweiterter EA).")
    return p


def _wege(args) -> list[str]:
    if args.alle_varianten or args.variante == "beide":
        return ["rueckwaerts", "matrix"]
    return [args.variante]


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    wege = _wege(args)
    aus: list[str] = []
    try:
        if args.mod is not None:
            aus.append(f"Annahme: Lösungsweg = {VARIANTEN_TEXT['beide' if len(wege) > 1 else wege[0]]}; "
                       f"Repräsentant = {args.rep}.")
            for i, w in enumerate(wege):
                if len(wege) > 1:
                    aus.append(f"\n=== Weg {i + 1}: {VARIANTEN_TEXT[w]} ===")
                aus.append(inverse_mit_weg(args.a, args.mod, args.rep, w).text)
        elif args.b is None:
            print("Hinweis: Bitte --b angeben (oder --mod für die Inverse).", file=sys.stderr)
            return 2
        elif args.c is not None:
            form = _form(args.form)
            aus.append(f"Annahme: Form a·x {'+' if form == '+' else '−'} b·y = c; "
                       f"Lösungsweg = {VARIANTEN_TEXT['beide' if len(wege) > 1 else wege[0]]}.")
            ergs = []
            for i, w in enumerate(wege):
                if len(wege) > 1:
                    aus.append(f"\n=== Weg {i + 1}: {VARIANTEN_TEXT[w]} ===")
                e = loese_diophantisch(args.a, args.b, args.c, form, w)
                ergs.append(e)
                aus.append(e.text)
                if args.erweitert and i == 0:
                    aus.append("\nEA als Tabelle:\n" + text_ea_tabelle(e.ea, erweitert=True))
                if not e.loesbar:
                    break
            if len(ergs) > 1 and ergs[0].loesbar and (ergs[0].x, ergs[0].y) != (ergs[1].x, ergs[1].y):
                aus.append(f"\nHinweis: Beide Wege sind richtig. Sie geben verschiedene spezielle Lösungen: "
                           f"({ergs[0].x}, {ergs[0].y}) und ({ergs[1].x}, {ergs[1].y}). "
                           f"Sie unterscheiden sich um ein Vielfaches der homogenen Lösung.")
        else:
            aus.append(f"Annahme: Lösungsweg = {VARIANTEN_TEXT['beide' if len(wege) > 1 else wege[0]]}.")
            for i, w in enumerate(wege):
                if len(wege) > 1:
                    aus.append(f"\n=== Weg {i + 1}: {VARIANTEN_TEXT[w]} ===")
                aus.append(bezout_mit_weg(args.a, args.b, w).text)
            aus.append("\nEA als Tabelle:\n" + text_ea_tabelle(euklid(args.a, args.b), erweitert=args.erweitert))
            if args.kgv:
                aus.append("")
                k = kgv_mit_weg(args.a, args.b)
                aus.append(k.text.split("\n\n", 1)[1])
    except KeinInversesFehler as f:
        print(f.text or str(f))
        print(f"Ergebnis: Kein Inverses. ggT({args.a % args.mod}, {args.mod}) = {f.ggt} ≠ 1.")
        return 1
    except ValueError as f:
        print(f"Hinweis: {f}")
        return 2
    print("\n".join(aus))
    return 0


if __name__ == "__main__":
    sys.exit(main())
