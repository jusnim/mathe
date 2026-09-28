"""Zyklische Codes über F_p (p Primzahl).

Zweck
-----
Das Skript rechnet die Aufgaben zu zyklischen Codes aus VL5, VL6 und Übung 6.
Ein zyklischer Code ist ein linearer Code. Jede zyklische Verschiebung eines
Codeworts ist wieder ein Codewort.

Das Skript zeigt:
- Kreisteilungsklassen von Z_n (Bahnen s, s·p, s·p², … modulo n),
- Zerlegung von X^n − 1 in irreduzible Polynome (durch Probedivision),
- alle Generatorpolynome g(X), also alle normierten Teiler von X^n − 1,
- zu einem g(X): zyklische Generatormatrix G, Kontrollpolynom h(X),
  reziprokes Polynom, Kontrollmatrix H, Mindestabstand d,
- Codieren c(X) = a(X)·g(X),
- Syndrom s(X) = Rest von y(X) : g(X) und Decodieren mit Klassenführern,
- alle irreduziblen Polynome vom Grad d (Ü6 6_1 a),
- das Produkt a * b in F_p[X]/(X^n − 1) (VL5, Beispiel F_5).

Jede Polynomdivision steht Schritt für Schritt in der Ausgabe.

Aufruf
------
    python -m skripte.zyklischer_code --p 2 --n 7                  (Ü6 6_1 b)
    python -m skripte.zyklischer_code --p 2 --irreduzibel 3        (Ü6 6_1 a)
    python -m skripte.zyklischer_code --p 2 --n 7 --k 4            (Ü6 6_2)
    python -m skripte.zyklischer_code --p 3 --n 4 --k 2 --d 3      (Ü6 6_6, 6_7)
    python -m skripte.zyklischer_code --p 2 --n 7 --g "X^3+X+1" --y "X^6+X+1"   (VL6)
    python -m skripte.zyklischer_code --p 5 --n 4 --produkt 0012 2314           (VL5)

Konventionen (wie VL5)
----------------------
- Stellen zählen ab 0: c_0 c_1 … c_{n−1}.
- Wort a_0 a_1 … a_{n−1} gehört zu a_0 + a_1·X + … + a_{n−1}·X^{n−1}.
  Der niedrigste Grad steht links. Beispiel: g = 1 + X + X³ ↔ 1101000.
- X * a(X) verschiebt das Wort zyklisch um eine Stelle nach rechts.

Varianten (--variante)
----------------------
- symmetrisch (Standard): Elemente von F_p als −(p−1)/2 … (p−1)/2.
  In F_3 steht also −1 statt 2. So schreiben VL5 (Golay G_11) und Ü6 6_6.
- standard: Elemente von F_p als 0 … p−1. In F_3 steht 2 statt −1.
Die Varianten ändern nur die Schreibweise. Die Codes sind gleich.
--alle-varianten gibt beide Versionen nacheinander aus.
"""

from __future__ import annotations

import argparse
import itertools
import math
import re
import sys
from dataclasses import dataclass, field

import sympy

VARIANTEN = {
    "symmetrisch": "Elemente von F_p als -(p-1)/2 … (p-1)/2 (in F_3: -1 statt 2), wie VL5 und Ü6 6_6; Kreisteilungsklassen mit ±, wie Ü6 6_1",
    "standard": "Elemente von F_p als 0 … p-1 (in F_3: 2 statt -1); Kreisteilungsklassen mit 0 … n-1",
}
STANDARD_VARIANTE = "symmetrisch"

# Grenzen, damit das Skript nicht zu lange rechnet.
MAX_KANDIDATEN = 300_000     # Probedivision: Kandidaten pro Grad
MAX_CODEWOERTER = 60_000     # Mindestabstand durch Aufzählen
MAX_TABELLE = 64             # Klassenführer-Tabelle nur bis zu so vielen Zeilen
MAX_MUSTER = 2_000_000       # Suche nach Klassenführern


class EingabeFehler(ValueError):
    """Die Eingabe passt nicht. Die Nachricht erklärt, was falsch ist."""


# ---------------------------------------------------------------------------
# Rechnen mit Polynomen über F_p
# Ein Polynom ist eine Liste [a_0, a_1, …, a_m] mit Werten 0 … p−1.
# Das Nullpolynom ist die leere Liste.
# ---------------------------------------------------------------------------

def norm(a: list[int], p: int) -> list[int]:
    """Rechnet alle Koeffizienten mod p und entfernt Nullen am Ende."""
    b = [x % p for x in a]
    while b and b[-1] == 0:
        b.pop()
    return b


def grad(a: list[int]) -> int:
    """Grad des Polynoms. Das Nullpolynom hat hier Grad −1."""
    return len(a) - 1


def p_add(a: list[int], b: list[int], p: int) -> list[int]:
    m = max(len(a), len(b))
    return norm([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(m)], p)


def p_sub(a: list[int], b: list[int], p: int) -> list[int]:
    return p_add(a, [-x for x in b], p)


def p_mul(a: list[int], b: list[int], p: int) -> list[int]:
    if not a or not b:
        return []
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return norm(c, p)


def monom(c: int, e: int, p: int) -> list[int]:
    """Das Polynom c·X^e."""
    return norm([0] * e + [c], p)


@dataclass
class DivSchritt:
    """Ein Schritt der Polynomdivision."""

    rest_vorher: list[int]
    term: list[int]          # c·X^e, der neue Summand des Quotienten
    produkt: list[int]       # term · Divisor
    rest_nachher: list[int]


@dataclass
class Division:
    dividend: list[int]
    divisor: list[int]
    quotient: list[int]
    rest: list[int]
    schritte: list[DivSchritt]


def p_divmod(a: list[int], b: list[int], p: int) -> Division:
    """Polynomdivision mit Rest: a = q·b + r mit grad r < grad b."""
    a = norm(a, p)
    b = norm(b, p)
    if not b:
        raise EingabeFehler("Division durch das Nullpolynom ist nicht erlaubt.")
    inv = pow(b[-1], -1, p)
    rest = a[:]
    q: list[int] = []
    schritte = []
    while rest and grad(rest) >= grad(b):
        e = grad(rest) - grad(b)
        c = rest[-1] * inv % p
        term = monom(c, e, p)
        prod = p_mul(term, b, p)
        neu = p_sub(rest, prod, p)
        schritte.append(DivSchritt(rest, term, prod, neu))
        q = p_add(q, term, p)
        rest = neu
    return Division(a, b, q, rest, schritte)


def x_hoch_n_minus_1(n: int, p: int) -> list[int]:
    return norm([-1] + [0] * (n - 1) + [1], p)


def ist_primzahl(p: int) -> bool:
    return p >= 2 and all(p % t for t in range(2, math.isqrt(p) + 1))


def reziprok(h: list[int]) -> list[int]:
    """Reziprokes Polynom X^k·h(1/X): Koeffizienten in umgekehrter Reihenfolge."""
    return list(reversed(h))


def verschiebungen(g: list[int], n: int, anzahl: int) -> list[list[int]]:
    """Zeilen g, X*g, …, X^(anzahl−1)*g als Wörter der Länge n."""
    zeilen = []
    for i in range(anzahl):
        w = [0] * n
        for j, c in enumerate(g):
            w[(i + j) % n] = c
        zeilen.append(w)
    return zeilen


def wort_zu_poly(w: list[int], p: int) -> list[int]:
    return norm(w, p)


def poly_zu_wort(a: list[int], n: int) -> list[int]:
    return a + [0] * (n - len(a))


def gewicht(w: list[int]) -> int:
    return sum(1 for x in w if x)


# ---------------------------------------------------------------------------
# Eingaben lesen
# ---------------------------------------------------------------------------

_X = sympy.Symbol("X")


def lies_polynom(text: str, p: int) -> list[int]:
    """Liest "1 + X + X^3" oder "(X+1)^2" oder "-1+X^2"."""
    t = text.replace("−", "-").replace("^", "**").replace("x", "X").replace("·", "*")
    t = re.sub(r"(\d)\s*X", r"\1*X", t)
    try:
        ausdruck = sympy.sympify(t, locals={"X": _X})
        koeff = sympy.Poly(ausdruck, _X).all_coeffs()
    except (sympy.SympifyError, sympy.PolynomialError, TypeError, SyntaxError) as fehler:
        raise EingabeFehler(f"Das Polynom '{text}' kann ich nicht lesen ({fehler}).") from None
    if any(not c.is_Integer for c in koeff):
        raise EingabeFehler(f"Das Polynom '{text}' hat Koeffizienten, die keine ganzen Zahlen sind.")
    return norm([int(c) for c in reversed(koeff)], p)


def lies_wort(text: str, p: int) -> list[int]:
    """Liest "1101000", "2 2 0 0", "-1 0 1 0" oder "2,2,0,0"."""
    t = text.replace("−", "-").strip()
    if re.search(r"[\s,;]", t):
        teile = [s for s in re.split(r"[\s,;]+", t) if s]
    elif p <= 10 and "-" not in t:
        teile = list(t)
    else:
        teile = [t]
    try:
        werte = [int(s) for s in teile]
    except ValueError:
        raise EingabeFehler(f"Das Wort '{text}' enthält Zeichen, die keine Zahlen sind.") from None
    return [x % p for x in werte]


def lies_poly_oder_wort(text: str, p: int) -> tuple[list[int], str]:
    """Gibt (Polynom, Art) zurück. Art ist "polynom" oder "wort"."""
    if re.search(r"[xX]", text):
        return lies_polynom(text, p), "polynom"
    return wort_zu_poly(lies_wort(text, p), p), "wort"


# ---------------------------------------------------------------------------
# Kreisteilungsklassen und Zerlegung von X^n − 1
# ---------------------------------------------------------------------------

def kreisteilungsklassen(n: int, p: int) -> list[list[int]] | None:
    """Bahnen B_s = {s, s·p, s·p², …} in Z_n. None, falls ggT(n, p) ≠ 1."""
    if math.gcd(n, p) != 1:
        return None
    gesehen: set[int] = set()
    klassen = []
    for s in range(n):
        if s in gesehen:
            continue
        bahn = [s]
        x = s * p % n
        while x != s:
            bahn.append(x)
            x = x * p % n
        gesehen.update(bahn)
        klassen.append(bahn)
    return klassen


def normierte_polynome(d: int, p: int, ohne_x_teiler: bool = True):
    """Alle normierten Polynome vom Grad d. Optional ohne Konstante 0."""
    for koeff in itertools.product(range(p), repeat=d):
        if ohne_x_teiler and d >= 1 and koeff[0] == 0:
            continue
        yield list(koeff) + [1]


def nullstellen(f: list[int], p: int) -> list[int]:
    return [a for a in range(p) if sum(c * pow(a, i, p) for i, c in enumerate(f)) % p == 0]


@dataclass
class IrrKandidat:
    poly: list[int]
    nullstellen: list[int]
    teiler: list[int] | None      # gefundener irreduzibler Teiler vom Grad ≥ 2
    irreduzibel: bool


@dataclass
class IrreduzibelErgebnis:
    p: int
    d: int
    kandidaten: list[IrrKandidat]
    irreduzibel: list[list[int]]


def irreduzible_polynome(d: int, p: int) -> IrreduzibelErgebnis:
    """Alle normierten irreduziblen Polynome vom Grad d über F_p."""
    if d < 1:
        raise EingabeFehler("Der Grad muss mindestens 1 sein.")
    if p ** d > MAX_KANDIDATEN:
        raise EingabeFehler(f"p^d = {p ** d} Kandidaten sind zu viele zum Probieren.")
    kleine: list[list[int]] = []
    for e in range(2, d // 2 + 1):
        kleine += irreduzible_polynome(e, p).irreduzibel
    kandidaten = []
    for f in normierte_polynome(d, p, ohne_x_teiler=(d > 1)):
        ns = nullstellen(f, p)
        teiler = None
        if d > 1 and not ns:
            for t in kleine:
                if not p_divmod(f, t, p).rest:
                    teiler = t
                    break
        irr = (d == 1) or (not ns and teiler is None)
        kandidaten.append(IrrKandidat(f, ns, teiler, irr))
    return IrreduzibelErgebnis(p, d, kandidaten, [k.poly for k in kandidaten if k.irreduzibel])


@dataclass
class ZerlegungsSchritt:
    faktor: list[int]
    division: Division


@dataclass
class Zerlegung:
    n: int
    p: int
    nullstellen_test: list[tuple[int, int]]       # (a, a^n mod p)
    schritte: list[ZerlegungsSchritt]
    faktoren: list[tuple[list[int], int]]          # (irreduzibles Polynom, Vielfachheit)
    getestet: dict[int, int]                       # Grad → Anzahl Kandidaten
    rest_irreduzibel: list[int] | None             # letzter Rest, ohne Test als irreduzibel erkannt
    sympy_benutzt: bool = False


def zerlege(n: int, p: int) -> Zerlegung:
    """Zerlegt X^n − 1 durch Probedivision mit normierten Polynomen steigenden Grades.

    Weil die Grade steigen, ist jeder gefundene Teiler automatisch irreduzibel.
    Ist der Rest vom Grad < 2d, dann ist er selbst irreduzibel.
    """
    rest = x_hoch_n_minus_1(n, p)
    ns_test = [(a, pow(a, n, p)) for a in range(p)]
    schritte: list[ZerlegungsSchritt] = []
    faktoren: dict[tuple[int, ...], int] = {}
    getestet: dict[int, int] = {}
    rest_irr = None
    sympy_benutzt = False
    d = 1
    while grad(rest) >= 1:
        if grad(rest) < 2 * d:
            rest_irr = rest
            faktoren[tuple(rest)] = faktoren.get(tuple(rest), 0) + 1
            schritte.append(ZerlegungsSchritt(rest, p_divmod(rest, rest, p)))
            rest = [1]
            break
        if p ** d > MAX_KANDIDATEN:
            # Zu viele Kandidaten: Rest mit sympy zerlegen und die Divisionen trotzdem zeigen.
            sympy_benutzt = True
            fl = sympy.factor_list(sympy.Poly(list(reversed(rest)), _X, modulus=p))
            for f, m in sorted(fl[1], key=lambda t: t[0].degree()):
                fk = norm([int(c) for c in reversed(f.all_coeffs())], p)
                for _ in range(m):
                    div = p_divmod(rest, fk, p)
                    schritte.append(ZerlegungsSchritt(fk, div))
                    rest = div.quotient
                faktoren[tuple(fk)] = faktoren.get(tuple(fk), 0) + m
            break
        anzahl = 0
        for f in normierte_polynome(d, p):
            anzahl += 1
            while True:
                div = p_divmod(rest, f, p)
                if div.rest:
                    break
                schritte.append(ZerlegungsSchritt(f, div))
                faktoren[tuple(f)] = faktoren.get(tuple(f), 0) + 1
                rest = div.quotient
            if grad(rest) < 2 * d:
                break
        getestet[d] = anzahl
        d += 1
    liste = sorted(((list(f), m) for f, m in faktoren.items()), key=lambda t: (len(t[0]), tuple(reversed(t[0]))))
    return Zerlegung(n, p, ns_test, schritte, liste, getestet, rest_irr, sympy_benutzt)


# ---------------------------------------------------------------------------
# Generatorpolynome und Codes
# ---------------------------------------------------------------------------

@dataclass
class Generator:
    g: list[int]
    exponenten: list[int]    # Exponent jedes irreduziblen Faktors in g
    r: int
    k: int
    d: int | None


def alle_generatorpolynome(zerl: Zerlegung) -> list[Generator]:
    """Alle normierten Teiler von X^n − 1."""
    n, p = zerl.n, zerl.p
    ergebnis = []
    for exps in itertools.product(*[range(m + 1) for _, m in zerl.faktoren]):
        g = [1]
        for (f, _), e in zip(zerl.faktoren, exps):
            for _ in range(e):
                g = p_mul(g, f, p)
        r = grad(g)
        ergebnis.append(Generator(g, list(exps), r, n - r, mindestabstand(g, n, p)))
    ergebnis.sort(key=lambda gen: (gen.r, tuple(reversed(gen.g))))
    return ergebnis


def mindestabstand(g: list[int], n: int, p: int) -> int | None:
    """Kleinstes Gewicht ≠ 0 aller Codewörter. None, wenn es zu viele Wörter sind."""
    k = n - grad(g)
    if k == 0:
        return None     # Nullcode: d ist nicht definiert
    if p ** k > MAX_CODEWOERTER:
        return None
    zeilen = verschiebungen(g, n, k)
    best = n + 1
    for a in itertools.product(range(p), repeat=k):
        if not any(a):
            continue
        c = [sum(a[i] * zeilen[i][j] for i in range(k)) % p for j in range(n)]
        best = min(best, gewicht(c))
    return best


@dataclass
class CodeDaten:
    n: int
    p: int
    g: list[int]
    r: int
    k: int
    G: list[list[int]]
    h_division: Division
    h: list[int]
    h_rez: list[int]
    H: list[list[int]]
    probe_GHt: list[list[int]]
    d: int | None


def code_daten(g: list[int], n: int, p: int) -> CodeDaten:
    """G, h(X), reziprokes h, H und d zu einem Generatorpolynom g."""
    div = p_divmod(x_hoch_n_minus_1(n, p), g, p)
    if div.rest:
        raise EingabeFehler("g(X) teilt X^n − 1 nicht.")
    r = grad(g)
    k = n - r
    h = div.quotient
    h_rez = reziprok(h)
    G = verschiebungen(g, n, k)
    H = verschiebungen(h_rez, n, r)
    probe = [[sum(G[i][t] * H[j][t] for t in range(n)) % p for j in range(r)] for i in range(k)]
    return CodeDaten(n, p, g, r, k, G, div, h, h_rez, H, probe, mindestabstand(g, n, p))


@dataclass
class Kodierung:
    a: list[int]            # Klartext als Polynom
    summanden: list[tuple[int, list[int]]]   # (i, a_i·X^i·g(X))
    c: list[int]


def codieren(a: list[int], g: list[int], n: int, p: int) -> Kodierung:
    """c(X) = a(X)·g(X) mit grad a < k."""
    k = n - grad(g)
    if grad(a) >= k:
        raise EingabeFehler(f"Das Klartextwort braucht grad a(X) < k = {k} (also höchstens {k} Stellen).")
    summanden = [(i, p_mul(monom(ai, i, p), g, p)) for i, ai in enumerate(a) if ai]
    return Kodierung(a, summanden, p_mul(a, g, p))


@dataclass
class ProduktModulo:
    a: list[int]
    b: list[int]
    produkt: list[int]       # in F_p[X]
    ergebnis: list[int]      # in F_p[X]/(X^n − 1)


def produkt_modulo(a: list[int], b: list[int], n: int, p: int) -> ProduktModulo:
    """a * b in F_p[X]/(X^n − 1): X^(n+j) wird zu X^j."""
    prod = p_mul(a, b, p)
    red = [0] * n
    for i, c in enumerate(prod):
        red[i % n] += c
    return ProduktModulo(a, b, prod, norm(red, p))


def syndrom_von(e: list[int], g: list[int], p: int) -> tuple[int, ...]:
    return tuple(p_divmod(e, g, p).rest)


def fehlermuster(n: int, p: int, w: int):
    """Alle Fehlermuster vom Gewicht w. Erst nach Fehlergröße, dann nach Stelle (wie Ü6 6_8)."""
    for werte in itertools.product(range(1, p), repeat=w):
        for stellen in itertools.combinations(range(n), w):
            e = [0] * n
            for s, v in zip(stellen, werte):
                e[s] = v
            yield norm(e, p)


def klassenfuehrer_tabelle(g: list[int], n: int, p: int) -> list[tuple[list[int], list[int]]] | None:
    """Klassenführer (kleinstes Gewicht) und Syndrom jeder Nebenklasse.

    None, wenn es mehr als MAX_TABELLE Nebenklassen gibt.
    """
    r = grad(g)
    anzahl = p ** r
    if anzahl > MAX_TABELLE:
        return None
    tabelle: dict[tuple[int, ...], list[int]] = {}
    for w in range(n + 1):
        for e in fehlermuster(n, p, w):
            s = syndrom_von(e, g, p)
            if s not in tabelle:
                tabelle[s] = e
        if len(tabelle) == anzahl:
            break
    return [(e, list(s)) for s, e in tabelle.items()]


@dataclass
class Decodierung:
    y: list[int]
    division: Division
    s: list[int]
    fuehrer: list[int] | None
    gleich_gute: int                # Anzahl Fehlermuster mit gleichem kleinsten Gewicht
    c: list[int] | None
    c_division: Division | None
    abstand: int | None
    abgebrochen: bool = False


def decodieren(y: list[int], g: list[int], n: int, p: int) -> Decodierung:
    """Syndrom s = y mod g. Dann Klassenführer e suchen und c = y − e."""
    div = p_divmod(y, g, p)
    s = div.rest
    if not s:
        return Decodierung(y, div, s, [], 1, y, div, 0)
    ziel = tuple(s)
    gezaehlt = 0
    for w in range(1, n + 1):
        gefunden = []
        for e in fehlermuster(n, p, w):
            gezaehlt += 1
            if gezaehlt > MAX_MUSTER:
                return Decodierung(y, div, s, None, 0, None, None, None, abgebrochen=True)
            if syndrom_von(e, g, p) == ziel:
                gefunden.append(e)
        if gefunden:
            e = gefunden[0]
            c = p_sub(y, e, p)
            return Decodierung(y, div, s, e, len(gefunden), c, p_divmod(c, g, p), gewicht(e))
    return Decodierung(y, div, s, None, 0, None, None, None)


# ---------------------------------------------------------------------------
# Gesamtlösung
# ---------------------------------------------------------------------------

@dataclass
class Fall:
    """Ein Code mit festem g(X), dazu Codieren und Decodieren."""

    code: CodeDaten
    kodierung: Kodierung | None = None
    decodierung: Decodierung | None = None


@dataclass
class Loesung:
    p: int
    n: int | None
    k: int | None
    d_soll: int | None
    hinweise: list[str] = field(default_factory=list)
    irreduzibel: IrreduzibelErgebnis | None = None
    produkt: ProduktModulo | None = None
    klassen: list[list[int]] | None = None
    zerlegung: Zerlegung | None = None
    generatoren: list[Generator] = field(default_factory=list)
    g_eingabe: list[int] | None = None
    g_pruefung: Division | None = None
    faelle: list[Fall] = field(default_factory=list)


def loese(p: int, n: int | None = None, g: str | None = None, k: int | None = None,
          a: str | None = None, y: str | None = None, d: int | None = None,
          irreduzibel: int | None = None, produkt: tuple[str, str] | None = None) -> Loesung:
    """Rechnet alles, was die Eingaben verlangen. Gibt nur Daten zurück."""
    if not ist_primzahl(p):
        raise EingabeFehler(f"p = {p} ist keine Primzahl. Das Skript rechnet nur über F_p mit p prim.")
    los = Loesung(p, n, k, d)
    if irreduzibel is not None:
        los.irreduzibel = irreduzible_polynome(irreduzibel, p)
    if n is None:
        if irreduzibel is None:
            raise EingabeFehler("Bitte die Wortlänge --n angeben (oder --irreduzibel D).")
        return los
    if n < 1:
        raise EingabeFehler("Die Wortlänge n muss mindestens 1 sein.")
    if k is not None and not 0 <= k <= n:
        raise EingabeFehler(f"k = {k} passt nicht. Es muss 0 ≤ k ≤ n = {n} gelten.")

    if produkt is not None:
        wa = lies_wort(produkt[0], p)
        wb = lies_wort(produkt[1], p)
        for w, t in ((wa, produkt[0]), (wb, produkt[1])):
            if len(w) != n:
                los.hinweise.append(f"Achtung: Das Wort '{t}' hat {len(w)} Stellen, n ist aber {n}.")
        los.produkt = produkt_modulo(norm(wa, p), norm(wb, p), n, p)

    if math.gcd(n, p) != 1:
        los.hinweise.append(
            f"ggT(n, p) = ggT({n}, {p}) = {math.gcd(n, p)} ≠ 1. Dann hat X^n − 1 mehrfache Nullstellen. "
            "Die Kreisteilungsklassen gelten nur für ggT(n, p) = 1 (VL5, Abschnitt 5.5). "
            "Die Zerlegung klappt trotzdem durch Probedivision.")
    los.klassen = kreisteilungsklassen(n, p)
    los.zerlegung = zerlege(n, p)
    los.generatoren = alle_generatorpolynome(los.zerlegung)

    gewaehlt: list[list[int]] = []
    if g is not None:
        gp = lies_polynom(g, p) if re.search(r"[xX]", g) else wort_zu_poly(lies_wort(g, p), p)
        if not gp:
            raise EingabeFehler("g(X) ist das Nullpolynom. Das ist kein Generatorpolynom.")
        if gp[-1] != 1:
            inv = pow(gp[-1], -1, p)
            los.hinweise.append(
                f"g(X) ist nicht normiert (höchster Koeffizient {gp[-1]}). "
                f"Das Skript multipliziert g(X) mit {inv} (Inverses mod {p}).")
            gp = norm([c * inv for c in gp], p)
        los.g_eingabe = gp
        los.g_pruefung = p_divmod(x_hoch_n_minus_1(n, p), gp, p)
        if los.g_pruefung.rest:
            los.hinweise.append(
                "g(X) teilt X^n − 1 nicht (Rest ≠ 0). Also ist g(X) kein Generatorpolynom "
                f"eines zyklischen Codes der Länge {n}.")
        else:
            gewaehlt = [gp]
            if k is not None and n - grad(gp) != k:
                los.hinweise.append(f"Achtung: grad g = {grad(gp)}, also k = n − r = {n - grad(gp)} und nicht {k}.")
    elif k is not None:
        gewaehlt = [gen.g for gen in los.generatoren if gen.k == k]
        if not gewaehlt:
            los.hinweise.append(f"Es gibt keinen zyklischen Code mit n = {n} und k = {k} über F_{p}.")

    if (a is not None or y is not None) and len(gewaehlt) != 1:
        if len(gewaehlt) > 1:
            los.hinweise.append(
                "Für Codieren/Decodieren gibt es mehrere passende g(X). "
                "Das Skript rechnet es für jeden Fall. Wähle g(X) mit --g, wenn die Aufgabe g vorgibt.")
        else:
            los.hinweise.append("Für Codieren/Decodieren fehlt ein gültiges g(X). Bitte --g oder --k angeben.")

    for gp in gewaehlt:
        fall = Fall(code_daten(gp, n, p))
        if a is not None:
            ap, art = lies_poly_oder_wort(a, p)
            if art == "wort" and len(lies_wort(a, p)) != fall.code.k:
                los.hinweise.append(f"Achtung: Das Klartextwort hat {len(lies_wort(a, p))} Stellen, k ist {fall.code.k}.")
            try:
                fall.kodierung = codieren(ap, gp, n, p)
            except EingabeFehler as fehler:
                los.hinweise.append(str(fehler))
        if y is not None:
            yp, art = lies_poly_oder_wort(y, p)
            if art == "wort" and len(lies_wort(y, p)) != n:
                los.hinweise.append(f"Achtung: Das empfangene Wort hat {len(lies_wort(y, p))} Stellen, n ist {n}.")
            if grad(yp) >= n:
                los.hinweise.append(f"Achtung: grad y(X) = {grad(yp)} ≥ n = {n}. Das ist kein Wort der Länge n.")
            fall.decodierung = decodieren(yp, gp, n, p)
        los.faelle.append(fall)
    return los


# ---------------------------------------------------------------------------
# Ausgabe
# ---------------------------------------------------------------------------

def zahl(c: int, p: int, variante: str) -> int:
    c %= p
    if variante == "symmetrisch" and c > p // 2:
        c -= p
    return c


def fmt_poly(a: list[int], p: int, variante: str, aufsteigend: bool = False) -> str:
    """Polynom als Text, z. B. "X^3 + X + 1" oder aufsteigend "1 + X + X^3"."""
    terme = [(e, zahl(c, p, variante)) for e, c in enumerate(a) if c % p]
    if not terme:
        return "0"
    if not aufsteigend:
        terme.reverse()
    teile = []
    for i, (e, c) in enumerate(terme):
        betrag = abs(c)
        x = "" if e == 0 else ("X" if e == 1 else f"X^{e}")
        if e == 0:
            t = str(betrag)
        else:
            t = x if betrag == 1 else f"{betrag}{x}"
        if i == 0:
            teile.append(("-" if c < 0 else "") + t)
        else:
            teile.append((" - " if c < 0 else " + ") + t)
    return "".join(teile)


def klammer(a: list[int], p: int, variante: str, aufsteigend: bool = False) -> str:
    s = fmt_poly(a, p, variante, aufsteigend)
    return f"({s})" if len([c for c in a if c % p]) > 1 else s


def fmt_wort(w: list[int], p: int, variante: str) -> str:
    werte = [str(zahl(x, p, variante)) for x in w]
    return "".join(werte) if p == 2 else " ".join(werte)


def fmt_matrix(m: list[list[int]], p: int, variante: str, einrueckung: str = "    ") -> list[str]:
    zellen = [[str(zahl(x, p, variante)) for x in zeile] for zeile in m]
    breite = max((len(z) for zeile in zellen for z in zeile), default=1)
    return [einrueckung + "( " + " ".join(z.rjust(breite) for z in zeile) + " )" for zeile in zellen]


def fmt_menge(xs: list[int], n: int, variante: str) -> str:
    werte = [x - n if variante == "symmetrisch" and x > n // 2 else x for x in xs]
    return "{" + ", ".join(str(v) for v in werte) + "}"


def text_division(div: Division, p: int, variante: str, name_a: str = "", name_b: str = "",
                  einrueckung: str = "  ") -> list[str]:
    """Polynomdivision Schritt für Schritt."""
    z = []
    fa = fmt_poly(div.dividend, p, variante)
    fb = fmt_poly(div.divisor, p, variante)
    z.append(f"{einrueckung}Polynomdivision in F_{p}[X]: ({fa}) : ({fb})")
    lead_b = monom(div.divisor[-1], grad(div.divisor), p)
    for i, s in enumerate(div.schritte, 1):
        lead_r = monom(s.rest_vorher[-1], grad(s.rest_vorher), p)
        z.append(f"{einrueckung}  Schritt {i}: {fmt_poly(lead_r, p, variante)} : {fmt_poly(lead_b, p, variante)}"
                 f" = {fmt_poly(s.term, p, variante)}")
        z.append(f"{einrueckung}{' ' * (13 + len(str(i)))}{fmt_poly(s.term, p, variante)} · ({fb}) = {fmt_poly(s.produkt, p, variante)}")
        z.append(f"{einrueckung}{' ' * (13 + len(str(i)))}({fmt_poly(s.rest_vorher, p, variante)}) - ({fmt_poly(s.produkt, p, variante)})"
                 f" = {fmt_poly(s.rest_nachher, p, variante)}")
    if not div.schritte:
        z.append(f"{einrueckung}  grad Dividend < grad Divisor: Quotient 0, Rest = Dividend.")
    z.append(f"{einrueckung}  Quotient q(X) = {fmt_poly(div.quotient, p, variante)}, Rest = {fmt_poly(div.rest, p, variante)}")
    rest = "" if not div.rest else f" + {klammer(div.rest, p, variante)}"
    links = name_a or fa
    z.append(f"{einrueckung}  ⇒ {links} = {klammer(div.quotient, p, variante)}·{klammer(div.divisor, p, variante)}{rest}")
    return z


def text_irreduzibel(ir: IrreduzibelErgebnis, variante: str) -> list[str]:
    p, d = ir.p, ir.d
    z = [f"Irreduzible normierte Polynome vom Grad {d} in F_{p}[X] (Ü6 6_1 a)", ""]
    z.append("Irreduzibel heißt: Das Polynom ist kein Produkt von zwei Polynomen kleineren Grades.")
    if d == 1:
        z.append("Jedes Polynom vom Grad 1 ist irreduzibel.")
    else:
        z.append("Die Konstante ist ≠ 0, sonst ist X ein Teiler.")
        if d <= 3:
            z.append(f"Bei Grad {d} ≤ 3 gilt: irreduzibel ⇔ keine Nullstelle in F_{p}.")
        else:
            z.append(f"Bei Grad {d} ≥ 4 reicht das nicht. Das Skript prüft zusätzlich die Teiler vom Grad 2 … {d // 2}.")
    z.append("")
    if len(ir.kandidaten) <= 64:
        z.append("Kandidat | Nullstellen | Ergebnis")
        for kd in ir.kandidaten:
            ns = ", ".join(f"X = {zahl(a, p, variante)}" for a in kd.nullstellen) or "keine"
            if kd.irreduzibel:
                erg = "irreduzibel"
            elif kd.nullstellen:
                erg = "zerlegbar (Nullstelle)"
            else:
                erg = f"zerlegbar (Teiler {fmt_poly(kd.teiler, p, variante)})"
            z.append(f"  {fmt_poly(kd.poly, p, variante)} | {ns} | {erg}")
    else:
        z.append(f"{len(ir.kandidaten)} Kandidaten geprüft (Tabelle zu lang).")
    z.append("")
    z.append(f"Ergebnis: {len(ir.irreduzibel)} irreduzible Polynome vom Grad {d}: "
             + ", ".join(fmt_poly(f, p, variante) for f in ir.irreduzibel))
    return z


def text_produkt(pr: ProduktModulo, n: int, p: int, variante: str) -> list[str]:
    z = [f"Produkt a * b in F_{p}[X]/(X^{n} - 1) (VL5, Abschnitt 5.2)", ""]
    z.append(f"  a = {fmt_wort(poly_zu_wort(pr.a, n), p, variante)} ≅ a(X) = {fmt_poly(pr.a, p, variante, True)}")
    z.append(f"  b = {fmt_wort(poly_zu_wort(pr.b, n), p, variante)} ≅ b(X) = {fmt_poly(pr.b, p, variante, True)}")
    for i, ai in enumerate(pr.a):
        if ai:
            teil = p_mul(monom(ai, i, p), pr.b, p)
            z.append(f"  {fmt_poly(monom(ai, i, p), p, variante)} · b(X) = {fmt_poly(teil, p, variante, True)}")
    z.append(f"  a(X)·b(X) = {fmt_poly(pr.produkt, p, variante, True)}   (in F_{p}[X])")
    z.append(f"  Modulo X^{n} - 1 gilt X^{n} = 1, also X^({n}+j) = X^j:")
    z.append(f"  a * b = {fmt_poly(pr.ergebnis, p, variante, True)}")
    z.append("")
    z.append(f"Ergebnis: a * b = {fmt_wort(poly_zu_wort(pr.ergebnis, n), p, variante)}")
    return z


def text_zerlegung(los: Loesung, variante: str) -> list[str]:
    n, p = los.n, los.p
    zerl = los.zerlegung
    z = []
    # Kreisteilungsklassen
    z.append(f"1) Kreisteilungsklassen von Z_{n} (Bahnen B_s = {{s, s·p, s·p^2, …}} mit p = {p})")
    if los.klassen is None:
        z.append(f"   ggT({n}, {p}) ≠ 1: Die Bahnen schließen sich nicht. Keine Kreisteilungsklassen.")
    else:
        for bahn in los.klassen:
            kette = " → ".join(str(x) for x in bahn) + f" → {bahn[-1]}·{p} = {bahn[-1] * p} ≡ {bahn[0]} (mod {n})"
            z.append(f"   B_{bahn[0]}: {kette}")
        z.append(f"   Z_{n} = " + " ∪ ".join(fmt_menge(b, n, variante) for b in los.klassen))
        laengen = sorted(len(b) for b in los.klassen)
        z.append(f"   Anzahl Klassen = {len(los.klassen)} = Anzahl irreduzibler Faktoren von X^{n} - 1.")
        z.append(f"   Längen der Klassen = Grade der Faktoren: {laengen}")
    z.append("")
    # Zerlegung
    z.append(f"2) Zerlegung von X^{n} - 1 in F_{p}[X] durch Probieren")
    z.append(f"   Nullstellen-Test (lineare Faktoren X - a): a^{n} = 1?")
    for a_, wert in zerl.nullstellen_test:
        z.append(f"     a = {zahl(a_, p, variante)}: a^{n} ≡ {zahl(wert, p, variante)}"
                 + ("  → Nullstelle" if wert == 1 else ""))
    for d_, anz in zerl.getestet.items():
        z.append(f"   Grad {d_}: {anz} normierte Kandidaten (Konstante ≠ 0) als Teiler getestet.")
    z.append("   Jeder gefundene Teiler ist irreduzibel. Denn alle Teiler kleineren Grades sind schon abgespalten.")
    if zerl.sympy_benutzt:
        z.append("   Hinweis: Für große Grade hat sympy die Faktoren gesucht. Die Divisionen stehen trotzdem unten.")
    z.append("")
    for i, s in enumerate(zerl.schritte, 1):
        if s.faktor is zerl.rest_irreduzibel and s.division.quotient == [1]:
            z.append(f"   Rest {fmt_poly(s.faktor, p, variante)}: Alle Teiler vom Grad ≤ {grad(s.faktor) // 2} sind geprüft."
                     " Er hat also keinen echten Teiler und ist irreduzibel.")
            continue
        z.append(f"   Teiler {i}: {fmt_poly(s.faktor, p, variante)}")
        z += text_division(s.division, p, variante, einrueckung="   ")
        z.append("")
    produkt = " · ".join(klammer(f, p, variante) + (f"^{m}" if m > 1 else "") for f, m in zerl.faktoren)
    z.append(f"   X^{n} - 1 = {produkt}")
    # Probe
    pr = [1]
    for f, m in zerl.faktoren:
        for _ in range(m):
            pr = p_mul(pr, f, p)
    ok = pr == x_hoch_n_minus_1(n, p)
    z.append(f"   Probe: Ausmultiplizieren ergibt {fmt_poly(pr, p, variante)} {'✓' if ok else '✗ FEHLER'}")
    if los.klassen is not None:
        grade = sorted(grad(f) for f, m in zerl.faktoren for _ in range(m))
        if grade == sorted(len(b) for b in los.klassen):
            z.append("   Die Grade passen zu den Längen der Kreisteilungsklassen ✓")
    if p == 2 and n == 15:
        z.append("   Quellenhinweis: VL5 (Z. 334) schreibt X^4 + X^3 + X^4 + X + 1. Richtig ist X^4 + X^3 + X^2 + X + 1.")
    return z


def text_generatoren(los: Loesung, variante: str) -> list[str]:
    n, p = los.n, los.p
    z = [f"3) Alle Generatorpolynome (normierte Teiler g(X) von X^{n} - 1, VL5 Abschnitt 5.3)",
         "   r = grad g(X), k = n - r (Dimension), d = Mindestabstand (durch Aufzählen aller Codewörter)"]
    namen = [klammer(f, p, variante) for f, _ in los.zerlegung.faktoren]
    for nr, gen in enumerate(los.generatoren, 1):
        teile = [nm + (f"^{e}" if e > 1 else "") for nm, e in zip(namen, gen.exponenten) if e]
        prod = "·".join(teile) + " = " if teile else ""
        mark = "  ←" if los.k is not None and gen.k == los.k else ""
        dtext = "-" if gen.d is None else str(gen.d)
        z.append(f"   {nr:>2}. r = {gen.r}, k = {gen.k}, d = {dtext}: g(X) = {prod}{fmt_poly(gen.g, p, variante, True)}"
                 f"   Wort {fmt_wort(poly_zu_wort(gen.g, n), p, variante)}{mark}")
    z.append(f"   Anzahl zyklischer Codes der Länge {n} über F_{p}: {len(los.generatoren)}")
    if any(gen.d is None and gen.k > 0 for gen in los.generatoren):
        z.append("   (d = -: zu viele Codewörter zum Aufzählen.)")
    return z


def text_fall(fall: Fall, nr: int, anzahl: int, variante: str) -> list[str]:
    c = fall.code
    n, p = c.n, c.p
    z = []
    kopf = f"{nr}. Fall: " if anzahl > 1 else ""
    z.append(f"{kopf}g(X) = {fmt_poly(c.g, p, variante, True)}   (r = grad g = {c.r}, k = n - r = {n} - {c.r} = {c.k})")
    z.append(f"  Wort g = {fmt_wort(poly_zu_wort(c.g, n), p, variante)}")
    name = BEKANNTE_CODES.get((p, n, tuple(c.g)))
    if name:
        z.append(f"  Das ist {name}.")
    z.append("  Zyklische Generatormatrix G: Zeilen g, X*g, …, X^(k-1)*g (jede Zeile eine Stelle nach rechts)")
    z += fmt_matrix(c.G, p, variante)
    if p == 2 and n == 7 and c.g == [1, 0, 1, 1]:
        z.append("  Quellenhinweis: Ü6 6_2 (2. Fall) hat als letzte Zeile 0001101. Das ist ein Tippfehler.")
        z.append("  Richtig ist X^3 * g = X^3 + X^5 + X^6 → 0001011.")
    z.append("")
    z.append("  Kontrollpolynom h(X) = (X^n - 1) : g(X)  (VL5 Abschnitt 5.5*)")
    z += text_division(c.h_division, p, variante, einrueckung="  ")
    z.append(f"  h(X) = {fmt_poly(c.h, p, variante, True)}   (grad h = k = {c.k})")
    z.append(f"  Reziprokes Polynom: h←(X) = X^k·h(1/X) = {fmt_poly(c.h_rez, p, variante, True)}"
             "   (Koeffizienten in umgekehrter Reihenfolge)")
    z.append("  Kontrollmatrix H: Zeilen h←, X*h←, …, X^(r-1)*h←")
    if c.r > 0:
        z += fmt_matrix(c.H, p, variante)
        ok = all(x == 0 for zeile in c.probe_GHt for x in zeile)
        z.append(f"  Probe: G·H^T = 0 {'✓' if ok else '✗ FEHLER'}")
    else:
        z.append("  r = 0: H hat keine Zeilen (C = F^n).")
    if c.d is not None:
        e = (c.d - 1) // 2
        z.append(f"  Mindestabstand d = {c.d} (kleinstes Gewicht aller {p}^{c.k} - 1 Codewörter ≠ 0), "
                 f"korrigiert e = ⌊(d-1)/2⌋ = {e} Fehler")
        z.append(f"  Parameter: [{n}, {c.k}, {c.d}]_{p}")
    else:
        z.append(f"  Parameter: [{n}, {c.k}]_{p} (d nicht berechnet)")
    if fall.kodierung is not None:
        z += [""] + text_kodierung(fall.kodierung, c, variante)
    if fall.decodierung is not None:
        z += [""] + text_decodierung(fall.decodierung, c, variante)
    return z


BEKANNTE_CODES = {
    (2, 23, (1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1)): "der binäre Golay-Code G_23 aus VL5 (Abschnitt 5.6*)",
    (3, 11, (2, 0, 1, 2, 1, 1)): "der ternäre Golay-Code G_11 aus VL5 (Abschnitt 5.7*)",
    (2, 7, (1, 1, 0, 1)): "der Hamming-Code H_2(3) mit g aus VL6 und Ü6 6_2 (1. Fall)",
    (2, 7, (1, 0, 1, 1)): "ein zu H_2(3) äquivalenter Code (Ü6 6_2, 2. Fall)",
}


def text_kodierung(ko: Kodierung, c: CodeDaten, variante: str) -> list[str]:
    n, p = c.n, c.p
    z = ["  Codieren: c(X) = a(X)·g(X)  (VL5, Algorithmus Kodieren)"]
    z.append(f"    a = {fmt_wort(poly_zu_wort(ko.a, c.k), p, variante)} ≅ a(X) = {fmt_poly(ko.a, p, variante, True)}")
    for i, s in ko.summanden:
        z.append(f"    {fmt_poly(monom(ko.a[i], i, p), p, variante)}·g(X) = {fmt_poly(s, p, variante, True)}")
    z.append(f"    c(X) = {fmt_poly(ko.c, p, variante, True)}")
    z.append(f"    c = {fmt_wort(poly_zu_wort(ko.c, n), p, variante)}  (gleich a·G)")
    return z


def text_decodierung(de: Decodierung, c: CodeDaten, variante: str) -> list[str]:
    n, p = c.n, c.p
    z = ["  Syndrom: y(X) mit Rest durch g(X) teilen (VL6, Syndromdecodierung bei zyklischen Codes)"]
    z.append(f"    y = {fmt_wort(poly_zu_wort(de.y, n), p, variante)} ≅ y(X) = {fmt_poly(de.y, p, variante)}")
    z += text_division(de.division, p, variante, name_a="y(X)", einrueckung="    ")
    z.append(f"    ⇒ s(X) = {fmt_poly(de.s, p, variante)}")
    if not de.s:
        z.append("    s(X) = 0: y ist ein Codewort. Keine Korrektur nötig.")
        return z
    tab = klassenfuehrer_tabelle(c.g, n, p)
    if tab is not None:
        z.append(f"    Klassenführer-Tabelle ({p}^{c.r} = {p ** c.r} Nebenklassen; Syndrom = Klassenführer mod g):")
        z.append("      Klassenführer | Syndrom")
        for e, s in tab:
            mark = "  ←" if s == de.s else ""
            z.append(f"      {fmt_poly(e, p, variante)} | {fmt_poly(s, p, variante)}{mark}")
    if de.abgebrochen:
        z.append("    Die Suche nach dem Klassenführer ist zu groß. Abbruch.")
        return z
    if de.fuehrer is None:
        z.append("    Kein Fehlermuster hat dieses Syndrom. Das Wort ist nicht decodierbar.")
        return z
    e_c = de.abstand
    z.append(f"    Klassenführer a_i(X) = {fmt_poly(de.fuehrer, p, variante)}   (Gewicht {e_c})")
    if de.gleich_gute > 1:
        z.append(f"    Achtung: {de.gleich_gute} Fehlermuster mit Gewicht {e_c} haben dieses Syndrom."
                 " Die Decodierung ist nicht eindeutig. Das Skript nimmt das erste.")
    if c.d is not None and e_c > (c.d - 1) // 2:
        z.append(f"    Achtung: Gewicht {e_c} > e = {(c.d - 1) // 2}. Der Code garantiert diese Korrektur nicht.")
    z.append(f"    ⇒ c(X) = y(X) - a_i(X) = {fmt_poly(de.c, p, variante)}")
    z.append(f"    Probe: c(X) = {klammer(de.c_division.quotient, p, variante)}·{klammer(c.g, p, variante)}"
             f" = q(X)·g(X), Rest {fmt_poly(de.c_division.rest, p, variante)} {'✓' if not de.c_division.rest else '✗'}")
    z.append(f"    y = {fmt_wort(poly_zu_wort(de.y, n), p, variante)}")
    z.append(f"    c = {fmt_wort(poly_zu_wort(de.c, n), p, variante)}     d(y, c) = {e_c}")
    return z


def loesungsweg(los: Loesung, variante: str = STANDARD_VARIANTE) -> str:
    p, n = los.p, los.n
    z = [f"Zyklische Codes über F_{p}" + (f", Wortlänge n = {n}" if n else ""),
         f"Variante '{variante}': {VARIANTEN[variante]}"]
    if n:
        z.append("Annahmen (VL5): Stellen ab 0; Wort a_0 a_1 … ↔ a_0 + a_1·X + … (niedrigster Grad links);"
                 " X * a(X) verschiebt nach rechts.")
    if p == 2:
        z.append("In F_2 gilt -1 = +1. Minus und Plus sind gleich.")
    z.append("")
    for h in los.hinweise:
        z.append(f"HINWEIS: {h}")
    if los.hinweise:
        z.append("")
    ergebnis: list[str] = []

    if los.irreduzibel is not None:
        teil = text_irreduzibel(los.irreduzibel, variante)
        ergebnis.append(teil[-1].removeprefix("Ergebnis: "))
        z += teil[:-1]
    if los.produkt is not None:
        teil = text_produkt(los.produkt, n, p, variante)
        ergebnis.append(teil[-1].removeprefix("Ergebnis: "))
        z += teil[:-1]
    if los.zerlegung is not None:
        z += text_zerlegung(los, variante) + [""]
        z += text_generatoren(los, variante) + [""]
        prod = " · ".join(klammer(f, p, variante) + (f"^{m}" if m > 1 else "") for f, m in los.zerlegung.faktoren)
        ergebnis.append(f"X^{n} - 1 = {prod}")
        if los.klassen is not None:
            ergebnis.append("Kreisteilungsklassen: " + ", ".join(fmt_menge(b, n, variante) for b in los.klassen))
    if los.g_pruefung is not None and los.g_pruefung.rest:
        z.append("Prüfung von g(X): Ist g(X) ein Teiler von X^n - 1?")
        z += text_division(los.g_pruefung, p, variante)
        z.append("")
        ergebnis.append("g(X) ist kein Generatorpolynom (teilt X^n - 1 nicht).")
    if los.faelle:
        z.append("4) Code(s) zum gewählten Generatorpolynom")
        for i, fall in enumerate(los.faelle, 1):
            z += text_fall(fall, i, len(los.faelle), variante) + [""]
            c = fall.code
            dtext = f", {c.d}" if c.d is not None else ""
            kopf = f"{i}. Fall: " if len(los.faelle) > 1 else ""
            ergebnis.append(f"{kopf}g(X) = {fmt_poly(c.g, p, variante, True)}, [{n}, {c.k}{dtext}]_{p}, "
                            f"h(X) = {fmt_poly(c.h, p, variante, True)}, G-Zeile 1 = {fmt_wort(c.G[0], p, variante)}")
            if fall.kodierung is not None:
                ergebnis.append(f"{kopf}c = {fmt_wort(poly_zu_wort(fall.kodierung.c, n), p, variante)}")
            de = fall.decodierung
            if de is not None:
                if not de.s:
                    ergebnis.append(f"{kopf}s(X) = 0, y ist Codewort")
                elif de.c is not None:
                    ergebnis.append(f"{kopf}s(X) = {fmt_poly(de.s, p, variante)}, "
                                    f"c(X) = {fmt_poly(de.c, p, variante)}, c = {fmt_wort(poly_zu_wort(de.c, n), p, variante)}")
                else:
                    ergebnis.append(f"{kopf}s(X) = {fmt_poly(de.s, p, variante)}, nicht decodierbar")
    if los.d_soll is not None and los.k is not None and los.zerlegung is not None:
        passend = [gen for gen in los.generatoren if gen.k == los.k]
        gute = [gen for gen in passend if gen.d is not None and gen.d >= los.d_soll]
        ds = ", ".join("-" if gen.d is None else str(gen.d) for gen in passend) or "keine Codes"
        if gute:
            antwort = (f"Ja. Es gibt einen zyklischen [{n}, {los.k}, ≥{los.d_soll}]_{p}-Code, "
                       f"z. B. g(X) = {fmt_poly(gute[0].g, p, variante, True)}.")
        else:
            antwort = (f"Nein. Die zyklischen [{n}, {los.k}]_{p}-Codes haben d = {ds}. "
                       f"Keiner erreicht d = {los.d_soll}.")
        z.append(f"5) Gibt es einen zyklischen Code mit n = {n}, k = {los.k}, d ≥ {los.d_soll}? (vgl. Ü6 6_7)")
        z.append(f"   {antwort}")
        z.append("")
        ergebnis.append(antwort)

    z.append("Ergebnis:")
    z += [f"  {e}" for e in ergebnis] if ergebnis else ["  (keine Rechnung angefordert)"]
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

BEISPIELE = """Beispiele aus den Quellen:
  X^7 - 1 über F_2 zerlegen (Ü6 6_1 b):
    python -m skripte.zyklischer_code --p 2 --n 7
  Irreduzible Polynome vom Grad 3 (Ü6 6_1 a):
    python -m skripte.zyklischer_code --p 2 --irreduzibel 3
  Hamming-Code H_2(3) als zyklischer Code, beide g (Ü6 6_2):
    python -m skripte.zyklischer_code --p 2 --n 7 --k 4
  Ternäre [4,2]-Codes und Frage nach d = 3 (Ü6 6_6, 6_7):
    python -m skripte.zyklischer_code --p 3 --n 4 --k 2 --d 3
  Syndrom und Decodieren (VL6, y(X) = X^6 + X + 1):
    python -m skripte.zyklischer_code --p 2 --n 7 --g "X^3+X+1" --y "X^6+X+1"
  Produkt in F_5[X]/(X^4 - 1) (VL5: 0012 * 2314 = 2102):
    python -m skripte.zyklischer_code --p 5 --n 4 --produkt 0012 2314 --variante standard
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m skripte.zyklischer_code",
        description="Zyklische Codes über F_p: Kreisteilungsklassen, X^n - 1 zerlegen, g(X), G, h(X), H, "
                    "Codieren, Syndrom und Decodieren. Der Lösungsweg steht Schritt für Schritt da.",
        epilog=BEISPIELE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--p", "--q", dest="p", type=int, required=True,
                        help="Primzahl p. Der Code hat Buchstaben aus F_p (zum Beispiel 2 oder 3).")
    parser.add_argument("--n", type=int, help="Wortlänge n. Das Skript zerlegt X^n - 1.")
    parser.add_argument("--k", type=int,
                        help="Dimension k. Das Skript zeigt alle g(X) vom Grad r = n - k mit G und H.")
    parser.add_argument("--g", help='Generatorpolynom, zum Beispiel "1+X+X^3" oder als Wort "1101".')
    parser.add_argument("--a", help='Klartextwort zum Codieren (k Stellen), zum Beispiel "1001" oder "1+X^3".')
    parser.add_argument("--y", help='Empfangenes Wort zum Decodieren, zum Beispiel "1100001" oder "X^6+X+1".')
    parser.add_argument("--d", type=int,
                        help="Gewünschter Mindestabstand. Frage: Gibt es einen zyklischen [n, k, d]-Code? (Ü6 6_7)")
    parser.add_argument("--irreduzibel", type=int, metavar="GRAD",
                        help="Alle irreduziblen normierten Polynome von diesem Grad auflisten (Ü6 6_1 a).")
    parser.add_argument("--produkt", nargs=2, metavar=("A", "B"),
                        help="Produkt a * b zweier Wörter in F_p[X]/(X^n - 1) berechnen (VL5).")
    parser.add_argument("--variante", choices=sorted(VARIANTEN), default=STANDARD_VARIANTE,
                        help="Schreibweise der Elemente von F_p: symmetrisch (-1 statt 2, Standard, wie VL5/Ü6 6_6) "
                             "oder standard (0 … p-1).")
    parser.add_argument("--alle-varianten", action="store_true",
                        help="Den Lösungsweg in allen Varianten nacheinander ausgeben.")
    args = parser.parse_args(argv)

    try:
        los = loese(args.p, n=args.n, g=args.g, k=args.k, a=args.a, y=args.y, d=args.d,
                    irreduzibel=args.irreduzibel, produkt=tuple(args.produkt) if args.produkt else None)
    except EingabeFehler as fehler:
        print(f"FEHLER in der Eingabe: {fehler}", file=sys.stderr)
        print(f"FEHLER in der Eingabe: {fehler}")
        return 2

    varianten = list(VARIANTEN) if args.alle_varianten else [args.variante]
    for i, v in enumerate(varianten):
        if i:
            print("\n" + "=" * 70 + "\n")
        print(loesungsweg(los, v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
