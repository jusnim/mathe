"""Lineare Codes über F_q: G und H, Codieren, Syndrom, Decodieren, Standardarray.

Zweck
-----
Das Skript rechnet Aufgaben zu linearen Codes wie in VL4 (lineare Codes,
Hamming-Codes), VL6 (Nebenklassen, Syndromdecodierung) und Übung 5/6.
q ist eine Primzahl. Alle Rechnungen laufen modulo q.

Begriffe:
- Erzeugermatrix G: Ihre Zeilen erzeugen den Code. Codieren: c = a·G.
- Kontrollmatrix H: Es gilt H·cᵀ = o genau für Codewörter c.
- Syndrom S(r) = H·rᵀ: Es zeigt, ob r ein Codewort ist.
- Klassenführer: Wort mit kleinstem Gewicht in einer Nebenklasse r + C.

Das Skript zeigt:
- Umrechnung G ↔ H in systematischer Form (mit Zeilenumformungen),
- Parameter [n, k, d]_q, d über den Dualitätssatz,
- Codieren c = a·G,
- Syndrom als Spalte H·rᵀ und als Zeile r·Hᵀ,
- Vergleich mit den Spalten von H: Fehlerstelle x, Fehlergröße y (= e_x),
- korrigiertes Codewort c mit Probe H·cᵀ = o,
- Tabelle Klassenführer/Syndrom und auf Wunsch das Standardarray.

Eingabe von Matrizen und Wörtern
--------------------------------
Zeilen trennt man mit ";". Für q ≤ 10 darf jede Ziffer ein Eintrag sein:
    --H "1011;0112"      --r 2200
Bei q > 10 (oder wenn man will) trennt man Einträge mit Leerzeichen/Komma:
    --H "1 1 1 1 1 1 1 1 1 1; 1 2 3 4 5 6 7 8 9 10"   --r "1,1,5,0,0,0,0,7,3,3"

Aufruf (Beispiele aus den Quellen)
----------------------------------
Ü6 6_8 (H₃(2), Klassenführer, Decodieren):
    python -m skripte.linearer_code --q 3 --typ hamming --m 2 --tabelle --r 2200 --r 0121
Ü5 5_4 (H₂(3) mit dem H der Aufgabe):
    python -m skripte.linearer_code --q 2 --H "1001101;0101011;0010111" --r 1101100 --r 1111111 --r 1111000
Ü5 5_2 (Paritätscode F₃⁵):
    python -m skripte.linearer_code --q 3 --typ paritaet --n 5 --a 2120 --r 22120 --r 11022 --r 11111
VL6 (Nebenklassen des [4,2,2]₂-Codes):
    python -m skripte.linearer_code --q 2 --G "1011;0101" --nebenklassen

Varianten (--variante, --alle-varianten)
----------------------------------------
Die systematische Form legt fest, wie man G und H ineinander umrechnet:
- A: H = (E | A), G = (−Aᵀ | E). Vorlesung VL4, Ü6 6_4, Ü6 6_8.
     Die Nachricht steht in den letzten k Stellen von c.
- B: G = (E | A), H = (−Aᵀ | E). Ü5 5_1, Ü5 5_2.
     Die Nachricht steht in den ersten k Stellen von c.
- C: H = (A | E), G = (E | −Aᵀ). VL6 (Abschnitt 6.3).
     Die Nachricht steht in den ersten k Stellen von c.
Standard: A. Nur bei --typ wiederholung/paritaet ist B Standard (wie Ü5).
Ist H oder G fest vorgegeben, dann rechnet das Skript immer mit dieser
Matrix. Die Variante ändert dann nur die berechnete zweite Matrix.

Weitere Wahl: --fuehrer erste|letzte. Das legt fest, welcher Vektor
Klassenführer wird, wenn mehrere das kleinste Gewicht haben
(VL6 wählt den ersten: 1000, 0100, 0010).
"""

from __future__ import annotations

import argparse
import itertools
import sys
from dataclasses import dataclass, field

from sympy import isprime

Matrix = list[list[int]]
Vektor = list[int]

VARIANTEN = {
    "A": "H = (E | A), G = (−Aᵀ | E)  (VL4, Ü6 6_4, Ü6 6_8)",
    "B": "G = (E | A), H = (−Aᵀ | E)  (Ü5 5_1, Ü5 5_2)",
    "C": "H = (A | E), G = (E | −Aᵀ)  (VL6, Abschnitt 6.3)",
}

# Grenzen, damit das Skript in der Klausur nicht hängen bleibt.
MAX_TEILMENGEN = 2_000_000
MAX_SYNDROME = 100_000
MAX_STANDARDARRAY = 3**7


class EingabeFehler(ValueError):
    """Fehler in der Eingabe. Der Text erklärt das Problem."""


# ---------------------------------------------------------------------------
# Einlesen
# ---------------------------------------------------------------------------

def lies_zeile(text: str, q: int) -> tuple[Vektor, list[str]]:
    """Liest ein Wort wie "2200" oder "1 0 10". Gibt Vektor und Hinweise zurück."""
    text = text.strip()
    if not text:
        raise EingabeFehler("Leeres Wort.")
    if any(z in text for z in " ,\t"):
        teile = [t for t in text.replace(",", " ").split() if t]
    else:
        teile = list(text)
    werte: Vektor = []
    hinweise: list[str] = []
    for t in teile:
        try:
            w = int(t)
        except ValueError as exc:
            raise EingabeFehler(f"'{t}' ist keine ganze Zahl (in '{text}').") from exc
        if not 0 <= w < q:
            hinweise.append(f"Eintrag {w} liegt nicht in F_{q}. Ich rechne mit {w} mod {q} = {w % q}.")
        werte.append(w % q)
    return werte, hinweise


def lies_matrix(text: str, q: int) -> tuple[Matrix, list[str]]:
    """Liest eine Matrix wie "1011;0112". Zeilen trennt ";"."""
    zeilen = [z for z in text.split(";") if z.strip()]
    if not zeilen:
        raise EingabeFehler("Leere Matrix.")
    matrix: Matrix = []
    hinweise: list[str] = []
    for z in zeilen:
        v, h = lies_zeile(z, q)
        matrix.append(v)
        hinweise += h
    laengen = {len(z) for z in matrix}
    if len(laengen) != 1:
        raise EingabeFehler(f"Die Zeilen der Matrix sind verschieden lang: {sorted(laengen)}.")
    return matrix, hinweise


# ---------------------------------------------------------------------------
# Lineare Algebra über F_q
# ---------------------------------------------------------------------------

def transponiert(M: Matrix) -> Matrix:
    return [list(s) for s in zip(*M)] if M else []


def spalte(M: Matrix, j: int) -> Vektor:
    return [z[j] for z in M]


def rang(M: Matrix, q: int) -> int:
    return len(zeilenstufenform(M, q)[1])


def zeilenstufenform(M: Matrix, q: int) -> tuple[Matrix, list[int]]:
    """Reduzierte Zeilenstufenform. Gibt Matrix und Pivotspalten zurück."""
    A = [z[:] for z in M]
    pivots: list[int] = []
    zeile = 0
    n = len(A[0]) if A else 0
    for j in range(n):
        p = next((i for i in range(zeile, len(A)) if A[i][j] % q), None)
        if p is None:
            continue
        A[zeile], A[p] = A[p], A[zeile]
        inv = pow(A[zeile][j], -1, q)
        A[zeile] = [x * inv % q for x in A[zeile]]
        for i in range(len(A)):
            if i != zeile and A[i][j]:
                f = A[i][j]
                A[i] = [(a - f * b) % q for a, b in zip(A[i], A[zeile])]
        pivots.append(j)
        zeile += 1
    return A, pivots


def zeilen_name(i: int) -> str:
    return f"Z{i + 1}"


def bringe_auf_einheit(M: Matrix, spalten: list[int], q: int) -> tuple[Matrix | None, list[str]]:
    """Formt M mit Zeilenumformungen so um, dass in `spalten` die Einheitsmatrix steht.

    Gibt (neue Matrix, Liste der Umformungen) zurück.
    Klappt es nicht, ist die Matrix None.
    """
    A = [z[:] for z in M]
    schritte: list[str] = []
    for i, j in enumerate(spalten):
        p = next((r for r in range(i, len(A)) if A[r][j] % q), None)
        if p is None:
            return None, schritte
        if p != i:
            A[i], A[p] = A[p], A[i]
            schritte.append(f"{zeilen_name(i)} ↔ {zeilen_name(p)}")
        if A[i][j] != 1:
            inv = pow(A[i][j], -1, q)
            A[i] = [x * inv % q for x in A[i]]
            schritte.append(f"{zeilen_name(i)} := {inv}·{zeilen_name(i)}")
        for r in range(len(A)):
            if r != i and A[r][j]:
                f = A[r][j]
                A[r] = [(a - f * b) % q for a, b in zip(A[r], A[i])]
                schritte.append(f"{zeilen_name(r)} := {zeilen_name(r)} − {f}·{zeilen_name(i)}")
    return A, schritte


def kern_basis(M: Matrix, n: int, q: int) -> Matrix:
    """Basis des Kerns {x | M·xᵀ = o} (Zeilen der Ergebnismatrix)."""
    if not M:
        return [[int(i == j) for j in range(n)] for i in range(n)]
    R, piv = zeilenstufenform(M, q)
    frei = [j for j in range(n) if j not in piv]
    basis: Matrix = []
    for f in frei:
        x = [0] * n
        x[f] = 1
        for zi, pj in enumerate(piv):
            x[pj] = (-R[zi][f]) % q
        basis.append(x)
    return basis


def mal_transponiert(H: Matrix, v: Vektor, q: int) -> Vektor:
    """H·vᵀ modulo q."""
    return [sum(h * x for h, x in zip(z, v)) % q for z in H]


def loese_aG(G: Matrix, c: Vektor, q: int) -> Vektor | None:
    """Findet a mit a·G = c (oder None)."""
    k = len(G)
    # Gᵀ·aᵀ = cᵀ: erweiterte Matrix (Gᵀ | cᵀ)
    erw = [row + [ci] for row, ci in zip(transponiert(G), c)]
    R, piv = zeilenstufenform(erw, q)
    if k in piv:
        return None
    a = [0] * k
    for zi, pj in enumerate(piv):
        a[pj] = R[zi][k]
    return a


# ---------------------------------------------------------------------------
# Standard-Codes
# ---------------------------------------------------------------------------

def hamming_spalten(q: int, m: int) -> list[Vektor]:
    """Alle Klassenführer (erste Nicht-Null = 1) in F_q^m, lexikografisch."""
    return [list(v) for v in itertools.product(range(q), repeat=m)
            if any(v) and v[next(i for i, x in enumerate(v) if x)] == 1]


def hamming_H(q: int, m: int, variante: str = "A") -> Matrix:
    """Kontrollmatrix von H_q(m). Variante A: E vorne, sonst E hinten.

    Die übrigen Spalten sind die Klassenführer in lexikografischer Reihenfolge
    (wie Ü5 5_3, Ü5 5_5, Ü6 6_4).
    """
    einheit = [[int(i == j) for i in range(m)] for j in range(m)]
    rest = [v for v in hamming_spalten(q, m) if v not in einheit]
    spalten = einheit + rest if variante == "A" else rest + einheit
    return transponiert(spalten)


def wiederholung_G(q: int, n: int) -> Matrix:
    """G = (1 1 … 1). Passt zu jeder Variante."""
    return [[1] * n]


def paritaet_G(q: int, n: int, variante: str = "B") -> Matrix:
    """Paritätscode: c = (a₁ … a_{n−1}, Σaᵢ). Bei Variante A steht die Prüfstelle vorne."""
    k = n - 1
    if variante == "A":
        return [[1] + [int(i == j) for j in range(k)] for i in range(k)]
    return [[int(i == j) for j in range(k)] + [1] for i in range(k)]


# ---------------------------------------------------------------------------
# Code aufbauen: G ↔ H
# ---------------------------------------------------------------------------

@dataclass
class Umrechnung:
    """Weg von der gegebenen Matrix zur anderen Matrix."""

    gegeben: str                    # "G" oder "H"
    variante: str
    matrix_ein: Matrix
    schritte: list[str]             # Zeilenumformungen
    matrix_sys: Matrix | None       # gegebene Matrix in systematischer Form
    A: Matrix | None                # Block A
    systematisch: bool
    hinweise: list[str] = field(default_factory=list)


@dataclass
class Code:
    q: int
    n: int
    k: int
    G: Matrix                       # zum Codieren
    H: Matrix                       # zum Syndrom
    variante: str
    umrechnung: Umrechnung
    info_stellen: list[int] | None  # 0-basiert, falls G systematisch
    beschreibung: str
    hinweise: list[str] = field(default_factory=list)


def _neg(M: Matrix, q: int) -> Matrix:
    return [[(-x) % q for x in z] for z in M]


def baue_code(q: int, G: Matrix | None = None, H: Matrix | None = None,
              variante: str = "A", beschreibung: str = "") -> Code:
    """Rechnet aus G die Matrix H (oder umgekehrt) in der gewählten systematischen Form."""
    if not isprime(q):
        raise EingabeFehler(f"q = {q} ist keine Primzahl. Das Skript rechnet nur in F_q mit Primzahl q.")
    if (G is None) == (H is None):
        raise EingabeFehler("Bitte genau eine Matrix angeben: G oder H.")
    hinweise: list[str] = []
    gegeben = "G" if G is not None else "H"
    M = G if G is not None else H
    assert M is not None
    n = len(M[0])
    r = rang(M, q)
    if r < len(M):
        hinweise.append(
            f"Die Zeilen von {gegeben} sind linear abhängig (Rang {r} < {len(M)} Zeilen). "
            f"Ich streiche überflüssige Zeilen.")
        R, _ = zeilenstufenform(M, q)
        M = R[:r]
    k = r if gegeben == "G" else n - r
    m = n - k

    # Wo muss die Einheitsmatrix stehen?
    if gegeben == "G":
        ziel = list(range(m, n)) if variante == "A" else list(range(k))
    else:
        ziel = list(range(m)) if variante == "A" else list(range(k, n))

    Msys, schritte = bringe_auf_einheit(M, ziel, q) if M else ([], [])
    info: list[int] | None
    if Msys is not None:
        if gegeben == "G":
            if variante == "A":        # G = (−Aᵀ | E), H = (E | A)
                minus_AT = [z[:m] for z in Msys]
                A = transponiert(_neg(minus_AT, q)) if minus_AT and m else [[] for _ in range(m)]
                Hn = [[int(i == j) for j in range(m)] + A[i] for i in range(m)]
            elif variante == "B":      # G = (E | A), H = (−Aᵀ | E)
                A = [z[k:] for z in Msys]
                Hn = [row + [int(i == j) for j in range(m)]
                      for i, row in enumerate(_neg(transponiert(A), q) or [[] for _ in range(m)])]
            else:                      # G = (E | −Aᵀ), H = (A | E)
                minus_AT = [z[k:] for z in Msys]
                A = transponiert(_neg(minus_AT, q)) or [[] for _ in range(m)]
                Hn = [A[i] + [int(i == j) for j in range(m)] for i in range(m)]
            G_code = M
            H_code = Hn
            info = list(range(m, n)) if variante == "A" else list(range(k))
            if M != Msys:
                info = None   # Codieren mit dem gegebenen G: Nachricht nicht direkt ablesbar
        else:
            if variante == "A":        # H = (E | A), G = (−Aᵀ | E)
                A = [z[m:] for z in Msys]
                Gn = [row + [int(i == j) for j in range(k)]
                      for i, row in enumerate(_neg(transponiert(A), q) or [[] for _ in range(k)])]
                info = list(range(m, n))
            elif variante == "B":      # H = (−Aᵀ | E), G = (E | A)
                minus_AT = [z[:k] for z in Msys]
                A = transponiert(_neg(minus_AT, q)) or [[] for _ in range(k)]
                Gn = [[int(i == j) for j in range(k)] + A[i] for i in range(k)]
                info = list(range(k))
            else:                      # H = (A | E), G = (E | −Aᵀ)
                A = [z[:k] for z in Msys]
                Gn = [[int(i == j) for j in range(k)] + row
                      for i, row in enumerate(_neg(transponiert(A), q) or [[] for _ in range(k)])]
                info = list(range(k))
            G_code = Gn
            H_code = M
        um = Umrechnung(gegeben, variante, M, schritte, Msys, A, True)
    else:
        um = Umrechnung(gegeben, variante, M, schritte, None, None, False)
        um.hinweise.append(
            f"{gegeben} lässt sich ohne Spaltentausch nicht in die Form der Variante {variante} bringen. "
            f"Ich berechne die andere Matrix allgemein als Basis des Kerns.")
        if gegeben == "G":
            G_code, H_code = M, kern_basis(M, n, q)
        else:
            H_code, G_code = M, kern_basis(M, n, q)
        info = None
    return Code(q, n, k, G_code, H_code, variante, um, info, beschreibung, hinweise)


# ---------------------------------------------------------------------------
# Mindestabstand über den Dualitätssatz
# ---------------------------------------------------------------------------

@dataclass
class Abstand:
    s: int                          # je s Spalten von H sind l.u.
    d: int
    la_spalten: list[int] | None    # Beispiel l.a. Spalten (0-basiert)
    min_zeilengewicht: int | None   # min |g_i|
    vollstaendig: bool


def mindestabstand(code: Code) -> Abstand:
    q, H, n = code.q, code.H, code.n
    m = len(H)
    min_g = min((sum(1 for x in z if x) for z in code.G), default=None)
    if m == 0:
        return Abstand(n, 1, None, min_g, True)
    anzahl = 0
    for t in range(1, m + 2):
        for teil in itertools.combinations(range(n), t):
            anzahl += 1
            if anzahl > MAX_TEILMENGEN:
                return Abstand(t - 1, t, None, min_g, False)
            if rang([[H[i][j] for j in teil] for i in range(m)], q) < t:
                return Abstand(t - 1, t, list(teil), min_g, True)
    return Abstand(n, n + 1, None, min_g, True)   # nur bei k = 0


# ---------------------------------------------------------------------------
# Codieren, Syndrom, Decodieren
# ---------------------------------------------------------------------------

@dataclass
class Codierung:
    a: Vektor
    roh: Vektor         # Σ aᵢ·gᵢ ohne mod
    c: Vektor
    probe: Vektor       # H·cᵀ


def codiere(code: Code, a: Vektor) -> Codierung:
    if len(a) != code.k:
        raise EingabeFehler(f"Die Nachricht a = {wort(a)} hat {len(a)} Stellen. Nötig sind k = {code.k}.")
    roh = [sum(ai * g[j] for ai, g in zip(a, code.G)) for j in range(code.n)]
    c = [x % code.q for x in roh]
    return Codierung(a, roh, c, mal_transponiert(code.H, c, code.q))


@dataclass
class Syndrom:
    r: Vektor
    roh: Vektor         # Summe ohne mod
    s: Vektor


def syndrom(code: Code, r: Vektor) -> Syndrom:
    roh = [sum(h * x for h, x in zip(z, r)) for z in code.H]
    return Syndrom(r, roh, [x % code.q for x in roh])


@dataclass
class Fuehrer:
    e: Vektor
    s: Vektor
    eindeutig: bool     # nur ein Vektor mit kleinstem Gewicht


def klassenfuehrer(code: Code, regel: str = "erste") -> list[Fuehrer] | None:
    """Tabelle Klassenführer/Syndrom. Reihenfolge: nach Gewicht, dann Wert, dann Stelle."""
    q, n, m = code.q, code.n, len(code.H)
    if q ** m > MAX_SYNDROME:
        return None
    gefunden: dict[tuple[int, ...], Fuehrer] = {}
    reihenfolge: list[tuple[int, ...]] = []
    for w in range(n + 1):
        kandidaten: dict[tuple[int, ...], list[Vektor]] = {}
        for werte in itertools.product(range(1, q), repeat=w):
            for stellen in itertools.combinations(range(n), w):
                e = [0] * n
                for st, wert in zip(stellen, werte):
                    e[st] = wert
                s = tuple(mal_transponiert(code.H, e, q))
                if s not in gefunden:
                    kandidaten.setdefault(s, []).append(e)
        for s, liste in kandidaten.items():
            wahl = liste[0] if regel == "erste" else liste[-1]
            gefunden[s] = Fuehrer(wahl, list(s), len(liste) == 1)
        # Reihenfolge wie Ü6 6_8: nach Position des Führers in der Aufzählung
        ordnung = []
        for werte in itertools.product(range(1, q), repeat=w):
            for stellen in itertools.combinations(range(n), w):
                e = [0] * n
                for st, wert in zip(stellen, werte):
                    e[st] = wert
                ordnung.append(e)
        index = {tuple(e): i for i, e in enumerate(ordnung)}
        neu = sorted(kandidaten, key=lambda s: index[tuple(gefunden[s].e)])
        reihenfolge += neu
        if len(gefunden) == q ** m:
            break
    return [gefunden[s] for s in reihenfolge]


@dataclass
class Decodierung:
    syn: Syndrom
    ist_codewort: bool
    spalten_treffer: list[tuple[int, int]]  # (x 0-basiert, y) mit S = y·h_x
    fuehrer: Fuehrer | None
    e: Vektor | None
    c: Vektor | None
    probe: Vektor | None
    a: Vektor | None
    t: int                                  # korrigierbare Fehler
    sicher: bool = True                     # Fehlergewicht ≤ t
    hinweise: list[str] = field(default_factory=list)


def decodiere(code: Code, r: Vektor, abstand: Abstand,
              tabelle: list[Fuehrer] | None) -> Decodierung:
    q, n = code.q, code.n
    if len(r) != n:
        raise EingabeFehler(f"Das Wort r = {wort(r)} hat {len(r)} Stellen. Nötig sind n = {n}.")
    syn = syndrom(code, r)
    t = (abstand.d - 1) // 2
    hinweise: list[str] = []
    if not any(syn.s):
        a = loese_aG(code.G, r, q)
        return Decodierung(syn, True, [], None, [0] * n, r[:], syn.s, a, t)
    treffer = []
    for x in range(n):
        h = spalte(code.H, x)
        for y in range(1, q):
            if [y * hi % q for hi in h] == syn.s:
                treffer.append((x, y))
    fuehrer = None
    if tabelle is not None:
        fuehrer = next(f for f in tabelle if f.s == syn.s)
    if len(treffer) == 1 and (fuehrer is None or sum(1 for v in fuehrer.e if v) == 1):
        x, y = treffer[0]
        e = [0] * n
        e[x] = y
    elif fuehrer is not None:
        e = fuehrer.e
        if len(treffer) > 1:
            hinweise.append("Das Syndrom passt zu mehreren Spalten von H. Ein einzelner Fehler ist nicht eindeutig.")
    else:
        hinweise.append("Das Syndrom ist kein Vielfaches einer Spalte von H. "
                        "Es gibt mehr als einen Fehler. Die Tabelle der Klassenführer ist zu groß.")
        return Decodierung(syn, False, treffer, None, None, None, None, None, t, False, hinweise)
    gewicht = sum(1 for v in e if v)
    if gewicht > t:
        hinweise.append(
            f"Der Fehler hat Gewicht {gewicht} > e = {t}. Der Code korrigiert nur {t} Fehler. "
            f"Das Ergebnis ist nicht sicher (Fehler erkannt, aber nicht eindeutig korrigierbar).")
    if fuehrer is not None and not fuehrer.eindeutig:
        hinweise.append("Mehrere Vektoren haben in dieser Nebenklasse das kleinste Gewicht. "
                        "Die Wahl des Klassenführers ist nicht eindeutig (siehe --fuehrer).")
    c = [(ri - ei) % q for ri, ei in zip(r, e)]
    a = loese_aG(code.G, c, q)
    return Decodierung(syn, False, treffer, fuehrer, e, c, mal_transponiert(code.H, c, q), a, t,
                       gewicht <= t, hinweise)


def codewoerter(code: Code) -> list[Vektor]:
    """Alle a·G, a₁ läuft am schnellsten (wie VL6: 0000, 1011, 0101, 1110)."""
    res = []
    for a in itertools.product(range(code.q), repeat=code.k):
        a = list(reversed(a))
        res.append([sum(ai * g[j] for ai, g in zip(a, code.G)) % code.q for j in range(code.n)])
    return res


def standardarray(code: Code, tabelle: list[Fuehrer]) -> list[list[Vektor]]:
    cw = codewoerter(code)
    return [[[(e + c) % code.q for e, c in zip(f.e, cv)] for cv in cw] for f in tabelle]


# Bekannte Fehler in den Quellen (PRIORISIERUNG.md, Abschnitt 5).
def quellen_hinweis(code: Code, r: Vektor, dec: Decodierung) -> str | None:
    if (code.q == 3 and code.H == hamming_H(3, 3, "A")
            and r == [0] * 10 + [1, 2, 2]):
        return ("Achtung Quelle: Ü5 5_5 d) schreibt für r₂ die Stelle x = 13. "
                "Richtig ist x = 11, denn (1 2 0)ᵀ ist die 11. Spalte von H. "
                "Das Ergebnis c₂ = 0…0 2 2 2 in der Quelle stimmt.")
    return None


# ---------------------------------------------------------------------------
# Ausgabe
# ---------------------------------------------------------------------------

def wort(v: Vektor) -> str:
    if any(x > 9 for x in v):
        return " ".join(str(x) for x in v)
    return "".join(str(x) for x in v)


def zahlen(v: Vektor) -> str:
    return " ".join(str(x).replace("-", "−") for x in v)


def spalte_txt(v: Vektor) -> str:
    return "(" + zahlen(v) + ")ᵀ"


def matrix_txt(M: Matrix, name: str, trenn: int | None = None) -> list[str]:
    """Matrix als Textblock. `trenn` setzt einen senkrechten Strich vor diese Spalte."""
    if not M or not M[0]:
        return [f"{name} = ( )  (leere Matrix)"]
    breite = max(len(str(x)) for z in M for x in z)
    zeilen = []
    for z in M:
        teile = [str(x).rjust(breite) for x in z]
        if trenn is not None and 0 < trenn < len(z):
            teile = teile[:trenn] + ["|"] + teile[trenn:]
        zeilen.append(" ".join(teile))
    vor = f"{name} = "
    leer = " " * len(vor)
    out = []
    for i, z in enumerate(zeilen):
        out.append((vor if i == len(zeilen) // 2 else leer) + "( " + z + " )")
    return out


def stelle(i: int) -> str:
    return str(i + 1)


def text_code(code: Code, abstand: Abstand) -> list[str]:
    q, n, k = code.q, code.n, code.k
    m = n - k
    um = code.umrechnung
    L: list[str] = []
    L.append(f"Annahme (Variante {code.variante}): {VARIANTEN[code.variante]}")
    if code.beschreibung:
        L.append(code.beschreibung)
    L.append("")
    L += [f"Hinweis: {h}" for h in code.hinweise + um.hinweise]
    L.append(f"Gegeben: {um.gegeben}")
    L += matrix_txt(um.matrix_ein, um.gegeben)
    if um.schritte:
        L.append("")
        L.append(f"Zeilenumformungen, bis {um.gegeben} systematisch ist:")
        L += ["  " + s for s in um.schritte]
    if um.systematisch and um.matrix_sys is not None:
        andere = "H" if um.gegeben == "G" else "G"
        v = code.variante
        # Position des Strichs in der systematischen Matrix
        if um.gegeben == "G":
            trenn = m if v == "A" else k
        else:
            trenn = m if v == "A" else k
        L.append("")
        L.append(f"Systematische Form von {um.gegeben}:")
        L += matrix_txt(um.matrix_sys, um.gegeben, trenn)
        L.append("")
        L += matrix_txt(um.A or [[]], "A")
        L.append("")
        if um.gegeben == "G":
            L.append({"A": "G = (−Aᵀ | E)  ⇒  H = (E | A):",
                      "B": "G = (E | A)  ⇒  H = (−Aᵀ | E):",
                      "C": "G = (E | −Aᵀ)  ⇒  H = (A | E):"}[v])
            trenn2 = m if v == "A" else k
            L += matrix_txt(code.H, "H", trenn2)
        else:
            L.append({"A": "H = (E | A)  ⇒  G = (−Aᵀ | E):",
                      "B": "H = (−Aᵀ | E)  ⇒  G = (E | A):",
                      "C": "H = (A | E)  ⇒  G = (E | −Aᵀ):"}[v])
            trenn2 = m if v == "A" else k
            L += matrix_txt(code.G, "G", trenn2)
        if um.gegeben == "G" and code.G != um.matrix_sys:
            L.append("Zum Codieren nehme ich das gegebene G.")
        if um.gegeben == "H" and code.H != um.matrix_sys:
            L.append("Für das Syndrom nehme ich das gegebene H.")
        del andere
    else:
        L.append("")
        L += matrix_txt(code.G, "G")
        L += matrix_txt(code.H, "H")
    # Probe H·Gᵀ = O
    HGt = [[sum(h * g for h, g in zip(hz, gz)) % q for gz in code.G] for hz in code.H]
    ok = all(x == 0 for z in HGt for x in z)
    L.append(f"Probe: H·Gᵀ = O  {'✓' if ok else '✗ (FEHLER!)'}")
    if code.info_stellen is not None and k:
        st = code.info_stellen
        L.append(f"Nachricht a steht in c an den Stellen {st[0] + 1} … {st[-1] + 1}.")
    L.append("")

    # Parameter
    L.append("Parameter:")
    L.append(f"  n = {n} (Wortlänge), k = Rang G = {k}, Rang H = n − k = {m}")
    L.append(f"  |C| = q^k = {q}^{k} = {q ** k}")
    L.append("  Dualitätssatz: Sind je s Spalten von H l.u., dann gilt s + 1 ≤ d_C ≤ |c| für alle c ≠ o in C.")
    if abstand.la_spalten is not None:
        sp = ", ".join(f"h{stelle(j)}" for j in abstand.la_spalten)
        L.append(f"  Je {abstand.s} Spalten von H sind l.u. ⇒ s = {abstand.s}.")
        L.append(f"  Die {len(abstand.la_spalten)} Spalten {sp} sind l.a. ⇒ d_C = s + 1 = {abstand.d}.")
    elif not abstand.vollstaendig:
        L.append(f"  Suche abgebrochen (zu viele Teilmengen). Sicher: d_C ≥ {abstand.d}.")
    if abstand.min_zeilengewicht is not None:
        L.append(f"  Kontrolle: {abstand.s} + 1 ≤ d_C ≤ min |g_i| = {abstand.min_zeilengewicht}")
    t = (abstand.d - 1) // 2
    L.append(f"  ⇒ [{n}, {k}, {abstand.d}]_{q}-Code, korrigiert e = ⌊(d − 1)/2⌋ = {t} Fehler, "
             f"erkennt d − 1 = {abstand.d - 1} Fehler.")
    L.append("")
    return L


def text_codierung(code: Code, cod: Codierung) -> list[str]:
    q = code.q
    L = [f"Codieren a = {wort(cod.a)}:"]
    summanden = [f"{ai}·g{i + 1}" for i, ai in enumerate(cod.a) if ai]
    L.append("  c = a·G = " + (" + ".join(summanden) if summanden else "o"))
    for i, ai in enumerate(cod.a):
        if ai:
            L.append(f"      g{i + 1} = {wort(code.G[i])}")
    if cod.roh != cod.c:
        L.append(f"    = ({' '.join(map(str, cod.roh))})  ≡ {wort(cod.c)} (mod {q})")
    L.append(f"  c = {wort(cod.c)}")
    ok = not any(cod.probe)
    L.append(f"  Probe: H·cᵀ = {spalte_txt(cod.probe)} {'= o ✓' if ok else '≠ o ✗'}")
    L.append(f"Ergebnis: a = {wort(cod.a)} wird zu c = {wort(cod.c)}.")
    L.append("")
    return L


def text_syndrom(code: Code, syn: Syndrom) -> list[str]:
    q = code.q
    teile = []
    for j, rj in enumerate(syn.r):
        if rj:
            h = spalte_txt(spalte(code.H, j))
            teile.append(h if rj == 1 else f"{h}·{rj}")
    zeile = "  S(r) = H·rᵀ = " + (" + ".join(teile) if teile else "o")
    L = [zeile]
    if syn.roh != syn.s:
        L.append(f"       = {spalte_txt(syn.roh)} ≡ {spalte_txt(syn.s)} (mod {q})")
    else:
        L.append(f"       = {spalte_txt(syn.s)}")
    L.append(f"  (als Zeile: r·Hᵀ = ({' '.join(map(str, syn.s))}))")
    return L


def text_decodierung(code: Code, dec: Decodierung, nr: int, quelle: str | None) -> list[str]:
    q = code.q
    r = dec.syn.r
    L = [f"Decodieren r{nr} = {wort(r)}:"]
    L += text_syndrom(code, dec.syn)
    if dec.ist_codewort:
        L.append("  S(r) = o ⇒ r ist ein Codewort, c = r.")
        if dec.a is not None:
            L.append(f"  Nachricht: a = {wort(dec.a)}")
        L.append(f"Ergebnis: r{nr} = {wort(r)} ∈ C, c = {wort(r)}.")
        L.append("")
        return L
    L.append("  S(r) ≠ o ⇒ r ist kein Codewort.")
    if dec.spalten_treffer:
        for x, y in dec.spalten_treffer:
            h = spalte_txt(spalte(code.H, x))
            L.append(f"  S(r) = {y}·{h} = {y}·h{stelle(x)}" if y != 1
                     else f"  S(r) = {h} = h{stelle(x)}")
    else:
        L.append("  S(r) ist kein Vielfaches einer Spalte von H ⇒ mehr als ein Fehler.")
    if dec.fuehrer is not None and (len(dec.spalten_treffer) != 1
                                    or sum(1 for v in dec.fuehrer.e if v) != 1):
        L.append(f"  Klassenführer zum Syndrom {' '.join(map(str, dec.syn.s))}: e = {wort(dec.fuehrer.e)}")
    L += [f"  Hinweis: {h}" for h in dec.hinweise]
    if dec.e is None or dec.c is None:
        L.append(f"Ergebnis: r{nr} ist nicht decodierbar (mehr als {dec.t} Fehler).")
        L.append("")
        return L
    fehler = [(j, v) for j, v in enumerate(dec.e) if v]
    if len(fehler) == 1:
        x, y = fehler[0]
        L.append(f"  ⇒ Fehlerstelle x = {x + 1}, Fehlergröße y = e_{x + 1} = {y}")
    L.append(f"  e = {wort(dec.e)}")
    roh = [ri - ei for ri, ei in zip(r, dec.e)]
    if roh != dec.c:
        L.append(f"  c = r − e = ({zahlen(roh)}) ≡ {wort(dec.c)} (mod {q})")
    else:
        L.append(f"  c = r − e = {wort(dec.c)}")
    ok = dec.probe is not None and not any(dec.probe)
    L.append(f"  Probe: H·cᵀ = {spalte_txt(dec.probe or [])} {'= o ✓' if ok else '≠ o ✗'}")
    if dec.a is not None:
        L.append(f"  Nachricht: a = {wort(dec.a)}")
    if quelle:
        L.append(f"  {quelle}")
    if dec.sicher:
        L.append(f"Ergebnis: r{nr} = {wort(r)} wird zu c = {wort(dec.c)} decodiert.")
    else:
        L.append(f"Ergebnis: r{nr} = {wort(r)} ist nicht eindeutig decodierbar (mehr als e = {dec.t} Fehler). "
                 f"Der Klassenführer liefert nur den Vorschlag c = {wort(dec.c)}.")
    L.append("")
    return L


def text_tabelle(code: Code, tabelle: list[Fuehrer] | None, regel: str) -> list[str]:
    if tabelle is None:
        return ["Tabelle der Klassenführer: zu groß (mehr als "
                f"{MAX_SYNDROME} Syndrome).", ""]
    L = [f"Klassenführer und Syndrome (q^(n−k) = {code.q}^{len(code.H)} = {len(tabelle)} Nebenklassen,"
         f" Gleichstand: '{regel}' gewählt):"]
    breite = max(len("Klassenführer"), len(wort(tabelle[0].e)) + 2)
    L.append(f"  {'Klassenführer'.ljust(breite)}  Syndrom")
    L.append("  " + "-" * (breite + 9))
    for f in tabelle:
        mark = "" if f.eindeutig else "   (nicht eindeutig)"
        L.append(f"  {' '.join(map(str, f.e)).ljust(breite)}  {' '.join(map(str, f.s))}{mark}")
    L.append("")
    return L


def text_standardarray(code: Code, tabelle: list[Fuehrer] | None) -> list[str]:
    if tabelle is None or code.q ** code.n > MAX_STANDARDARRAY:
        return ["Standardarray: zu groß für die Ausgabe.", ""]
    arr = standardarray(code, tabelle)
    L = ["Nebenklassen / Standardarray (erste Spalte: Klassenführer, rechts: Syndrom):"]
    for f, zeile in zip(tabelle, arr):
        L.append(f"  {wort(f.e)} + C = {{ {', '.join(wort(v) for v in zeile)} }}   "
                 f"Syn = {' '.join(map(str, f.s))}")
    L.append("")
    return L


# ---------------------------------------------------------------------------
# Gesamtablauf
# ---------------------------------------------------------------------------

@dataclass
class Auftrag:
    q: int
    G: Matrix | None = None
    H: Matrix | None = None
    typ: str | None = None
    n: int | None = None
    m: int | None = None
    nachrichten: list[Vektor] = field(default_factory=list)
    woerter: list[Vektor] = field(default_factory=list)
    tabelle: bool = False
    nebenklassen: bool = False
    fuehrer: str = "erste"


def standard_variante(auftrag: Auftrag) -> str:
    return "B" if auftrag.typ in ("wiederholung", "paritaet") else "A"


def code_aus_auftrag(a: Auftrag, variante: str) -> Code:
    q = a.q
    if not isprime(q):
        raise EingabeFehler(f"q = {q} ist keine Primzahl. Das Skript rechnet nur in F_q mit Primzahl q.")
    if a.typ == "hamming":
        if not a.m or a.m < 2:
            raise EingabeFehler("Für den Hamming-Code bitte --m angeben (m ≥ 2).")
        n = (q ** a.m - 1) // (q - 1)
        text = (f"Hamming-Code H_{q}({a.m}): n = (q^m − 1)/(q − 1) = ({q}^{a.m} − 1)/{q - 1} = {n}, "
                f"k = n − m = {n - a.m}. Spalten von H: Klassenführer (erste Nicht-Null = 1), "
                f"Einheitsvektoren {'vorne' if variante == 'A' else 'hinten'}, sonst lexikografisch.")
        return baue_code(q, H=hamming_H(q, a.m, variante), variante=variante, beschreibung=text)
    if a.typ == "wiederholung":
        if not a.n or a.n < 2:
            raise EingabeFehler("Für den Wiederholungscode bitte --n angeben (n ≥ 2).")
        return baue_code(q, G=wiederholung_G(q, a.n), variante=variante,
                         beschreibung=f"Wiederholungscode in F_{q}^{a.n}: a → (a, a, …, a), G = (1 1 … 1).")
    if a.typ == "paritaet":
        if not a.n or a.n < 2:
            raise EingabeFehler("Für den Paritätscode bitte --n angeben (n ≥ 2).")
        wo = "vorne" if variante == "A" else "hinten"
        return baue_code(q, G=paritaet_G(q, a.n, variante), variante=variante,
                         beschreibung=f"Paritätscode in F_{q}^{a.n}: Prüfstelle Σ aᵢ steht {wo}.")
    return baue_code(q, G=a.G, H=a.H, variante=variante)


def loesungsweg(auftrag: Auftrag, variante: str | None = None) -> str:
    """Erzeugt den ganzen Lösungsweg als Text."""
    v = variante or standard_variante(auftrag)
    code = code_aus_auftrag(auftrag, v)
    abstand = mindestabstand(code)
    L = text_code(code, abstand)
    tab = None
    if auftrag.tabelle or auftrag.nebenklassen or auftrag.woerter:
        tab = klassenfuehrer(code, auftrag.fuehrer)
    if auftrag.tabelle:
        L += text_tabelle(code, tab, auftrag.fuehrer)
    if auftrag.nebenklassen:
        L += text_standardarray(code, tab)
    for a in auftrag.nachrichten:
        L += text_codierung(code, codiere(code, a))
    for i, r in enumerate(auftrag.woerter, 1):
        dec = decodiere(code, r, abstand, tab)
        L += text_decodierung(code, dec, i, quellen_hinweis(code, r, dec))
    if not auftrag.nachrichten and not auftrag.woerter:
        L.append(f"Ergebnis: [{code.n}, {code.k}, {abstand.d}]_{code.q}-Code (G und H siehe oben).")
    return "\n".join(L)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m skripte.linearer_code",
        description="Lineare Codes über F_q: G ↔ H, [n, k, d], Codieren, Syndrom, Decodieren.",
        epilog=(
            "Beispiele aus den Quellen:\n"
            "  Ü6 6_8:  --q 3 --typ hamming --m 2 --tabelle --r 2200 --r 0121\n"
            "  Ü5 5_4:  --q 2 --H \"1001101;0101011;0010111\" --r 1101100 --r 1111111 --r 1111000\n"
            "  Ü5 5_1:  --q 3 --typ wiederholung --n 5 --a 2 --r 22120 --r 11022 --r 11111\n"
            "  VL6:     --q 2 --G \"1011;0101\" --nebenklassen\n"
            "Matrizen: Zeilen mit ';' trennen. Bei q > 10 Einträge mit Leerzeichen oder Komma trennen."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--q", type=int, default=2, help="Primzahl q, Rechnen in F_q (Standard: 2).")
    parser.add_argument("--G", help='Erzeugermatrix G, z. B. "1011;0101".')
    parser.add_argument("--H", help='Kontrollmatrix H, z. B. "1011;0112".')
    parser.add_argument("--typ", choices=["hamming", "wiederholung", "paritaet"],
                        help="Fester Code-Typ statt G oder H.")
    parser.add_argument("--n", type=int, help="Wortlänge n (für Wiederholungs- und Paritätscode).")
    parser.add_argument("--m", type=int, help="m für den Hamming-Code H_q(m).")
    parser.add_argument("--a", action="append", default=[],
                        help="Nachricht a zum Codieren (c = a·G). Mehrfach möglich.")
    parser.add_argument("--r", action="append", default=[],
                        help="Empfangenes Wort r zum Prüfen/Decodieren. Mehrfach möglich.")
    parser.add_argument("--tabelle", action="store_true",
                        help="Tabelle Klassenführer/Syndrom ausgeben.")
    parser.add_argument("--nebenklassen", action="store_true",
                        help="Alle Nebenklassen (Standardarray) ausgeben. Nur für kleine Codes.")
    parser.add_argument("--fuehrer", choices=["erste", "letzte"], default="erste",
                        help="Wahl des Klassenführers bei gleichem Gewicht (Standard: erste, wie VL6).")
    parser.add_argument("--variante", choices=list(VARIANTEN),
                        help="Systematische Form: A = H(E|A), G(−Aᵀ|E) (Standard); "
                             "B = G(E|A), H(−Aᵀ|E) (Standard bei wiederholung/paritaet); "
                             "C = H(A|E), G(E|−Aᵀ).")
    parser.add_argument("--alle-varianten", action="store_true",
                        help="Alle Varianten A, B, C nacheinander ausgeben.")
    args = parser.parse_args(argv)

    try:
        q = args.q
        if not isprime(q):
            raise EingabeFehler(f"q = {q} ist keine Primzahl. Das Skript rechnet nur in F_q mit Primzahl q.")
        anzahl = sum(x is not None for x in (args.G, args.H, args.typ))
        if anzahl != 1:
            raise EingabeFehler("Bitte genau eines angeben: --G, --H oder --typ.")
        hinweise: list[str] = []
        auftrag = Auftrag(q=q, typ=args.typ, n=args.n, m=args.m,
                          tabelle=args.tabelle, nebenklassen=args.nebenklassen, fuehrer=args.fuehrer)
        if args.G:
            auftrag.G, h = lies_matrix(args.G, q)
            hinweise += h
        if args.H:
            auftrag.H, h = lies_matrix(args.H, q)
            hinweise += h
        for t in args.a:
            v, h = lies_zeile(t, q)
            auftrag.nachrichten.append(v)
            hinweise += h
        for t in args.r:
            v, h = lies_zeile(t, q)
            auftrag.woerter.append(v)
            hinweise += h
        for h in hinweise:
            print(f"Hinweis: {h}")
        if args.alle_varianten:
            for v in VARIANTEN:
                print("=" * 70)
                print(f"Variante {v}")
                print("=" * 70)
                print(loesungsweg(auftrag, v))
        else:
            print(loesungsweg(auftrag, args.variante))
    except EingabeFehler as fehler:
        print(f"FEHLER in der Eingabe: {fehler}", file=sys.stderr)
        print(f"FEHLER in der Eingabe: {fehler}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
