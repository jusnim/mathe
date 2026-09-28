"""[10, 8, 3]₁₁-Code (Modulo-11-Code): Kodieren und Decodieren mit Lösungsweg.

Zweck
-----
Das Skript rechnet Aufgaben zum [10, 8, 3]₁₁-Code wie in VL6 (Abschnitt 6.3),
Übung 6_9 und Klausur A5. Es gibt jeden Rechenschritt als Text aus.
Den Text kann man auf Papier abschreiben.

Begriffe:
- Klartextwort a: 8 Ziffern aus ℤ₁₁ = {0, …, 10}.
- Codewort c: 10 Ziffern. c₁ … c₈ = a, dazu die Kontrollziffern c₉ und c₁₀.
- Kontrollmatrix H: Für jedes Codewort gilt H·cᵀ = 0.
- Syndrom S(r) = H·rᵀ: Es zeigt, ob und wo das empfangene Wort r falsch ist.

Aufruf
------
    python -m skripte.mod11_code --a 11500005 --r 1150000511 1150000733

Das ist Klausur A5: H, c₉, c₁₀, c, Syndrome, Decodierung und Probe.
Die Ziffer 10 schreibt man als X (oder man trennt alle Ziffern mit Kommas).

Varianten (--variante)
----------------------
- grundform (Standard, VL6 und Ü6 6_9): H = (1 1 … 1 ; 1 2 … 10).
- grundform-getauscht: Zeile 1 2 … 10 oben, Einsenzeile unten.
- systematisch-rechts: H = (A | E), Einheitsmatrix rechts (VL6 Z. 224).
- systematisch-links: H = (E | A), Einheitsmatrix links (VL4, Ü6 6_4).
Alle vier Matrizen beschreiben denselben Code. Die Syndrome sind aber
verschieden. --alle-varianten gibt alle vier Rechenwege aus.

Weitere Optionen: --repraesentanten (0 … 10 oder −5 … 5) und
--namen (s₁, s₂, eₓ wie VL6 oder s₀, s₁, y wie Ü6).
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field

P = 11
N = 10
K = 8

VARIANTEN = {
    "grundform": "H in Grundform, Einsenzeile oben (VL6, Ü6 6_9). Standard.",
    "grundform-getauscht": "H in Grundform, Zeile 1 2 … 10 oben, Einsenzeile unten.",
    "systematisch-rechts": "H = (A | E) mit Einheitsmatrix rechts (VL6 Z. 224).",
    "systematisch-links": "H = (E | A) mit Einheitsmatrix links (VL4, Ü6 6_4).",
}

NAMEN = {
    # Name: (erstes Syndrom, zweites Syndrom, Fehlergröße)
    "vl": ("s₁", "s₂", "eₓ"),
    "uebung": ("s₀", "s₁", "y"),
}

_TIEF = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def tief(i: int) -> str:
    """Schreibt eine Zahl als Index (tiefgestellt)."""
    return str(i).translate(_TIEF)


# ---------------------------------------------------------------------------
# Rechnen in ℤ₁₁
# ---------------------------------------------------------------------------

def inverse_mod(a: int, m: int = P) -> int:
    """Gibt a⁻¹ mod m zurück. Rechnet mit dem erweiterten Euklidischen Algorithmus.

    Das Inverse existiert nur, wenn ggT(a, m) = 1 ist.
    """
    r0, r1 = a % m, m
    x0, x1 = 1, 0
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        x0, x1 = x1, x0 - q * x1
    if r0 != 1:
        raise ValueError(f"{a} hat kein Inverses modulo {m}, denn ggT = {r0}.")
    return x0 % m


def inversentabelle() -> dict[int, int]:
    """Tabelle a → a⁻¹ mod 11 für a = 1 … 10."""
    return {a: inverse_mod(a) for a in range(1, P)}


def rep(v: int, symmetrisch: bool = False) -> int:
    """Wählt den Vertreter der Restklasse: 0 … 10 oder −5 … 5."""
    v %= P
    if symmetrisch and v > P // 2:
        v -= P
    return v


# ---------------------------------------------------------------------------
# Kontrollmatrix H
# ---------------------------------------------------------------------------

Matrix = list[list[int]]


@dataclass
class Kontrollmatrix:
    variante: str
    H: Matrix
    grundform: Matrix
    schritte: list[tuple[str, Matrix]] = field(default_factory=list)
    G: Matrix | None = None
    info_stellen: list[int] = field(default_factory=list)  # Stellen 1 … 10 der Klartextziffern


def h_grundform() -> Matrix:
    return [[1] * N, list(range(1, N + 1))]


def _gauss_jordan(M: Matrix, pivots: list[int]) -> tuple[Matrix, list[tuple[str, Matrix]]]:
    """Formt M mit Zeilenumformungen so um, dass die Spalten `pivots` die Einheitsmatrix bilden."""
    M = [row[:] for row in M]
    schritte: list[tuple[str, Matrix]] = []
    for k, col in enumerate(pivots):
        r = next(i for i in range(k, len(M)) if M[i][col] % P)
        if r != k:
            M[k], M[r] = M[r], M[k]
            schritte.append((f"Z{tief(k + 1)} ↔ Z{tief(r + 1)}", [row[:] for row in M]))
        f = M[k][col] % P
        if f != 1:
            inv = inverse_mod(f)
            M[k] = [(inv * v) % P for v in M[k]]
            schritte.append((f"Z{tief(k + 1)} := {inv}·Z{tief(k + 1)}  (denn {f}⁻¹ = {inv})",
                             [row[:] for row in M]))
        for i in range(len(M)):
            if i != k and M[i][col] % P:
                g = M[i][col] % P
                M[i] = [(a - g * b) % P for a, b in zip(M[i], M[k])]
                schritte.append((f"Z{tief(i + 1)} := Z{tief(i + 1)} − {g}·Z{tief(k + 1)}",
                                 [row[:] for row in M]))
    return M, schritte


def kontrollmatrix(variante: str = "grundform") -> Kontrollmatrix:
    """Baut H für die gewählte Variante. Bei systematischer Form auch G und die Umformungen."""
    if variante not in VARIANTEN:
        raise ValueError(f"Unbekannte Variante: {variante}")
    g = h_grundform()
    if variante == "grundform":
        return Kontrollmatrix(variante, g, g, info_stellen=list(range(1, K + 1)))
    if variante == "grundform-getauscht":
        return Kontrollmatrix(variante, [g[1][:], g[0][:]], g, info_stellen=list(range(1, K + 1)))
    if variante == "systematisch-rechts":
        H, schritte = _gauss_jordan(g, [N - 2, N - 1])
        A = [row[:K] for row in H]
        G = [[1 if j == i else 0 for j in range(K)] + [(-A[0][i]) % P, (-A[1][i]) % P]
             for i in range(K)]
        return Kontrollmatrix(variante, H, g, schritte, G, list(range(1, K + 1)))
    # systematisch-links
    H, schritte = _gauss_jordan(g, [0, 1])
    A = [row[2:] for row in H]
    G = [[(-A[0][i]) % P, (-A[1][i]) % P] + [1 if j == i else 0 for j in range(K)]
         for i in range(K)]
    return Kontrollmatrix(variante, H, g, schritte, G, list(range(3, N + 1)))


# ---------------------------------------------------------------------------
# Kodieren
# ---------------------------------------------------------------------------

@dataclass
class Kodierung:
    a: list[int]
    c9_terme: list[tuple[int, int]]   # (Gewicht, Ziffer)
    c9_summe: int
    c9: int
    c10_terme: list[tuple[int, int]]
    c10_summe: int
    c10: int
    c: list[int]
    c_ueber_G: list[int]
    G: Matrix


def kodieren(a: list[int]) -> Kodierung:
    """Berechnet c₉, c₁₀ und c = a·G (VL6 Z. 246: c₉ ≡ Σ(1+i)cᵢ, c₁₀ ≡ Σ(9−i)cᵢ)."""
    if len(a) != K:
        raise ValueError(f"Das Klartextwort a braucht {K} Ziffern, es hat {len(a)}.")
    t9 = [(1 + i, a[i - 1]) for i in range(1, K + 1)]
    t10 = [(9 - i, a[i - 1]) for i in range(1, K + 1)]
    s9 = sum(w * z for w, z in t9)
    s10 = sum(w * z for w, z in t10)
    c = a + [s9 % P, s10 % P]
    G = kontrollmatrix("systematisch-rechts").G
    cG = [sum(a[i] * G[i][j] for i in range(K)) % P for j in range(N)]
    return Kodierung(a, t9, s9, s9 % P, t10, s10, s10 % P, c, cG, G)


# ---------------------------------------------------------------------------
# Syndrom und Decodieren
# ---------------------------------------------------------------------------

@dataclass
class Decodierung:
    r: list[int]
    variante: str
    H: Matrix
    zeilen_terme: list[list[tuple[int, int]]]  # je Zeile von H: (Eintrag, Ziffer), ohne Nullen
    summen: list[int]
    S: tuple[int, int]
    status: str                  # "codewort", "ein_fehler", "nicht_decodierbar"
    grund: str
    x: int | None = None         # Fehlerstelle 1 … 10
    e: int | None = None         # Fehlergröße
    c: list[int] | None = None
    probe_summen: list[int] | None = None
    hinweis_quelle: str = ""


def syndrom(r: list[int], H: Matrix) -> tuple[list[list[tuple[int, int]]], list[int], tuple[int, int]]:
    terme = [[(h, z) for h, z in zip(zeile, r) if h % P and z % P] for zeile in H]
    summen = [sum(h * z for h, z in zeile) for zeile in terme]
    return terme, summen, (summen[0] % P, summen[1] % P)


def _spalte_suchen(S: tuple[int, int], H: Matrix) -> tuple[int, int] | None:
    """Sucht die Spalte hⱼ von H mit S = y·hⱼ (y ≠ 0). Gibt (j, y) zurück, j ab 1."""
    for j in range(N):
        h = (H[0][j] % P, H[1][j] % P)
        for y in range(1, P):
            if ((y * h[0]) % P, (y * h[1]) % P) == S:
                return j + 1, y
    return None


def decodieren(r: list[int], variante: str = "grundform") -> Decodierung:
    """Berechnet S(r) = H·rᵀ und decodiert r, falls höchstens ein Fehler vorliegt."""
    if len(r) != N:
        raise ValueError(f"Das empfangene Wort braucht {N} Ziffern, es hat {len(r)}.")
    km = kontrollmatrix(variante)
    H = km.H
    terme, summen, S = syndrom(r, H)
    d = Decodierung(r, variante, H, terme, summen, S, "", "")

    if S == (0, 0):
        d.status = "codewort"
        d.grund = "S(r) = 0. Also ist r ein Codewort: c = r."
        d.c = r[:]
    elif variante in ("grundform", "grundform-getauscht"):
        # e = Eintrag zur Einsenzeile, x·e = Eintrag zur Zeile 1 … 10
        e_idx, xe_idx = (0, 1) if variante == "grundform" else (1, 0)
        e, xe = S[e_idx], S[xe_idx]
        if e == 0:
            d.status = "nicht_decodierbar"
            d.grund = ("Die Fehlergröße wäre 0, aber S ≠ 0. "
                       "Ein einzelner Fehler hat Fehlergröße ≠ 0. "
                       "Also hat r mehr als einen Fehler (VL6 Z. 279).")
        else:
            x = (xe * inverse_mod(e)) % P
            if x == 0:
                d.status = "nicht_decodierbar"
                d.grund = (f"Die Fehlerstelle wäre x = {xe}/{e} ≡ 0. "
                           "Die Stellen sind aber 1 … 10. x = 0 ist keine Stelle. "
                           "Also hat r mindestens 2 Fehler.")
                d.hinweis_quelle = ("Die Vorlesung (VL6 Z. 279) nennt nur den Fall "
                                    "„Fehlergröße = 0, zweiter Eintrag ≠ 0“. Hier liegt der andere "
                                    "Fall vor: zweiter Eintrag = 0, Fehlergröße ≠ 0. "
                                    "Auch er bedeutet: nicht decodierbar. Das muss man selbst begründen.")
                d.x, d.e = 0, e
            else:
                d.status = "ein_fehler"
                d.x, d.e = x, e
    else:
        treffer = _spalte_suchen(S, H)
        if treffer is None:
            d.status = "nicht_decodierbar"
            d.grund = ("S ist kein Vielfaches einer Spalte von H. "
                       "Ein einzelner Fehler an Stelle x gibt S = eₓ·hₓ. "
                       "Also hat r mindestens 2 Fehler.")
        else:
            d.status = "ein_fehler"
            d.x, d.e = treffer

    if d.status == "ein_fehler":
        c = r[:]
        c[d.x - 1] = (r[d.x - 1] - d.e) % P
        d.c = c
        d.grund = f"Ein Fehler an Stelle x = {d.x} mit Fehlergröße {d.e}."
    if d.c is not None:
        d.probe_summen = [sum(h * z for h, z in zip(zeile, d.c)) for zeile in H]
    return d


# ---------------------------------------------------------------------------
# Eingabe
# ---------------------------------------------------------------------------

def ziffern_lesen(text: str) -> list[int]:
    """Liest ein Wort. '1150000546' oder '1,1,5,…' oder mit X für 10."""
    t = text.strip()
    if "," in t or " " in t:
        teile = [s for s in t.replace(",", " ").split() if s]
    else:
        teile = list(t)
    werte = []
    for s in teile:
        if s.upper() in ("X", "A"):
            werte.append(10)
        elif s.lstrip("-").isdigit():
            werte.append(int(s) % P)
        else:
            raise ValueError(f"„{s}“ ist keine Ziffer in ℤ₁₁.")
    return werte


def null_einfuegen(r: list[int], stelle: int) -> list[int]:
    """Fügt eine 0 so ein, dass sie danach an Stelle `stelle` (1 … 10) steht."""
    return r[:stelle - 1] + [0] + r[stelle - 1:]


def stelle_in_nullblock(r: list[int]) -> int:
    """Stelle für die fehlende 0: am Ende des längsten Blocks aus Nullen."""
    beste, laenge, start, i = 1, 0, None, 0
    for i in range(len(r) + 1):
        if i < len(r) and r[i] == 0:
            if start is None:
                start = i
        else:
            if start is not None and i - start > laenge:
                laenge, beste = i - start, i + 1
            start = None
    return beste


# ---------------------------------------------------------------------------
# Ausgabe
# ---------------------------------------------------------------------------

def wort(v: list[int], sym: bool = False) -> str:
    """Schreibt ein Wort. Die Ziffern bleiben immer 0 … 10 (wie in den Aufgaben)."""
    return " ".join(str(z) for z in v)


def kl(v: int) -> str:
    """Setzt negative Zahlen in Klammern."""
    return f"({v})" if v < 0 else str(v)


def matrix_text(M: Matrix, sym: bool = False, einzug: str = "    ") -> str:
    werte = [[str(rep(v, sym)) for v in row] for row in M]
    breite = max(len(w) for row in werte for w in row)
    zeilen = []
    for i, row in enumerate(werte):
        links, rechts = ("⎛", "⎞") if i == 0 else (("⎝", "⎠") if i == len(werte) - 1 else ("⎜", "⎟"))
        if len(werte) == 1:
            links, rechts = "(", ")"
        zeilen.append(einzug + links + " " + " ".join(w.rjust(breite) for w in row) + " " + rechts)
    return "\n".join(zeilen)


def summe_text(terme: list[tuple[int, int]], sym: bool = False) -> str:
    """Schreibt Σ hᵢ·rᵢ aus und hängt den Wert der Summe an: '1 + 2·1 + … = 77'."""
    if not terme:
        return "0 = 0"
    teile, summe = [], 0
    for h, z in terme:
        h = rep(h, sym)
        summe += h * z
        teile.append(str(z) if h == 1 else f"{kl(h)}·{z}")
    return " + ".join(teile) + f" = {summe}"


def text_inversentabelle(sym: bool) -> str:
    tab = inversentabelle()
    reihenfolge = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    a = [str(rep(v, sym)) for v in reihenfolge]
    inv = [str(rep(tab[v], sym)) for v in reihenfolge]
    b = max(len(s) for s in a + inv)
    return ("    a    | " + " ".join(s.rjust(b) for s in a) + "\n"
            + "    " + "-" * (7 + (b + 1) * len(a)) + "\n"
            + "    a⁻¹  | " + " ".join(s.rjust(b) for s in inv))


def text_kontrollmatrix(km: Kontrollmatrix, sym: bool) -> str:
    z = [f"Kontrollmatrix H ({VARIANTEN[km.variante]})", ""]
    if km.variante != "grundform":
        z.append("Grundform:")
        z.append(matrix_text(km.grundform, sym))
    if km.variante == "grundform-getauscht":
        z.append("Zeilen getauscht (Zeile 1 2 … 10 oben):")
        z.append("H =")
        z.append(matrix_text(km.H, sym))
    elif km.variante.startswith("systematisch"):
        seite = "rechts" if km.variante.endswith("rechts") else "links"
        z.append(f"Zeilenumformungen, bis die Einheitsmatrix E {seite} steht:")
        for text, M in km.schritte:
            z.append(f"  {text}:")
            z.append(matrix_text(M, sym))
        if seite == "rechts":
            z.append("H = (A | E) =")
            z.append(matrix_text(km.H, sym))
            z.append("Generatormatrix G = (E | −Aᵀ) =")
        else:
            z.append("H = (E | A) =")
            z.append(matrix_text(km.H, sym))
            z.append("Generatormatrix G = (−Aᵀ | E) =")
            z.append("  (Achtung: Mit diesem G stehen die Klartextziffern an den Stellen 3 … 10.)")
        z.append(matrix_text(km.G, sym))
    else:
        z.append("H =")
        z.append(matrix_text(km.H, sym))
    z.append("Alle Formen von H beschreiben denselben Code. Die Syndrome hängen aber von H ab.")
    return "\n".join(z)


def text_kodierung(k: Kodierung, sym: bool) -> str:
    z = ["Kodierung", ""]
    z.append(f"a = {wort(k.a, sym)}.  Es gilt cᵢ = aᵢ für i = 1 … 8.")
    z.append("Formeln (VL6): c₉ ≡ Σ (1+i)·cᵢ,  c₁₀ ≡ Σ (9−i)·cᵢ  (i = 1 … 8, mod 11)")
    for name, terme, summe, wert in (("c₉", k.c9_terme, k.c9_summe, k.c9),
                                     ("c₁₀", k.c10_terme, k.c10_summe, k.c10)):
        ausgeschrieben = " + ".join(f"{w}·{rep(zf, sym)}" for w, zf in terme)
        z.append(f"{name} ≡ {ausgeschrieben}")
        z.append(f"    = {summe} ≡ {rep(wert, sym)} (mod 11)")
    z.append(f"c = {wort(k.c, sym)}")
    z.append("")
    z.append("Kontrolle mit c = a·G, G = (E | −Aᵀ) aus H = (A | E):")
    z.append(matrix_text(k.G, sym))
    gleich = "stimmt überein" if k.c == k.c_ueber_G else "WEICHT AB"
    z.append(f"a·G = {wort(k.c_ueber_G, sym)}  ({gleich})")
    return "\n".join(z)


def text_decodierung(d: Decodierung, name: str, namen: str, sym: bool) -> str:
    n1, n2, ne = NAMEN[namen]
    z = [f"Syndrom S({name}) = H·{name}ᵀ", ""]
    z.append(f"{name} = {wort(d.r, sym)}")
    s_namen = (n1, n2)
    for i in range(2):
        z.append(f"  {s_namen[i]} = {summe_text(d.zeilen_terme[i], sym)}"
                 f" ≡ {rep(d.S[i], sym)} (mod 11)")
    z.append(f"S({name}) = ({rep(d.S[0], sym)}, {rep(d.S[1], sym)})ᵀ")
    z.append("")

    grund = d.variante in ("grundform", "grundform-getauscht")
    if d.status == "codewort":
        z.append(d.grund)
    elif grund:
        oben = d.variante == "grundform"
        s_e, s_xe = (n1, n2) if oben else (n2, n1)
        e_wert = d.S[0] if oben else d.S[1]
        xe_wert = d.S[1] if oben else d.S[0]
        if oben:
            z.append(f"Ansatz: S = (1, x)ᵀ·{ne}. Also {ne} = {n1} und x = {n2}/{n1}.")
        else:
            z.append(f"Ansatz: S = (x, 1)ᵀ·{ne}. Also {ne} = {n2} und x = {n1}/{n2}.")
        if e_wert == 0:
            z.append(f"{s_e} = 0, aber {s_xe} = {rep(xe_wert, sym)} ≠ 0.")
        else:
            inv = inverse_mod(e_wert)
            x = (xe_wert * inv) % P
            a_, b_, i_ = rep(xe_wert, sym), rep(e_wert, sym), rep(inv, sym)
            z.append(f"x = {s_xe}/{s_e} = {a_}/{kl(b_)} = {a_}·{kl(b_)}⁻¹"
                     f" = {a_}·{kl(i_)} = {a_ * i_} ≡ {x} (mod 11)")
            z.append(f"  (Inverse aus der Tabelle: {b_}·{kl(i_)} = {b_ * i_} ≡ 1 mod 11)")
        if d.status == "nicht_decodierbar":
            z.append(f"=> {name} ist NICHT decodierbar. {d.grund}")
        else:
            z.append(f"=> Fehlerstelle x = {d.x}, Fehlergröße {ne} = {rep(d.e, sym)}.")
    else:
        z.append("Ansatz: Ein Fehler an Stelle x gibt S = eₓ·hₓ (hₓ = Spalte x von H).")
        if d.status == "nicht_decodierbar":
            z.append(f"=> {name} ist NICHT decodierbar. {d.grund}")
        else:
            h = (d.H[0][d.x - 1], d.H[1][d.x - 1])
            z.append(f"Spalte h{tief(d.x)} = ({rep(h[0], sym)}, {rep(h[1], sym)})ᵀ,"
                     f"  {rep(d.e, sym)}·h{tief(d.x)} = ({rep(d.e * h[0], sym)}, {rep(d.e * h[1], sym)})ᵀ = S")
            z.append(f"=> Fehlerstelle x = {d.x}, Fehlergröße {ne} = {rep(d.e, sym)}.")
    if d.hinweis_quelle:
        z.append(f"Hinweis zur Quelle: {d.hinweis_quelle}")

    if d.status == "ein_fehler":
        x = d.x
        z.append(f"c{tief(x)} = r{tief(x)} − {ne} = {d.r[x - 1]} − {kl(rep(d.e, sym))}"
                 f" = {d.r[x - 1] - rep(d.e, sym)} ≡ {d.c[x - 1]} (mod 11)")
        z.append(f"c = {name} − e = {wort(d.c, sym)}")
    if d.c is not None:
        z.append("")
        z.append("Probe H·cᵀ = 0:")
        pterme = [[(h, zf) for h, zf in zip(zeile, d.c) if h % P and zf % P] for zeile in d.H]
        for i in range(2):
            z.append(f"  Zeile {i + 1}: {summe_text(pterme[i], sym)}"
                     f" ≡ 0 (mod 11)" if d.probe_summen[i] % P == 0 else
                     f"  Zeile {i + 1}: {summe_text(pterme[i], sym)} ≢ 0 (mod 11)")
        ok = all(s % P == 0 for s in d.probe_summen)
        z.append("  => H·cᵀ = 0. c ist ein Codewort." if ok else "  => FEHLER: c ist kein Codewort.")
    return "\n".join(z)


def loesungsweg(a: list[int] | None, woerter: list[tuple[str, list[int]]], variante: str,
                namen: str = "vl", sym: bool = False) -> str:
    """Erzeugt den ganzen Lösungsweg als Text für eine Variante."""
    teile = [f"=== Variante: {variante} ===",
             f"Annahme: {VARIANTEN[variante]}",
             "Stellen werden ab 1 gezählt (c₁ … c₁₀), wie in VL6, Ü6 6_9 und der Klausur.",
             f"Vertreter der Restklassen: {'−5 … 5' if sym else '0 … 10'}.",
             ""]
    km = kontrollmatrix(variante)
    teile.append(text_kontrollmatrix(km, sym))
    if a is not None:
        teile += ["", text_kodierung(kodieren(a), sym)]
    ergebnisse = []
    if woerter:
        teile += ["", "Inverse modulo 11 (für die Division):", text_inversentabelle(sym)]
        for name, r in woerter:
            d = decodieren(r, variante)
            teile += ["", text_decodierung(d, name, namen, sym)]
            ergebnisse.append((name, d))
    teile += ["", "Ergebnis:"]
    if a is not None:
        k = kodieren(a)
        teile.append(f"  c₉ = {rep(k.c9, sym)}, c₁₀ = {rep(k.c10, sym)}, c = {wort(k.c, sym)}")
    for name, d in ergebnisse:
        S = f"S({name}) = ({rep(d.S[0], sym)}, {rep(d.S[1], sym)})ᵀ"
        if d.status == "codewort":
            teile.append(f"  {S}: Codewort, c = {wort(d.c, sym)}")
        elif d.status == "ein_fehler":
            teile.append(f"  {S}: x = {d.x}, Fehlergröße {rep(d.e, sym)}, c = {wort(d.c, sym)}")
        else:
            teile.append(f"  {S}: nicht decodierbar (mindestens 2 Fehler)")
    return "\n".join(teile)


def text_luecke(name: str, r9: list[int], stelle: int) -> str:
    """Hinweis bei 9 Ziffern: welche 0 ergänzt wird und was andere Stellen ergeben."""
    z = [f"WARNUNG: {name} = {''.join(map(str, r9))} hat nur 9 Ziffern. Der Code hat Länge 10.",
         "Wahrscheinlich fehlt eine 0 (GP A4, vergleiche Klausur A5).",
         f"Annahme: Die 0 steht an Stelle {stelle} (im längsten Block aus Nullen).",
         "Andere mögliche Stellen für die 0 (Grundform von H):"]
    gruppen: dict[tuple[int, ...], list[int]] = {}
    for s in range(1, N + 1):
        gruppen.setdefault(tuple(null_einfuegen(r9, s)), []).append(s)
    for r, stellen in gruppen.items():
        d = decodieren(list(r), "grundform")
        if d.status == "codewort":
            info = "Codewort"
        elif d.status == "ein_fehler":
            info = f"x = {d.x}, Fehlergröße {d.e}"
        else:
            info = "nicht decodierbar"
        st = str(stellen[0]) if len(stellen) == 1 else f"{stellen[0]} … {stellen[-1]}"
        markierung = "  <- Annahme" if stelle in stellen else ""
        z.append(f"  0 an Stelle {st:>6}: {''.join(str(v) if v < 10 else 'X' for v in r)}"
                 f"  S = ({d.S[0]}, {d.S[1]})ᵀ  {info}{markierung}")
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m skripte.mod11_code",
        description="[10, 8, 3]₁₁-Code: H aufstellen, c₉ und c₁₀ berechnen, "
                    "Syndrome berechnen, decodieren, Probe.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Beispiel (Klausur A5):\n"
               "  python -m skripte.mod11_code --a 11500005 --r 1150000511 1150000733\n"
               "Beispiel (Ü6 6_9, Schreibweise der Übung):\n"
               "  python -m skripte.mod11_code --r 1151111103 1156111103 1111111103 "
               "--namen uebung --repraesentanten symmetrisch\n\n"
               "Varianten:\n" + "\n".join(f"  {k}: {v}" for k, v in VARIANTEN.items()))
    parser.add_argument("--a", help="Klartextwort mit 8 Ziffern, z. B. 11500005. "
                                    "Die Ziffer 10 als X schreiben.")
    parser.add_argument("--r", nargs="+", action="extend", default=[],
                        help="Ein oder mehrere empfangene Wörter mit 10 Ziffern, z. B. 1150000511.")
    parser.add_argument("--variante", choices=list(VARIANTEN), default="grundform",
                        help="Form der Kontrollmatrix H. Standard: grundform (wie VL6 und Ü6).")
    parser.add_argument("--alle-varianten", action="store_true",
                        help="Gibt den Lösungsweg für alle Formen von H nacheinander aus.")
    parser.add_argument("--repraesentanten", choices=["0-10", "symmetrisch"], default="0-10",
                        help="Vertreter der Restklassen: 0 … 10 (Standard) oder −5 … 5 (wie Ü6 6_9).")
    parser.add_argument("--namen", choices=list(NAMEN), default="vl",
                        help="Namen im Syndrom: vl = s₁, s₂, eₓ (VL6, Standard); "
                             "uebung = s₀, s₁, y (Ü6 6_9).")
    parser.add_argument("--luecke", type=int, metavar="STELLE",
                        help="Nur bei Wörtern mit 9 Ziffern: Stelle (1 … 10), an der die fehlende 0 "
                             "ergänzt wird. Ohne Angabe: im längsten Block aus Nullen.")
    args = parser.parse_args(argv)

    if args.a is None and not args.r:
        parser.print_help()
        print("\nHinweis: Bitte --a oder --r angeben.")
        return 1

    sym = args.repraesentanten == "symmetrisch"
    a = None
    if args.a is not None:
        try:
            a = ziffern_lesen(args.a)
        except ValueError as fehler:
            print(f"FEHLER bei --a: {fehler}")
            return 1
        if len(a) != K:
            print(f"FEHLER: Das Klartextwort a braucht {K} Ziffern. Eingegeben: {len(a)} Ziffern.")
            return 1

    woerter: list[tuple[str, list[int]]] = []
    for i, text in enumerate(args.r, start=1):
        name = f"r{tief(i)}" if len(args.r) > 1 else "r"
        try:
            r = ziffern_lesen(text)
        except ValueError as fehler:
            print(f"FEHLER bei {name}: {fehler} Das Wort wird übersprungen.\n")
            continue
        if len(r) == N - 1:
            stelle = args.luecke or stelle_in_nullblock(r)
            if not 1 <= stelle <= N:
                print(f"FEHLER: --luecke muss zwischen 1 und {N} liegen.")
                return 1
            print(text_luecke(name, r, stelle) + "\n")
            r = null_einfuegen(r, stelle)
        elif len(r) != N:
            print(f"FEHLER: {name} hat {len(r)} Ziffern, der Code braucht {N}. "
                  "Das Wort wird übersprungen.\n")
            continue
        woerter.append((name, r))

    varianten = list(VARIANTEN) if args.alle_varianten else [args.variante]
    for v in varianten:
        print(loesungsweg(a, woerter, v, args.namen, sym))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
