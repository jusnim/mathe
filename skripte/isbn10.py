"""ISBN-10: Prüfziffer berechnen, Nummer prüfen, Zahlendreher erkennen.

Zweck
-----
Das Skript zeigt den Rechenweg Schritt für Schritt.
Eine ISBN-10 ist ein Wort c = c₁ c₂ … c₁₀ mit Buchstaben aus ℤ₁₁.
Die Prüfziffer c₁₀ erfüllt die Prüfsumme

    Σ i·cᵢ ≡ 0 (mod 11),   i = 1 … 10.

Daraus folgt c₁₀ ≡ Σ_{i=1}^{9} i·cᵢ (mod 11). Für 10 schreibt man X.

Aufruf
------
Prüfziffer ergänzen (9 Ziffern, Bindestriche sind erlaubt):
    python -m skripte.isbn10 --nummer 3-05-501517
Nummer prüfen (10 Zeichen):
    python -m skripte.isbn10 --nummer 3-528-06580-X
Zahlendreher an den Stellen 3 und 7 zeigen:
    python -m skripte.isbn10 --nummer 0-387-98999-4 --dreher 3 7
Allgemeiner Beweis „ISBN-10 erkennt Zahlendreher“:
    python -m skripte.isbn10 --dreher-beweis

Varianten (--variante)
----------------------
trick   (Standard) Gewichte 6 … 10 werden als −5 … −1 geschrieben.
        Dann fasst man Paare zusammen: 2·(c₂ − c₉) usw. Das spart
        Rechenarbeit.
direkt  Summe 1·c₁ + 2·c₂ + … + 9·c₉ ohne Trick.
praxis  Gewichte 10, 9, …, 1 wie auf echten Büchern. Die Prüfziffer ist
        gleich, die Zwischensummen sind anders.
Mit --alle-varianten zeigt das Skript alle drei Wege nacheinander.

--rep legt fest, wie die Summanden im Trick-Weg reduziert werden:
symmetrisch (−5 … 5, Standard) oder standard (0 … 10).
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field

M = 11
VARIANTEN = ("trick", "direkt", "praxis")
VARIANTEN_TEXT = {
    "trick": "Gewichte 1 … 10, dabei 6 … 10 als −5 … −1 geschrieben "
             "(Rechentrick mit Paaren)",
    "direkt": "Gewichte 1 … 10, Summe Σ i·cᵢ direkt",
    "praxis": "Gewichte 10, 9, …, 1 (Praxis). Gleiche Prüfziffer, "
              "andere Zwischensummen",
}


class EingabeFehler(ValueError):
    """Die Eingabe ist keine gültige ISBN-Eingabe."""


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def zeichen(wert: int) -> str:
    """Schreibt einen Wert aus 0 … 10 als Ziffer. 10 wird zu X."""
    return "X" if wert == 10 else str(wert)


def sym(wert: int, m: int = M) -> int:
    """Symmetrischer Repräsentant: Zahl aus −5 … 5 (für m = 11)."""
    r = wert % m
    return r - m if r > m // 2 else r


def reduziere(wert: int, rep: str) -> int:
    return sym(wert) if rep == "symmetrisch" else wert % M


def lies_nummer(text: str) -> list[int]:
    """Liest 9 oder 10 Zeichen. Bindestriche und Leerzeichen fallen weg.

    X (oder x) ist nur an Stelle 10 erlaubt und bedeutet 10.
    """
    roh = text.replace("-", "").replace(" ", "")
    if len(roh) not in (9, 10):
        raise EingabeFehler(
            f"Die Nummer hat {len(roh)} Zeichen. Nötig sind 9 Zeichen "
            "(Prüfziffer berechnen) oder 10 Zeichen (Nummer prüfen)."
        )
    werte = []
    for pos, z in enumerate(roh, start=1):
        if z.isdigit():
            werte.append(int(z))
        elif z in "Xx":
            if pos != 10:
                raise EingabeFehler(
                    f"X steht an Stelle {pos}. X (= 10) ist nur als "
                    "Prüfziffer an Stelle 10 erlaubt."
                )
            werte.append(10)
        else:
            raise EingabeFehler(f"Das Zeichen „{z}“ an Stelle {pos} ist keine Ziffer.")
    return werte


def formatiere(c: list[int], vorlage: str = "") -> str:
    """Schreibt die Nummer mit Bindestrichen, z. B. 3-05-501517-7.

    Die Gruppen sind bei jeder ISBN verschieden lang. Deshalb übernimmt
    das Skript die Bindestriche aus der Eingabe (vorlage). Ohne Vorlage
    nutzt es die Form c₁-c₂c₃c₄-c₅…c₉-c₁₀.
    """
    s = "".join(zeichen(x) for x in c)
    gruppen = [len(g) for g in vorlage.replace(" ", "-").split("-") if g]
    if sum(gruppen) not in (9, 10) or len(gruppen) < 2:
        gruppen = [1, 3, 5, 1]
    teile, pos = [], 0
    for laenge in gruppen:
        if pos >= len(s):
            break
        teile.append(s[pos:pos + laenge])
        pos += laenge
    if pos < len(s):
        teile.append(s[pos:])
    return "-".join(teile)


def zahl(v: int) -> str:
    """Schreibt eine ganze Zahl mit echtem Minuszeichen."""
    return str(v) if v >= 0 else f"−{-v}"


def tief(i: int) -> str:
    """Schreibt einen Index tiefgestellt, z. B. 10 → ₁₀."""
    return str(i).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))


def gewicht(i: int, variante: str) -> int:
    """Gewicht der Stelle i (1 … 10) in der gewählten Variante."""
    if variante == "praxis":
        return 11 - i
    if variante == "trick" and i >= 6:
        return i - 11
    return i


# ---------------------------------------------------------------------------
# Rechnen
# ---------------------------------------------------------------------------

@dataclass
class Pruefziffer:
    ziffern: list[int]            # c₁ … c₉
    variante: str
    rep: str
    gewichte: list[int]
    summanden: list[int]          # gewicht · cᵢ
    summe: int                    # Summe der Summanden (ganze Zahl)
    paare: list[tuple[int, int, int]] = field(default_factory=list)  # (i, cᵢ, c₁₁₋ᵢ)
    paar_werte: list[int] = field(default_factory=list)   # c₁ und i·(cᵢ − c₁₁₋ᵢ)
    paar_reduziert: list[int] = field(default_factory=list)
    summe_reduziert: int = 0
    c10: int = 0
    codewort: list[int] = field(default_factory=list)
    probe_summe: int = 0          # Σ_{i=1}^{10} i·cᵢ


def berechne_pruefziffer(ziffern: list[int], variante: str = "trick",
                         rep: str = "symmetrisch") -> Pruefziffer:
    """Berechnet c₁₀ zu c₁ … c₉ im gewählten Rechenweg."""
    if len(ziffern) != 9:
        raise EingabeFehler("Für die Prüfziffer braucht man genau 9 Ziffern.")
    if any(not 0 <= z <= 9 for z in ziffern):
        raise EingabeFehler("c₁ … c₉ müssen Ziffern 0 … 9 sein.")
    gew = [gewicht(i, variante) for i in range(1, 10)]
    summanden = [g * c for g, c in zip(gew, ziffern)]
    summe = sum(summanden)
    erg = Pruefziffer(ziffern, variante, rep, gew, summanden, summe)
    if variante == "trick":
        # c₁ + 2(c₂ − c₉) + 3(c₃ − c₈) + 4(c₄ − c₇) + 5(c₅ − c₆)
        erg.paar_werte.append(ziffern[0])
        for i in range(2, 6):
            ci, cj = ziffern[i - 1], ziffern[11 - i - 1]
            erg.paare.append((i, ci, cj))
            erg.paar_werte.append(i * (ci - cj))
        erg.paar_reduziert = [reduziere(w, rep) for w in erg.paar_werte]
        erg.summe_reduziert = sum(erg.paar_reduziert)
        erg.c10 = erg.summe_reduziert % M
    elif variante == "praxis":
        # Σ(11 − i)cᵢ + 1·c₁₀ ≡ 0  →  c₁₀ ≡ −Σ(11 − i)cᵢ
        erg.c10 = (-summe) % M
    else:
        erg.c10 = summe % M
    erg.codewort = ziffern + [erg.c10]
    erg.probe_summe = sum(i * c for i, c in enumerate(erg.codewort, start=1))
    return erg


@dataclass
class Pruefung:
    c: list[int]
    variante: str
    gewichte: list[int]
    summanden: list[int]
    summe: int
    rest: int                 # summe mod 11 (0 … 10)
    gueltig: bool
    richtige_c10: int


def pruefe_nummer(c: list[int], variante: str = "trick") -> Pruefung:
    """Prüft, ob h·cᵀ ≡ 0 (mod 11) gilt."""
    if len(c) != 10:
        raise EingabeFehler("Zum Prüfen braucht man genau 10 Zeichen.")
    gew = [gewicht(i, variante) for i in range(1, 11)]
    summanden = [g * x for g, x in zip(gew, c)]
    summe = sum(summanden)
    rest = summe % M
    richtige = berechne_pruefziffer(c[:9], "direkt").c10 if all(
        0 <= x <= 9 for x in c[:9]) else -1
    return Pruefung(c, variante, gew, summanden, summe, rest, rest == 0, richtige)


@dataclass
class Dreher:
    c: list[int]
    x: int
    y: int
    r: list[int]
    e: list[int]
    cx: int
    cy: int
    teil_x: int               # x·(c_y − c_x)
    teil_y: int               # y·(c_x − c_y)
    he: int                   # (x − y)(c_y − c_x)
    he_mod: int
    hr: int                   # Σ i·rᵢ direkt
    hr_mod: int
    erkannt: bool
    kein_fehler: bool         # c_x = c_y: r = c


def zahlendreher(c: list[int], x: int, y: int) -> Dreher:
    """Vertauscht c_x und c_y und berechnet h·rᵀ = h·eᵀ."""
    if len(c) != 10:
        raise EingabeFehler("Für den Zahlendreher braucht man eine 10-stellige Nummer.")
    for s in (x, y):
        if not 1 <= s <= 10:
            raise EingabeFehler(f"Die Stelle {s} liegt nicht in 1 … 10.")
    if x == y:
        raise EingabeFehler("x und y müssen verschiedene Stellen sein.")
    r = list(c)
    r[x - 1], r[y - 1] = c[y - 1], c[x - 1]
    e = [a - b for a, b in zip(r, c)]
    cx, cy = c[x - 1], c[y - 1]
    he = (x - y) * (cy - cx)
    hr = sum(i * v for i, v in enumerate(r, start=1))
    return Dreher(c, x, y, r, e, cx, cy, x * (cy - cx), y * (cx - cy),
                  he, he % M, hr, hr % M, he % M != 0, cx == cy)


# ---------------------------------------------------------------------------
# Ausgabe
# ---------------------------------------------------------------------------

def _p(v: int, erstes: bool = False) -> str:
    """Schreibt ein Vorzeichen mit Leerzeichen: + a bzw. − a."""
    if erstes:
        return str(v) if v >= 0 else f"−{-v}"
    return f"+ {v}" if v >= 0 else f"− {-v}"


def _kette(werte: list[int]) -> str:
    return " ".join(_p(v, k == 0) for k, v in enumerate(werte))


def _mal_kette(gewichte: list[int], ziffern: list[int]) -> str:
    teile = []
    for k, (g, c) in enumerate(zip(gewichte, ziffern)):
        if k == 0:
            teile.append(f"{g}·{zeichen(c)}" if g >= 0 else f"−{-g}·{zeichen(c)}")
        else:
            teile.append(f"+ {g}·{zeichen(c)}" if g >= 0 else f"− {-g}·{zeichen(c)}")
    return " ".join(teile)


def text_pruefziffer(p: Pruefziffer, vorlage: str = "") -> str:
    z = []
    z.append(f"Variante: {p.variante} – {VARIANTEN_TEXT[p.variante]}")
    z.append(f"Gegeben: {formatiere(p.ziffern, vorlage)}-c₁₀")
    z.append("c = " + ",  ".join(f"c{tief(i)} = {zeichen(v)}" for i, v in enumerate(p.ziffern, 1)))
    z.append("")
    if p.variante == "praxis":
        z.append("Prüfsumme (Praxis): Σ_{i=1}^{10} (11 − i)·cᵢ ≡ 0 (mod 11)")
        z.append("c₁₀ hat das Gewicht 1. Also gilt: c₁₀ ≡ −Σ_{i=1}^{9} (11 − i)·cᵢ (mod 11)")
        z.append("")
        z.append("Σ (11 − i)·cᵢ = " + _mal_kette(p.gewichte, p.ziffern))
        z.append(f"             = {_kette(p.summanden)}")
        z.append(f"             = {p.summe}")
        z.append(f"{p.summe} = {p.summe // M}·11 + {p.summe % M}  →  "
                 f"Σ ≡ {p.summe % M} (mod 11)")
        z.append(f"c₁₀ ≡ −{p.summe % M} ≡ {p.c10} (mod 11)")
    else:
        z.append("Prüfsumme: Σ_{i=1}^{10} i·cᵢ ≡ 0 (mod 11)")
        z.append("Wegen 10 ≡ −1 gilt: Σ_{i=1}^{9} i·cᵢ − c₁₀ ≡ 0")
        z.append("  →  c₁₀ ≡ Σ_{i=1}^{9} i·cᵢ (mod 11)")
        z.append("")
        if p.variante == "trick":
            z.append("Rechentrick: 6 ≡ −5, 7 ≡ −4, 8 ≡ −3, 9 ≡ −2 (mod 11).")
            z.append("")
            z.append("c₁₀ = " + _mal_kette(p.gewichte, p.ziffern))
            paar_txt = [f"{p.ziffern[0]}"] + [
                f"{i}·({ci} − {cj})" for i, ci, cj in p.paare]
            z.append("    = " + " + ".join(paar_txt))
            z.append("    = " + _kette(p.paar_werte)
                     + f"      (ungekürzt, Summe {zahl(sum(p.paar_werte))})")
            bereich = "−5 … 5" if p.rep == "symmetrisch" else "0 … 10"
            z.append("    = " + _kette(p.paar_reduziert)
                     + f"      (jeder Summand mod 11, Bereich {bereich})")
            z.append(f"    = {zahl(p.summe_reduziert)}")
            if p.summe_reduziert != p.c10:
                z.append(f"    ≡ {p.c10} (mod 11)")
        else:
            z.append("c₁₀ ≡ " + _mal_kette(p.gewichte, p.ziffern))
            z.append(f"    = {_kette(p.summanden)}")
            z.append(f"    = {p.summe}")
            z.append(f"{p.summe} = {p.summe // M}·11 + {p.summe % M}")
            z.append(f"c₁₀ ≡ {p.c10} (mod 11)")
    if p.c10 == 10:
        z.append("Der Wert 10 wird als X geschrieben.")
    z.append("")
    z.append(f"Ergebnis: c₁₀ = {zeichen(p.c10)},  ISBN = {formatiere(p.codewort, vorlage)}")
    z.append("")
    z.append("Probe: h·cᵀ = Σ_{i=1}^{10} i·cᵢ = "
             + _mal_kette(list(range(1, 11)), p.codewort)
             + f" = {p.probe_summe}")
    z.append(f"       {p.probe_summe} = {p.probe_summe // M}·11 + {p.probe_summe % M}"
             f"  →  h·cᵀ ≡ {p.probe_summe % M} (mod 11)"
             + ("  ✓" if p.probe_summe % M == 0 else "  ✗ FEHLER"))
    return "\n".join(z)


def text_pruefung(p: Pruefung, vorlage: str = "") -> str:
    z = []
    z.append(f"Variante: {p.variante} – {VARIANTEN_TEXT[p.variante]}")
    z.append(f"Gegeben: r = {formatiere(p.c, vorlage)}")
    if p.variante == "praxis":
        z.append("Prüfe: Σ_{i=1}^{10} (11 − i)·rᵢ ≡ 0 (mod 11)?")
    else:
        z.append("Prüfe: h·rᵀ = Σ_{i=1}^{10} i·rᵢ ≡ 0 (mod 11)?  mit h = 1 2 3 4 5 6 7 8 9 10")
        if p.variante == "trick":
            z.append("Rechentrick: 6 … 10 ≡ −5 … −1 (mod 11).")
    z.append("")
    z.append("Σ = " + _mal_kette(p.gewichte, p.c))
    z.append(f"  = {_kette(p.summanden)}")
    z.append(f"  = {zahl(p.summe)}")
    z.append(f"{zahl(p.summe)} ≡ {p.rest} (mod 11)")
    z.append("")
    if p.gueltig:
        z.append(f"Ergebnis: Die Prüfsumme ist ≡ 0. {formatiere(p.c, vorlage)} ist eine gültige ISBN-10.")
    else:
        z.append(f"Ergebnis: Die Prüfsumme ist ≡ {p.rest} ≠ 0. "
                 f"{formatiere(p.c, vorlage)} ist KEINE gültige ISBN-10.")
        z.append("Hinweis: Es gibt mindestens einen Fehler. Der Code erkennt ihn, "
                 "kann ihn aber nicht korrigieren.")
        if p.richtige_c10 >= 0:
            z.append(f"Zu c₁ … c₉ passt die Prüfziffer c₁₀ = {zeichen(p.richtige_c10)} "
                     "(falls der Fehler in c₁₀ liegt).")
    return "\n".join(z)


def text_dreher(d: Dreher) -> str:
    z = []
    x, y = d.x, d.y
    z.append("Zahlendreher: Die Zeichen an den Stellen x und y sind vertauscht.")
    z.append(f"x = {x}, y = {y}")
    z.append(f"c = {' '.join(zeichen(v) for v in d.c)}")
    z.append(f"r = {' '.join(zeichen(v) for v in d.r)}")
    z.append("e = r − c = " + " ".join(zahl(v) for v in d.e)
             + f"   (an Stelle {x}: c_y − c_x, an Stelle {y}: c_x − c_y)")
    z.append("")
    z.append("h·rᵀ = h·cᵀ + h·eᵀ ≡ h·eᵀ, denn h·cᵀ ≡ 0 (c ist Codewort).")
    z.append("h·eᵀ = x·(c_y − c_x) + y·(c_x − c_y) = (x − y)·(c_y − c_x)")
    z.append(f"     = {x}·({d.cy} − {d.cx}) + {y}·({d.cx} − {d.cy})"
             f" = {d.teil_x} {_p(d.teil_y)}")
    z.append(f"     = ({x} − {y})·({d.cy} − {d.cx}) = {zahl(x - y)}·{zahl(d.cy - d.cx)} = {zahl(d.he)}")
    z.append(f"     ≡ {d.he_mod} (mod 11)")
    z.append("")
    if d.kein_fehler:
        z.append("Hinweis: c_x = c_y. Das Vertauschen ändert nichts. Es gibt keinen Fehler.")
        z.append("Ergebnis: r = c, h·rᵀ ≡ 0. Kein Fehler vorhanden.")
    elif d.erkannt:
        z.append(f"Ergebnis: h·rᵀ ≡ {d.he_mod} ≠ 0 (mod 11). Der Zahlendreher wird erkannt.")
    else:
        z.append("Ergebnis: h·rᵀ ≡ 0. Der Fehler wird NICHT erkannt. "
                 "(Das darf bei ISBN-10 nicht vorkommen. Eingabe prüfen.)")
    z.append("")
    z.append(f"Probe: Σ i·rᵢ direkt = {d.hr} ≡ {d.hr_mod} (mod 11)"
             + ("  ✓" if d.hr_mod == d.he_mod else "  ✗ (c ist kein Codewort)"))
    return "\n".join(z)


def text_dreher_beweis() -> str:
    return "\n".join([
        "Behauptung: ISBN-10 erkennt jeden Zahlendreher.",
        "",
        "c = … c_x … c_y …   (Codewort, h·cᵀ ≡ 0)",
        "r = … c_y … c_x …   (Stellen x und y vertauscht, c_x ≠ c_y)",
        "h = 1 2 3 4 5 6 7 8 9 10",
        "",
        "e = r − c = 0 … (c_y − c_x) … (c_x − c_y) … 0",
        "            (Stelle x)      (Stelle y)",
        "",
        "⇒ h·rᵀ = h·cᵀ + h·eᵀ ≡ h·eᵀ",
        "       = x·(c_y − c_x) + y·(c_x − c_y)",
        "       = (x − y)·(c_y − c_x)",
        "",
        "Begründung, warum das ≠ 0 in ℤ₁₁ ist:",
        "- x ≠ y und 1 ≤ x, y ≤ 10. Also 0 < |x − y| ≤ 9. 11 teilt x − y nicht.",
        "- c_x ≠ c_y und 0 ≤ c_x, c_y ≤ 10. Also 0 < |c_y − c_x| ≤ 10. "
        "11 teilt c_y − c_x nicht.",
        "- ℤ₁₁ ist ein Körper (11 ist Primzahl). Ein Körper hat keine Nullteiler.",
        "  Nullteiler: Zahlen a, b ≠ 0 mit a·b = 0.",
        "",
        "Ergebnis: h·rᵀ = (x − y)·(c_y − c_x) ≠ 0 in ℤ₁₁. "
        "ISBN-10 erkennt den Zahlendreher.",
    ])


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

BEISPIEL = """Beispiele:
  python -m skripte.isbn10 --nummer 3-05-501517  (Prüfziffer berechnen: c₁₀ = 7)
  python -m skripte.isbn10 --nummer 0-471-82819  (Prüfziffer berechnen: c₁₀ = X)
  python -m skripte.isbn10 --nummer 3-528-06580-X  (Nummer prüfen)
  python -m skripte.isbn10 --nummer 0-387-98999-4 --dreher 3 7  (Zahlendreher an Stelle 3 und 7)
  python -m skripte.isbn10 --dreher-beweis  (allgemeiner Beweis für Zahlendreher)
"""


def baue_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m skripte.isbn10",
        description="ISBN-10: Prüfziffer berechnen, Nummer prüfen, Zahlendreher zeigen.",
        epilog=BEISPIEL,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--nummer",
                    help="ISBN mit 9 Zeichen (Prüfziffer wird berechnet) oder "
                         "10 Zeichen (Nummer wird geprüft). Bindestriche sind erlaubt. "
                         "X steht für 10.")
    ap.add_argument("--dreher", nargs=2, type=int, metavar=("X", "Y"),
                    help="Zahlendreher: Vertauscht die Stellen X und Y (1 … 10) "
                         "der Nummer und zeigt h·rᵀ = (x − y)(c_y − c_x).")
    ap.add_argument("--dreher-beweis", action="store_true",
                    help="Zeigt den allgemeinen Beweis, dass ISBN-10 jeden "
                         "Zahlendreher erkennt.")
    ap.add_argument("--variante", choices=VARIANTEN, default="trick",
                    help="Rechenweg: trick (Standard, Gewichte 6 … 9 "
                         "als −5 … −2), direkt (Σ i·cᵢ), praxis (Gewichte 10 … 1).")
    ap.add_argument("--alle-varianten", action="store_true",
                    help="Zeigt alle Rechenwege nacheinander.")
    ap.add_argument("--rep", choices=("symmetrisch", "standard"), default="symmetrisch",
                    help="Reste der Summanden im Trick-Weg: symmetrisch (−5 … 5, "
                         "Standard) oder standard (0 … 10).")
    return ap


def loese(args: argparse.Namespace) -> str:
    """Erzeugt den gesamten Text für die gegebenen Optionen."""
    teile: list[str] = []
    if args.dreher_beweis:
        teile.append(text_dreher_beweis())
    if args.nummer is not None:
        werte = lies_nummer(args.nummer)
        if args.dreher:
            if len(werte) == 9:
                p = berechne_pruefziffer(werte, "direkt")
                teile.append(f"Hinweis: Nur 9 Zeichen. Ich ergänze zuerst "
                             f"c₁₀ = {zeichen(p.c10)}.")
                werte = p.codewort
            elif not pruefe_nummer(werte, "direkt").gueltig:
                teile.append("Hinweis: Die Nummer ist selbst kein Codewort "
                             "(h·cᵀ ≢ 0). Die Rechnung zeigt nur h·eᵀ.")
            teile.append(text_dreher(zahlendreher(werte, *args.dreher)))
        else:
            varianten = VARIANTEN if args.alle_varianten else (args.variante,)
            for k, v in enumerate(varianten):
                if len(varianten) > 1:
                    teile.append(f"===== Variante {k + 1} von {len(varianten)}: {v} =====")
                if len(werte) == 9:
                    teile.append(text_pruefziffer(berechne_pruefziffer(werte, v, args.rep), args.nummer))
                else:
                    teile.append(text_pruefung(pruefe_nummer(werte, v), args.nummer))
    elif args.dreher:
        raise EingabeFehler("Für --dreher braucht man --nummer. "
                            "Für den allgemeinen Beweis: --dreher-beweis.")
    if not teile:
        raise EingabeFehler("Bitte --nummer oder --dreher-beweis angeben. Hilfe: --help")
    return "\n\n".join(teile)


def main(argv: list[str] | None = None) -> int:
    ap = baue_parser()
    args = ap.parse_args(argv)
    try:
        print(loese(args))
    except EingabeFehler as fehler:
        print(f"Hinweis zur Eingabe: {fehler}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
