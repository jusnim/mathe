"""Primfaktorzerlegung und Eulersche φ-Funktion mit Lösungsweg.

Zweck
-----
Das Skript zerlegt eine Zahl n in Primfaktoren. Eine Primzahl ist eine Zahl > 1,
die nur durch 1 und sich selbst teilbar ist. Danach berechnet das Skript φ(n).
φ(n) (Eulersche φ-Funktion) ist die Anzahl der Zahlen a von 1 bis n mit ggT(a, n) = 1.
Die Ausgabe zeigt jeden Schritt.

Aufruf
------
    python -m skripte.primfaktor_phi --n 360
    python -m skripte.primfaktor_phi --n 2451 --faktoren 43 57
    python -m skripte.primfaktor_phi --n 2803 --phi-phi
    python -m skripte.primfaktor_phi --n 1000 --alle-varianten

Varianten (--variante)
----------------------
- produktformel (Standard):
  φ(n) = n · (1 − 1/p₁) · … · (1 − 1/pᵣ)
- multiplikativ:
  φ(n) = φ(p₁^k₁) · … · φ(pᵣ^kᵣ) mit φ(pᵏ) = pᵏ − p^(k−1)
Mit --alle-varianten erscheinen beide Wege nacheinander.

Angabe „m = p · q“
------------------
Das Skript zerlegt immer vollständig. Beispiel: m = 43 · 57, aber 57 = 3 · 19.
Dann gilt φ(m) ≠ (43 − 1)(57 − 1). Mit --faktoren warnt das Skript davor.

Nutzung als Modul
-----------------
    from skripte.primfaktor_phi import ist_prim, primfaktorzerlegung, phi, phi_mit_weg
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd, prod

VARIANTEN = ("produktformel", "multiplikativ")
STANDARD_VARIANTE = "produktformel"

# Bis zu dieser Größe zählt die Probe alle a mit ggT(a, n) = 1 direkt.
PROBE_GRENZE = 200_000

_HOCH = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


# ---------------------------------------------------------------------------
# Rechnen
# ---------------------------------------------------------------------------

def ist_prim(n: int) -> bool:
    """Gibt True zurück, wenn n eine Primzahl ist."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    t = 3
    while t * t <= n:
        if n % t == 0:
            return False
        t += 2
    return True


@dataclass
class Divisionsschritt:
    """Ein Schritt der Probedivision: rest : teiler."""

    rest: int
    teiler: int
    geht: bool
    quotient: int | None = None


@dataclass
class Zerlegung:
    """Ergebnis der Probedivision."""

    n: int
    faktoren: dict[int, int]
    schritte: list[Divisionsschritt] = field(default_factory=list)
    schluss: str = ""


def probedivision(n: int) -> Zerlegung:
    """Zerlegt n durch Probedivision und merkt sich jeden Schritt.

    Probedivision heißt: Man teilt n nacheinander durch 2, 3, 5, 7, …
    Man hört auf, wenn teiler² > Rest ist. Dann ist der Rest 1 oder eine Primzahl.
    """
    if n < 1:
        raise ValueError("n muss eine natürliche Zahl ≥ 1 sein.")
    faktoren: dict[int, int] = {}
    schritte: list[Divisionsschritt] = []
    rest = n
    teiler = 2
    while teiler * teiler <= rest:
        if rest % teiler == 0:
            q = rest // teiler
            schritte.append(Divisionsschritt(rest, teiler, True, q))
            faktoren[teiler] = faktoren.get(teiler, 0) + 1
            rest = q
        else:
            schritte.append(Divisionsschritt(rest, teiler, False))
            # Nur Primzahlen als Teiler probieren (9, 15, … überspringen).
            teiler = 3 if teiler == 2 else teiler + 2
            while not ist_prim(teiler):
                teiler += 2
    if rest > 1:
        faktoren[rest] = faktoren.get(rest, 0) + 1
        if teiler * teiler > rest and schritte:
            schluss = (f"{teiler}² = {teiler * teiler} > {rest}. "
                       f"Also ist {rest} eine Primzahl.")
        else:
            schluss = f"{rest} ist eine Primzahl."
    else:
        schluss = "Der Rest ist 1. Die Zerlegung ist fertig."
    return Zerlegung(n, dict(sorted(faktoren.items())), schritte, schluss)


def primfaktorzerlegung(n: int) -> dict[int, int]:
    """Gibt die Zerlegung als {Primzahl: Exponent} zurück. Beispiel: 360 → {2: 3, 3: 2, 5: 1}."""
    return probedivision(n).faktoren


def phi_primpotenz(p: int, k: int) -> int:
    """φ(pᵏ) = pᵏ − p^(k−1)."""
    return p ** k - p ** (k - 1)


def phi(n: int) -> int:
    """Eulersche φ-Funktion: Anzahl der a von 1 bis n mit ggT(a, n) = 1."""
    if n < 1:
        raise ValueError("n muss eine natürliche Zahl ≥ 1 sein.")
    return prod(phi_primpotenz(p, k) for p, k in primfaktorzerlegung(n).items())


def phi_direkt(n: int) -> int:
    """Zählt direkt alle a von 1 bis n mit ggT(a, n) = 1 (nur für die Probe)."""
    return sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


@dataclass
class FaktorPruefung:
    """Prüfung der vorgegebenen Faktoren (z. B. m = 43 · 57)."""

    faktoren: list[int]
    produkt: int
    produkt_stimmt: bool
    keine_primzahl: dict[int, dict[int, int]]  # Faktor → seine Zerlegung
    falsches_phi: int | None  # ∏(f − 1), falls man die Faktoren für prim hält


def pruefe_faktoren(n: int, faktoren: list[int]) -> FaktorPruefung:
    """Prüft, ob die vorgegebenen Faktoren n ergeben und ob sie Primzahlen sind."""
    produkt = prod(faktoren)
    keine = {f: primfaktorzerlegung(f) for f in faktoren if not ist_prim(f)}
    falsch = prod(f - 1 for f in faktoren) if keine else None
    return FaktorPruefung(list(faktoren), produkt, produkt == n, keine, falsch)


@dataclass
class PhiErgebnis:
    """Alle Werte für den Lösungsweg zu φ(n)."""

    n: int
    zerlegung: Zerlegung
    phi: int
    # Produktformel
    brueche: list[Fraction]
    bruchprodukt: Fraction
    # Multiplikativ
    teil_phis: list[tuple[int, int, int]]  # (p, k, φ(pᵏ))
    probe_direkt: int | None
    faktor_pruefung: FaktorPruefung | None = None
    phi_phi: "PhiErgebnis | None" = None


def phi_mit_weg(n: int, faktoren: list[int] | None = None,
                phi_phi: bool = False) -> PhiErgebnis:
    """Berechnet Zerlegung und φ(n) und sammelt alle Zwischenwerte.

    faktoren: vorgegebene Faktoren (werden geprüft).
    phi_phi: True berechnet zusätzlich φ(φ(n)) (für RSA, Weg über Euler).
    """
    z = probedivision(n)
    brueche = [1 - Fraction(1, p) for p in z.faktoren]
    bruchprodukt = prod(brueche, start=Fraction(1))
    teil = [(p, k, phi_primpotenz(p, k)) for p, k in z.faktoren.items()]
    wert = prod(t[2] for t in teil)
    assert n * bruchprodukt == wert
    probe = phi_direkt(n) if n <= PROBE_GRENZE else None
    pruef = pruefe_faktoren(n, faktoren) if faktoren else None
    ergebnis = PhiErgebnis(n, z, wert, brueche, bruchprodukt, teil, probe, pruef)
    if phi_phi:
        ergebnis.phi_phi = phi_mit_weg(wert)
    return ergebnis


# ---------------------------------------------------------------------------
# Ausgeben
# ---------------------------------------------------------------------------

def hoch(k: int) -> str:
    """Schreibt k als hochgestellte Zahl."""
    return str(k).translate(_HOCH)


def potenz(p: int, k: int) -> str:
    return str(p) if k == 1 else f"{p}{hoch(k)}"


def zerlegung_text(faktoren: dict[int, int]) -> str:
    """Beispiel: {2: 3, 3: 2, 5: 1} → '2³ · 3² · 5'."""
    if not faktoren:
        return "1"
    return " · ".join(potenz(p, k) for p, k in faktoren.items())


def _bruch(b: Fraction) -> str:
    return str(b.numerator) if b.denominator == 1 else f"{b.numerator}/{b.denominator}"


def _weg_produktformel(e: PhiErgebnis) -> list[str]:
    n = e.n
    ps = list(e.zerlegung.faktoren)
    zeilen = ["φ(n) = n · ∏ (1 − 1/p)   (Produkt über alle Primteiler p von n)"]
    if not ps:
        return zeilen + ["φ(1) = 1 (n = 1 hat keine Primteiler)."]
    klammern = " · ".join(f"(1 − 1/{p})" for p in ps)
    brueche = " · ".join(_bruch(b) for b in e.brueche)
    zeilen.append(f"φ({n}) = {n} · {klammern}")
    zeilen.append(f"       = {n} · {brueche}")
    if len(ps) > 1:
        zeilen.append(f"       = {n} · {_bruch(e.bruchprodukt)}")
    # Erst n durch den Nenner teilen, dann mit dem Zähler malen.
    nenner = e.bruchprodukt.denominator
    zaehler = e.bruchprodukt.numerator
    if nenner > 1 and zaehler > 1:
        zeilen.append(f"       = {n // nenner} · {zaehler}")
    zeilen.append(f"       = {e.phi}")
    return zeilen


def _weg_multiplikativ(e: PhiErgebnis) -> list[str]:
    n = e.n
    zeilen = ["φ(m · n) = φ(m) · φ(n), wenn ggT(m, n) = 1.",
              "φ(p) = p − 1 und φ(pᵏ) = pᵏ − p^(k−1) für eine Primzahl p."]
    if not e.teil_phis:
        return zeilen + ["φ(1) = 1."]
    teile = " · ".join(f"φ({potenz(p, k)})" for p, k, _ in e.teil_phis)
    formeln = []
    for p, k, _ in e.teil_phis:
        if k == 1:
            formeln.append(f"({p} − 1)")
        else:
            formeln.append(f"({potenz(p, k)} − {potenz(p, k - 1)})")
    werte = " · ".join(str(w) for _, _, w in e.teil_phis)
    zeilen.append(f"φ({n}) = {teile}")
    zeilen.append(f"       = {' · '.join(formeln)}")
    if len(e.teil_phis) > 1:
        zeilen.append(f"       = {werte}")
    zeilen.append(f"       = {e.phi}")
    return zeilen


def _zerlegung_zeilen(e: PhiErgebnis) -> list[str]:
    z = e.zerlegung
    zeilen = [f"Primfaktorzerlegung von {z.n} (Probedivision):",
              "Teile durch die Primzahlen 2, 3, 5, 7, 11, … Stopp, wenn teiler² > Rest."]
    if z.n == 1:
        return zeilen + ["n = 1 hat keine Primfaktoren."]
    for s in z.schritte:
        if s.geht:
            zeilen.append(f"  {s.rest} : {s.teiler} = {s.quotient}")
        else:
            zeilen.append(f"  {s.rest} : {s.teiler} geht nicht")
    zeilen.append(f"  {z.schluss}")
    zeilen.append(f"{z.n} = {zerlegung_text(z.faktoren)}")
    einzeln = " · ".join(" · ".join([str(p)] * k) for p, k in z.faktoren.items())
    zeilen.append(f"Probe: {einzeln} = {prod(p ** k for p, k in z.faktoren.items())} ✓")
    if ist_prim(z.n):
        zeilen.append(f"Also ist {z.n} eine Primzahl. Dann gilt φ({z.n}) = {z.n} − 1.")
    return zeilen


def _faktor_zeilen(e: PhiErgebnis) -> list[str]:
    f = e.faktor_pruefung
    if f is None:
        return []
    angabe = " · ".join(str(x) for x in f.faktoren)
    zeilen = [f"Prüfung der Angabe n = {angabe}:"]
    if not f.produkt_stimmt:
        zeilen.append(f"  ACHTUNG: {angabe} = {f.produkt} ≠ {e.n}. "
                      "Die Angabe passt nicht zu n. Das Skript rechnet mit n.")
    else:
        zeilen.append(f"  {angabe} = {f.produkt} ✓")
    if not f.keine_primzahl:
        zeilen.append("  Alle angegebenen Faktoren sind Primzahlen.")
        return zeilen
    for x, zx in f.keine_primzahl.items():
        zeilen.append(f"  ACHTUNG: {x} ist keine Primzahl. {x} = {zerlegung_text(zx)}.")
    zeilen.append("  Die Formel φ = (p − 1)(q − 1) gilt nur für Primzahlen p, q.")
    formel = " · ".join(f"({x} − 1)" for x in f.faktoren)
    zeilen.append(f"  Falsch wäre: {formel} = {f.falsches_phi}.")
    zeilen.append("  Deshalb: n vollständig in Primfaktoren zerlegen (siehe unten).")
    return zeilen


def loesungsweg(e: PhiErgebnis, varianten: list[str] | None = None) -> str:
    """Erzeugt den Lösungsweg als Text aus einem PhiErgebnis."""
    varianten = varianten or [STANDARD_VARIANTE]
    zeilen: list[str] = []
    fz = _faktor_zeilen(e)
    if fz:
        zeilen += fz + [""]
    zeilen += _zerlegung_zeilen(e)
    namen = {"produktformel": "Produktformel n · ∏(1 − 1/p)",
             "multiplikativ": "Multiplikativität φ = ∏ φ(pᵏ)"}
    for v in varianten:
        zeilen += ["", f"Variante „{v}“: {namen[v]}"]
        zeilen += _weg_produktformel(e) if v == "produktformel" else _weg_multiplikativ(e)
    if e.probe_direkt is not None:
        ok = "✓" if e.probe_direkt == e.phi else "FEHLER"
        zeilen += ["", f"Probe: Zählen aller a von 1 bis {e.n} mit ggT(a, {e.n}) = 1 "
                       f"ergibt {e.probe_direkt} {ok}"]
    zeilen += ["", f"Ergebnis: {e.n} = {zerlegung_text(e.zerlegung.faktoren)}, "
                   f"φ({e.n}) = {e.phi}"]
    if e.phi_phi is not None:
        pp = e.phi_phi
        zeilen += ["", "=" * 60, f"Zusatz: φ(φ({e.n})) = φ({e.phi})", "=" * 60]
        zeilen += [loesungsweg(pp, varianten)]
        zeilen += ["", f"Ergebnis: φ(φ({e.n})) = φ({e.phi}) = {pp.phi}"]
    return "\n".join(zeilen)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m skripte.primfaktor_phi",
        description="Zerlegt n in Primfaktoren und berechnet φ(n) mit Lösungsweg.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Beispiele:\n"
            "  python -m skripte.primfaktor_phi --n 360 --alle-varianten   (beide Rechenwege)\n"
            "  python -m skripte.primfaktor_phi --n 2451 --faktoren 43 57   (Faktoren prüfen)\n"
            "  python -m skripte.primfaktor_phi --n 2803 --phi-phi   (zusätzlich φ(φ(n)))"
        ),
    )
    p.add_argument("--n", type=int, required=True,
                   help="Die Zahl n, die zerlegt wird (natürliche Zahl ≥ 1).")
    p.add_argument("--faktoren", type=int, nargs="+", metavar="F",
                   help="Vorgegebene Faktoren, z. B. 43 57 bei m = 43 · 57. "
                        "Das Skript prüft, ob sie Primzahlen sind.")
    p.add_argument("--phi-phi", action="store_true",
                   help="Berechnet zusätzlich φ(φ(n)). Das braucht man bei RSA für d über Euler.")
    p.add_argument("--variante", choices=VARIANTEN, default=STANDARD_VARIANTE,
                   help="Rechenweg für φ(n): produktformel = n · ∏(1 − 1/p) (Standard), "
                        "multiplikativ = ∏(pᵏ − p^(k−1)).")
    p.add_argument("--alle-varianten", action="store_true",
                   help="Zeigt beide Rechenwege nacheinander.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.n < 1:
        print(f"Hinweis: n = {args.n} ist keine natürliche Zahl ≥ 1. "
              "φ(n) ist nur für n ≥ 1 definiert.")
        return 1
    if args.faktoren and any(f < 1 for f in args.faktoren):
        print("Hinweis: Alle Faktoren müssen natürliche Zahlen ≥ 1 sein.")
        return 1
    varianten = list(VARIANTEN) if args.alle_varianten else [args.variante]
    e = phi_mit_weg(args.n, args.faktoren, args.phi_phi)
    print("Annahme: n wird immer vollständig in Primfaktoren zerlegt.")
    print("Rechenweg: " + ", ".join(varianten)
          + (" (Standard)" if varianten == [STANDARD_VARIANTE] else ""))
    print()
    print(loesungsweg(e, varianten))
    return 0


if __name__ == "__main__":
    sys.exit(main())
