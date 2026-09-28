"""Schnelles Potenzieren a^k mod m, Satz von Fermat/Euler, Ordnung und primitive Elemente.

Zweck
-----
Das Skript rechnet a^k mod m wie in der Vorlesung (VL3, Abschnitt 3.1).
Es zeigt jeden Schritt:
1. den Exponenten k binär (höchstes Bit links) und als Summe von Zweierpotenzen,
2. die Quadrattabelle a, a^2, a^4, a^8, ... mod m,
3. die Auswahl der Potenzen mit Bit 1,
4. das Produkt Schritt für Schritt mod m.
Am Ende steht "Ergebnis: ..." und eine Probe.

Zusätzlich kann das Skript die Ordnung von Elementen in Z_m* bestimmen.
Die Ordnung von a ist die kleinste Zahl k > 0 mit a^k = 1.
Ein Element heißt primitiv, wenn seine Ordnung gleich phi(m) = |Z_m*| ist.

Aufruf (Beispiele aus den Quellen)
----------------------------------
  python -m skripte.schnell_potenzieren --a 12 --k 100 --m 34          (VL3: 30)
  python -m skripte.schnell_potenzieren --a 3 --k 1000 --m 7 --variante euler   (Ü4 A7: 4)
  python -m skripte.schnell_potenzieren --ordnung --a 2 3 4 5 --m 17   (Ü4 A8)
  python -m skripte.schnell_potenzieren --primitiv --m 17

Varianten
---------
--variante binaer  (Standard, VL3 Beispiel 12^100): nur schnelles Potenzieren.
--variante euler   (VL3 Bemerkung, Ü4 A7): zuerst k mod phi(m) rechnen.
                   Das geht nur bei ggT(a, m) = 1. Sonst warnt das Skript.
--reste normal       Reste 0 ... m-1 (Ü4 A9). Standard beim Potenzieren.
--reste symmetrisch  Reste zwischen -m/2 und m/2 (VL3: -4, -16; Ü4 A8).
                     Standard bei --ordnung.
--alle-varianten     gibt alle Kombinationen nacheinander aus.

Nutzung als Modul (zum Beispiel in rsa.py)
------------------------------------------
  from skripte.schnell_potenzieren import schnell_potenzieren, potenz_text
  erg = schnell_potenzieren(715, 113, 2803)   # erg.ergebnis == 708
  print(potenz_text(erg, name="c"))
Weitere Funktionen: binaer_zerlegung, binaer_text, quadrat_tabelle,
potenz_mit_reduktion, reduktion_text, ordnung, ordnung_text,
primitive_elemente, primitiv_text, phi, phi_text, rest.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field

from sympy import factorint, isprime

VARIANTEN = ("binaer", "euler")
RESTE = ("normal", "symmetrisch")

# Bekannte Schreibfehler in den Musterlösungen: (a, m, 2^i) -> Hinweis.
QUELLEN_HINWEISE = {
    (715, 2803, 16): "Ü4 A9 schreibt 715¹⁶ ≡ 163. Richtig ist 1653 (eine Ziffer fehlt). "
                     "Mit 163 käme nicht c = 708 heraus.",
}

HOCH = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def hoch(n: int) -> str:
    """Schreibt eine Zahl als Hochzahl, zum Beispiel 16 -> ¹⁶."""
    return str(n).translate(HOCH)


def z(x: int) -> str:
    """Schreibt negative Zahlen mit echtem Minuszeichen."""
    return str(x).replace("-", "−")


def zk(x: int) -> str:
    """Wie z(), aber negative Zahlen stehen in Klammern: (−16)."""
    return f"({z(x)})" if x < 0 else z(x)


# ---------------------------------------------------------------------------
# Rechnen
# ---------------------------------------------------------------------------

def rest(x: int, m: int, reste: str = "normal") -> int:
    """Rest von x modulo m.

    normal: 0 ... m-1. symmetrisch: der Rest mit dem kleinsten Betrag
    (-m/2 < r <= m/2).
    """
    r = x % m
    if reste == "symmetrisch" and r > m // 2:
        r -= m
    return r


@dataclass
class Binaer:
    """Binärdarstellung von k."""
    k: int
    bits: str                 # höchstes Bit links, zum Beispiel "1110001"
    exponenten: list[int]     # i mit Bit k_i = 1, aufsteigend


def binaer_zerlegung(k: int) -> Binaer:
    """Zerlegt k >= 0 in Zweierpotenzen."""
    if k < 0:
        raise ValueError("k muss >= 0 sein.")
    bits = bin(k)[2:]
    exps = [i for i in range(k.bit_length()) if (k >> i) & 1]
    return Binaer(k, bits, exps)


def quadrat_tabelle(a: int, m: int, anzahl: int, reste: str = "normal") -> list[tuple[int, int, int]]:
    """Tabelle a^(2^i) mod m für i = 0 ... anzahl-1.

    Jede Zeile ist (Exponent 2^i, Quadrat des Vorgängers ohne Reduktion,
    reduzierter Wert). In Zeile 0 ist das "Quadrat" einfach a.
    """
    zeilen = []
    wert = rest(a, m, reste)
    zeilen.append((1, a, wert))
    for i in range(1, anzahl):
        quadrat = wert * wert
        wert = rest(quadrat, m, reste)
        zeilen.append((2 ** i, quadrat, wert))
    return zeilen


@dataclass
class PotenzErgebnis:
    """Alle Daten des Lösungswegs für a^k mod m."""
    a: int
    k: int
    m: int
    reste: str
    binaer: Binaer
    tabelle: list[tuple[int, int, int]]      # (2^i, Quadrat, Rest)
    faktoren: list[tuple[int, int]]          # (2^i, Rest) mit Bit 1
    schritte: list[tuple[int, int, int, int]]  # (links, rechts, Produkt, Rest)
    ergebnis: int                            # immer 0 ... m-1


def schnell_potenzieren(a: int, k: int, m: int, reste: str = "normal") -> PotenzErgebnis:
    """Berechnet a^k mod m mit Square-and-Multiply und speichert jeden Schritt."""
    if m < 1:
        raise ValueError("Der Modul m muss >= 1 sein.")
    if k < 0:
        raise ValueError("Der Exponent k muss >= 0 sein.")
    b = binaer_zerlegung(k)
    tab = quadrat_tabelle(a, m, max(k.bit_length(), 1), reste)
    faktoren = [(tab[i][0], tab[i][2]) for i in b.exponenten]
    schritte = []
    if faktoren:
        akt = faktoren[0][1]
        for _, f in faktoren[1:]:
            produkt = akt * f
            neu = rest(produkt, m, reste)
            schritte.append((akt, f, produkt, neu))
            akt = neu
        ergebnis = akt % m
    else:
        ergebnis = 1 % m
    return PotenzErgebnis(a, k, m, reste, b, tab, faktoren, schritte, ergebnis)


def phi(m: int) -> int:
    """Eulersche phi-Funktion: Anzahl der Zahlen 1 ... m mit ggT(x, m) = 1."""
    ergebnis = m
    for p in factorint(m):
        ergebnis = ergebnis // p * (p - 1)
    return ergebnis


@dataclass
class ReduktionErgebnis:
    """Daten für den Weg "zuerst k mod phi(m)" (Euler/Fermat)."""
    a: int
    k: int
    m: int
    ggt: int
    phi_m: int
    q: int | None               # k = q * phi(m) + r
    r: int | None
    moeglich: bool              # nur bei ggT(a, m) = 1
    potenz: PotenzErgebnis      # Rechnung für a^r (oder a^k, wenn nicht möglich)
    ergebnis: int


def potenz_mit_reduktion(a: int, k: int, m: int, reste: str = "normal") -> ReduktionErgebnis:
    """a^k mod m. Falls ggT(a, m) = 1: zuerst k durch k mod phi(m) ersetzen."""
    g = math.gcd(a, m)
    ph = phi(m)
    if g == 1 and m > 1:
        q, r = divmod(k, ph)
        pot = schnell_potenzieren(a, r, m, reste)
        return ReduktionErgebnis(a, k, m, g, ph, q, r, True, pot, pot.ergebnis)
    pot = schnell_potenzieren(a, k, m, reste)
    return ReduktionErgebnis(a, k, m, g, ph, None, None, False, pot, pot.ergebnis)


@dataclass
class OrdnungErgebnis:
    """Daten für die Ordnung von a in Z_m*."""
    a: int
    m: int
    reste: str
    einheit: bool               # ggT(a, m) = 1?
    potenzen: list[int]         # a, a^2, ..., a^ord (= 1), als Reste
    ordnung: int | None
    phi_m: int
    primitiv: bool
    inverse: int | None         # a^(ord-1), 0 ... m-1


def ordnung(a: int, m: int, reste: str = "symmetrisch") -> OrdnungErgebnis:
    """Bestimmt ord(a) in Z_m* durch Aufschreiben aller Potenzen bis 1."""
    if m < 2:
        raise ValueError("Der Modul m muss >= 2 sein.")
    ph = phi(m)
    if math.gcd(a, m) != 1:
        return OrdnungErgebnis(a, m, reste, False, [], None, ph, False, None)
    potenzen = []
    wert = 1
    for _ in range(ph):
        wert = wert * a % m
        potenzen.append(rest(wert, m, reste))
        if wert == 1:
            break
    o = len(potenzen)
    inverse = pow(a, o - 1, m)
    return OrdnungErgebnis(a, m, reste, True, potenzen, o, ph, o == ph, inverse)


def primitive_elemente(m: int) -> list[int]:
    """Alle primitiven Elemente von Z_m* (leer, wenn Z_m* nicht zyklisch ist)."""
    if m < 2:
        return []
    ph = phi(m)
    teiler = list(factorint(ph)) if ph > 1 else []
    ergebnis = []
    for g in range(1, m):
        if math.gcd(g, m) != 1:
            continue
        # g ist primitiv, wenn g^(ph/p) != 1 für jeden Primteiler p von ph.
        if all(pow(g, ph // p, m) != 1 for p in teiler):
            ergebnis.append(g)
    return ergebnis


# ---------------------------------------------------------------------------
# Ausgeben
# ---------------------------------------------------------------------------

def binaer_text(k: int) -> str:
    """Zeile wie in den Quellen: k binär und als Summe von Zweierpotenzen."""
    b = binaer_zerlegung(k)
    if k == 0:
        return "0 = 0₂ (keine Zweierpotenz)"
    absteigend = sorted(b.exponenten, reverse=True)
    summe_hoch = " + ".join(f"2{hoch(i)}" for i in absteigend)
    summe_zahl = " + ".join(str(2 ** i) for i in absteigend)
    teile = [f"{k} binär: {b.bits}₂ (höchstes Bit links)", f"{k} = {summe_hoch}"]
    if len(absteigend) > 1:
        teile.append(f"{k} = {summe_zahl}")
    return "\n".join(teile)


def _zerlegt(x: int, m: int, r: int) -> str:
    """Schreibt x = q·m + r (r darf negativ sein)."""
    q = (x - r) // m
    if q == 0:
        return ""
    if r == 0:
        return f" = {q} · {m}"
    zeichen = "+" if r > 0 else "−"
    return f" = {q} · {m} {zeichen} {abs(r)}"


def potenz_text(erg: PotenzErgebnis, name: str | None = None) -> str:
    """Lösungsweg für a^k mod m als Text (Stil VL3 und Ü4 A9)."""
    a, k, m = erg.a, erg.k, erg.m
    basis = zk(a)
    zeilen = [f"Gesucht: {basis}{hoch(k)} mod {m}"]
    if erg.reste == "symmetrisch":
        zeilen.append("Annahme: symmetrische Reste (Betrag ≤ m/2), wie VL3 (−4, −16).")
    else:
        zeilen.append("Annahme: Reste 0 … m−1, wie Ü4 A9.")
    if k == 0:
        zeilen.append(f"k = 0, also {basis}⁰ = 1.")
        zeilen.append(f"Ergebnis: {basis}⁰ ≡ {erg.ergebnis} (mod {m})")
        return "\n".join(zeilen)

    zeilen.append("")
    zeilen.append("1. Exponent binär schreiben")
    zeilen.append(binaer_text(k))
    produkt_potenzen = " · ".join(f"{basis}{hoch(e)}" for e, _ in erg.faktoren)
    zeilen.append(f"⇒ {basis}{hoch(k)} = {produkt_potenzen}")

    zeilen.append("")
    zeilen.append(f"2. Quadrattabelle (alles mod {m}, jede Zeile ist das Quadrat der Zeile davor)")
    breite = len(f"{basis}{hoch(erg.tabelle[-1][0])}")
    vorher = None
    for e, quadrat, wert in erg.tabelle:
        links = f"{basis}{hoch(e)}".ljust(breite)
        markiert = "  ← Bit 1" if e in [f[0] for f in erg.faktoren] else ""
        if vorher is None:
            if wert == a:
                zeilen.append(f"{links} = {z(a)}{markiert}")
            else:
                zeilen.append(f"{links} = {z(a)} ≡ {z(wert)}{markiert}")
        else:
            zeile = f"{links} ≡ {zk(vorher)}² = {z(quadrat)}{_zerlegt(quadrat, m, wert)}"
            if quadrat != wert:
                zeile += f" ≡ {z(wert)}"
            zeilen.append(zeile + markiert)
        vorher = wert

    for e, _, _ in erg.tabelle:
        hinweis = QUELLEN_HINWEISE.get((a, m, e))
        if hinweis:
            zeilen.append(f"HINWEIS zur Quelle: {hinweis}")

    zeilen.append("")
    zeilen.append("3. Potenzen mit Bit 1 auswählen")
    auswahl = " · ".join(zk(w) for _, w in erg.faktoren)
    zeilen.append(f"{basis}{hoch(k)} = {produkt_potenzen} ≡ {auswahl} (mod {m})")

    zeilen.append("")
    zeilen.append(f"4. Produkt Schritt für Schritt mod {m}")
    if not erg.schritte:
        zeilen.append("Nur ein Faktor, nichts zu multiplizieren.")
    for links, rechts, produkt, neu in erg.schritte:
        zeile = f"{zk(links)} · {zk(rechts)} = {z(produkt)}"
        if produkt != neu:
            zeile += f"{_zerlegt(produkt, m, neu)} ≡ {z(neu)}"
        zeilen.append(zeile + f" (mod {m})")
    schluss = erg.faktoren[0][1] if not erg.schritte else erg.schritte[-1][3]
    if schluss != erg.ergebnis:
        zeilen.append(f"{z(schluss)} ≡ {erg.ergebnis} (mod {m})")

    zeilen.append("")
    vorn = f"{name} = " if name else ""
    zeilen.append(f"Ergebnis: {vorn}{basis}{hoch(k)} ≡ {erg.ergebnis} (mod {m})")
    kontrolle = pow(a, k, m)
    ok = "stimmt" if kontrolle == erg.ergebnis else "FEHLER"
    zeilen.append(f"Probe (Kontrolle mit Python pow): {kontrolle} → {ok}")
    return "\n".join(zeilen)


def phi_text(m: int) -> str:
    """Kurze Rechnung für phi(m)."""
    if isprime(m):
        return f"{m} ist eine Primzahl, also φ({m}) = {m} − 1 = {m - 1}."
    fak = factorint(m)
    zerl = " · ".join(f"{p}{hoch(e) if e > 1 else ''}" for p, e in fak.items())
    teile = []
    for p, e in fak.items():
        teile.append(f"{p ** (e - 1) * (p - 1)}")
    einzel = " · ".join(
        f"φ({p}{hoch(e) if e > 1 else ''})" for p, e in fak.items()
    )
    rechnung = f"{einzel} = {' · '.join(teile)}"
    if len(teile) > 1:
        rechnung += f" = {phi(m)}"
    return f"{m} = {zerl}, also φ({m}) = {rechnung}."


def reduktion_text(erg: ReduktionErgebnis) -> str:
    """Lösungsweg für den Weg über Fermat/Euler (Stil Ü4 A7)."""
    a, k, m = erg.a, erg.k, erg.m
    basis = zk(a)
    zeilen = [f"Gesucht: {basis}{hoch(k)} mod {m}",
              "Annahme: zuerst den Exponenten mod φ(m) verkleinern (VL3 Bemerkung, Ü4 A7)."]
    zeilen.append(phi_text(m))
    if not erg.moeglich:
        zeilen.append(f"HINWEIS: ggT({a}, {m}) = {erg.ggt} ≠ 1.")
        zeilen.append("Der Satz von Euler gilt nur für ggT(a, m) = 1.")
        zeilen.append(f"Man darf {k} also NICHT durch {k} mod {erg.phi_m} ersetzen.")
        zeilen.append("(Warnbeispiel VL3: 12¹⁶ ≡ −16 ≢ 1 mod 34, denn ggT(12, 34) = 2.)")
        zeilen.append("Das Skript rechnet deshalb nur mit schnellem Potenzieren.")
        zeilen.append("")
        zeilen.append(potenz_text(erg.potenz))
        return "\n".join(zeilen)
    satz = "Fermat" if isprime(m) else "Euler"
    zeilen.append(f"ggT({a}, {m}) = 1. Nach {satz} ist {basis}{hoch(erg.phi_m)} ≡ 1 (mod {m}).")
    zeilen.append(f"{k} = {erg.q} · {erg.phi_m} + {erg.r}")
    zeilen.append(f"Also {basis}{hoch(k)} = ({basis}{hoch(erg.phi_m)}){hoch(erg.q)} · {basis}{hoch(erg.r)}"
                  f" ≡ 1{hoch(erg.q)} · {basis}{hoch(erg.r)} = {basis}{hoch(erg.r)} (mod {m})")
    zeilen.append("")
    if erg.r == 0:
        zeilen.append(f"Ergebnis: {basis}{hoch(k)} ≡ {erg.ergebnis} (mod {m})")
    else:
        zeilen.append(f"Jetzt {basis}{hoch(erg.r)} mod {m} schnell berechnen:")
        zeilen.append(potenz_text(erg.potenz))
        zeilen.append(f"Also Ergebnis: {basis}{hoch(k)} ≡ {erg.ergebnis} (mod {m})")
    kontrolle = pow(a, k, m)
    zeilen.append(f"Probe (direkt {basis}{hoch(k)} mit Python pow): {kontrolle} → "
                  f"{'stimmt' if kontrolle == erg.ergebnis else 'FEHLER'}")
    return "\n".join(zeilen)


def ordnung_text(erg: OrdnungErgebnis) -> str:
    """Lösungsweg für die Ordnung (Stil Ü4 A8)."""
    a, m = erg.a, erg.m
    if not erg.einheit:
        return (f"HINWEIS: ggT({a}, {m}) = {math.gcd(a, m)} ≠ 1. {a} ist keine Einheit in ℤ_{m}*.\n"
                f"Keine Potenz von {a} ist ≡ 1. {a} hat keine Ordnung und keine Inverse.")
    liste = ", ".join(z(x) for x in erg.potenzen)
    zeilen = [f"Potenzen {a}, {a}², {a}³, … mod {m}:",
              f"{liste}  ⇒  ord({a}) = {erg.ordnung}"]
    if erg.ordnung % 2 == 0 and erg.ordnung >= 2:
        halb = erg.ordnung // 2
        if erg.potenzen[halb - 1] % m == m - 1:
            zeilen.append(f"(Abkürzung: {a}{hoch(halb)} ≡ −1, also ist {a}{hoch(erg.ordnung)} = (−1)² = 1"
                          f" die erste 1.)")
    zeilen.append(f"ord({a}) teilt φ({m}) = {erg.phi_m}: {erg.phi_m} = {erg.phi_m // erg.ordnung} · {erg.ordnung}.")
    if erg.primitiv:
        zeilen.append(f"ord({a}) = φ({m}) = {erg.phi_m} ⇒ {a} ist primitiv (erzeugt ℤ_{m}*).")
    else:
        zeilen.append(f"ord({a}) = {erg.ordnung} ≠ φ({m}) = {erg.phi_m} ⇒ {a} ist nicht primitiv.")
    inv_liste = erg.potenzen[erg.ordnung - 2] if erg.ordnung >= 2 else 1
    inv_zeile = f"Inverse: {a}⁻¹ = {a}{hoch(erg.ordnung - 1)} ≡ {z(inv_liste)}"
    if inv_liste != erg.inverse:
        inv_zeile += f" ≡ {erg.inverse}"
    zeilen.append(inv_zeile + f" (mod {m})")
    zeilen.append(f"Probe: {a} · {erg.inverse} = {a * erg.inverse} ≡ {a * erg.inverse % m} (mod {m})")
    zeilen.append(f"Ergebnis: ord({a}) = {erg.ordnung}, primitiv: {'ja' if erg.primitiv else 'nein'}, "
                  f"{a}⁻¹ = {erg.inverse}")
    return "\n".join(zeilen)


def primitiv_text(m: int) -> str:
    """Liste der primitiven Elemente von Z_m* mit kurzer Begründung."""
    ph = phi(m)
    prim = primitive_elemente(m)
    zeilen = [phi_text(m) if m > 1 else "m muss ≥ 2 sein.",
              f"g ist primitiv, wenn ord(g) = φ({m}) = {ph} ist.",
              "Prüfregel: g ist primitiv, wenn g^(φ/p) ≢ 1 für jeden Primteiler p von φ."]
    if ph > 1:
        zeilen.append(f"Primteiler von {ph}: {', '.join(str(p) for p in factorint(ph))}")
    if prim:
        zeilen.append(f"Ergebnis: primitive Elemente von ℤ_{m}*: {', '.join(map(str, prim))} "
                      f"(Anzahl {len(prim)} = φ(φ({m})) = {phi(ph)})")
    else:
        zeilen.append(f"Ergebnis: ℤ_{m}* hat kein primitives Element. Die Gruppe ist nicht zyklisch.")
    return "\n".join(zeilen)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m skripte.schnell_potenzieren",
        description="Schnelles Potenzieren a^k mod m mit Lösungsweg. "
                    "Außerdem: Ordnung und primitive Elemente in Z_m*.",
        epilog="Beispiele aus den Quellen:\n"
               "  --a 12 --k 100 --m 34                 VL3: 12^100 ≡ 30 (mod 34)\n"
               "  --a 12 --k 100 --m 34 --reste symmetrisch   wie VL3 mit −4, −16\n"
               "  --a 3 --k 1000 --m 7 --variante euler  Ü4 A7: 3^1000 ≡ 4 (mod 7)\n"
               "  --a 715 --k 113 --m 2803              Ü4 A9: c = 708\n"
               "  --ordnung --a 2 3 4 5 --m 17          Ü4 A8: Ordnungen 8, 16, 4, 16\n"
               "  --primitiv --m 17                     alle primitiven Elemente von Z_17*",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--a", type=int, nargs="+", help="Basis a (bei --ordnung auch mehrere Zahlen).")
    p.add_argument("--k", type=int, help="Exponent k (ganze Zahl ≥ 0).")
    p.add_argument("--m", type=int, required=True, help="Modul m (ganze Zahl ≥ 2).")
    p.add_argument("--variante", choices=VARIANTEN, default="binaer",
                   help="binaer: nur schnelles Potenzieren (Standard). "
                        "euler: zuerst k mod φ(m) (nur bei ggT(a, m) = 1).")
    p.add_argument("--reste", choices=RESTE, default=None,
                   help="normal: Reste 0 … m−1 (Standard beim Potenzieren). "
                        "symmetrisch: Reste mit kleinstem Betrag (Standard bei --ordnung).")
    p.add_argument("--alle-varianten", action="store_true",
                   help="Alle Varianten und Rest-Schreibweisen nacheinander ausgeben.")
    p.add_argument("--ordnung", action="store_true",
                   help="Ordnung von a in Z_m* bestimmen, dazu primitiv ja/nein und Inverse.")
    p.add_argument("--primitiv", action="store_true",
                   help="Alle primitiven Elemente von Z_m* auflisten.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    m = args.m
    if m < 2:
        print(f"HINWEIS: Der Modul m = {m} ist zu klein. Bitte m ≥ 2 wählen.")
        return 1
    ausgaben: list[str] = []

    if args.primitiv:
        ausgaben.append(primitiv_text(m))

    if args.ordnung:
        if not args.a:
            print("HINWEIS: Für --ordnung bitte --a angeben, zum Beispiel --a 2 3 4 5.")
            return 1
        reste_liste = RESTE if args.alle_varianten else [args.reste or "symmetrisch"]
        for reste in reste_liste:
            ausgaben.append(f"=== Ordnungen in ℤ_{m}*, Reste: {reste} ===")
            ausgaben.append(phi_text(m))
            for a in args.a:
                ausgaben.append("")
                ausgaben.append(ordnung_text(ordnung(a, m, reste)))

    if args.k is not None or not (args.ordnung or args.primitiv):
        if not args.a or args.k is None:
            print("HINWEIS: Zum Potenzieren bitte --a und --k angeben, zum Beispiel --a 12 --k 100 --m 34.")
            return 1
        if args.k < 0:
            print(f"HINWEIS: Der Exponent k = {args.k} ist negativ. Bitte k ≥ 0 wählen.")
            print("Tipp: a^(−k) = (a⁻¹)^k. Berechne zuerst a⁻¹ (zum Beispiel mit --ordnung).")
            return 1
        if args.alle_varianten:
            kombis = [(v, r) for v in VARIANTEN for r in RESTE]
        else:
            kombis = [(args.variante, args.reste or "normal")]
        for a in args.a:
            for variante, reste in kombis:
                ausgaben.append(f"=== Variante: {variante}, Reste: {reste} ===")
                if variante == "euler":
                    ausgaben.append(reduktion_text(potenz_mit_reduktion(a, args.k, m, reste)))
                else:
                    ausgaben.append(potenz_text(schnell_potenzieren(a, args.k, m, reste)))
                ausgaben.append("")

    print("\n".join(ausgaben).rstrip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
