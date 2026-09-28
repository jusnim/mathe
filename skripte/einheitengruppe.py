"""Einheitengruppe ℤ_m*: Einheiten, Nullteiler, Inverse, Gruppentafel, Ordnungen.

Zweck
-----
Das Skript berechnet die Einheitengruppe ℤ_m* und zeigt jeden Rechenschritt.

Begriffe
--------
- Einheit: ein Element k von ℤ_m mit einem Inversen. Das gilt genau bei ggT(k, m) = 1.
- ℤ_m*: die Menge aller Einheiten. Sie ist eine Gruppe mit |ℤ_m*| = φ(m).
- Nullteiler: ein Element a ≠ 0 mit a·b ≡ 0 für ein b ≠ 0. Jede Nicht-Einheit a ≠ 0 ist ein Nullteiler.
- Gruppentafel: Tabelle mit allen Produkten a·b in ℤ_m*.
- Ordnung ord(a): die kleinste Zahl k > 0 mit a^k ≡ 1.
- primitiv: a ist primitiv, wenn ord(a) = φ(m) ist. Dann erzeugt a die ganze Gruppe.

Aufruf (Beispiele)
------------------
    python -m skripte.einheitengruppe --m 12                          (alles zu ℤ_12*)
    python -m skripte.einheitengruppe --m 18 --modus tafel            (Gruppentafel von ℤ_18*)
    python -m skripte.einheitengruppe --m 17 --modus ordnungen --elemente 2 3 4 5   (Ordnungen in ℤ_17*)
    python -m skripte.einheitengruppe --m 11 --modus inverse --reihenfolge zahlen   (Inverse mod 11)
    python -m skripte.einheitengruppe --m 6 --modus tafel --mit-nullteilern --variante standard  (Tafel von ℤ_6)
    python -m skripte.einheitengruppe --m 7 --modus ordnungen --variante standard   (Ordnungen in ℤ_7*)
    python -m skripte.einheitengruppe --m 18 --modus inverse --inverse-weg 7        (Inverse mit EA)

Modi (--modus)
--------------
einheiten, inverse, tafel, ordnungen oder alles (Standard).

Varianten (--variante)
----------------------
- symmetrisch (Standard): Repräsentanten −m/2 … m/2, zum Beispiel ℤ_12* = {±1, ±5}.
- standard: Repräsentanten 0 … m−1, zum Beispiel ℤ_12* = {1, 5, 7, 11}.
Beide Varianten sind richtig. Sie beschreiben dieselben Restklassen.
--alle-varianten gibt beide Varianten nacheinander aus.

Reihenfolge (--reihenfolge)
---------------------------
- paare (Standard): 1, −1, 5, −5, …
- zahlen: 1, 2, 3, …, m−1
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass

from skripte.euklid import KeinInversesFehler, ggt, inverse, inverse_mit_weg
from skripte.primfaktor_phi import phi
from skripte.schnell_potenzieren import ordnung, primitive_elemente

VARIANTEN = ("symmetrisch", "standard")
REIHENFOLGEN = ("paare", "zahlen")
MODI = ("einheiten", "inverse", "tafel", "ordnungen", "alles")
MAX_TAFEL = 40   # größere Tafeln passen nicht auf Papier

ANNAHME = {
    "symmetrisch": "Repräsentanten symmetrisch (−m/2 … m/2).",
    "standard": "Repräsentanten 0 … m−1.",
}


# ---------------------------------------------------------------------------
# Schreibweise
# ---------------------------------------------------------------------------

def z(x: int) -> str:
    """Zahl mit echtem Minuszeichen."""
    return str(x).replace("-", "−")


def zk(x: int) -> str:
    """Negative Zahl in Klammern: (−7)."""
    return f"({z(x)})" if x < 0 else z(x)


def hoch(n: int) -> str:
    return str(n).translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))


def rep(x: int, m: int, variante: str) -> int:
    """Repräsentant von x mod m nach der Variante."""
    r = x % m
    if variante == "symmetrisch" and r > m // 2:
        r -= m
    return r


def menge_text(werte: list[int], m: int, variante: str) -> str:
    """Menge als Text. Bei symmetrisch werden a und −a zu ±a zusammengefasst."""
    if variante != "symmetrisch":
        return "{" + ", ".join(z(w) for w in werte) + "}"
    teile, gesehen = [], set()
    for w in werte:
        if w in gesehen:
            continue
        if w > 0 and -w in werte and w != -w % m:
            teile.append(f"±{w}")
            gesehen.update({w, -w})
        else:
            teile.append(z(w))
            gesehen.add(w)
    return "{" + ", ".join(teile) + "}"


# ---------------------------------------------------------------------------
# Rechnen
# ---------------------------------------------------------------------------

def pruefe_modul(m: int) -> None:
    if m < 2:
        raise ValueError(f"Der Modul m = {m} ist zu klein. Bitte m ≥ 2 wählen.")


def ordne(reste: list[int], m: int, variante: str, reihenfolge: str) -> list[int]:
    """Reste (0 … m−1) als Repräsentanten in der gewünschten Reihenfolge."""
    werte = [rep(r, m, variante) for r in sorted(reste)]
    if reihenfolge == "paare" and variante == "symmetrisch":
        werte.sort(key=lambda w: (abs(w), w < 0))
    return werte


@dataclass
class EinheitenErgebnis:
    """Einheiten und Nullteiler von ℤ_m."""
    m: int
    phi_m: int
    ggts: list[tuple[int, int]]            # (k, ggT(k, m)) für k = 1 … m−1
    einheiten: list[int]                   # Reste 0 … m−1
    nullteiler: list[tuple[int, int]]      # (a, b) mit a·b ≡ 0, b ≠ 0


def einheiten(m: int) -> EinheitenErgebnis:
    """Bestimmt ℤ_m* mit dem ggT-Test. Zu jedem Nullteiler a gibt es b = m/ggT(a, m)."""
    pruefe_modul(m)
    ggts = [(k, ggt(k, m)) for k in range(1, m)]
    ein = [k for k, g in ggts if g == 1]
    null = [(k, m // g) for k, g in ggts if g != 1]
    return EinheitenErgebnis(m, phi(m), ggts, ein, null)


@dataclass
class InversZeile:
    a: int          # Repräsentant
    inv: int        # Repräsentant von a⁻¹
    produkt: int    # a·inv als ganze Zahl


def inversen(m: int, variante: str = "symmetrisch", reihenfolge: str = "paare",
             elemente: list[int] | None = None) -> list[InversZeile]:
    """Inversentabelle von ℤ_m*. Wirft KeinInversesFehler für eine Nicht-Einheit in elemente."""
    pruefe_modul(m)
    if elemente is None:
        reste = einheiten(m).einheiten
        werte = ordne(reste, m, variante, reihenfolge)
    else:
        werte = [rep(a, m, variante) for a in elemente]
    zeilen = []
    for a in werte:
        b = rep(inverse(a, m), m, variante)
        zeilen.append(InversZeile(a, b, a * b))
    return zeilen


@dataclass
class Tafel:
    m: int
    elemente: list[int]          # Kopfzeile und Kopfspalte (Repräsentanten)
    werte: list[list[int]]       # werte[i][j] = elemente[i]·elemente[j] als Repräsentant
    mit_nullteilern: bool


def gruppentafel(m: int, variante: str = "symmetrisch", reihenfolge: str = "paare",
                 mit_nullteilern: bool = False) -> Tafel:
    """Gruppentafel von ℤ_m* (oder von ℤ_m ∖ {0}, wenn mit_nullteilern)."""
    pruefe_modul(m)
    reste = list(range(1, m)) if mit_nullteilern else einheiten(m).einheiten
    el = ordne(reste, m, variante, reihenfolge)
    werte = [[rep(a * b, m, variante) for b in el] for a in el]
    return Tafel(m, el, werte, mit_nullteilern)


@dataclass
class OrdnungZeile:
    a: int
    einheit: bool
    potenzen: list[int]      # a, a², …, a^ord = 1 (Repräsentanten)
    ordnung: int | None
    inv: int | None          # a^(ord−1) als Repräsentant
    primitiv: bool


@dataclass
class OrdnungenErgebnis:
    m: int
    phi_m: int
    zeilen: list[OrdnungZeile]
    primitive: list[int]     # alle primitiven Elemente (Repräsentanten 0 … m−1)


def ordnungen(m: int, variante: str = "symmetrisch", reihenfolge: str = "paare",
              elemente: list[int] | None = None) -> OrdnungenErgebnis:
    """Ordnung jedes Elements (oder der gewählten Elemente) und alle primitiven Elemente."""
    pruefe_modul(m)
    if elemente is None:
        elemente = ordne(einheiten(m).einheiten, m, variante, reihenfolge)
    reste = "symmetrisch" if variante == "symmetrisch" else "normal"
    zeilen = []
    for a in elemente:
        o = ordnung(a % m, m, reste)
        if not o.einheit:
            zeilen.append(OrdnungZeile(a, False, [], None, None, False))
            continue
        zeilen.append(OrdnungZeile(a, True, o.potenzen, o.ordnung,
                                   rep(o.inverse, m, variante), o.primitiv))
    return OrdnungenErgebnis(m, phi(m), zeilen, primitive_elemente(m))


# ---------------------------------------------------------------------------
# Ausgeben
# ---------------------------------------------------------------------------

def text_einheiten(e: EinheitenErgebnis, variante: str = "symmetrisch",
                   reihenfolge: str = "paare") -> str:
    m = e.m
    z_ = [f"Einheiten von ℤ_{m}",
          f"Regel: k ist Einheit ⇔ ggT(k, {m}) = 1. Sonst ist k ≠ 0 ein Nullteiler.",
          "Nullteiler: Es gibt b ≠ 0 mit k·b ≡ 0. Wähle b = m / ggT(k, m).",
          "",
          f"{'k':>4} | ggT(k, {m}) | Ergebnis"]
    z_.append("-" * 44)
    null = dict(e.nullteiler)
    for k, g in e.ggts:
        if g == 1:
            erg = "Einheit"
        else:
            b = null[k]
            erg = f"Nullteiler: {k}·{b} = {k * b} ≡ 0"
        z_.append(f"{k:>4} | {g:>{len(f'ggT(k, {m})')}} | {erg}")
    z_.append("")
    werte = ordne(e.einheiten, m, variante, reihenfolge)
    z_.append(f"ℤ_{m}* = {menge_text(werte, m, variante)}")
    if variante == "symmetrisch":
        z_.append(f"      = {menge_text(e.einheiten, m, 'standard')}   (in 0 … {m - 1})")
    z_.append(f"Probe: |ℤ_{m}*| = {len(e.einheiten)} = φ({m}) = {e.phi_m} "
              + ("✓" if len(e.einheiten) == e.phi_m else "✗ FEHLER"))
    if e.nullteiler:
        z_.append(f"Nullteiler: {', '.join(str(k) for k, _ in e.nullteiler)}")
        z_.append(f"Ergebnis: ℤ_{m}* = {menge_text(werte, m, variante)}. "
                  f"Alle anderen Elemente ≠ 0 sind Nullteiler. ℤ_{m} ist kein Körper.")
    else:
        z_.append(f"Ergebnis: ℤ_{m}* = ℤ_{m} ∖ {{0}}. Es gibt keine Nullteiler. "
                  f"{m} ist eine Primzahl, also ist ℤ_{m} ein Körper.")
    return "\n".join(z_)


def text_inversen(m: int, zeilen: list[InversZeile]) -> str:
    z_ = [f"Inverse in ℤ_{m}*",
          f"Gesucht ist zu jedem a ein b mit a·b ≡ 1 (mod {m}). Dann ist b = a⁻¹.",
          ""]
    for zl in zeilen:
        p = zl.produkt
        q, r = divmod(p, m)
        z_.append(f"{zk(zl.a)}·{zk(zl.inv)} = {z(p)} = {zk(q)}·{m} + {r} ≡ 1  ⇒  {zk(zl.a)}⁻¹ = {z(zl.inv)}")
    z_.append("")
    breite = max(len(z(x)) for zl in zeilen for x in (zl.a, zl.inv))
    breite = max(breite, 3)
    z_.append(f"{'a':>4} | " + " ".join(f"{z(zl.a):>{breite}}" for zl in zeilen))
    z_.append("-" * (7 + (breite + 1) * len(zeilen)))
    z_.append(f"{'a⁻¹':>4} | " + " ".join(f"{z(zl.inv):>{breite}}" for zl in zeilen))
    selbst = [z(zl.a) for zl in zeilen if zl.a == zl.inv]
    if selbst:
        z_.append(f"Selbstinvers (a⁻¹ = a): {', '.join(selbst)}")
    z_.append("Ergebnis: " + ", ".join(f"{zk(zl.a)}⁻¹ = {z(zl.inv)}" for zl in zeilen))
    return "\n".join(z_)


def text_tafel(t: Tafel) -> str:
    m = t.m
    name = f"ℤ_{m} ∖ {{0}}" if t.mit_nullteilern else f"ℤ_{m}*"
    b = max(max(len(z(x)) for x in t.elemente), 3)
    kopf_b = max(len(name), b)
    z_ = [f"Gruppentafel von ({name}, ·)" if not t.mit_nullteilern
          else f"Multiplikationstafel von ({name}, ·)",
          f"Eintrag in Zeile a, Spalte b: a·b mod {m}.", ""]
    z_.append(f"{name:<{kopf_b}} | " + " ".join(f"{z(x):>{b}}" for x in t.elemente))
    z_.append("-" * (kopf_b + 3 + (b + 1) * len(t.elemente)))
    for a, zeile in zip(t.elemente, t.werte):
        z_.append(f"{z(a):<{kopf_b}} | " + " ".join(f"{z(x):>{b}}" for x in zeile))
    z_.append("")
    if t.mit_nullteilern:
        null = sorted({a for a, zeile in zip(t.elemente, t.werte) if 0 in zeile}, key=lambda x: x % m)
        ein = [a for a in t.elemente if 1 in [w % m for w in t.werte[t.elemente.index(a)]]]
        z_.append("Eine 0 in der Zeile zeigt einen Nullteiler. Eine 1 zeigt eine Einheit.")
        z_.append(f"Ergebnis: Einheiten {', '.join(z(x) for x in ein)}; "
                  f"Nullteiler {', '.join(z(x) for x in null) or 'keine'}. "
                  + ("Keine Gruppe (Axiom G4 fehlt)." if null else "Das ist eine Gruppe."))
    else:
        z_.append("Kontrolle: In jeder Zeile und Spalte steht jedes Element genau einmal.")
        ok = all(sorted(zeile) == sorted(t.elemente) for zeile in t.werte)
        z_.append(f"Ergebnis: Gruppentafel von ℤ_{m}* mit {len(t.elemente)} Elementen. "
                  + ("Kontrolle erfüllt ✓" if ok else "Kontrolle NICHT erfüllt ✗"))
    return "\n".join(z_)


def text_ordnungen(o: OrdnungenErgebnis, variante: str = "symmetrisch") -> str:
    m = o.m
    z_ = [f"Ordnungen in ℤ_{m}*",
          f"|ℤ_{m}*| = φ({m}) = {o.phi_m}. Jede Ordnung teilt {o.phi_m} (Satz von Euler/Lagrange).",
          f"Weg: Potenzen a, a², a³, … mod {m} aufschreiben, bis 1 kommt.",
          "Inverse: a⁻¹ = a^(ord(a)−1), das ist die vorletzte Zahl der Liste.",
          ""]
    for zl in o.zeilen:
        a = zl.a
        if not zl.einheit:
            z_.append(f"HINWEIS: ggT({a}, {m}) = {ggt(a, m)} ≠ 1. {a} ist keine Einheit in ℤ_{m}*. "
                      f"{a} hat keine Ordnung und kein Inverses.")
            continue
        liste = ", ".join(z(x) for x in zl.potenzen)
        z_.append(f"{zk(a)}: {liste}  ⇒  ord({z(a)}) = {zl.ordnung}  ⇒  "
                  f"{zk(a)}⁻¹ = {zk(a)}{hoch(zl.ordnung - 1)} = {z(zl.inv)}"
                  + ("   (primitiv)" if zl.primitiv else ""))
    z_.append("")
    gueltig = [zl for zl in o.zeilen if zl.einheit]
    for zl in gueltig:
        p = zl.a * zl.inv
        q, r = divmod(p, m)
        z_.append(f"Probe: {zk(zl.a)}·{zk(zl.inv)} = {z(p)} = {zk(q)}·{m} + {r} ≡ 1 (mod {m})")
    prim_gew = [z(zl.a) for zl in gueltig if zl.primitiv]
    z_.append("")
    z_.append(f"Primitiv heißt: ord(a) = φ({m}) = {o.phi_m}.")
    z_.append(f"Primitiv unter den Elementen oben: {', '.join(prim_gew) if prim_gew else 'keines'}")
    if o.primitive:
        alle = [rep(g, m, variante) for g in o.primitive]
        z_.append(f"Alle primitiven Elemente von ℤ_{m}*: {', '.join(z(g) for g in alle)} "
                  f"(Anzahl φ(φ({m})) = φ({o.phi_m}) = {len(o.primitive)}). ℤ_{m}* ist zyklisch.")
    else:
        z_.append(f"ℤ_{m}* hat kein primitives Element. Die Gruppe ist nicht zyklisch.")
    z_.append("Ergebnis: " + "; ".join(
        f"ord({z(zl.a)}) = {zl.ordnung}, {zk(zl.a)}⁻¹ = {z(zl.inv)}" for zl in gueltig)
        + f"; primitiv: {', '.join(prim_gew) if prim_gew else 'keines'}")
    return "\n".join(z_)


def loesungsweg(m: int, modus: str = "alles", variante: str = "symmetrisch",
                reihenfolge: str = "paare", elemente: list[int] | None = None,
                mit_nullteilern: bool = False, inverse_weg: list[int] | None = None) -> str:
    """Der ganze Lösungsweg als Text für einen Modus und eine Variante."""
    teile = [f"=== ℤ_{m}*, Variante: {variante} ===", f"Annahme: {ANNAHME[variante]}"]
    if elemente:
        fremd = [a for a in elemente if ggt(a, m) != 1]
        for a in fremd:
            teile.append(f"HINWEIS: ggT({a}, {m}) = {ggt(a, m)} ≠ 1. {a} ist keine Einheit. "
                         f"{a} fehlt deshalb bei Inversen und Ordnungen.")
    einheits_el = [a for a in elemente if ggt(a, m) == 1] if elemente else None
    if elemente and not einheits_el and modus in ("inverse", "ordnungen"):
        return "\n".join(teile)
    if modus in ("einheiten", "alles"):
        teile += ["", text_einheiten(einheiten(m), variante, reihenfolge)]
    if modus in ("inverse", "alles"):
        teile += ["", text_inversen(m, inversen(m, variante, reihenfolge, einheits_el))]
    if inverse_weg:
        for a in inverse_weg:
            teile.append("")
            try:
                teile.append(inverse_mit_weg(a, m, "symmetrisch" if variante == "symmetrisch"
                                             else "standard").text)
            except KeinInversesFehler as f:
                teile.append(f.text or str(f))
    if modus in ("tafel", "alles"):
        t = gruppentafel(m, variante, reihenfolge, mit_nullteilern)
        if len(t.elemente) > MAX_TAFEL:
            teile += ["", f"HINWEIS: Die Tafel hätte {len(t.elemente)}×{len(t.elemente)} Einträge. "
                          f"Das ist zu groß für Papier. Die Tafel wird nicht ausgegeben."]
        else:
            teile += ["", text_tafel(t)]
    if modus in ("ordnungen", "alles"):
        teile += ["", text_ordnungen(ordnungen(m, variante, reihenfolge, einheits_el), variante)]
    return "\n".join(teile)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m skripte.einheitengruppe",
        description="Einheitengruppe ℤ_m*: Einheiten, Nullteiler, Inverse, Gruppentafel, "
                    "Ordnungen und primitive Elemente mit Lösungsweg.",
        epilog="Beispiele:\n"
               "  python -m skripte.einheitengruppe --m 12  (alles zu ℤ_12*)\n"
               "  python -m skripte.einheitengruppe --m 18 --modus tafel  (Gruppentafel von ℤ_18*)\n"
               "  python -m skripte.einheitengruppe --m 17 --modus ordnungen --elemente 2 3 4 5  (Ordnungen in ℤ_17*)\n"
               "  python -m skripte.einheitengruppe --m 11 --modus inverse --reihenfolge zahlen  (Inverse mod 11)\n"
               "  python -m skripte.einheitengruppe --m 6 --modus tafel --mit-nullteilern --variante standard  (Tafel von ℤ_6 mit Nullteilern)\n"
               "  python -m skripte.einheitengruppe --m 7 --modus ordnungen --variante standard  (Ordnungen in ℤ_7*)\n"
               "  python -m skripte.einheitengruppe --m 18 --modus inverse --inverse-weg 7  (Inverse mit dem EA)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--m", type=int, required=True, help="Modul m (ganze Zahl ≥ 2).")
    p.add_argument("--modus", choices=MODI, default="alles",
                   help="Was soll berechnet werden? einheiten, inverse, tafel, ordnungen "
                        "oder alles (Standard).")
    p.add_argument("--variante", choices=VARIANTEN, default="symmetrisch",
                   help="symmetrisch: Repräsentanten −m/2 … m/2 (Standard). "
                        "standard: Repräsentanten 0 … m−1.")
    p.add_argument("--alle-varianten", action="store_true",
                   help="Beide Varianten nacheinander ausgeben.")
    p.add_argument("--reihenfolge", choices=REIHENFOLGEN, default="paare",
                   help="paare: 1, −1, 5, −5, … (Standard). "
                        "zahlen: 1, 2, …, m−1.")
    p.add_argument("--elemente", type=int, nargs="+",
                   help="Nur diese Elemente bei Inversen und Ordnungen, zum Beispiel 2 3 4 5.")
    p.add_argument("--mit-nullteilern", action="store_true",
                   help="Tafel für ℤ_m ∖ {0} mit allen Elementen, auch den Nullteilern.")
    p.add_argument("--inverse-weg", type=int, nargs="+", metavar="A",
                   help="Zusätzlich a⁻¹ mit dem Euklidischen Algorithmus (EA rückwärts) zeigen.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.m < 2:
        print(f"HINWEIS: Der Modul m = {args.m} ist zu klein. Bitte m ≥ 2 wählen.")
        return 1
    varianten = VARIANTEN if args.alle_varianten else (args.variante,)
    ausgaben = [loesungsweg(args.m, args.modus, v, args.reihenfolge, args.elemente,
                            args.mit_nullteilern, args.inverse_weg) for v in varianten]
    print("\n\n".join(ausgaben))
    return 0


if __name__ == "__main__":
    sys.exit(main())
