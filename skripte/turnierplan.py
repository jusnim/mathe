"""Turnierplan für ein Rundenturnier.

Zweck
-----
Ein Rundenturnier ist ein Turnier, in dem jede Mannschaft gegen jede andere
genau einmal spielt. Bei 2m Mannschaften braucht man 2m − 1 Runden
(Spieltage). Das Skript stellt den Plan nach dieser Regel auf
(Algorithmus Rundenturnier):

    Runde r: x spielt gegen y, wenn x + y ≡ r (mod 2m − 1), 1 ≤ x, y ≤ 2m − 1.
    Falls x + x ≡ r (mod 2m − 1), spielt x gegen 2m.

Es berechnet auch den Gegner von Mannschaft 2m in Runde r:

    2x ≡ r (mod 2m − 1) und 2⁻¹ ≡ m  ⇒  x ≡ r · m (mod 2m − 1).

Aufruf
------
    python -m skripte.turnierplan --mannschaften 8
    python -m skripte.turnierplan --mannschaften 8 --gegner
    python -m skripte.turnierplan --mannschaften 8 --runde 3
    python -m skripte.turnierplan --mannschaften 8 --alle-varianten

Varianten (--variante)
----------------------
standard   (Standard) Regel x + y ≡ r (mod 2m − 1). Die Restklasse 0 wird
           durch 2m − 1 dargestellt. Die Runden heißen 1 … 2m − 1.
rest-null  Gleiche Regel. Die Restklasse 0 wird durch 0 dargestellt.
           Die Runden heißen 0 … 2m − 2. Runde 0 ist Runde 2m − 1 der
           Standard-Variante. Der Rest 0 steht für Mannschaft 2m − 1.
summe-2m   Die Regel „x + y = 2m“. Sie ist ungenau. Wörtlich genommen
           liefert sie nur Runde 1. Das Skript zeigt das und stellt dann
           den richtigen Plan auf.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field

VARIANTEN = ("standard", "rest-null", "summe-2m")


# ---------------------------------------------------------------------------
# Rechnen
# ---------------------------------------------------------------------------


def rest(a: int, modul: int, variante: str = "standard") -> int:
    """Gibt den Repräsentanten von a mod modul zurück.

    standard/summe-2m: Werte 1 … modul (0 wird zu modul).
    rest-null: Werte 0 … modul − 1.
    """
    r = a % modul
    if variante != "rest-null" and r == 0:
        return modul
    return r


@dataclass
class Spiel:
    """Ein Spiel: zwei Mannschaften, dazu der Rechenweg."""

    a: int
    b: int
    summe: int  # x + y bzw. x + x
    sonderfall: bool  # True: x + x ≡ r, also Spiel gegen 2m


@dataclass
class Runde:
    nummer: int  # Name der Runde in der gewählten Variante
    r: int  # Wert von r in 1 … 2m − 1 (Standard)
    spiele: list[Spiel] = field(default_factory=list)


@dataclass
class Plan:
    n: int  # Anzahl der Mannschaften (gerade)
    m: int
    modul: int  # 2m − 1
    variante: str
    runden: list[Runde]
    n_eingabe: int  # ursprüngliche Eingabe (kann ungerade sein)


def runde_berechnen(n: int, r: int, variante: str = "standard") -> Runde:
    """Berechnet die Spiele von Runde r (r in 1 … n − 1) bei n Mannschaften."""
    modul = n - 1
    spiele: list[Spiel] = []
    for x in range(1, n):
        y = rest(r - x, modul)  # y in 1 … modul
        if y == x:
            spiele.append(Spiel(x, n, 2 * x, True))
        elif x < y:
            spiele.append(Spiel(x, y, x + y, False))
    spiele.sort(key=lambda s: s.a)
    nummer = r % modul if variante == "rest-null" else r
    return Runde(nummer, r, spiele)


def turnierplan(n: int, variante: str = "standard") -> Plan:
    """Stellt den ganzen Plan für n Mannschaften auf.

    Ist n ungerade, kommt eine Mannschaft n + 1 dazu ("spielfrei").
    """
    if variante not in VARIANTEN:
        raise ValueError(f"Unbekannte Variante: {variante}")
    if n < 2:
        raise ValueError("Es braucht mindestens 2 Mannschaften.")
    n_eingabe = n
    if n % 2:
        n += 1
    modul = n - 1
    runden = [runde_berechnen(n, r, variante) for r in range(1, modul + 1)]
    if variante == "rest-null":
        runden.sort(key=lambda rd: rd.nummer)
    return Plan(n, n // 2, modul, variante, runden, n_eingabe)


def paare(plan: Plan) -> dict[int, list[tuple[int, int]]]:
    """Gibt die Spiele als {Rundenname: [(a, b), ...]} zurück."""
    return {rd.nummer: [(s.a, s.b) for s in rd.spiele] for rd in plan.runden}


@dataclass
class GegnerSchritt:
    r: int  # Runde (Standard, 1 … 2m − 1)
    produkt: int  # r · m
    x: int  # Repräsentant in der gewählten Variante
    mannschaft: int  # tatsächliche Mannschaftsnummer (1 … 2m − 1)


def gegner_von_2m(n: int, r: int, variante: str = "standard") -> GegnerSchritt:
    """Gegner von Mannschaft n = 2m in Runde r: x ≡ r · m (mod 2m − 1)."""
    if n % 2 or n < 2:
        raise ValueError("Die Anzahl der Mannschaften muss gerade und ≥ 2 sein.")
    m = n // 2
    modul = n - 1
    if not 1 <= r <= modul and not (variante == "rest-null" and r == 0):
        raise ValueError(f"Runde r muss zwischen 1 und {modul} liegen.")
    if r == 0:
        r = modul
    produkt = r * m
    x = rest(produkt, modul, variante)
    return GegnerSchritt(r, produkt, x, rest(produkt, modul))


def probe(plan: Plan) -> list[str]:
    """Prüft den Plan. Gibt eine Liste von Fehlern zurück (leer = alles gut)."""
    fehler: list[str] = []
    gesehen: set[tuple[int, int]] = set()
    for rd in plan.runden:
        teams = [t for s in rd.spiele for t in (s.a, s.b)]
        if sorted(teams) != list(range(1, plan.n + 1)):
            fehler.append(f"Runde {rd.nummer}: nicht jede Mannschaft spielt genau einmal.")
        for s in rd.spiele:
            p = (min(s.a, s.b), max(s.a, s.b))
            if p in gesehen:
                fehler.append(f"Spiel {p} kommt doppelt vor.")
            gesehen.add(p)
    if len(gesehen) != plan.n * (plan.n - 1) // 2:
        fehler.append("Nicht alle Paare spielen gegeneinander.")
    return fehler


# ---------------------------------------------------------------------------
# Ausgeben
# ---------------------------------------------------------------------------


def _annahme(variante: str, modul: int) -> str:
    if variante == "rest-null":
        return (f"Annahme (Variante rest-null): Reste liegen in 0 … {modul - 1}. "
                f"Runde 0 entspricht Runde {modul} der Standard-Variante. "
                f"Rest 0 steht für Mannschaft {modul}.")
    if variante == "summe-2m":
        return ("Annahme (Variante summe-2m): Regel „x + y = 2m“. "
                "Diese Regel ist ungenau. Unten steht der richtige Weg.")
    return (f"Annahme (Variante standard): Reste liegen in 1 … {modul}. "
            f"Der Rest 0 wird als {modul} geschrieben.")


def _tabelle(plan: Plan) -> list[str]:
    breite = max(7, len(f"r={plan.runden[-1].nummer}") + 2, len(str(plan.n)) * 2 + 3)
    kopf = "".join(f"r={rd.nummer}".center(breite) for rd in plan.runden)
    zeilen = [kopf, "-" * len(kopf)]
    for i in range(plan.m):
        zeilen.append("".join(
            f"{rd.spiele[i].a} {rd.spiele[i].b}".center(breite) for rd in plan.runden))
    return zeilen


def text_plan(plan: Plan) -> str:
    n, m, modul, v = plan.n, plan.m, plan.modul, plan.variante
    z: list[str] = []
    z.append(f"Turnierplan für 2m = {n} Mannschaften an 2m − 1 = {modul} Spieltagen")
    z.append("=" * len(z[-1]))
    z.append(_annahme(v, modul))
    if plan.n_eingabe != n:
        z.append(f"Hinweis: {plan.n_eingabe} ist ungerade. Ein Rundenturnier braucht "
                 f"eine gerade Anzahl. Ich füge Mannschaft {n} hinzu. "
                 f"Wer gegen {n} spielt, hat in dieser Runde spielfrei.")
    z.append("")

    if v == "summe-2m":
        z.append("Regel „x + y = 2m“ wörtlich: x + y = 2m = %d für x ≠ y, sonst x gegen 2m." % n)
        erste = [f"{x} {n - x}" for x in range(1, m) ] + [f"{m} {n}"]
        z.append("  Das ergibt nur eine Runde: " + ", ".join(erste))
        z.append(f"  Grund: 2m = {n} ≡ 1 (mod {modul}). Die Regel ist also nur Runde r = 1.")
        z.append("  Ein Plan braucht aber alle Runden. Sie folgen aus der Regel")
        z.append("  x + y ≡ r (mod 2m − 1). Diese Regel nutze ich jetzt.")
        z.append("")

    z.append("Regel (Algorithmus Rundenturnier):")
    z.append(f"  Runde r: x spielt gegen y, wenn x + y ≡ r (mod {modul}) für 1 ≤ x, y ≤ {modul}.")
    z.append(f"  Falls x + x ≡ r (mod {modul}), spielt x gegen 2m = {n}.")
    z.append("")
    z.append("Rechenweg je Runde:")
    for rd in plan.runden:
        teile = []
        for s in rd.spiele:
            if s.sonderfall:
                teile.append(f"{s.a}+{s.a}={s.summe}≡{rest(s.summe, modul, v)} ⇒ {s.a} gegen {n}")
            else:
                teile.append(f"{s.a}+{s.b}={s.summe}≡{rest(s.summe, modul, v)}")
        z.append(f"  r = {rd.nummer}: " + ", ".join(teile))
    z.append("")
    z.append("Tabelle:")
    z.extend("  " + zeile for zeile in _tabelle(plan))
    z.append("")
    fehler = probe(plan)
    anzahl = n * (n - 1) // 2
    if fehler:
        z.append("Probe: FEHLER")
        z.extend("  " + f for f in fehler)
    else:
        z.append(f"Probe: In jeder Runde spielt jede der {n} Mannschaften genau einmal.")
        z.append(f"       Es gibt {modul} Runden mit je {m} Spielen: {modul} · {m} = {anzahl}"
                 f" = ({n} über 2). Jedes Paar spielt genau einmal.")
    z.append("")
    z.append(f"Ergebnis: Turnierplan für {n} Mannschaften, siehe Tabelle oben.")
    return "\n".join(z)


def text_gegner(n: int, runden: list[int], variante: str = "standard") -> str:
    m = n // 2
    modul = n - 1
    z: list[str] = []
    z.append(f"Gegner von Mannschaft 2m = {n} in Runde r")
    z.append("=" * len(z[-1]))
    z.append(_annahme("standard" if variante == "summe-2m" else variante, modul))
    z.append(f"  x + x = 2 · x ≡ r        (mod {modul})")
    z.append(f"  2m = {n} ≡ 1              (mod {modul})")
    z.append(f"  ⇒ 2⁻¹ ≡ m = {m}           (mod {modul})")
    z.append(f"  ⇒ x ≡ r · m = r · {m}     (mod {modul})")
    z.append("")
    ergebnisse = []
    for r in runden:
        g = gegner_von_2m(n, r, variante)
        name = g.r % modul if variante == "rest-null" else g.r
        produkt = name * m
        zeile = f"  r = {name}:  x = {name} · {m} = {produkt}"
        if produkt != g.x:
            zeile += f" ≡ {g.x} (mod {modul})"
        if g.x != g.mannschaft:
            zeile += f"  → Mannschaft {g.mannschaft}"
        z.append(zeile)
        ergebnisse.append((name, g.mannschaft))
    z.append("")
    z.append("Probe: 2 · x ≡ r?")
    for r in runden:
        g = gegner_von_2m(n, r, variante)
        z.append(f"  2 · {g.mannschaft} = {2 * g.mannschaft} ≡ "
                 f"{rest(2 * g.mannschaft, modul, variante)} = r (mod {modul})  ✓")
    z.append("")
    z.append("Ergebnis: " + "; ".join(f"Runde {r}: {n} gegen {x}" for r, x in ergebnisse))
    return "\n".join(z)


def loesung(n: int, variante: str, gegner: bool, runde: int | None) -> str:
    plan = turnierplan(n, variante)
    teile = [text_plan(plan)]
    if runde is not None:
        teile.append(text_gegner(plan.n, [runde], variante))
    elif gegner:
        runden = [rd.nummer for rd in plan.runden] if variante == "rest-null" else \
            list(range(1, plan.modul + 1))
        teile.append(text_gegner(plan.n, runden, variante))
    return "\n\n".join(teile)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m skripte.turnierplan",
        description="Stellt einen Turnierplan (Rundenturnier) für 2m Mannschaften auf. "
                    "Zeigt auch den Gegner von Mannschaft 2m in Runde r.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Beispiele:\n"
               "  python -m skripte.turnierplan --mannschaften 8   (Plan für 8 Mannschaften)\n"
               "  python -m skripte.turnierplan --mannschaften 8 --gegner   (mit Gegner von Mannschaft 8)\n"
               "  python -m skripte.turnierplan --mannschaften 8 --runde 3   (Gegner von Mannschaft 8 in Runde 3)\n"
               "  python -m skripte.turnierplan --mannschaften 7   (ungerade Anzahl, mit spielfrei)\n"
               "  python -m skripte.turnierplan --mannschaften 8 --alle-varianten   (alle Varianten)",
    )
    parser.add_argument("--mannschaften", "-n", type=int, required=True,
                        help="Anzahl 2m der Mannschaften, zum Beispiel 8.")
    parser.add_argument("--gegner", action="store_true",
                        help="Zeigt für jede Runde den Gegner von Mannschaft 2m.")
    parser.add_argument("--runde", "-r", type=int,
                        help="Zeigt den Gegner von Mannschaft 2m nur in dieser Runde r.")
    parser.add_argument("--variante", choices=VARIANTEN, default="standard",
                        help="Konvention: standard (Standard, Reste 1 … 2m − 1), "
                             "rest-null (Reste 0 … 2m − 2), "
                             "summe-2m (Regel „x + y = 2m“, mit Korrektur).")
    parser.add_argument("--alle-varianten", action="store_true",
                        help="Gibt den Lösungsweg für alle Varianten nacheinander aus.")
    args = parser.parse_args(argv)

    n = args.mannschaften
    if n < 2:
        print(f"Hinweis: {n} Mannschaften sind zu wenig. Es braucht mindestens 2.")
        return
    n_gerade = n + (n % 2)
    if args.runde is not None:
        lo = 0 if args.variante == "rest-null" and not args.alle_varianten else 1
        if not lo <= args.runde <= n_gerade - 1:
            print(f"Hinweis: Runde {args.runde} gibt es nicht. "
                  f"Erlaubt ist r = {lo} … {n_gerade - 1}.")
            return
    varianten = VARIANTEN if args.alle_varianten else (args.variante,)
    for i, v in enumerate(varianten):
        if len(varianten) > 1:
            print(("\n" if i else "") + f"##### Variante: {v} #####\n")
        print(loesung(n, v, args.gegner, args.runde))


if __name__ == "__main__":
    main()
