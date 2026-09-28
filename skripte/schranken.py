"""Hamming-Schranke und Singleton-Schranke für lineare [n, k, d]_q-Codes.

Zweck
-----
Das Skript prüft für einen linearen Code zwei obere Schranken:
die Hamming-Schranke und die Singleton-Schranke.

Ein [n, k, d]_q-Code hat Wortlänge n, Dimension k und Mindestabstand d.
Die Buchstaben kommen aus dem Körper F_q mit q Elementen.

- Hamming-Schranke: |C| · Σ_{i=0}^{e} C(n, i)·(q − 1)^i ≤ q^n.
  Gilt "=", dann heißt der Code perfekt.
- Singleton-Schranke: k + d ≤ n + 1.
  Gilt "=", dann heißt der Code MDS-Code. Man nennt ihn auch "optimal".

Aufruf
------
Einen Code prüfen:
    python -m skripte.schranken --n 10 --k 8 --d 3 --q 11

Kleinstes n bestimmen, dazu --n weglassen:
    python -m skripte.schranken --k 5 --d 5 --q 2

Varianten
---------
Es gibt keine Konvention, die zu anderen Zahlen führt. Deshalb gibt es keine
Option --variante. Die Ausgabe nennt trotzdem die Annahmen:
- e kommt aus 2e + 1 = d. Bei geradem d hat diese Gleichung
  keine ganze Lösung. Dann gilt e = ⌊(d − 1)/2⌋.
- "optimal" heißt hier: Singleton-Schranke mit "=" (MDS).
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from math import comb, prod

from sympy import factorint


# ---------------------------------------------------------------------------
# Rechnen
# ---------------------------------------------------------------------------

@dataclass
class Summand:
    """Ein Summand C(n, i)·(q − 1)^i der Kugelgröße |B(e)|."""

    i: int
    binom: int      # C(n, i)
    faktor: int     # (q − 1)^i
    wert: int       # C(n, i)·(q − 1)^i


@dataclass
class HammingErgebnis:
    n: int
    k: int
    d: int
    q: int
    e: int
    d_gerade: bool
    summanden: list[Summand]
    kugel: int          # |B(e)| = Σ C(n, i)(q − 1)^i
    code_groesse: int   # |C| = q^k
    rest: int           # q^(n − k)
    erfuellt: bool      # |B(e)| ≤ q^(n − k)
    perfekt: bool       # |B(e)| = q^(n − k)


@dataclass
class SingletonErgebnis:
    n: int
    k: int
    d: int
    links: int          # k + d
    rechts: int         # n + 1
    erfuellt: bool
    mds: bool


@dataclass
class KleinstesNZeile:
    n: int
    kugel: int
    rest: int           # q^(n − k)
    erfuellt: bool


@dataclass
class KleinstesNErgebnis:
    k: int
    d: int
    q: int
    e: int
    n_singleton: int
    zeilen: list[KleinstesNZeile] = field(default_factory=list)
    n_hamming: int = 0
    n_min: int = 0


def e_aus_d(d: int) -> int:
    """Gibt e = ⌊(d − 1)/2⌋ zurück. Für ungerades d gilt 2e + 1 = d."""
    return (d - 1) // 2


def ist_primzahlpotenz(q: int) -> bool:
    """Prüft, ob es einen Körper F_q gibt, also ob q = p^r mit Primzahl p."""
    return q >= 2 and len(factorint(q)) == 1


def kugel_summanden(n: int, e: int, q: int) -> list[Summand]:
    """Summanden von |B(e)| = Σ_{i=0}^{e} C(n, i)·(q − 1)^i."""
    return [
        Summand(i, comb(n, i), (q - 1) ** i, comb(n, i) * (q - 1) ** i)
        for i in range(e + 1)
    ]


def kugel_groesse(n: int, e: int, q: int) -> int:
    """|B(e)|: Anzahl der Wörter mit Abstand höchstens e zu einem Codewort."""
    return sum(s.wert for s in kugel_summanden(n, e, q))


def hamming(n: int, k: int, d: int, q: int) -> HammingErgebnis:
    """Prüft die Hamming-Schranke q^k · |B(e)| ≤ q^n."""
    e = e_aus_d(d)
    summanden = kugel_summanden(n, e, q)
    kugel = sum(s.wert for s in summanden)
    rest = q ** (n - k)
    return HammingErgebnis(
        n=n, k=k, d=d, q=q, e=e, d_gerade=(d % 2 == 0),
        summanden=summanden, kugel=kugel, code_groesse=q ** k, rest=rest,
        erfuellt=kugel <= rest, perfekt=kugel == rest,
    )


def singleton(n: int, k: int, d: int) -> SingletonErgebnis:
    """Prüft die Singleton-Schranke k + d ≤ n + 1."""
    return SingletonErgebnis(
        n=n, k=k, d=d, links=k + d, rechts=n + 1,
        erfuellt=k + d <= n + 1, mds=k + d == n + 1,
    )


def kleinstes_n(k: int, d: int, q: int) -> KleinstesNErgebnis:
    """Kleinstes n, für das beide Schranken einen [n, k, d]_q-Code erlauben.

    Start ist n = k + d − 1 aus der Singleton-Schranke. Dann steigt n,
    bis die Hamming-Schranke |B(e)| ≤ q^(n − k) gilt.
    """
    e = e_aus_d(d)
    n_s = k + d - 1
    erg = KleinstesNErgebnis(k=k, d=d, q=q, e=e, n_singleton=n_s)
    n = n_s
    while True:
        kugel = kugel_groesse(n, e, q)
        rest = q ** (n - k)
        ok = kugel <= rest
        erg.zeilen.append(KleinstesNZeile(n, kugel, rest, ok))
        if ok:
            break
        n += 1
    erg.n_hamming = n
    erg.n_min = max(n, n_s)
    return erg


def eingabe_fehler(n: int | None, k: int, d: int, q: int) -> list[str]:
    """Sammelt Hinweise zu ungültigen Eingaben. Leere Liste heißt: alles gut."""
    fehler = []
    if not ist_primzahlpotenz(q):
        fehler.append(
            f"q = {q} ist keine Primzahlpotenz. Einen Körper F_{q} gibt es nicht."
        )
    if k < 1:
        fehler.append(f"k = {k} ist zu klein. Die Dimension k muss mindestens 1 sein.")
    if d < 1:
        fehler.append(f"d = {d} ist zu klein. Der Mindestabstand d muss mindestens 1 sein.")
    if n is not None:
        if n < 1:
            fehler.append(f"n = {n} ist zu klein.")
        if k > n:
            fehler.append(f"k = {k} ist größer als n = {n}. Das geht nicht.")
        if d > n:
            fehler.append(f"d = {d} ist größer als n = {n}. Das geht nicht.")
    return fehler


# ---------------------------------------------------------------------------
# Ausgeben
# ---------------------------------------------------------------------------

def _potenz(x: int, q: int) -> str | None:
    """Gibt "q^j" zurück, wenn x = q^j ist. Sonst None."""
    j, y = 0, 1
    while y < x:
        y *= q
        j += 1
    return f"{q}^{j}" if y == x else None


def _summand_formel(s: Summand, n: int, q: int) -> str:
    """Summand als Formel, z. B. 23·22/(1·2) oder 11·10/(1·2)·2^2."""
    if s.i == 0:
        return "1"
    oben = "·".join(str(n - j) for j in range(s.i))
    if s.i == 1:
        teil = oben
    else:
        unten = "·".join(str(j) for j in range(1, s.i + 1))
        teil = f"{oben}/({unten})"
    if q == 2:
        return teil
    if s.i == 1:
        return f"{teil}·{q - 1}"
    return f"{teil}·{q - 1}^{s.i}"


def _summand_zahl(s: Summand, q: int) -> str:
    """Summand C(n, i)·(q − 1)^i als Zahlen, z. B. 45·100."""
    if s.i == 0:
        return "1"
    if q == 2:
        return str(s.binom)
    return f"{s.binom}·{s.faktor}"


def text_e(d: int, e: int) -> list[str]:
    if d % 2 == 1:
        return [f"e = {e}   (aus 2e + 1 = d = {d})"]
    return [
        f"Hinweis: d = {d} ist gerade. 2e + 1 = d hat keine ganze Lösung.",
        f"Annahme: e = ⌊(d − 1)/2⌋ = ⌊{d - 1}/2⌋ = {e}.",
        f"e = {e}",
    ]


def text_kugel(n: int, e: int, q: int, summanden: list[Summand], kugel: int) -> list[str]:
    zeilen = [f"Σ_{{i=0}}^{{{e}}} C({n}, i)·({q} − 1)^i"]
    formel = " + ".join(_summand_formel(s, n, q) for s in summanden)
    zeilen.append("  = " + formel)
    zahl = " + ".join(_summand_zahl(s, q) for s in summanden)
    werte = " + ".join(str(s.wert) for s in summanden)
    if zahl not in (formel, werte):
        zeilen.append("  = " + zahl)
    if len(summanden) > 1:
        zeilen.append("  = " + werte)
    ende = f"  = {kugel}"
    p = _potenz(kugel, q)
    if p and p != str(kugel) and len(summanden) > 1:
        ende += f" = {p}"
    zeilen.append(ende)
    return zeilen


def text_hamming(h: HammingErgebnis) -> str:
    n, k, q = h.n, h.k, h.q
    z = [
        "Hamming-Schranke:  |C| · Σ_{i=0}^{e} C(n, i)·(q − 1)^i ≤ |F_q^n| = q^n",
        "(|B(e)| = Σ ... ist die Anzahl der Wörter mit Abstand ≤ e zu einem Codewort.)",
        "",
    ]
    z += text_e(h.d, h.e)
    z.append(f"|C| = q^k = {q}^{k} = {h.code_groesse}")
    z += text_kugel(n, h.e, q, h.summanden, h.kugel)
    z.append("")
    z.append(f"|F_{q}^{n}| = {q}^{n} = {q}^{k}·{q}^{n - k} = {q}^{k}·{h.rest}")
    z.append(f"|C|·|B({h.e})| = {q}^{k}·{h.kugel}")
    if h.perfekt:
        z.append(f"Vergleich: {q}^{k}·{h.kugel} = {q}^{k}·{h.rest} = {q}^{n}")
        z.append("Die Hamming-Schranke gilt mit \"=\". Der Code ist perfekt.")
        z.append(f"Probe: {q}^{k}·{h.kugel} = {q ** k * h.kugel} = {q}^{n} = {q ** n}")
    elif h.erfuellt:
        diff = h.rest - h.kugel
        z.append(f"Vergleich: {q}^{k}·{h.kugel} ≤ {q}^{k}·{h.rest}, also {h.kugel} ≤ {h.rest}")
        z.append("Die Hamming-Schranke ist erfüllt, aber nicht mit \"=\".")
        z.append("Der Code ist nicht perfekt.")
        z.append(
            f"Es gibt {q}^{k}·{h.rest} − {q}^{k}·{h.kugel} = {q}^{k}·{diff} Wörter"
            f" in F_{q}^{n}, deren Abstand zu einem Codewort mindestens {h.e + 1} beträgt."
        )
    else:
        z.append(f"Vergleich: {q}^{k}·{h.kugel} > {q}^{k}·{h.rest}, also {h.kugel} > {h.rest}")
        z.append("Die Hamming-Schranke ist verletzt.")
        z.append(f"Einen linearen [{n}, {k}, {h.d}]_{q}-Code gibt es nicht.")
    return "\n".join(z)


def text_singleton(s: SingletonErgebnis, q: int) -> str:
    z = [
        "Singleton-Schranke:  k + d ≤ n + 1",
        f"k = {s.k}, d = {s.d}, n = {s.n}",
    ]
    if s.mds:
        z.append(f"{s.k} + {s.d} = {s.links} = {s.n} + 1")
        z.append("Die Singleton-Schranke gilt mit \"=\".")
        z.append(f"D. h. der Code mit den Parametern [{s.n}, {s.k}, {s.d}]_{q} ist optimal"
                 " (MDS-Code).")
    elif s.erfuellt:
        z.append(f"{s.k} + {s.d} = {s.links} < {s.rechts} = {s.n} + 1")
        z.append("Die Singleton-Schranke ist erfüllt, aber nicht mit \"=\".")
        z.append("Der Code ist nicht optimal (kein MDS-Code).")
    else:
        z.append(f"{s.k} + {s.d} = {s.links} > {s.rechts} = {s.n} + 1")
        z.append("Die Singleton-Schranke ist verletzt.")
        z.append(f"Einen linearen [{s.n}, {s.k}, {s.d}]_{q}-Code gibt es nicht.")
    return "\n".join(z)


def text_code(n: int, k: int, d: int, q: int) -> str:
    h = hamming(n, k, d, q)
    s = singleton(n, k, d)
    teile = [
        f"[{n}, {k}, {d}]_{q}",
        "Annahme: \"optimal\" heißt Singleton-Schranke mit \"=\" (MDS).",
        "",
        text_hamming(h),
        "",
        text_singleton(s, q),
        "",
    ]
    if h.erfuellt and s.erfuellt:
        teile.append(
            f"Ergebnis: [{n}, {k}, {d}]_{q} erfüllt beide Schranken; "
            f"perfekt: {'ja' if h.perfekt else 'nein'}; "
            f"optimal (MDS): {'ja' if s.mds else 'nein'}."
        )
    else:
        verletzt = [name for name, ok in (("Hamming", h.erfuellt), ("Singleton", s.erfuellt))
                    if not ok]
        teile.append(
            f"Ergebnis: [{n}, {k}, {d}]_{q} verletzt die {' und '.join(verletzt)}-Schranke. "
            "Einen solchen linearen Code gibt es nicht."
        )
    return "\n".join(teile)


def text_kleinstes_n(r: KleinstesNErgebnis) -> str:
    k, d, q, e = r.k, r.d, r.q, r.e
    z = [f"[n, {k}, {d}]_{q}: Wie groß muss n mindestens sein?", ""]
    z.append("Singleton-Schranke:  k + d ≤ n + 1")
    z.append(f"  ⇒ {k} + {d} ≤ n + 1  ⇒  {r.n_singleton} ≤ n")
    z.append("")
    z.append("Hamming-Schranke:  |C| · Σ_{i=0}^{e} C(n, i)·(q − 1)^i ≤ |F_q^n|,  2e + 1 = d")
    z += text_e(d, e)
    # Allgemeine Form der Summe mit n als Variable.
    teile = []
    for i in range(e + 1):
        if i == 0:
            teile.append("1")
            continue
        oben = "n" if i == 1 else "n·" + "·".join(f"(n−{j})" for j in range(1, i))
        t = oben if i == 1 else f"{oben}/{prod(range(1, i + 1))}"
        if q != 2:
            t += f"·{q - 1}" if i == 1 else f"·{q - 1}^{i}"
        teile.append(t)
    summe = " + ".join(teile)
    z.append(f"  ⇒ {q}^k·( {summe} ) ≤ {q}^n,   q = {q}, k = {k}")
    z.append(f"  ⇒ {summe} ≤ {q}^(n−{k})")
    z.append("")
    breite = max(len(str(zeile.n)) for zeile in r.zeilen)
    for zeile in r.zeilen:
        summanden = " + ".join(str(s.wert) for s in kugel_summanden(zeile.n, e, q))
        zeichen = "≤" if zeile.erfuellt else ">"
        status = "erfüllt" if zeile.erfuellt else "verletzt"
        z.append(
            f"n = {zeile.n:>{breite}}  ⇒  {summanden} = {zeile.kugel}"
            f"  {zeichen}  {q}^{zeile.n - k} = {zeile.rest}   ({status})"
        )
    z.append("")
    z.append(f"Die Hamming-Schranke gilt erst ab n = {r.n_hamming}.")
    if r.n_hamming > r.n_singleton:
        z.append(f"Das ist mehr als n ≥ {r.n_singleton} aus der Singleton-Schranke.")
    z.append(f"Ergebnis: n ≥ {r.n_min}")
    # Probe: n_min erfüllt beide Schranken, n_min − 1 nicht.
    h = hamming(r.n_min, k, d, q)
    z.append(
        f"Probe: n = {r.n_min}: {h.kugel} ≤ {h.rest} und {k} + {d} ≤ {r.n_min} + 1."
    )
    if r.n_min - 1 >= 1:
        m = r.n_min - 1
        hk, hr = kugel_groesse(m, e, q), q ** (m - k) if m >= k else 0
        grund = []
        if k + d > m + 1:
            grund.append(f"{k} + {d} > {m} + 1")
        if hk > hr:
            grund.append(f"{hk} > {q}^{m - k} = {hr}" if m >= k else "n < k")
        z.append(f"       n = {m}: " + " und ".join(grund) + " (verletzt).")
    z.append(
        "Hinweis: Die Schranken sind nur notwendige Bedingungen. Sie zeigen, dass"
        f" n < {r.n_min} unmöglich ist. Ob ein Code mit n = {r.n_min} wirklich"
        " existiert, sagen sie nicht."
    )
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m skripte.schranken",
        description=(
            "Prüft Hamming- und Singleton-Schranke für einen linearen [n, k, d]_q-Code. "
            "Sagt, ob der Code perfekt oder optimal (MDS) ist. "
            "Ohne --n bestimmt das Skript das kleinste mögliche n."
        ),
        epilog=(
            "Beispiele:\n"
            "  python -m skripte.schranken --n 10 --k 8 --d 3 --q 11  (MDS-Code prüfen)\n"
            "  python -m skripte.schranken --n 23 --k 12 --d 7 --q 2  (perfekter Golay-Code)\n"
            "  python -m skripte.schranken --k 5 --d 5 --q 2  (kleinstes n bestimmen)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--n", type=int, default=None,
                        help="Wortlänge n. Weglassen, um das kleinste n zu bestimmen.")
    parser.add_argument("--k", type=int, required=True, help="Dimension k des Codes.")
    parser.add_argument("--d", type=int, required=True, help="Mindestabstand d des Codes.")
    parser.add_argument("--q", type=int, required=True,
                        help="Anzahl q der Buchstaben (Körper F_q, q = Primzahlpotenz).")
    args = parser.parse_args(argv)

    fehler = eingabe_fehler(args.n, args.k, args.d, args.q)
    if fehler:
        print("Achtung, die Eingabe passt nicht:")
        for f in fehler:
            print("  - " + f)
        return 1

    if args.n is None:
        print(text_kleinstes_n(kleinstes_n(args.k, args.d, args.q)))
    else:
        print(text_code(args.n, args.k, args.d, args.q))
    return 0


if __name__ == "__main__":
    sys.exit(main())
