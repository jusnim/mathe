"""RSA-Verfahren komplett: verschlüsseln, entschlüsseln, φ(m), privater Schlüssel d.

Zweck
-----
Das Skript rechnet RSA-Aufgaben Schritt für Schritt durch.
RSA ist ein Verschlüsselungsverfahren mit zwei Schlüsseln:
- öffentlicher Schlüssel (e, m): Damit verschlüsselt Alice. Formel: c = wᵉ mod m.
- privater Schlüssel d: Damit entschlüsselt Bob. Formel: w = cᵈ mod m.
Dabei gilt e·d ≡ 1 (mod φ(m)). φ(m) (Eulersche φ-Funktion) ist die Anzahl
der Zahlen von 1 bis m, die mit m keinen gemeinsamen Teiler haben.

Die Ausgabe zeigt jeden Schritt:
 a) e binär,  b) c = wᵉ mod m mit schnellem Potenzieren,  c) Formel w = cᵈ mod m,
 d) Primfaktorzerlegung von m und φ(m),
 e) d auf zwei Wegen: Weg 1 mit EA und Matrix Q (Lemma von Bezout),
    Weg 2 mit dem Satz von Euler: d ≡ e^(φ(φ(m)) − 1) mod φ(m),
 dann die Probe cᵈ mod m = w und den Text „Was braucht Oscar?“.

Aufruf (Beispiele)
------------------
    python -m skripte.rsa --p 43 --q 57 --e 221 --w 1511      (57 ist keine Primzahl: c = 896, d = 821)
    python -m skripte.rsa --m 2803 --e 113 --a 715            (m prim: c = 708, d = 1463)
    python -m skripte.rsa --m 101 --e 13                      (nur d berechnen: d = 77)
    python -m skripte.rsa --m 2803 --schluessel 113 --w 715 --alle-varianten   (Schlüssel als e und als d)
    python -m skripte.rsa --m 2803 --d 1463 --c 708           (nur entschlüsseln)

Varianten (--variante)
----------------------
Die Variante sagt, welcher Schlüssel mit --schluessel gegeben ist.
- oeffentlich (Standard): Die Zahl ist e. Das Skript berechnet d.
- privat: Die Zahl ist d. Das Skript berechnet e.
  Mit m = 2803, d = 113, w = 715 folgt e = 1463 und c = 265.
--alle-varianten gibt beide Deutungen nacheinander aus.
Mit --e oder --d ist die Deutung fest.

Wichtige Prüfungen
------------------
- Das Skript zerlegt m immer vollständig in Primfaktoren.
  Beispiel: m = 43·57, aber 57 = 3·19 ist keine Primzahl.
  Falsch wäre φ = (43 − 1)(57 − 1) = 2352 und d = 149. Richtig ist φ = 1512 und d = 821.
  Das Skript zeigt beide Rechnungen und erklärt den Fehler.
- ggT(e, φ(m)) ≠ 1: Dann gibt es kein d. Das Skript sagt das deutlich.

Weitere Optionen
----------------
- --ea-weg matrix (Standard), rueckwaerts oder beide:
  Form von Weg 1 (Matrix Q oder „EA rückwärts einsetzen“).
- --rep standard (Reste 0 … m−1, Standard) oder symmetrisch (Reste −m/2 … m/2)
  für die Quadrattabellen. d wird am Ende immer positiv angegeben (z. B. −1339 ≡ 1463).

Nutzung als Modul
-----------------
    from skripte.rsa import rsa_berechnen, loesungsweg
    erg = rsa_berechnen(m=2451, schluessel=221, w=1511, p=43, q=57)
    erg.c, erg.phi_m, erg.d      # 896, 1512, 821
    print(loesungsweg(erg))
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from math import gcd

if __package__ in (None, ""):  # direkt gestartet: python3 skripte/<name>.py
    import os as _os
    import sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from skripte.euklid import KeinInversesFehler, euklid, inverse_mit_weg, text_ea
from skripte.primfaktor_phi import (
    PhiErgebnis,
    ist_prim,
    loesungsweg as phi_loesungsweg,
    phi_mit_weg,
    primfaktorzerlegung,
    zerlegung_text,
)
from skripte.schnell_potenzieren import (
    PotenzErgebnis,
    binaer_text,
    binaer_zerlegung,
    potenz_text,
    schnell_potenzieren,
)

VARIANTEN = ("oeffentlich", "privat")
STANDARD_VARIANTE = "oeffentlich"
EA_WEGE = ("matrix", "rueckwaerts", "beide")
REPS = ("standard", "symmetrisch")

_HOCH = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def _hoch(k: int) -> str:
    return str(k).translate(_HOCH)


# ---------------------------------------------------------------------------
# Rechnen
# ---------------------------------------------------------------------------

@dataclass
class FalscherWeg:
    """Rechnung mit φ = (p − 1)(q − 1), obwohl p oder q keine Primzahl ist."""

    phi_falsch: int
    partner_falsch: int | None      # falsches d (bzw. e); None, wenn ggT ≠ 1
    ggt_falsch: int
    klartext_falsch: int | None     # Ergebnis der Entschlüsselung mit falschem d


@dataclass
class RSAErgebnis:
    """Alle Werte einer RSA-Aufgabe. Das Feld text enthält den Lösungsweg."""

    m: int
    p: int | None
    q: int | None
    variante: str                   # "oeffentlich": Schlüssel ist e; "privat": Schlüssel ist d
    schluessel: int
    klartext_name: str              # Name des Klartexts: "w" oder "a"
    ea_weg: str
    reste: str                      # "normal" oder "symmetrisch" (für schnell_potenzieren)
    phi: PhiErgebnis
    phi_m: int
    ggt_schluessel_phi: int
    e: int | None = None
    d: int | None = None
    partner_roh: int | None = None  # Wert direkt aus Matrix Q / EA (z. B. −1339)
    inverse_wege: list = field(default_factory=list)   # InverseErgebnis je EA-Weg
    ea_fehler_text: str = ""        # EA-Text, wenn ggT ≠ 1
    phi_phi: PhiErgebnis | None = None
    phi_phi_m: int | None = None
    euler: PotenzErgebnis | None = None    # schluessel^(φ(φ(m)) − 1) mod φ(m)
    w: int | None = None
    c: int | None = None
    c_gegeben: int | None = None
    verschluesselung: PotenzErgebnis | None = None   # wᵉ mod m
    entschluesselung: PotenzErgebnis | None = None   # cᵈ mod m (Probe oder Entschlüsselung)
    falscher_weg: FalscherWeg | None = None
    hinweise: list[str] = field(default_factory=list)
    text: str = ""


def rsa_berechnen(m: int | None = None, schluessel: int | None = None, w: int | None = None,
                  c: int | None = None, p: int | None = None, q: int | None = None,
                  variante: str = STANDARD_VARIANTE, klartext_name: str = "w",
                  ea_weg: str = "matrix", rep: str = "standard") -> RSAErgebnis:
    """Rechnet eine RSA-Aufgabe durch und gibt alle Zwischenwerte zurück.

    m: Modul. Oder p und q angeben, dann ist m = p·q.
    schluessel: e (variante "oeffentlich") oder d (variante "privat"). Darf fehlen.
    w: Klartext. c: Codewort (wird entschlüsselt, wenn w fehlt).
    Wirft ValueError bei unbrauchbaren Eingaben (zum Beispiel m < 3).
    """
    if variante not in VARIANTEN:
        raise ValueError(f"Unbekannte Variante {variante!r}. Erlaubt: {', '.join(VARIANTEN)}.")
    if ea_weg not in EA_WEGE:
        raise ValueError(f"Unbekannter EA-Weg {ea_weg!r}. Erlaubt: {', '.join(EA_WEGE)}.")
    hinweise: list[str] = []
    if m is None:
        if p is None or q is None:
            raise ValueError("Bitte m angeben oder p und q.")
        m = p * q
    elif p is not None and q is not None and p * q != m:
        hinweise.append(f"ACHTUNG: p·q = {p}·{q} = {p * q} ≠ m = {m}. Das Skript rechnet mit m = {m}.")
    if m < 3:
        raise ValueError(f"Der Modul m = {m} ist zu klein. Bitte m ≥ 3 wählen.")

    faktoren = [f for f in (p, q) if f is not None] or None
    phi_erg = phi_mit_weg(m, faktoren=faktoren)
    phi_m = phi_erg.phi
    reste = "normal" if rep == "standard" else "symmetrisch"

    erg = RSAErgebnis(m=m, p=p, q=q, variante=variante, schluessel=schluessel or 0,
                      klartext_name=klartext_name, ea_weg=ea_weg, reste=reste,
                      phi=phi_erg, phi_m=phi_m, ggt_schluessel_phi=0, hinweise=hinweise)

    # Hinweise zu m
    zerl = phi_erg.zerlegung.faktoren
    for f in (p, q):
        if f is not None and not ist_prim(f):
            zf = zerlegung_text(primfaktorzerlegung(f)) if f > 1 else str(f)
            hinweise.append(f"ACHTUNG: {f} ist keine Primzahl ({f} = {zf}). "
                            "Dann gilt φ(m) ≠ (p − 1)(q − 1). Richtig: m vollständig zerlegen (Teil d).")
    if ist_prim(m):
        hinweise.append(f"m = {m} ist eine Primzahl. Dann gilt φ(m) = m − 1 = {m - 1}. "
                        "Oscar kann d sofort berechnen. RSA ist so unsicher.")
    elif any(k > 1 for k in zerl.values()):
        hinweise.append(f"m = {zerlegung_text(zerl)} enthält eine Primzahl mehrfach. "
                        "Dann klappt die Entschlüsselung nur sicher für Klartexte mit ggT(w, m) = 1.")
    elif len(zerl) != 2:
        hinweise.append(f"m = {zerlegung_text(zerl)} hat {len(zerl)} Primfaktoren, nicht 2. "
                        "Die Formel (p − 1)(q − 1) passt dann nicht. Richtig ist φ(m) aus der vollen Zerlegung.")

    # Schlüssel und Partner (d zu e oder e zu d)
    if schluessel is not None:
        if schluessel < 1:
            raise ValueError(f"Der Schlüssel muss positiv sein (gegeben: {schluessel}).")
        g = gcd(schluessel, phi_m)
        erg.ggt_schluessel_phi = g
        if variante == "oeffentlich":
            erg.e = schluessel
        else:
            erg.d = schluessel
        if g != 1:
            erg.ea_fehler_text = text_ea(euklid(phi_m, schluessel % phi_m or phi_m))
            hinweise.append(f"ACHTUNG: ggT({schluessel}, φ(m)) = ggT({schluessel}, {phi_m}) = {g} ≠ 1. "
                            f"Dann gibt es kein passendes {'d' if variante == 'oeffentlich' else 'e'}. "
                            "Der Schlüssel muss in ℤ_φ(m)* liegen.")
        else:
            wege = ["matrix", "rueckwaerts"] if ea_weg == "beide" else [ea_weg]
            for weg in wege:
                erg.inverse_wege.append(inverse_mit_weg(schluessel, phi_m, "standard", weg))
            partner = erg.inverse_wege[0].inverse
            erg.partner_roh = erg.inverse_wege[0].roh
            if variante == "oeffentlich":
                erg.d = partner
            else:
                erg.e = partner
            # Weg 2: Satz von Euler
            erg.phi_phi = phi_mit_weg(phi_m)
            erg.phi_phi_m = erg.phi_phi.phi
            erg.euler = schnell_potenzieren(schluessel, erg.phi_phi_m - 1, phi_m, reste)
            assert erg.euler.ergebnis == partner

    # Klartext und Codewort
    if w is not None:
        erg.w = w
        if not 0 <= w < m:
            hinweise.append(f"ACHTUNG: Der Klartext {klartext_name} = {w} liegt nicht in 0 … m−1 = 0 … {m - 1}. "
                            "RSA kann nur Zahlen kleiner als m verschlüsseln. "
                            f"Das Skript rechnet mit {w} mod {m} = {w % m}.")
        elif gcd(w, m) != 1:
            hinweise.append(f"Hinweis: ggT({klartext_name}, m) = ggT({w}, {m}) = {gcd(w, m)} ≠ 1. "
                            f"{klartext_name} liegt nicht in ℤ_m*. Das kommt praktisch nie vor.")
    if c is not None:
        erg.c_gegeben = c
    if w is not None and erg.e is not None:
        erg.verschluesselung = schnell_potenzieren(w, erg.e, m, reste)
        erg.c = erg.verschluesselung.ergebnis
        if c is not None and c % m != erg.c:
            hinweise.append(f"ACHTUNG: Das gegebene Codewort c = {c} passt nicht. "
                            f"Aus {klartext_name}ᵉ mod m folgt c = {erg.c}.")
    elif c is not None:
        erg.c = c % m
    if erg.c is not None and erg.d is not None:
        erg.entschluesselung = schnell_potenzieren(erg.c, erg.d, m, reste)
        if w is None:
            erg.w = erg.entschluesselung.ergebnis

    # Falscher Weg φ = (p − 1)(q − 1), wenn p oder q keine Primzahl ist
    fp = phi_erg.faktor_pruefung
    if fp is not None and fp.keine_primzahl and schluessel is not None:
        pf = fp.falsches_phi
        gf = gcd(schluessel, pf)
        partner_f = pow(schluessel, -1, pf) if gf == 1 else None
        klar_f = None
        if partner_f is not None and erg.w is not None and erg.c is not None:
            d_f = partner_f if variante == "oeffentlich" else schluessel
            klar_f = pow(erg.c, d_f, m) if variante == "oeffentlich" else pow(erg.w, partner_f, m)
        erg.falscher_weg = FalscherWeg(pf, partner_f, gf, klar_f)

    erg.text = loesungsweg(erg)
    return erg


# ---------------------------------------------------------------------------
# Ausgeben
# ---------------------------------------------------------------------------

def _titel(t: str) -> list[str]:
    return ["", t, "-" * len(t)]


def _namen(erg: RSAErgebnis) -> tuple[str, str]:
    """(Name des gegebenen Schlüssels, Name des gesuchten Schlüssels)."""
    return ("e", "d") if erg.variante == "oeffentlich" else ("d", "e")


def _teil_binaer(erg: RSAErgebnis) -> list[str]:
    z = _titel("a) e im Binärsystem")
    if erg.e is None:
        return z + ["e ist nicht bekannt."]
    b = binaer_zerlegung(erg.e)
    z.append(binaer_text(erg.e))
    z.append(f"Ergebnis a): e = {erg.e} = {b.bits}₂")
    return z


def _teil_verschluesseln(erg: RSAErgebnis) -> list[str]:
    w = erg.klartext_name
    z = _titel(f"b) Verschlüsseln: c = {w}ᵉ mod m")
    z.append(f"Formel (Alice): c ≡ {w}ᵉ (mod m)")
    if erg.verschluesselung is None:
        if erg.w is None:
            z.append(f"Kein Klartext {w} gegeben. Keine Rechnung.")
        else:
            z.append("e ist nicht bekannt. Keine Rechnung.")
        return z
    z.append(f"Hier: c ≡ {erg.w}{_hoch(erg.e)} (mod {erg.m})")
    z.append("")
    z.append(potenz_text(erg.verschluesselung, name="c"))
    return z


def _teil_entschluesseln_formel(erg: RSAErgebnis) -> list[str]:
    w = erg.klartext_name
    z = _titel(f"c) Entschlüsseln: {w} = cᵈ mod m")
    z.append(f"Formel (Bob): {w} ≡ cᵈ (mod m), denn cᵈ ≡ ({w}ᵉ)ᵈ = {w}^(e·d) ≡ {w} (mod m),")
    z.append("weil e·d ≡ 1 (mod φ(m)) gilt (Satz von Euler: a^φ(m) ≡ 1).")
    if erg.c is not None and erg.d is not None:
        z.append(f"Hier: {w} ≡ {erg.c}{_hoch(erg.d)} (mod {erg.m})   (Rechnung: siehe Probe)")
    return z


def _teil_phi(erg: RSAErgebnis) -> list[str]:
    z = _titel("d) Primfaktorzerlegung von m und φ(m)")
    z.append(phi_loesungsweg(erg.phi, ["produktformel", "multiplikativ"]))
    return z


def _teil_falsch(erg: RSAErgebnis) -> list[str]:
    f = erg.falscher_weg
    if f is None:
        return []
    fp = erg.phi.faktor_pruefung
    geg, ges = _namen(erg)
    w = erg.klartext_name
    z = _titel("ACHTUNG: Falle „φ(m) = (p − 1)(q − 1)“")
    for x, zx in fp.keine_primzahl.items():
        z.append(f"{x} ist keine Primzahl: {x} = {zerlegung_text(zx)}.")
    z.append("Die Formel φ(p·q) = (p − 1)(q − 1) gilt nur für Primzahlen p ≠ q.")
    formel = " · ".join(f"({x} − 1)" for x in fp.faktoren)
    z.append(f"Falsch: φ = {formel} = {f.phi_falsch}.")
    if f.partner_falsch is None:
        z.append(f"Mit diesem falschen φ wäre ggT({erg.schluessel}, {f.phi_falsch}) = {f.ggt_falsch} ≠ 1.")
    else:
        z.append(f"Daraus folgt das falsche {ges} = {erg.schluessel}⁻¹ mod {f.phi_falsch} = {f.partner_falsch}.")
        if f.klartext_falsch is not None:
            if erg.variante == "oeffentlich":
                z.append(f"Test: c^{ges} = {erg.c}{_hoch(f.partner_falsch)} mod {erg.m} = {f.klartext_falsch}"
                         f" ≠ {w} = {erg.w}. Die Entschlüsselung geht schief.")
            else:
                z.append(f"Test: {w}^{ges} = {erg.w}{_hoch(f.partner_falsch)} mod {erg.m} = {f.klartext_falsch}"
                         f" ≠ c = {erg.c}. Das passt nicht.")
    z.append(f"Richtig: φ({erg.m}) = {erg.phi_m} aus der vollen Zerlegung "
             f"{erg.m} = {zerlegung_text(erg.phi.zerlegung.faktoren)}. "
             f"Damit ist {ges} = {erg.d if ges == 'd' else erg.e}.")
    return z


def _teil_schluessel(erg: RSAErgebnis) -> list[str]:
    geg, ges = _namen(erg)
    k, pm = erg.schluessel, erg.phi_m
    if erg.variante == "oeffentlich":
        z = _titel("e) Oscar berechnet den privaten Schlüssel d")
        z.append("Oscar kennt m und e. Er rechnet: m zerlegen → φ(m) → d.")
    else:
        z = _titel("Öffentlichen Schlüssel e aus d berechnen")
        z.append(f"Annahme: {k} ist der private Schlüssel d. Aus e·d ≡ 1 (mod φ(m)) folgt e = d⁻¹ mod φ(m).")
    z.append(f"Bedingung: e·d ≡ 1 (mod φ(m)), also {geg}·{ges} ≡ 1 (mod {pm}).")
    z.append(f"{ges} existiert nur, wenn ggT({geg}, φ(m)) = 1.")
    if erg.ggt_schluessel_phi != 1:
        z.append("")
        z.append(erg.ea_fehler_text)
        z.append(f"ggT({k}, {pm}) = {erg.ggt_schluessel_phi} ≠ 1 ⇒ Es gibt kein {ges}. "
                 f"Wählen Sie ein {geg} ohne gemeinsamen Teiler mit {pm}.")
        return z

    # Weg 1
    z.append("")
    z.append("1. Lösungsweg: Lemma von Bezout, EA")
    z.append(f"  {geg}·{ges} − φ(m)·y = 1  ⇔  {geg}·{ges} ≡ 1 (mod φ(m))")
    z.append(f"  {k}·{ges} ≡ 1 (mod {pm}).  EA mit {pm}, {k % pm}")
    for inv in erg.inverse_wege:
        z.append("")
        z.append(f"  [Form: {'Matrix Q' if inv.weg == 'matrix' else 'EA rückwärts einsetzen'}]")
        z.append(inv.text)
        if inv.roh != inv.inverse:
            z.append(f"⇒ {ges} = {inv.roh} ≡ {inv.inverse} (mod {pm})")
        else:
            z.append(f"⇒ {ges} = {inv.inverse}")

    # Weg 2
    pp = erg.phi_phi
    ex = erg.phi_phi_m - 1
    z.append("")
    z.append("2. Lösungsweg: Satz von Euler, schnelles Potenzieren")
    z.append(f"  Euler: a^φ(n) ≡ 1 (mod n) für ggT(a, n) = 1. Mit a = {geg} und n = φ(m):")
    z.append(f"  {geg}^φ(φ(m)) ≡ 1 ⇒ {geg}·{geg}^(φ(φ(m)) − 1) ≡ 1 ⇒ {ges} ≡ {geg}^(φ(φ(m)) − 1) (mod φ(m))")
    z.append(f"  Gesucht also φ(φ(m)) = φ({pm}).")
    z.append("")
    z.append(phi_loesungsweg(pp, ["multiplikativ"]))
    z.append(f"⇒ φ(φ(m)) = φ({pm}) = {erg.phi_phi_m}")
    z.append(f"⇒ φ(φ(m)) − 1 = {ex}")
    z.append(f"⇒ {ges} ≡ {k}{_hoch(ex)} (mod {pm})")
    z.append("")
    z.append(potenz_text(erg.euler, name=ges))
    z.append("")
    wert = erg.d if ges == "d" else erg.e
    z.append(f"Beide Wege geben {ges} = {wert}. ✓")
    z.append(f"Kontrolle: e·d = {erg.e}·{erg.d} = {erg.e * erg.d} = "
             f"{erg.e * erg.d // pm}·{pm} + {erg.e * erg.d % pm} ≡ 1 (mod {pm}) ✓")
    z.append(f"Ergebnis: {ges} = {wert}")
    return z


def _teil_probe(erg: RSAErgebnis) -> list[str]:
    if erg.entschluesselung is None:
        return []
    w = erg.klartext_name
    if erg.verschluesselung is not None:
        z = _titel(f"Probe: cᵈ mod m = {w}")
    else:
        z = _titel(f"Entschlüsseln: {w} = cᵈ mod m")
    z.append(f"{w} ≡ cᵈ = {erg.c}{_hoch(erg.d)} (mod {erg.m})")
    z.append("")
    z.append(potenz_text(erg.entschluesselung, name=w))
    if erg.verschluesselung is not None:
        ok = erg.entschluesselung.ergebnis == erg.w % erg.m
        z.append(f"Probe: {erg.entschluesselung.ergebnis} = {w} = {erg.w} "
                 + ("✓ Die Entschlüsselung liefert den Klartext." if ok else
                    "✗ Der Klartext kommt NICHT heraus (siehe Hinweise oben)."))
    return z


def _teil_oscar(erg: RSAErgebnis) -> list[str]:
    z = _titel("Was braucht Oscar?")
    z.append("Öffentlich sind e und m. Privat sind d und φ(m).")
    z.append("Oscar muss m in Primfaktoren zerlegen. Dann kennt er φ(m).")
    z.append("Dann rechnet er d aus e·d ≡ 1 (mod φ(m)): mit dem EA (Weg 1) oder mit")
    z.append("d ≡ e^(φ(φ(m)) − 1) mod φ(m) (Weg 2). Danach entschlüsselt er: w = cᵈ mod m.")
    z.append("Fall m Primzahl: φ(m) = m − 1. Oscar braucht keine Zerlegung.")
    z.append("Fall m = p·q mit großen Primzahlen p, q: Die Zerlegung ist praktisch unmöglich.")
    z.append("Nur dann ist RSA sicher.")
    return z


def loesungsweg(erg: RSAErgebnis) -> str:
    """Erzeugt den ganzen Lösungsweg als Text aus einem RSAErgebnis."""
    w = erg.klartext_name
    z = [f"RSA mit m = {erg.m}"
         + (f" = {erg.p}·{erg.q}" if erg.p is not None and erg.q is not None else "")]
    if erg.schluessel:
        geg, _ = _namen(erg)
        z.append(f"Variante „{erg.variante}“: Annahme, die gegebene Zahl {erg.schluessel} ist "
                 + ("der öffentliche Schlüssel e." if geg == "e" else "der private Schlüssel d."))
    if erg.w is not None and erg.verschluesselung is not None:
        z.append(f"Klartext {w} = {erg.w}")
    if erg.c_gegeben is not None:
        z.append(f"Codewort c = {erg.c_gegeben}")
    z.append(f"Reste in den Tabellen: {'0 … m−1' if erg.reste == 'normal' else 'symmetrisch (−m/2 … m/2)'}")
    if erg.hinweise:
        z.append("")
        z += erg.hinweise

    if erg.variante == "oeffentlich":
        teile = [_teil_binaer, _teil_verschluesseln, _teil_entschluesseln_formel, _teil_phi,
                 _teil_falsch, _teil_schluessel, _teil_probe, _teil_oscar]
    else:
        # Erst φ(m) und e bestimmen, dann verschlüsseln.
        teile = [_teil_phi, _teil_falsch, _teil_schluessel, _teil_binaer, _teil_verschluesseln,
                 _teil_entschluesseln_formel, _teil_probe, _teil_oscar]
    for t in teile:
        z += t(erg)

    teile_erg = []
    if erg.e is not None:
        teile_erg.append(f"e = {erg.e} = {binaer_zerlegung(erg.e).bits}₂")
    if erg.verschluesselung is not None:
        teile_erg.append(f"c = {erg.c}")
    elif erg.entschluesselung is not None:
        teile_erg.append(f"{w} = {erg.w}")
    teile_erg.append(f"φ(m) = {erg.phi_m}")
    if erg.phi_phi_m is not None:
        teile_erg.append(f"φ(φ(m)) = {erg.phi_phi_m}")
    if erg.d is not None:
        teile_erg.append(f"d = {erg.d}")
    elif erg.ggt_schluessel_phi not in (0, 1):
        teile_erg.append("kein passender Schlüssel (ggT ≠ 1)")
    z += ["", "=" * 60, "Ergebnis: " + ", ".join(teile_erg)]
    return "\n".join(z)


# ---------------------------------------------------------------------------
# Kommandozeile
# ---------------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m skripte.rsa",
        description="RSA: e binär, c = wᵉ mod m, φ(m), d auf zwei Wegen, Probe.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Beispiele:\n"
            "  python -m skripte.rsa --p 43 --q 57 --e 221 --w 1511   (57 ist keine Primzahl)\n"
            "  python -m skripte.rsa --m 2803 --e 113 --a 715   (verschlüsseln und d berechnen)\n"
            "  python -m skripte.rsa --m 101 --e 13   (nur d berechnen)\n"
            "  python -m skripte.rsa --m 2803 --schluessel 113 --w 715 --alle-varianten   (Schlüssel als e und als d)\n"
            "  python -m skripte.rsa --m 2803 --d 1463 --c 708   (nur entschlüsseln)\n"
        ),
    )
    ap.add_argument("--m", type=int, help="Modul m. Alternativ --p und --q angeben.")
    ap.add_argument("--p", type=int, help="Erster Faktor von m (wird auf Primzahl geprüft).")
    ap.add_argument("--q", type=int, help="Zweiter Faktor von m (wird auf Primzahl geprüft).")
    sg = ap.add_mutually_exclusive_group()
    sg.add_argument("--e", type=int, help="Öffentlicher Schlüssel e. Das Skript berechnet d.")
    sg.add_argument("--d", type=int, help="Privater Schlüssel d. Das Skript berechnet e.")
    sg.add_argument("--schluessel", type=int,
                    help="Schlüssel, dessen Rolle unklar ist. Rolle mit --variante wählen.")
    kg = ap.add_mutually_exclusive_group()
    kg.add_argument("--w", type=int, help="Klartext w.")
    kg.add_argument("--a", type=int, help="Klartext a (gleiche Rolle wie --w, nur anderer Name).")
    ap.add_argument("--c", type=int, help="Codewort c. Ohne Klartext wird c entschlüsselt.")
    ap.add_argument("--variante", choices=VARIANTEN, default=None,
                    help="Rolle von --schluessel: oeffentlich (= e, Standard) oder privat (= d).")
    ap.add_argument("--alle-varianten", action="store_true",
                    help="Gibt den gegebenen Schlüssel einmal als e und einmal als d aus.")
    ap.add_argument("--ea-weg", choices=EA_WEGE, default="matrix",
                    help="Form von Weg 1: matrix (Matrix Q, Standard), rueckwaerts (EA rückwärts) oder beide.")
    ap.add_argument("--rep", choices=REPS, default="standard",
                    help="Reste in den Tabellen: standard (0 … m−1) oder symmetrisch (−m/2 … m/2).")
    return ap


def main(argv: list[str] | None = None) -> int:
    ap = _parser()
    args = ap.parse_args(argv)
    if args.m is None and (args.p is None or args.q is None):
        print("HINWEIS: Bitte --m angeben oder --p und --q. Beispiel: --p 43 --q 57 --e 221 --w 1511")
        return 1
    if args.e is not None:
        schluessel, feste = args.e, "oeffentlich"
    elif args.d is not None:
        schluessel, feste = args.d, "privat"
    else:
        schluessel, feste = args.schluessel, None
    if feste is not None and args.variante is not None and args.variante != feste:
        print(f"HINWEIS: --{'e' if feste == 'oeffentlich' else 'd'} passt nicht zu --variante {args.variante}. "
              "Nutzen Sie --schluessel mit --variante.")
        return 1
    if args.alle_varianten and schluessel is not None:
        varianten = list(VARIANTEN)
    else:
        varianten = [feste or args.variante or STANDARD_VARIANTE]
    klartext = args.w if args.w is not None else args.a
    name = "a" if args.a is not None else "w"

    ausgaben = []
    for v in varianten:
        try:
            erg = rsa_berechnen(m=args.m, schluessel=schluessel, w=klartext, c=args.c, p=args.p, q=args.q,
                                variante=v, klartext_name=name, ea_weg=args.ea_weg, rep=args.rep)
        except (ValueError, KeinInversesFehler) as fehler:
            print(f"HINWEIS: {fehler}")
            return 1
        if len(varianten) > 1:
            ausgaben.append("#" * 60)
            ausgaben.append(f"# Variante: {v}")
            ausgaben.append("#" * 60)
        ausgaben.append(erg.text)
        ausgaben.append("")
    print("\n".join(ausgaben).rstrip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
