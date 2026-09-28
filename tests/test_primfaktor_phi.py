"""Tests für skripte/primfaktor_phi.py mit den Aufgaben aus PRIORISIERUNG.md."""

import pytest

from skripte.primfaktor_phi import (
    VARIANTEN,
    ist_prim,
    loesungsweg,
    main,
    phi,
    phi_mit_weg,
    primfaktorzerlegung,
    pruefe_faktoren,
)

# (Quelle, n, Zerlegung, φ(n))
FAELLE = [
    ("KV A4", 360, {2: 3, 3: 2, 5: 1}, 96),
    ("Ü3 A2", 1000, {2: 3, 5: 3}, 400),
    # Quelle schreibt "φ(1000)" statt "φ(1001)" (Schreibfehler, Abschnitt 5.1); Wert 720 stimmt.
    ("Ü3 A2", 1001, {7: 1, 11: 1, 13: 1}, 720),
    ("Ü3 A2", 1002, {2: 1, 3: 1, 167: 1}, 332),
    ("VL2 Z. 256", 60, {2: 2, 3: 1, 5: 1}, 16),
    ("Ü4 A9", 2802, {2: 1, 3: 1, 467: 1}, 932),
    ("Ü4 A9", 2803, {2803: 1}, 2802),
    ("K A1 d", 2451, {3: 1, 19: 1, 43: 1}, 1512),
    ("VL1 Z. 658", 6, {2: 1, 3: 1}, 2),
]


@pytest.mark.parametrize("quelle,n,zerlegung,wert", FAELLE)
def test_zerlegung_und_phi(quelle, n, zerlegung, wert):
    assert primfaktorzerlegung(n) == zerlegung
    assert phi(n) == wert


@pytest.mark.parametrize("variante", VARIANTEN)
@pytest.mark.parametrize("quelle,n,zerlegung,wert", FAELLE)
def test_jede_variante_zeigt_ergebnis(variante, quelle, n, zerlegung, wert):
    text = loesungsweg(phi_mit_weg(n), [variante])
    assert f"Ergebnis: {n} = " in text
    assert f"φ({n}) = {wert}" in text
    assert f"       = {wert}" in text


def test_kv_a4_rechenweg_wie_musterloesung():
    text = loesungsweg(phi_mit_weg(360), list(VARIANTEN))
    assert "360 = 2³ · 3² · 5" in text
    assert "360 · 1/2 · 2/3 · 4/5" in text
    assert "360 · 4/15" in text
    # KV A4 schreibt "= 12 · 8"; das Skript kürzt 360/15 = 24 und schreibt "24 · 4". Beides = 96.
    assert "= 24 · 4" in text
    assert "φ(2³) · φ(3²) · φ(5)" in text
    assert "4 · 6 · 4" in text


def test_ue3_multiplikativ_wie_musterloesung():
    text = loesungsweg(phi_mit_weg(1000), ["multiplikativ"])
    assert "(2³ − 2²) · (5³ − 5²)" in text
    assert "4 · 100" in text


def test_ue4_a9_phi_phi():
    e = phi_mit_weg(2803, phi_phi=True)
    assert e.phi == 2802
    assert e.phi_phi.phi == 932
    assert "φ(φ(2803)) = φ(2802) = 932" in loesungsweg(e)


def test_k_a1_warnung_57_keine_primzahl():
    # Quelle K A1: m = 43 · 57, aber 57 = 3 · 19. (43 − 1)(57 − 1) = 2352 wäre falsch.
    f = pruefe_faktoren(2451, [43, 57])
    assert f.produkt_stimmt
    assert f.keine_primzahl == {57: {3: 1, 19: 1}}
    assert f.falsches_phi == 2352
    text = loesungsweg(phi_mit_weg(2451, [43, 57]))
    assert "57 ist keine Primzahl" in text
    assert "Falsch wäre: (43 − 1) · (57 − 1) = 2352" in text
    assert "φ(2451) = 1512" in text


def test_probedivision_nur_primteiler():
    from skripte.primfaktor_phi import probedivision
    teiler = {s.teiler for s in probedivision(2451).schritte}
    assert teiler == {2, 3, 5, 7, 11, 13, 17, 19}


def test_faktoren_passen_nicht():
    assert not pruefe_faktoren(100, [3, 7]).produkt_stimmt


@pytest.mark.parametrize("n,prim", [(1, False), (2, True), (3, True), (4, False),
                                    (57, False), (167, True), (467, True), (2803, True)])
def test_ist_prim(n, prim):
    assert ist_prim(n) is prim


def test_phi_1_und_fehler():
    assert phi(1) == 1
    with pytest.raises(ValueError):
        phi(0)


def test_phi_gegen_direktes_zaehlen():
    from math import gcd
    for n in range(1, 300):
        assert phi(n) == sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


def test_cli(capsys):
    assert main(["--n", "360", "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert "Ergebnis: 360 = 2³ · 3² · 5, φ(360) = 96" in out
    assert main(["--n", "0"]) == 1
    assert "Hinweis" in capsys.readouterr().out
