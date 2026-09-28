"""Tests für skripte/crt.py mit festen Beispielen und Zufallssystemen."""

import itertools
import random

import pytest

from skripte.crt import kleinste_positive, loese_paar, loese_system, main

ALLE = [(v, w) for v in ("erste", "zweite") for w in ("rueckwaerts", "matrix")]


def _pruefe(erg, kongruenzen, modul, kp):
    assert erg.loesbar
    assert erg.modul == modul
    assert erg.kleinste_positive == kp
    for a, m in kongruenzen:
        assert (erg.z0 - a) % m == 0
    assert "Ergebnis:" in erg.text


# --- z ≡ 8 mod 17, z ≡ 5 mod 15 → 110 mod 255 -----------------------

@pytest.mark.parametrize("variante,weg", ALLE)
def test_crt_17_15_alle_varianten(variante, weg):
    k = [(8, 17), (5, 15)]
    _pruefe(loese_system(k, variante, weg), k, 255, 110)


def test_crt_17_15_standard_z0():
    # EA rückwärts: 1 = 8·15 − 7·17, also 17·(−7) − 15·(−8) = 1; ·(−3) ⇒ x = 21, z = 365.
    erg = loese_system([(8, 17), (5, 15)])
    assert erg.z0 == 365
    assert (erg.schritte[0].x, erg.schritte[0].y) == (21, 24)


def test_crt_17_15_mit_vorgabe():
    # Bekannte Bezout-Lösung: 17·8 − 15·9 = 1.
    erg = loese_system([(8, 17), (5, 15)], vorgabe=(8, 9))
    assert erg.schritte[0].weg == "vorgabe"
    assert erg.z0 == -400
    assert erg.kleinste_positive == 110
    # Andere Bezout-Lösung: 17·(−7) − 15·(−8) = 1 ⇒ z = 365.
    assert loese_system([(8, 17), (5, 15)], vorgabe=(-7, -8)).z0 == 365


def test_falsche_vorgabe_faellt_auf_ea_zurueck():
    erg = loese_system([(8, 17), (5, 15)], vorgabe=(1, 1))
    assert erg.schritte[0].weg == "rueckwaerts"
    assert "Hinweis" in erg.text
    assert erg.kleinste_positive == 110


# --- z ≡ 5 mod 13, z ≡ 4 mod 15 → −86 ≡ 109 mod 195 -----------------

@pytest.mark.parametrize("variante,weg", ALLE)
def test_crt_13_15_alle_varianten(variante, weg):
    k = [(5, 13), (4, 15)]
    erg = loese_system(k, variante, weg)
    _pruefe(erg, k, 195, 109)
    assert erg.z0 == -86


# --- z ≡ 3 mod 5, z ≡ 4 mod 7 → 18 mod 35 ------------------------

@pytest.mark.parametrize("variante,weg", ALLE)
def test_crt_5_7(variante, weg):
    k = [(3, 5), (4, 7)]
    erg = loese_system(k, variante, weg)
    _pruefe(erg, k, 35, 18)


def test_crt_5_7_zwischenwerte():
    # Variante 1: 1 = 3·5 − 2·7 ⇒ z = 3·5 + 3 = 2·7 + 4 = 18.
    p = loese_system([(3, 5), (4, 7)]).schritte[0]
    assert (p.x, p.y, p.z0) == (3, 2, 18)
    # Variante 2: Q = (7 3; 5 2), Det Q = −1 ⇒ z = 7·2 + 4 = 5·3 + 3 = 18.
    p = loese_system([(3, 5), (4, 7)], weg="matrix").schritte[0]
    assert p.matrix.Q == ((7, 3), (5, 2)) and p.matrix.det == -1
    assert p.z0 == 18


# --- z ≡ −1 mod 12, z ≡ 2 mod 5 → −73 ≡ 47 mod 60 -----------------------

@pytest.mark.parametrize("variante,weg", ALLE)
def test_crt_12_5_negativer_rest(variante, weg):
    k = [(-1, 12), (2, 5)]
    _pruefe(loese_system(k, variante, weg), k, 60, 47)


def test_crt_12_5_z0():
    p = loese_system([(-1, 12), (2, 5)]).schritte[0]
    assert p.d == 3
    assert p.z0 == -73


# --- Drei Kongruenzen → −74 ≡ 1246 mod 1320 ------------------------------

@pytest.mark.parametrize("variante,weg", ALLE)
def test_drei_kongruenzen(variante, weg):
    k = [(3, 11), (6, 8), (1, 15)]
    _pruefe(loese_system(k, variante, weg), k, 1320, 1246)


def test_drei_kongruenzen_zwischenschritte():
    erg = loese_system([(3, 11), (6, 8), (1, 15)])
    s1, s2 = erg.schritte
    assert (s1.x, s1.y, s1.z0, s1.modul) == (1, 1, 14, 88)     # 11 = 1·8 + 3
    assert (s2.a, s2.d) == (14, -13)
    assert (s2.x, s2.y, s2.z0) == (-1, -5, -74)                 # −88 + 14 = −5·15 + 1
    assert erg.z0 == -74 and erg.kleinste_positive == 1246


# --- Nicht teilerfremde Moduln --------------------------------------------

def test_nicht_teilerfremd_loesbar():
    erg = loese_system([(2, 4), (4, 6)])
    assert erg.schritte[0].ggt == 2
    _pruefe(erg, [(2, 4), (4, 6)], 12, 10)
    assert "nicht teilerfremd" in erg.text


def test_nicht_teilerfremd_unloesbar():
    erg = loese_system([(1, 4), (2, 6)])
    assert not erg.loesbar
    assert erg.z0 is None
    assert "keine Lösung" in erg.text


def test_unloesbar_erst_im_zweiten_schritt():
    erg = loese_system([(1, 3), (2, 5), (0, 6)])
    assert not erg.loesbar
    assert len(erg.schritte) == 2


@pytest.mark.parametrize("variante,weg", ALLE)
def test_vergleich_mit_ausprobieren(variante, weg):
    """Zufällige Systeme: Ergebnis stimmt mit Durchprobieren überein."""
    rnd = random.Random(4711)
    for _ in range(150):
        anz = rnd.randint(2, 4)
        k = [(rnd.randint(-20, 20), rnd.randint(1, 18)) for _ in range(anz)]
        M = 1
        for _, m in k:
            M = M * m // __import__("math").gcd(M, m)
        treffer = [z for z in range(1, M + 1) if all((z - a) % m == 0 for a, m in k)]
        erg = loese_system(k, variante, weg)
        if treffer:
            assert erg.loesbar and erg.modul == M and erg.kleinste_positive == treffer[0]
            assert len(treffer) == 1
        else:
            assert not erg.loesbar


# --- Kleinigkeiten ---------------------------------------------------------

def test_kleinste_positive():
    assert kleinste_positive(-86, 195) == 109
    assert kleinste_positive(0, 7) == 7
    assert kleinste_positive(365, 255) == 110


def test_symmetrischer_repraesentant():
    erg = loese_system([(8, 17), (5, 15)], rep="symmetrisch")
    assert erg.repraesentant == 110
    assert loese_system([(5, 13), (4, 15)], rep="symmetrisch").repraesentant == -86


def test_gleiche_reste():
    p = loese_paar(3, 5, 3, 7)
    assert (p.x, p.y, p.z0) == (0, 0, 3)


@pytest.mark.parametrize("k", [[(1, 5)], [(1, 0), (2, 3)], [(1, -4), (2, 3)]])
def test_eingabefehler(k):
    with pytest.raises(ValueError):
        loese_system(k)


def test_main(capsys):
    assert main(["-k", "5", "13", "-k", "4", "15"]) == 0
    out = capsys.readouterr().out
    assert "Ergebnis: z ≡ -86 ≡ 109 (mod 195)" in out
    assert main(["-k", "8", "17", "-k", "5", "15", "--vorgabe", "8", "9", "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert "-400" in out and "365" in out and "Vergleich" in out
    assert main(["-k", "1", "4", "-k", "2", "6"]) == 1
    assert main(["-k", "1", "0", "-k", "2", "6"]) == 2
    assert "Hinweis" in capsys.readouterr().out
