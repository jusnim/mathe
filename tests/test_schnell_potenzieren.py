"""Tests für skripte/schnell_potenzieren.py mit Aufgaben aus den Quellen."""

import pytest

from skripte.schnell_potenzieren import (
    binaer_zerlegung,
    main,
    ordnung,
    ordnung_text,
    phi,
    potenz_mit_reduktion,
    potenz_text,
    primitive_elemente,
    quadrat_tabelle,
    reduktion_text,
    schnell_potenzieren,
)

RESTE = ["normal", "symmetrisch"]

# (a, k, m, Ergebnis, Quelle)
POTENZEN = [
    (12, 100, 34, 30, "VL3 Beispiel 12^100"),
    (11, 16, 34, 1, "VL3 Beispiel 11^16"),
    (8, 13, 101, 18, "VL3 RSA-Beispiel"),
    (3, 1000, 7, 4, "Ü4 A7"),
    (715, 113, 2803, 708, "Ü4 A9 a"),
    (113, 931, 2802, 1463, "Ü4 A9 c, 2. Weg"),
    (708, 1463, 2803, 715, "Ü4 A9 Probe"),
    (13, 39, 100, 77, "KV A7, 1. Weg"),
    (1511, 221, 2451, 896, "K A1 b"),
    (715, 1463, 2803, 265, "GP A1, falls d = 113 privat wäre (PRIORISIERUNG Abschn. 4)"),
]


@pytest.mark.parametrize("reste", RESTE)
@pytest.mark.parametrize("a,k,m,erwartet,quelle", POTENZEN)
def test_variante_binaer(a, k, m, erwartet, quelle, reste):
    erg = schnell_potenzieren(a, k, m, reste)
    assert erg.ergebnis == erwartet, quelle
    assert f"Ergebnis: {a}" in potenz_text(erg)


@pytest.mark.parametrize("reste", RESTE)
@pytest.mark.parametrize("a,k,m,erwartet,quelle", POTENZEN)
def test_variante_euler(a, k, m, erwartet, quelle, reste):
    erg = potenz_mit_reduktion(a, k, m, reste)
    assert erg.ergebnis == erwartet, quelle
    reduktion_text(erg)


def test_ue4_a7_reduktion():
    erg = potenz_mit_reduktion(3, 1000, 7)
    assert (erg.phi_m, erg.q, erg.r) == (6, 166, 4)
    assert "1000 = 166 · 6 + 4" in reduktion_text(erg)


def test_vl3_warnung_ggt():
    # VL3: ggT(12, 34) = 2, Euler gilt nicht. 12^16 ≡ −16 ≠ 1.
    erg = potenz_mit_reduktion(12, 100, 34)
    assert not erg.moeglich
    assert "HINWEIS" in reduktion_text(erg)
    assert schnell_potenzieren(12, 16, 34, "symmetrisch").tabelle[4][2] == -16


def test_binaer():
    assert binaer_zerlegung(113).bits == "1110001"         # Ü4 A9
    assert binaer_zerlegung(221).bits == "11011101"        # K A1 a
    assert binaer_zerlegung(100).exponenten == [2, 5, 6]   # VL3
    assert binaer_zerlegung(39).exponenten == [0, 1, 2, 5]  # KV A7: 1+2+4+32
    assert binaer_zerlegung(931).exponenten == [0, 1, 5, 7, 8, 9]  # Ü4 A9
    assert binaer_zerlegung(13).exponenten == [0, 2, 3]    # VL3: 1 + 2^2 + 2^3


def test_tabelle_vl3_symmetrisch():
    werte = [w for _, _, w in quadrat_tabelle(12, 34, 7, "symmetrisch")]
    assert werte == [12, 8, -4, 16, -16, -16, -16]


def test_tabelle_vl3_11():
    werte = [w for _, _, w in quadrat_tabelle(11, 34, 5)]
    assert werte == [11, 19, 21, 33, 1]


def test_tabelle_ue4_a9_715():
    werte = [w for _, _, w in quadrat_tabelle(715, 2803, 7)]
    # Ü4 A9 schreibt 715^16 ≡ 163. Richtig ist 1653 (2557² = 2332·2803 + 1653).
    # Mit 163 käme nicht 708 heraus; mit 1653 schon. Es fehlt eine Ziffer in der Quelle.
    assert werte[4:] == [1653, 2287, 2774]
    erg = schnell_potenzieren(715, 113, 2803)
    assert [w for _, w in erg.faktoren] == [715, 1653, 2287, 2774]


def test_tabelle_ue4_a9_113():
    werte = [w for _, _, w in quadrat_tabelle(113, 2802, 10)]
    assert werte == [113, 1561, 1783, 1621, 2167, 2539, 1921, 7, 49, 2401]


def test_tabelle_ue4_a9_probe_708():
    werte = [w for _, _, w in quadrat_tabelle(708, 2803, 11)]
    assert werte == [708, 2330, 2292, 442, 1957, 951, 1835, 822, 161, 694, 2323]


def test_tabelle_kv_a7():
    werte = [w for _, _, w in quadrat_tabelle(13, 100, 6)]
    assert werte[:3] == [13, 69, 61] and werte[5] == 81


def test_tabelle_vl3_8():
    werte = [w for _, _, w in quadrat_tabelle(8, 101, 4)]
    assert werte == [8, 64, 56, 5]


@pytest.mark.parametrize("reste", RESTE)
def test_ue4_a8_ordnungen(reste):
    erwartet = {2: (8, False, 9), 3: (16, True, 6), 4: (4, False, 13), 5: (16, True, 7)}
    for a, (o, prim, inv) in erwartet.items():
        erg = ordnung(a, 17, reste)
        assert (erg.ordnung, erg.primitiv, erg.inverse) == (o, prim, inv)
        ordnung_text(erg)


def test_ue4_a8_potenzliste_symmetrisch():
    assert ordnung(2, 17).potenzen == [2, 4, 8, -1, -2, -4, -8, 1]
    assert ordnung(3, 17).potenzen[:8] == [3, -8, -7, -4, 5, -2, -6, -1]
    assert ordnung(4, 17).potenzen == [4, -1, -4, 1]
    assert ordnung(5, 17).potenzen[:8] == [5, 8, 6, -4, -3, 2, -7, -1]


def test_vl3_z7():
    # VL3: ord(2) = 3, ord(3) = 6, 3 ist primitiv in Z_7*.
    assert ordnung(2, 7).ordnung == 3
    assert ordnung(3, 7).potenzen == [3, 2, -1, -3, -2, 1]
    assert 3 in primitive_elemente(7)


def test_primitive_elemente():
    assert primitive_elemente(17) == [3, 5, 6, 7, 10, 11, 12, 14]
    assert primitive_elemente(8) == []


def test_keine_einheit():
    erg = ordnung(12, 34)
    assert erg.ordnung is None
    assert "HINWEIS" in ordnung_text(erg)


def test_phi():
    assert phi(34) == 16 and phi(2802) == 932 and phi(100) == 40 and phi(2451) == 1512


def test_randfaelle():
    assert schnell_potenzieren(5, 0, 7).ergebnis == 1
    assert schnell_potenzieren(5, 1, 7).ergebnis == 5
    assert potenz_mit_reduktion(3, 12, 7).ergebnis == 1  # r = 0


def test_cli(capsys):
    assert main(["--a", "12", "--k", "100", "--m", "34", "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert out.count("Ergebnis: 12¹⁰⁰ ≡ 30 (mod 34)") == 4
    assert main(["--ordnung", "--a", "2", "3", "4", "5", "--m", "17"]) == 0
    assert "ord(3) = 16" in capsys.readouterr().out
    assert main(["--a", "3", "--k", "-1", "--m", "7"]) == 1
