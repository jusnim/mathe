"""Tests für skripte/schnell_potenzieren.py mit Rechenbeispielen."""

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

# (a, k, m, Ergebnis, Beschreibung)
POTENZEN = [
    (12, 100, 34, 30, "12^100 mod 34"),
    (11, 16, 34, 1, "11^16 mod 34"),
    (8, 13, 101, 18, "kleines RSA-Beispiel"),
    (3, 1000, 7, 4, "3^1000 mod 7"),
    (715, 113, 2803, 708, "RSA verschlüsseln, n = 2803"),
    (113, 931, 2802, 1463, "Exponent mod phi, n = 2802"),
    (708, 1463, 2803, 715, "RSA entschlüsseln, n = 2803"),
    (13, 39, 100, 77, "13^39 mod 100"),
    (1511, 221, 2451, 896, "RSA, n = 2451"),
    (715, 1463, 2803, 265, "RSA mit vertauschten Schlüsseln, n = 2803"),
]


@pytest.mark.parametrize("reste", RESTE)
@pytest.mark.parametrize("a,k,m,erwartet,beschreibung", POTENZEN)
def test_variante_binaer(a, k, m, erwartet, beschreibung, reste):
    erg = schnell_potenzieren(a, k, m, reste)
    assert erg.ergebnis == erwartet, beschreibung
    assert f"Ergebnis: {a}" in potenz_text(erg)


@pytest.mark.parametrize("reste", RESTE)
@pytest.mark.parametrize("a,k,m,erwartet,beschreibung", POTENZEN)
def test_variante_euler(a, k, m, erwartet, beschreibung, reste):
    erg = potenz_mit_reduktion(a, k, m, reste)
    assert erg.ergebnis == erwartet, beschreibung
    reduktion_text(erg)


def test_reduktion_3_hoch_1000():
    erg = potenz_mit_reduktion(3, 1000, 7)
    assert (erg.phi_m, erg.q, erg.r) == (6, 166, 4)
    assert "1000 = 166 · 6 + 4" in reduktion_text(erg)


def test_warnung_ggt_ungleich_1():
    # ggT(12, 34) = 2, Euler gilt nicht. 12^16 ≡ −16 ≠ 1.
    erg = potenz_mit_reduktion(12, 100, 34)
    assert not erg.moeglich
    assert "HINWEIS" in reduktion_text(erg)
    assert schnell_potenzieren(12, 16, 34, "symmetrisch").tabelle[4][2] == -16


def test_binaer():
    assert binaer_zerlegung(113).bits == "1110001"
    assert binaer_zerlegung(221).bits == "11011101"
    assert binaer_zerlegung(100).exponenten == [2, 5, 6]
    assert binaer_zerlegung(39).exponenten == [0, 1, 2, 5]  # 1+2+4+32
    assert binaer_zerlegung(931).exponenten == [0, 1, 5, 7, 8, 9]
    assert binaer_zerlegung(13).exponenten == [0, 2, 3]    # 1 + 2^2 + 2^3


def test_tabelle_12_mod_34_symmetrisch():
    werte = [w for _, _, w in quadrat_tabelle(12, 34, 7, "symmetrisch")]
    assert werte == [12, 8, -4, 16, -16, -16, -16]


def test_tabelle_11_mod_34():
    werte = [w for _, _, w in quadrat_tabelle(11, 34, 5)]
    assert werte == [11, 19, 21, 33, 1]


def test_tabelle_715_mod_2803():
    werte = [w for _, _, w in quadrat_tabelle(715, 2803, 7)]
    # 715^16 ≡ 1653 (2557² = 2332·2803 + 1653).
    assert werte[4:] == [1653, 2287, 2774]
    erg = schnell_potenzieren(715, 113, 2803)
    assert [w for _, w in erg.faktoren] == [715, 1653, 2287, 2774]


def test_tabelle_113_mod_2802():
    werte = [w for _, _, w in quadrat_tabelle(113, 2802, 10)]
    assert werte == [113, 1561, 1783, 1621, 2167, 2539, 1921, 7, 49, 2401]


def test_tabelle_708_mod_2803():
    werte = [w for _, _, w in quadrat_tabelle(708, 2803, 11)]
    assert werte == [708, 2330, 2292, 442, 1957, 951, 1835, 822, 161, 694, 2323]


def test_tabelle_13_mod_100():
    werte = [w for _, _, w in quadrat_tabelle(13, 100, 6)]
    assert werte[:3] == [13, 69, 61] and werte[5] == 81


def test_tabelle_8_mod_101():
    werte = [w for _, _, w in quadrat_tabelle(8, 101, 4)]
    assert werte == [8, 64, 56, 5]


@pytest.mark.parametrize("reste", RESTE)
def test_ordnungen_mod_17(reste):
    erwartet = {2: (8, False, 9), 3: (16, True, 6), 4: (4, False, 13), 5: (16, True, 7)}
    for a, (o, prim, inv) in erwartet.items():
        erg = ordnung(a, 17, reste)
        assert (erg.ordnung, erg.primitiv, erg.inverse) == (o, prim, inv)
        ordnung_text(erg)


def test_potenzliste_mod_17_symmetrisch():
    assert ordnung(2, 17).potenzen == [2, 4, 8, -1, -2, -4, -8, 1]
    assert ordnung(3, 17).potenzen[:8] == [3, -8, -7, -4, 5, -2, -6, -1]
    assert ordnung(4, 17).potenzen == [4, -1, -4, 1]
    assert ordnung(5, 17).potenzen[:8] == [5, 8, 6, -4, -3, 2, -7, -1]


def test_ordnung_z7():
    # ord(2) = 3, ord(3) = 6, 3 ist primitiv in Z_7*.
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
