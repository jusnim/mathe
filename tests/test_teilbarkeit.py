"""Tests für skripte/teilbarkeit.py mit Aufgaben aus Ü1 und Ü4."""

import pytest

from skripte.teilbarkeit import (
    VARIANTEN,
    ausfuehren,
    bloecke,
    main,
    mersenne_daten,
    parser_bauen,
    regel_1001,
    regel_11,
    regel_9,
    teiler_daten,
    zahl_lesen,
)

X_LICHT = 299792          # Ü1 A1, A2
X_GROSS = 10**19 + 1      # Ü1 A2


# --- Ü1 A1: Regel für 9 -----------------------------------------------------

def test_ue1_a1_neun():
    r = regel_9(X_LICHT)
    assert r.ziffern == [2, 9, 9, 7, 9, 2]
    assert r.quersumme == 38
    assert r.rest == 2          # Quelle: ≡ 2 + 7 + 2 ≡ 2 (mod 9)
    assert r.ohne_neunen == [2, 7, 2]
    assert not r.teilbar


# --- Ü1 A1: Regel für 11, beide Varianten ------------------------------------

def test_ue1_a1_elf_standard():
    r = regel_11(X_LICHT)
    assert r.variante == "einer"
    assert r.terme == [-2, 9, -9, 7, -9, 2]   # Quelle: −2 + 9 − 9 + 7 − 9 + 2
    assert r.summe == -2
    assert r.reste_x[11] == 9
    assert not r.teilbar[11]


def test_ue1_a1_elf_links():
    r = regel_11(X_LICHT, "links")
    assert r.terme == [2, -9, 9, -7, 9, -2]
    assert r.summe == 2
    assert r.faktor == -1                      # Summe ≡ −x
    assert r.reste_x[11] == 9
    assert not r.teilbar[11]


@pytest.mark.parametrize("variante", list(VARIANTEN))
@pytest.mark.parametrize("x", [X_LICHT, X_GROSS, 121, 918082, 1, 0, 1001, 7 * 10**12])
def test_elf_und_1001_alle_varianten(x, variante):
    r11 = regel_11(x, variante)
    assert (r11.summe - r11.faktor * x) % 11 == 0
    assert r11.teilbar[11] == (x % 11 == 0)
    r7 = regel_1001(x, variante)
    for p in (7, 11, 13):
        assert (r7.summe - r7.faktor * x) % p == 0
        assert r7.teilbar[p] == (x % p == 0)
        assert r7.reste_x[p] == x % p


# --- Ü1 A2: Regel für 7 mit 1001 ---------------------------------------------

def test_ue1_a2_sieben_299792():
    r = regel_1001(X_LICHT)
    assert r.teile == [299, 792]
    assert r.terme == [-299, 792]              # Quelle: −299 + 792
    assert r.summe == 493                      # Quelle: 493 = 7 · 70 + 3
    assert r.reste_summe[7] == 3
    assert not r.teilbar[7]


def test_ue1_a2_sieben_10hoch19_plus_1():
    r = regel_1001(X_GROSS)
    assert r.teile == [10, 0, 0, 0, 0, 0, 1]
    assert r.summe == 11                       # Quelle: ≡ 11 (mod 7)
    assert r.reste_x[7] == 4
    assert not r.teilbar[7]


@pytest.mark.parametrize("variante", list(VARIANTEN))
def test_ue1_a2_varianten(variante):
    r = regel_1001(X_LICHT, variante)
    assert abs(r.summe) == 493
    assert r.reste_x[7] == 3
    r = regel_1001(X_GROSS, variante)
    assert r.summe == 11                       # m = 6 gerade: beide Varianten gleich


def test_bloecke():
    assert bloecke(299792) == [299, 792]
    assert bloecke(1001) == [1, 1]
    assert bloecke(5) == [5]


# --- Ü4 A3: d(n) und σ(n) ----------------------------------------------------

@pytest.mark.parametrize("n, d, sigma", [
    (1, 1, 1), (6, 4, 12), (28, 6, 56), (496, 10, 992), (12, 6, 28), (60, 12, 168), (13, 2, 14),
])
def test_ue4_a3_teilerfunktionen(n, d, sigma):
    t = teiler_daten(n)
    assert t.d == d
    assert t.sigma == sigma


def test_ue4_a3_formel_28():
    t = teiler_daten(28)
    assert t.faktoren == {2: 2, 7: 1}
    assert t.d_formel_faktoren == [3, 2]
    assert t.sigma_formel_faktoren == [7, 8]
    assert t.teiler == [1, 2, 4, 7, 14, 28]


# --- Ü4 A4, A6 b: vollkommene Zahlen -----------------------------------------

@pytest.mark.parametrize("n", [6, 28, 496, 8128])
def test_ue4_a4_vollkommen(n):
    assert teiler_daten(n).art == "vollkommen"


@pytest.mark.parametrize("n, art", [(12, "abundant"), (8, "defizient"), (1, "defizient")])
def test_nicht_vollkommen(n, art):
    assert teiler_daten(n).art == art


# --- Ü4 A5, A6: Mersenne -----------------------------------------------------

@pytest.mark.parametrize("k, m, vollkommen", [(2, 3, 6), (3, 7, 28), (5, 31, 496), (7, 127, 8128)])
def test_ue4_a5_a6_mersenne_prim(k, m, vollkommen):
    md = mersenne_daten(k)
    assert md.m == m and md.ist_prim
    assert md.vollkommen == vollkommen
    assert teiler_daten(vollkommen).sigma == 2 * vollkommen


def test_ue4_a5_2hoch11():
    md = mersenne_daten(11)
    assert md.m == 2047
    assert not md.ist_prim
    assert md.faktoren == {23: 1, 89: 1}      # Quelle: 2047 = 23 · 89


def test_ue4_a6c_k_zusammengesetzt():
    md = mersenne_daten(6)
    assert md.k_zerlegung == (2, 3)
    assert md.m % (2**2 - 1) == 0


# --- Eingabe und Kommandozeile -----------------------------------------------

@pytest.mark.parametrize("text, wert", [
    ("299792", 299792), ("299 792", 299792), ("10^19 + 1", X_GROSS),
    ("10**19+1", X_GROSS), ("2^11 - 1", 2047), ("1000 0000 0000 0000 0001", X_GROSS),
])
def test_zahl_lesen(text, wert):
    assert zahl_lesen(text) == wert


def test_zahl_lesen_fehler():
    with pytest.raises(ValueError):
        zahl_lesen("abc")


def test_ausgabe_regeln():
    args = parser_bauen().parse_args(["regeln", "299792"])
    text = ausfuehren(args)
    assert "−2 + 9 − 9 + 7 − 9 + 2" in text
    assert "−299 + 792" in text
    assert "493 = 7 · 70 + 3" in text
    assert text.count("Ergebnis:") == 3


def test_ausgabe_alle_varianten():
    args = parser_bauen().parse_args(["elf", "299792", "--alle-varianten"])
    text = ausfuehren(args)
    assert "Variante „einer“" in text and "Variante „links“" in text


def test_hinweis_negativ_und_null():
    text = ausfuehren(parser_bauen().parse_args(["teiler", "0"]))
    assert "Hinweis" in text


def test_main_fehlercode(capsys):
    assert main(["neun", "abc"]) == 2
    assert "Fehler" in capsys.readouterr().err
