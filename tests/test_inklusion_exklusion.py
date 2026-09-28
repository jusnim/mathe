"""Tests für skripte/inklusion_exklusion.py mit Beispielrechnungen."""

from fractions import Fraction

import pytest

from skripte.inklusion_exklusion import (
    berechne_derangement,
    berechne_mengen,
    berechne_teiler,
    genau_in_mengen,
    kgv,
    main,
    text_derangement,
    text_mengen,
    text_teiler,
)


# --- Modus Teiler -----------------------------------------------------------

@pytest.mark.parametrize("variante", ["kgv", "produkt"])
def test_teiler_720_3_4_5(variante):
    # 720 − 564 + 144 − 12 = 288. 3, 4, 5 sind paarweise teilerfremd,
    # deshalb geben beide Varianten dasselbe.
    e = berechne_teiler(720, [3, 4, 5], variante)
    assert e.alphas == [564, 144, 12]
    assert e.ergebnis == 288
    assert e.probe == 288
    assert e.produktformel == 288 and e.produktformel_exakt


@pytest.mark.parametrize("variante", ["kgv", "produkt"])
def test_teiler_840_2_3_5(variante):
    e = berechne_teiler(840, [2, 3, 5], variante)
    assert e.ergebnis == 224
    assert e.probe == 224
    assert e.produktformel_exakt


@pytest.mark.parametrize("variante", ["kgv", "produkt"])
def test_teiler_1000_2_3_5(variante):
    e = berechne_teiler(1000, [2, 3, 5], variante)
    assert [s.anzahl for s in e.stufen[0]] == [500, 333, 200]
    assert [s.anzahl for s in e.stufen[1]] == [166, 100, 66]
    assert e.alphas == [1033, 332, 33]
    assert e.ergebnis == 266
    # Produktformel nur näherungsweise: 800/3 = 266,66…
    assert e.produktformel == Fraction(800, 3)
    assert not e.produktformel_exakt


def test_nicht_teilerfremd_varianten_unterscheiden_sich():
    # Bei 2 und 4 ist nur kgV richtig: Zahlen ≤ 100, die nicht durch 2 oder 4 teilbar sind = 50.
    richtig = berechne_teiler(100, [2, 4], "kgv")
    falsch = berechne_teiler(100, [2, 4], "produkt")
    assert richtig.ergebnis == 50 == richtig.probe
    assert falsch.ergebnis == 100 - 75 + 12  # 37, falsch
    assert falsch.richtiges_ergebnis == 50
    assert "FALSCH" in text_teiler(falsch)


def test_text_teiler_enthaelt_rechenweg():
    t = text_teiler(berechne_teiler(720, [3, 4, 5]))
    assert "720 − 564 + 144 − 12 = 288" in t
    assert "Ergebnis: 288" in t
    assert "A₃ ∩ A₄" in t
    t2 = text_teiler(berechne_teiler(720, [3, 4, 5]), indizes="nummer")
    assert "A₁ ∩ A₂" in t2


def test_eingabe_fehler():
    with pytest.raises(ValueError):
        berechne_teiler(0, [2])
    with pytest.raises(ValueError):
        berechne_teiler(10, [])
    e = berechne_teiler(30, [2, 2, 3])
    assert e.teiler == [2, 3]
    assert any("doppelt" in h for h in e.hinweise)


def test_kgv():
    assert kgv(4, 6) == 12
    assert kgv(3, 5) == 15


# --- Modus Mengen -----------------------------------------------------------

def test_mengen_sportverein_55():
    e = berechne_mengen(55, [[35, 27, 12], [13, 7, 5], [2]])
    assert e.alphas == [74, 25, 2]
    assert e.ergebnis == 4
    assert sum(genau_in_mengen(e)) == 55
    assert "Ergebnis: 4" in text_mengen(e)


def test_mengen_musikschule_73():
    e = berechne_mengen(73, [[20, 25, 52], [7, 12, 17], [1]], ["Flöte", "Geige", "Klavier"])
    assert e.alphas == [97, 36, 1]
    assert e.ergebnis == 11
    assert "73 − (20 + 25 + 52) + (7 + 12 + 17) − 1" in text_mengen(e)


def test_mengen_falsche_anzahl():
    with pytest.raises(ValueError):
        berechne_mengen(55, [[35, 27, 12], [13, 7]])


def test_mengen_fehlende_stufe_wird_null():
    e = berechne_mengen(10, [[3, 4]])
    assert e.ergebnis == 3
    assert e.hinweise


# --- Modus Derangement ------------------------------------------------------

def test_derangement_4():
    e = berechne_derangement(4)
    assert e.d_n == 9
    assert e.wahrscheinlichkeit == Fraction(3, 8)
    assert len(e.liste) == 9
    assert (2, 3, 4, 1) in e.liste and (4, 3, 2, 1) in e.liste
    assert "Ergebnis: d₄ = 9" in text_derangement(e)


@pytest.mark.parametrize("n, d", [(1, 0), (2, 1), (3, 2), (5, 44), (6, 265)])
def test_derangement_werte(n, d):
    assert berechne_derangement(n).d_n == d


# --- Kommandozeile ----------------------------------------------------------

def test_main_alle_varianten(capsys):
    assert main(["--n", "720", "--teiler", "3", "4", "5", "--alle-varianten"]) == 0
    aus = capsys.readouterr().out
    assert "Variante: kgv" in aus and "Variante: produkt" in aus
    assert aus.count("Ergebnis: 288") == 2


def test_main_fehler(capsys):
    assert main(["--n", "0", "--teiler", "2"]) == 2
    assert "Fehler" in capsys.readouterr().err
