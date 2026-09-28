"""Tests für skripte/isbn10.py."""

import itertools

import pytest

from skripte.isbn10 import (
    VARIANTEN, EingabeFehler, baue_parser, berechne_pruefziffer, lies_nummer,
    loese, pruefe_nummer, zahlendreher, zeichen,
)

# Bekannte Nummern und ihre Prüfziffern
FAELLE = [
    ("3-05-501517", 7),
    ("0-471-82819", 10),  # → X
    ("4-263-79215", 7),
    ("0-387-98999", 4),   # 0-387-98999-4
    ("3-528-06580", 10),  # 3-528-06580-X
]


@pytest.mark.parametrize("variante", VARIANTEN)
@pytest.mark.parametrize("nummer, c10", FAELLE)
def test_pruefziffer_alle_varianten(nummer, c10, variante):
    for rep in ("symmetrisch", "standard"):
        p = berechne_pruefziffer(lies_nummer(nummer), variante, rep)
        assert p.c10 == c10
        assert p.probe_summe % 11 == 0


def test_trick_zwischenschritte_3055015177():
    # 3 + 2.(0-7) + 3.(5-1) + 4.(5-5) + 5.(0-1) = 3 - 3 + 1 - 5 = -4 = 7
    p = berechne_pruefziffer(lies_nummer("3-05-501517"))
    assert p.gewichte == [1, 2, 3, 4, 5, -5, -4, -3, -2]
    assert p.paar_werte == [3, -14, 12, 0, -5]
    assert p.paar_reduziert == [3, -3, 1, 0, -5]
    assert p.summe_reduziert == -4


def test_trick_zwischenschritte_047182819x():
    # 2.(4-9) + 3.(7-1) + 4.(1-8) + 5.(8-2) = 1 + 7 + 5 - 3 = 10 = X
    p = berechne_pruefziffer(lies_nummer("0-471-82819"))
    assert p.paar_werte == [0, -10, 18, -28, 30]
    # Bereich 0 … 10 ergibt 1, 7, 5, 8.
    # Symmetrisch (−5 … 5) ergibt sich 1, −4, 5, −3; die Summe ist ≡ gleich.
    assert berechne_pruefziffer(lies_nummer("0-471-82819"), rep="standard").paar_reduziert \
        == [0, 1, 7, 5, 8]
    assert p.c10 == 10 and zeichen(p.c10) == "X"


def test_trick_zwischenschritte_4263792157():
    # 4 + 2.(2-5) + 3.(6-1) + 4.(3-2) + 5.(7-9) = 4 - 6 + 4 + 4 + 1 = 7
    p = berechne_pruefziffer(lies_nummer("4-263-79215"))
    assert p.paar_werte == [4, -6, 15, 4, -10]
    assert sum(p.paar_werte) % 11 == 7


@pytest.mark.parametrize("variante", VARIANTEN)
@pytest.mark.parametrize("nummer", ["0-387-98999-4", "3-528-06580-X",
                                    "3-05-501517-7", "0-471-82819-X", "4-263-79215-7"])
def test_pruefen_gueltig(nummer, variante):
    assert pruefe_nummer(lies_nummer(nummer), variante).gueltig


def test_pruefen_ungueltig():
    p = pruefe_nummer(lies_nummer("3-05-501517-8"))
    assert not p.gueltig
    assert p.richtige_c10 == 7


def test_dreher_beispiel():
    d = zahlendreher(lies_nummer("0-387-98999-4"), 3, 7)
    # c₃ = 8, c₇ = 9: (3 − 7)(9 − 8) = −4 ≡ 7
    assert d.he == -4 and d.he_mod == 7
    assert d.hr_mod == d.he_mod
    assert d.erkannt


def test_dreher_allgemein():
    # Jeder echte Zahlendreher wird erkannt.
    for nummer in ["0-387-98999-4", "3-528-06580-X", "4-263-79215-7", "0-471-82819-X"]:
        c = lies_nummer(nummer)
        for x, y in itertools.permutations(range(1, 11), 2):
            d = zahlendreher(c, x, y)
            assert d.he == (x - y) * (c[y - 1] - c[x - 1])
            assert d.hr_mod == d.he_mod
            if c[x - 1] != c[y - 1]:
                assert d.erkannt
            else:
                assert d.kein_fehler and d.he_mod == 0


@pytest.mark.parametrize("text", ["12345678", "12345678901", "X23456789", "12a456789"])
def test_eingabefehler(text):
    with pytest.raises(EingabeFehler):
        lies_nummer(text)


def test_dreher_gleiche_stelle():
    with pytest.raises(EingabeFehler):
        zahlendreher(lies_nummer("0-387-98999-4"), 3, 3)


def test_cli_text():
    ap = baue_parser()
    out = loese(ap.parse_args(["--nummer", "3-05-501517", "--alle-varianten"]))
    assert out.count("Ergebnis: c₁₀ = 7") == 3
    out = loese(ap.parse_args(["--nummer", "0-471-82819"]))
    assert "Ergebnis: c₁₀ = X" in out
    out = loese(ap.parse_args(["--dreher-beweis"]))
    assert "(x − y)·(c_y − c_x)" in out
