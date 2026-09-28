"""Tests für skripte/mod11_code.py mit Aufgaben aus PRIORISIERUNG.md, Abschnitt 3."""

import pytest

from skripte.mod11_code import (
    VARIANTEN,
    decodieren,
    inverse_mod,
    inversentabelle,
    kodieren,
    kontrollmatrix,
    loesungsweg,
    main,
    null_einfuegen,
    rep,
    stelle_in_nullblock,
    ziffern_lesen,
)

W = ziffern_lesen

# Klausur A5
K_A = W("11500005")
K_R1 = W("1150000511")
K_R2 = W("1150000733")
K_C2 = W("1150020733")


def test_inversentabelle_ue6_6_9():
    # Ü6 6_9 / wissen/06_uebung.tex Z. 28, in Vertretern −5 … 5
    erwartet = {1: 1, 2: -5, 3: 4, 4: 3, 5: -2, -5: 2, -4: -3, -3: -4, -2: 5, -1: -1}
    tab = inversentabelle()
    for a, inv in erwartet.items():
        assert rep(tab[a % 11], True) == inv
    assert inverse_mod(6, 11) == 2
    with pytest.raises(ValueError):
        inverse_mod(0, 11)


def test_kodierung_klausur_a5():
    k = kodieren(K_A)
    assert (k.c9_summe, k.c9) == (70, 4)
    assert (k.c10_summe, k.c10) == (50, 6)
    assert k.c == W("1150000546")
    assert k.c_ueber_G == k.c  # Weg über G = (E | −Aᵀ) gibt dasselbe


def test_systematische_form_vl6():
    # VL6 Z. 224–229: (A | E)
    km = kontrollmatrix("systematisch-rechts")
    assert km.H == [[9, 8, 7, 6, 5, 4, 3, 2, 1, 0], [3, 4, 5, 6, 7, 8, 9, 10, 0, 1]]
    assert [row[8:] for row in km.G] == [[2, 8], [3, 7], [4, 6], [5, 5], [6, 4], [7, 3], [8, 2], [9, 1]]


# Syndrome von K r₁ und r₂ je Variante (PRIORISIERUNG.md, Abschnitt 4, Punkt 5)
SYNDROME_K = {
    "grundform": ((3, 0), (9, 10)),
    "grundform-getauscht": ((0, 3), (10, 9)),
    "systematisch-rechts": ((8, 6), (3, 6)),
    "systematisch-links": ((6, 8), (8, 1)),
}


@pytest.mark.parametrize("variante", list(VARIANTEN))
def test_klausur_a5_decodierung(variante):
    s1, s2 = SYNDROME_K[variante]
    d1 = decodieren(K_R1, variante)
    assert d1.S == s1
    assert d1.status == "nicht_decodierbar"
    d2 = decodieren(K_R2, variante)
    assert d2.S == s2
    assert d2.status == "ein_fehler"
    assert (d2.x, d2.e) == (6, 9)
    assert d2.c == K_C2
    assert all(s % 11 == 0 for s in d2.probe_summen)


def test_klausur_r1_hinweis_x_null():
    # Abschnitt 5.2: S = (3, 0) gibt x = 0. Die VL nennt diesen Fall nicht.
    d = decodieren(K_R1)
    assert d.x == 0 and d.hinweis_quelle


@pytest.mark.parametrize("variante", list(VARIANTEN))
def test_ue6_6_9(variante):
    d1 = decodieren(W("1151111103"), variante)
    assert (d1.x, d1.e) == (3, 4)
    assert d1.c == W("1111111103")
    d2 = decodieren(W("1156111103"), variante)
    assert (d2.x, rep(d2.e, True)) == (6, -2)
    assert d2.c == W("1156131103")
    d3 = decodieren(W("1111111103"), variante)
    assert d3.status == "codewort" and d3.S == (0, 0)


def test_ue6_6_9_syndrome_grundform():
    assert decodieren(W("1151111103")).S == (4, 1)          # (15, 4·3) ≡ (4, 1)
    d2 = decodieren(W("1156111103"))
    assert (rep(d2.S[0], True), rep(d2.S[1], True)) == (-2, -1)


@pytest.mark.parametrize("variante", list(VARIANTEN))
def test_vl6_beispiele(variante):
    # wissen/06_decodierung.tex Z. 262 und Z. 270
    d = decodieren(W("0000000050"), variante)
    assert (d.x, d.e) == (9, 5) and d.c == [0] * 10
    d = decodieren(W("3500000000"), variante)
    assert (d.x, d.e) == (3, 8) and d.c == W("3530000000")


def test_vl6_syndrome_grundform():
    assert decodieren(W("0000000050")).S == (5, 1)   # (5, 45)
    assert decodieren(W("3500000000")).S == (8, 2)   # (8, 13)


def test_nicht_decodierbar_s1_null():
    # VL6 Z. 279: s₁ = 0, s₂ ≠ 0 (zwei Fehler, Summe der Fehler 0)
    r = W("1500000000")
    r[1] = 10  # r = 1 10 0 … → s₁ = 11 ≡ 0, s₂ = 21 ≡ 10
    d = decodieren(r)
    assert d.S == (0, 10) and d.status == "nicht_decodierbar"


def test_gp_a4_fehlende_null():
    # GP A4: Wörter mit 9 Ziffern (Abschnitt 5.1). Mit ergänzter 0 im Nullblock:
    r1 = W("115000711")
    r2 = W("115000733")
    s = stelle_in_nullblock(r1)
    assert null_einfuegen(r1, s) == W("1150000711")
    assert null_einfuegen(r2, stelle_in_nullblock(r2)) == K_R2
    d2 = decodieren(null_einfuegen(r2, s))
    assert d2.c == K_C2
    # Hinweis: Mit dieser Annahme ist auch r₁ decodierbar (x = 1, Fehlergröße 5).
    # Das Gedächtnisprotokoll ist hier also wohl ungenau erinnert (vgl. K A5: r₁ = …0511).
    d1 = decodieren(null_einfuegen(r1, s))
    assert (d1.x, d1.e) == (1, 5) and d1.c == W("7150000711")


def test_ausgabe_und_kommandozeile(capsys):
    text = loesungsweg(K_A, [("r₁", K_R1), ("r₂", K_R2)], "grundform")
    assert "Ergebnis:" in text and "NICHT decodierbar" in text and "1 1 5 0 0 2 0 7 3 3" in text
    assert main(["--a", "11500005", "--r", "1150000511", "1150000733", "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert out.count("=== Variante:") == 4
    assert main(["--r", "115000711", "--r", "115000733"]) == 0
    assert "WARNUNG" in capsys.readouterr().out
    assert main(["--a", "1150000"]) == 1
    assert "FEHLER" in capsys.readouterr().out
    assert main(["--r", "12345"]) == 0
    assert "übersprungen" in capsys.readouterr().out


def test_eingabe_x_fuer_zehn():
    assert W("123456789X") == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    assert W("1,1,10,0") == [1, 1, 10, 0]
