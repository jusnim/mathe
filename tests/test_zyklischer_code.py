"""Tests für skripte/zyklischer_code.py."""

import pytest

from skripte.zyklischer_code import (
    VARIANTEN,
    EingabeFehler,
    code_daten,
    codieren,
    decodieren,
    fmt_menge,
    fmt_poly,
    irreduzible_polynome,
    kreisteilungsklassen,
    lies_polynom,
    loese,
    loesungsweg,
    main,
    p_divmod,
    produkt_modulo,
    x_hoch_n_minus_1,
    zerlege,
)


def faktoren(n, p):
    """Irreduzible Faktoren als Menge von Tupeln (a_0, a_1, …) mit Vielfachheit."""
    return sorted((tuple(f), m) for f, m in zerlege(n, p).faktoren)


def P(text, p):
    return tuple(lies_polynom(text, p))


# --- Irreduzible Polynome und Zerlegung ----------------------------------------

def test_irreduzibel_grad_3_ueber_f2():
    erg = irreduzible_polynome(3, 2)
    assert sorted(map(tuple, erg.irreduzibel)) == sorted([P("X^3+X+1", 2), P("X^3+X^2+1", 2)])


def test_irreduzibel_grad_2_und_4_ueber_f2():
    assert [tuple(f) for f in irreduzible_polynome(2, 2).irreduzibel] == [P("X^2+X+1", 2)]
    assert sorted(map(tuple, irreduzible_polynome(4, 2).irreduzibel)) == sorted(
        [P("X^4+X+1", 2), P("X^4+X^3+1", 2), P("X^4+X^3+X^2+X+1", 2)])


def test_zerlegung_x7_ueber_f2():
    assert faktoren(7, 2) == sorted([(P("X+1", 2), 1), (P("X^3+X+1", 2), 1), (P("X^3+X^2+1", 2), 1)])
    klassen = kreisteilungsklassen(7, 2)
    assert sorted(map(frozenset, klassen), key=min) == [frozenset({0}), frozenset({1, 2, 4}), frozenset({3, 5, 6})]
    # Symmetrische Schreibweise: {1, 2, -3}
    assert fmt_menge([1, 2, 4], 7, "symmetrisch") == "{1, 2, -3}"


# --- Hamming- und Paritätscode -----------------------------------------------

def test_hamming_7_4_zwei_generatoren():
    los = loese(2, n=7, k=4)
    gs = [tuple(f.code.g) for f in los.faelle]
    assert gs == [P("1+X+X^3", 2), P("1+X^2+X^3", 2)]
    G1 = los.faelle[0].code.G
    assert G1 == [[1, 1, 0, 1, 0, 0, 0], [0, 1, 1, 0, 1, 0, 0], [0, 0, 1, 1, 0, 1, 0], [0, 0, 0, 1, 1, 0, 1]]
    G2 = los.faelle[1].code.G
    # X^3 * g mit g = 1 + X^2 + X^3 ist 0001011.
    assert G2 == [[1, 0, 1, 1, 0, 0, 0], [0, 1, 0, 1, 1, 0, 0], [0, 0, 1, 0, 1, 1, 0], [0, 0, 0, 1, 0, 1, 1]]
    assert all(f.code.d == 3 for f in los.faelle)


def test_paritaetscode_7_6():
    los = loese(2, n=7, k=6)
    assert len(los.faelle) == 1
    c = los.faelle[0].code
    assert tuple(c.g) == P("X+1", 2)
    assert c.G == [[1 if j in (i, i + 1) else 0 for j in range(7)] for i in range(6)]
    assert c.d == 2


# --- Ternäre Codes der Länge 4 ----------------------------------------------

def test_zerlegung_x4_ueber_f3():
    assert faktoren(4, 3) == sorted([(P("X-1", 3), 1), (P("X+1", 3), 1), (P("X^2+1", 3), 1)])
    klassen = kreisteilungsklassen(4, 3)
    assert sorted(map(frozenset, klassen), key=min) == [frozenset({0}), frozenset({1, 3}), frozenset({2})]


def test_ternaere_4_2_codes():
    los = loese(3, n=4, k=2)
    assert [tuple(f.code.g) for f in los.faelle] == [P("X^2+1", 3), P("X^2-1", 3)]
    assert los.faelle[0].code.G == [[1, 0, 1, 0], [0, 1, 0, 1]]
    assert los.faelle[1].code.G == [[2, 0, 1, 0], [0, 2, 0, 1]]      # -1 = 2 in F_3
    assert [f.code.d for f in los.faelle] == [2, 2]
    # Kontrollpolynome h(X)
    assert tuple(los.faelle[0].code.h) == P("X^2-1", 3)
    assert tuple(los.faelle[1].code.h) == P("X^2+1", 3)


@pytest.mark.parametrize("variante", list(VARIANTEN))
def test_ternaere_4_2_codes_varianten_schreibweise(variante):
    los = loese(3, n=4, k=2)
    text = loesungsweg(los, variante)
    if variante == "symmetrisch":
        assert "g(X) = -1 + X^2" in text
        assert "( -1  0  1  0 )" in text
    else:
        assert "g(X) = 2 + X^2" in text
        assert "( 2 0 1 0 )" in text
    assert f"Variante '{variante}'" in text


def test_kein_zyklischer_4_2_3_code_ueber_f3():
    los = loese(3, n=4, k=2, d=3)
    text = loesungsweg(los)
    assert "Nein." in text
    assert not any(f.code.d >= 3 for f in los.faelle)


# --- Weitere Zerlegungen, Golay-Codes, Produkt --------------------------------

def test_n4_ueber_f2():
    assert faktoren(4, 2) == [(P("X+1", 2), 4)]
    assert kreisteilungsklassen(4, 2) is None      # ggT(4, 2) ≠ 1
    los = loese(2, n=4)
    assert [g.k for g in los.generatoren] == [4, 3, 2, 1, 0]
    # C_2 = {0000, 1010, 0101, 1111}: G = H
    c2 = code_daten(lies_polynom("(1+X)^2", 2), 4, 2)
    assert c2.G == c2.H == [[1, 0, 1, 0], [0, 1, 0, 1]]


def test_n5_n9_n15_ueber_f2():
    assert faktoren(5, 2) == sorted([(P("X+1", 2), 1), (P("X^4+X^3+X^2+X+1", 2), 1)])
    assert faktoren(9, 2) == sorted([(P("X+1", 2), 1), (P("X^2+X+1", 2), 1), (P("X^6+X^3+1", 2), 1)])
    assert faktoren(15, 2) == sorted([(P(t, 2), 1) for t in
                                      ["X+1", "X^2+X+1", "X^4+X+1", "X^4+X^3+1", "X^4+X^3+X^2+X+1"]])
    assert sorted(len(b) for b in kreisteilungsklassen(15, 2)) == [1, 2, 4, 4, 4]


def test_golay_g23():
    g1 = P("1+X^2+X^4+X^5+X^6+X^10+X^11", 2)
    g2 = P("1+X+X^5+X^6+X^7+X^9+X^11", 2)
    assert faktoren(23, 2) == sorted([(P("X+1", 2), 1), (g1, 1), (g2, 1)])
    c = code_daten(list(g1), 23, 2)
    assert (c.k, c.d) == (12, 7)


def test_golay_g11():
    g1 = P("-1+X^2-X^3+X^4+X^5", 3)
    g2 = P("-1-X+X^2-X^3+X^5", 3)
    assert faktoren(11, 3) == sorted([(P("X-1", 3), 1), (g1, 1), (g2, 1)])
    c = code_daten(list(g1), 11, 3)
    assert (c.k, c.d) == (6, 5)


def test_produkt_f5():
    # 0012 * 2314 = 2102 in F_5[X]/(X^4 - 1)
    pr = produkt_modulo([0, 0, 1, 2], [2, 3, 1, 4], 4, 5)
    assert pr.ergebnis + [0] * (4 - len(pr.ergebnis)) == [2, 1, 0, 2]


def test_paritaetscode_n4_kontrollmatrix():
    c = code_daten([1, 1], 4, 2)
    assert c.G == [[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]]
    assert c.H == [[1, 1, 1, 1]]


# --- Syndromdecodierung ------------------------------------------------------

def test_syndrom_und_decodierung_hamming():
    g = lies_polynom("X^3+X+1", 2)
    de = decodieren(lies_polynom("X^6+X+1", 2), g, 7, 2)
    assert tuple(de.division.quotient) == P("X^3+X+1", 2)
    assert tuple(de.s) == P("X^2+X", 2)
    assert tuple(de.fuehrer) == P("X^4", 2)
    assert tuple(de.c) == P("X^6+X^4+X+1", 2)
    assert tuple(de.c_division.quotient) == P("X^3+1", 2)
    assert de.abstand == 1


def test_codewort_hat_syndrom_null():
    g = lies_polynom("X^3+X+1", 2)
    ko = codieren([1, 0, 0, 1], g, 7, 2)
    de = decodieren(ko.c, g, 7, 2)
    assert de.s == [] and de.c == ko.c


# --- Polynomdivision, Eingabeprüfung, Kommandozeile ----------------------------

def test_division_probe():
    for n, p in [(7, 2), (11, 3), (8, 5)]:
        for f, _ in zerlege(n, p).faktoren:
            assert p_divmod(x_hoch_n_minus_1(n, p), f, p).rest == []


def test_keine_primzahl():
    with pytest.raises(EingabeFehler):
        loese(4, n=7)


def test_g_kein_teiler_gibt_hinweis():
    los = loese(2, n=7, g="X^2+1")
    assert los.faelle == []
    assert any("teilt X^n − 1 nicht" in h for h in los.hinweise)


def test_fmt_poly():
    assert fmt_poly([2, 0, 1], 3, "symmetrisch") == "X^2 - 1"
    assert fmt_poly([2, 0, 1], 3, "standard") == "X^2 + 2"
    assert fmt_poly([1, 1, 0, 1], 2, "symmetrisch", aufsteigend=True) == "1 + X + X^3"


def test_main_laeuft(capsys):
    assert main(["--p", "2", "--n", "7", "--g", "X^3+X+1", "--y", "1100001", "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert "Ergebnis:" in out and "c = 1100101" in out
    assert main(["--p", "9", "--n", "4"]) == 2
