"""Tests für skripte/euklid.py mit den Aufgaben aus den Quellen (PRIORISIERUNG.md, Abschnitt 3)."""

import pytest

from skripte.euklid import (
    KeinInversesFehler, bezout, bezout_mit_weg, ea_matrix, ea_rueckwaerts, erweiterter_ea,
    euklid, ggt, inverse, inverse_mit_weg, kgv, kgv_mit_weg, loese_diophantisch, main,
    repraesentant,
)


# --- EA und ggT -----------------------------------------------------------

def test_vl1_1965_225_ea():
    ea = euklid(1965, 225)
    assert [(z.dividend, z.q, z.divisor, z.rest) for z in ea.zeilen] == [
        (1965, 8, 225, 165), (225, 1, 165, 60), (165, 2, 60, 45), (60, 1, 45, 15), (45, 3, 15, 0)]
    assert ea.ggt == 15


def test_kv_a3_ggt():
    ea = euklid(2406, 654)
    assert ea.quotienten == [3, 1, 2, 8, 1, 3]
    assert ea.ggt == 6


def test_klausur_a3_ggt():
    assert ggt(17, 15) == 1


def test_ue2_a4_tabelle():
    ea, tab = erweiterter_ea(1001, 840)
    assert [t.q for t in tab] == [None, 1, 5, 4, 1, 1, 2, None]
    assert [t.r for t in tab] == [1001, 840, 161, 35, 21, 14, 7, 0]
    for t in tab[:-1]:
        assert t.x * 1001 + t.y * 840 == t.r


def test_ggt_reihenfolge_und_vorzeichen():
    assert ggt(654, 2406) == 6
    assert ggt(-12, 18) == 6
    assert ggt(0, 5) == 5
    with pytest.raises(ValueError):
        ggt(0, 0)


def test_kgv():
    assert kgv(12, 18) == 36
    assert kgv(1965, 225) == 29475
    assert "Ergebnis: kgV(12, 18) = 36" in kgv_mit_weg(12, 18).text


# --- Bezout ---------------------------------------------------------------

def test_vl1_matrix_q():
    mx = ea_matrix(1965, 225)
    assert mx.Q == ((131, 35), (15, 4))
    assert mx.det == -1
    # VL1 Z. 501: 1965·(−4) + 225·35 = 15
    assert bezout(1965, 225, "matrix") == (15, -4, 35)


def test_kv_a3_matrix_q():
    mx = ea_matrix(2406, 654)
    assert mx.Q == ((401, 103), (109, 28))
    assert mx.det == 1


def test_kv_a7_matrix_q():
    mx = ea_matrix(100, 13)
    assert mx.Q == ((100, 23), (13, 3))
    assert mx.det == 1


@pytest.mark.parametrize("weg", ["rueckwaerts", "matrix"])
@pytest.mark.parametrize("a,b", [(1965, 225), (2406, 654), (17, 15), (6, 11), (840, 1001), (-12, 18), (7, 7)])
def test_bezout_gleichung(a, b, weg):
    g, x, y = bezout(a, b, weg)
    assert a * x + b * y == g == ggt(a, b)
    assert "Ergebnis:" in bezout_mit_weg(a, b, weg).text


def test_ue3_a1_rueckwaerts():
    # Ü3 A1 a): 1 = 5 − 2·2 = 5 − 2·(7 − 1·5) = 3·5 − 2·7
    rw = ea_rueckwaerts(7, 5)
    assert (rw.koef_gross, rw.koef_klein) == (-2, 3)


# --- Diophantische Gleichungen --------------------------------------------

def test_kv_a3_rueckwaerts():
    e = loese_diophantisch(2406, 654, 24, "-", "rueckwaerts")
    assert (e.x, e.y) == (3, 11)
    assert (e.dx, e.dy) == (109, 401)   # KV A3 c): x = 3 + 109t, y = 11 + 401t
    assert "24 = 3·2406 − 11·654" in e.text


def test_kv_a3_matrix():
    e = loese_diophantisch(2406, 654, 24, "-", "matrix")
    assert (e.x, e.y) == (112, 412)
    assert (e.dx, e.dy) == (109, 401)
    assert "2406·112 − 654·412 = 24" in e.text


def test_ue2_a4_i_und_iii():
    e = loese_diophantisch(1001, 840, 98, "+", "matrix")
    assert (e.x, e.y) == (658, -784)
    assert (e.dx, e.dy) == (-120, 143)
    assert (e.x_klein, e.y_klein) == (58, -69)   # „z. B. x = (58, −69)“


def test_ue2_a4_ii_homogen():
    e = loese_diophantisch(1001, 840, 0, "+", "matrix")
    assert 1001 * e.x + 840 * e.y == 0
    assert 1001 * e.dx + 840 * e.dy == 0
    assert (e.dx, e.dy) == (-120, 143)


@pytest.mark.parametrize("weg", ["rueckwaerts", "matrix"])
def test_klausur_a3(weg):
    e = loese_diophantisch(17, 15, 1, "-", weg)
    assert 17 * e.x - 15 * e.y == 1
    assert (e.x, e.y) == (-7, -8)
    assert (e.x_klein, e.y_klein) == (8, 9)


@pytest.mark.parametrize("weg", ["rueckwaerts", "matrix"])
@pytest.mark.parametrize("a,b", [(84, 54), (78, 48), (18, 12)])
def test_gp_a3_form_rechte_seite_30(a, b, weg):
    # GP A3: a·x − b·y = 30. Die Werte sind nicht überliefert. Wir wählen selbst.
    e = loese_diophantisch(a, b, 30, "-", weg)
    assert e.loesbar
    assert a * e.x - b * e.y == 30
    x1, y1 = e.loesung(5)
    assert a * x1 - b * y1 == 30


def test_nicht_loesbar():
    e = loese_diophantisch(84, 54, 31, "-")
    assert not e.loesbar
    assert "keine ganzzahlige Lösung" in e.text


# --- Inverse --------------------------------------------------------------

@pytest.mark.parametrize("a,m,erwartet", [(6, 11, 2), (6, 17, 3), (3, 10, 7), (5, 12, 5)])
@pytest.mark.parametrize("weg", ["rueckwaerts", "matrix"])
def test_ue2_a5_inverse(a, m, erwartet, weg):
    assert inverse_mit_weg(a, m, "standard", weg).inverse == erwartet


def test_ue2_a5_c_symmetrisch():
    # Ü2 A5 c): 3⁻¹ = −3 = 7 in ℤ_10
    e = inverse_mit_weg(3, 10, "symmetrisch")
    assert e.roh == -3 and e.inverse == -3
    assert inverse(3, 10) == 7


def test_kv_a7_d():
    # KV A7: d = −23 ≡ 77 (mod 100)
    e = inverse_mit_weg(13, 100, "standard", "matrix")
    assert e.roh == -23 and e.inverse == 77


def test_ue6_inversentabelle_mod_11_symmetrisch():
    # Ü6 6_9: Inverse mod 11 mit Repräsentanten −5 … 5
    for a in range(1, 11):
        inv = inverse(a, 11, "symmetrisch")
        assert -5 <= inv <= 5
        assert (a * inv) % 11 == 1


def test_inverse_fehler():
    with pytest.raises(KeinInversesFehler) as f:
        inverse(6, 12)
    assert f.value.ggt == 6
    with pytest.raises(KeinInversesFehler):
        inverse(0, 7)


def test_repraesentant():
    assert repraesentant(-1339, 2802) == 1463   # Ü4 A9: d = −1339 ≡ 1463
    assert repraesentant(7, 10, "symmetrisch") == -3


# --- Kommandozeile --------------------------------------------------------

def test_cli_laeuft(capsys):
    assert main(["--a", "2406", "--b", "654", "--c", "24", "--form", "-"]) == 0
    out = capsys.readouterr().out
    assert "x = 3, y = 11" in out and "x = 112, y = 412" in out
    assert main(["--a", "6", "--mod", "12"]) == 1
    assert "Kein Inverses" in capsys.readouterr().out
    assert main(["--a", "12", "--b", "18", "--kgv", "--erweitert"]) == 0
