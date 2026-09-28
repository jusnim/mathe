"""Tests für skripte/schranken.py (Hamming- und Singleton-Schranke)."""

import pytest

from skripte.schranken import (
    e_aus_d,
    eingabe_fehler,
    hamming,
    kleinstes_n,
    main,
    singleton,
    text_code,
    text_kleinstes_n,
)


# Golay-Codes sind perfekt
@pytest.mark.parametrize(
    "n,k,d,q,e,kugel,rest",
    [
        (23, 12, 7, 2, 3, 2048, 2**11),   # 1 + 23 + 253 + 1771 = 2048
        (11, 6, 5, 3, 2, 243, 3**5),      # 1 + 22 + 220 = 243
    ],
)
def test_golay_perfekt(n, k, d, q, e, kugel, rest):
    h = hamming(n, k, d, q)
    assert h.e == e
    assert h.kugel == kugel
    assert h.rest == rest
    assert h.perfekt
    assert h.code_groesse * h.kugel == q**n


def test_golay_summanden_23():
    assert [s.wert for s in hamming(23, 12, 7, 2).summanden] == [1, 23, 253, 1771]


# H_2(m) ist ein [2^m − 1, 2^m − 1 − m, 3]_2-Code und perfekt
@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_hamming_code_perfekt(m):
    n = 2**m - 1
    h = hamming(n, n - m, 3, 2)
    assert h.kugel == 1 + n == 2**m
    assert h.perfekt


# [10, 8, 3]_11: Hamming erfüllt, MDS
def test_10_8_3_11():
    h = hamming(10, 8, 3, 11)
    assert (h.e, h.kugel, h.rest) == (1, 101, 121)
    assert h.erfuellt and not h.perfekt
    assert h.rest - h.kugel == 20
    s = singleton(10, 8, 3)
    assert s.erfuellt and s.mds


def test_10_6_5_11():
    h = hamming(10, 6, 5, 11)
    assert (h.e, h.kugel, h.rest) == (2, 4601, 14641)
    assert h.rest - h.kugel == 10040
    assert h.erfuellt and not h.perfekt
    # [10, 6, 5]_11 ist optimal, denn 6 + 5 = 10 + 1.
    s = singleton(10, 6, 5)
    assert s.mds
    assert "[10, 6, 5]_11 ist optimal" in text_code(10, 6, 5, 11)


# Tabelle: n = 11: 67 > 64, n = 12: 79 ≤ 128. Also n ≥ 12.
def test_kleinstes_n_5_5_2():
    r = kleinstes_n(5, 5, 2)
    assert r.n_singleton == 9
    assert [(z.n, z.kugel, z.rest, z.erfuellt) for z in r.zeilen] == [
        (9, 46, 16, False),
        (10, 56, 32, False),
        (11, 67, 64, False),
        (12, 79, 128, True),
    ]
    assert r.n_min == 12
    assert "Ergebnis: n ≥ 12" in text_kleinstes_n(r)


def test_e_aus_d():
    assert [e_aus_d(d) for d in (1, 2, 3, 4, 5, 7)] == [0, 0, 1, 1, 2, 3]


def test_verletzte_schranken():
    assert not singleton(10, 8, 4).erfuellt
    assert not hamming(11, 5, 5, 2).erfuellt
    assert "gibt es nicht" in text_code(11, 5, 5, 2)


def test_eingabe_pruefung(capsys):
    assert eingabe_fehler(10, 8, 3, 11) == []
    assert any("Primzahlpotenz" in f for f in eingabe_fehler(10, 8, 3, 6))
    assert any("größer als n" in f for f in eingabe_fehler(5, 8, 3, 2))
    assert main(["--n", "10", "--k", "8", "--d", "3", "--q", "6"]) == 1
    assert "Primzahlpotenz" in capsys.readouterr().out


def test_ausgabe_ergebnis(capsys):
    assert main(["--n", "23", "--k", "12", "--d", "7", "--q", "2"]) == 0
    out = capsys.readouterr().out
    assert "Ergebnis:" in out and "perfekt: ja" in out


def test_epilog_beispiele(capsys):
    with pytest.raises(SystemExit):
        main(["--help"])
    out = capsys.readouterr().out
    assert "Beispiele:" in out
    assert "python -m skripte.schranken --k 5 --d 5 --q 2  (kleinstes n bestimmen)" in out
