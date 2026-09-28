"""Tests für skripte/rsa.py mit den RSA-Aufgaben aus PRIORISIERUNG.md (Abschnitt 3 und 4)."""

import pytest

from skripte.rsa import VARIANTEN, main, rsa_berechnen
from skripte.schnell_potenzieren import binaer_zerlegung


def test_klausur_a1():
    # K A1: m = 43·57 = 2451, e = 221, w = 1511.
    # 57 = 3·19 ist keine Primzahl (PRIORISIERUNG 5.2). Richtig: φ = 1512, d = 821.
    erg = rsa_berechnen(p=43, q=57, schluessel=221, w=1511)
    assert erg.m == 2451
    assert binaer_zerlegung(erg.e).bits == "11011101"
    assert erg.c == 896
    assert erg.phi.zerlegung.faktoren == {3: 1, 19: 1, 43: 1}
    assert erg.phi_m == 1512
    assert erg.phi_phi_m == 432
    assert erg.d == 821
    assert erg.euler.ergebnis == 821
    assert erg.entschluesselung.ergebnis == 1511
    # Falscher Weg (p − 1)(q − 1) = 2352 → d = 149, entschlüsselt falsch.
    f = erg.falscher_weg
    assert f.phi_falsch == 2352
    assert f.partner_falsch == 149
    assert f.klartext_falsch != 1511
    assert "57 ist keine Primzahl" in erg.text
    assert "Ergebnis: e = 221 = 11011101₂, c = 896, φ(m) = 1512" in erg.text


def test_uebung4_a9():
    # Ü4 A9: m = 2803 (prim), e = 113, a = 715 → c = 708, d = −1339 ≡ 1463, φ(φ(m)) = 932.
    erg = rsa_berechnen(m=2803, schluessel=113, w=715, klartext_name="a")
    assert binaer_zerlegung(113).bits == "1110001"
    assert erg.c == 708
    assert erg.phi_m == 2802
    assert erg.phi_phi_m == 932
    assert erg.partner_roh == -1339
    assert erg.d == 1463
    assert erg.euler.k == 931
    assert erg.entschluesselung.ergebnis == 715
    assert erg.falscher_weg is None
    # Quadrattabelle Weg 2 wie in der Musterlösung (e^(2^i) mod 2802).
    assert [r for _, _, r in erg.euler.tabelle] == [113, 1561, 1783, 1621, 2167, 2539, 1921, 7, 49, 2401]
    # Probe-Tabelle c^(2^i) mod 2803 wie in der Musterlösung.
    assert [r for _, _, r in erg.entschluesselung.tabelle] == [
        708, 2330, 2292, 442, 1957, 951, 1835, 822, 161, 694, 2323]


def test_kv_a7_beide_wege():
    # KV A7: m = 101, e = 13 → d = −23 ≡ 77; Euler: 13^39, Faktoren 13·69·61·81.
    erg = rsa_berechnen(m=101, schluessel=13, ea_weg="beide")
    assert erg.phi_m == 100 and erg.phi_phi_m == 40
    assert erg.partner_roh == -23
    assert erg.d == 77
    assert erg.euler.k == 39
    assert [f for _, f in erg.euler.faktoren] == [13, 69, 61, 81]
    assert [inv.weg for inv in erg.inverse_wege] == ["matrix", "rueckwaerts"]
    assert all(inv.inverse == 77 for inv in erg.inverse_wege)
    assert erg.inverse_wege[0].bezout.matrix.Q == ((100, 23), (13, 3))


def test_vl3_beispiel():
    # VL3: m = 101, e = 13, a = 8 → c = 18, d = 77.
    erg = rsa_berechnen(m=101, schluessel=13, w=8, klartext_name="a")
    assert erg.c == 18
    assert erg.d == 77
    assert erg.entschluesselung.ergebnis == 8


def test_vl3_bob():
    # VL3 Bob: p = 101, q = 113, e = 3533 → d = 6597; a = 9726 → c = 5761.
    erg = rsa_berechnen(p=101, q=113, schluessel=3533, w=9726, klartext_name="a")
    assert erg.m == 11413
    assert erg.phi_m == 11200
    assert erg.d == 6597
    assert erg.c == 5761
    assert erg.entschluesselung.ergebnis == 9726
    assert erg.falscher_weg is None


@pytest.mark.parametrize("variante, e, d, c", [
    ("oeffentlich", 113, 1463, 708),   # wie Ü4 A9 (wahrscheinlich gemeint)
    ("privat", 1463, 113, 265),        # Wortlaut GP A1 „d = 113“
])
def test_gedaechtnisprotokoll_a1(variante, e, d, c):
    erg = rsa_berechnen(m=2803, schluessel=113, w=715, variante=variante)
    assert (erg.e, erg.d, erg.c) == (e, d, c)
    assert erg.euler.ergebnis == (d if variante == "oeffentlich" else e)
    assert erg.entschluesselung.ergebnis == 715


def test_nur_entschluesseln():
    erg = rsa_berechnen(m=2803, schluessel=1463, c=708, variante="privat")
    assert erg.e == 113
    assert erg.w == 715


def test_ggt_nicht_eins():
    erg = rsa_berechnen(p=11, q=13, schluessel=12, w=5)
    assert erg.ggt_schluessel_phi == 12
    assert erg.d is None
    assert "kein passendes d" in erg.text


def test_symmetrische_reste_gleiches_ergebnis():
    erg = rsa_berechnen(m=2803, schluessel=113, w=715, rep="symmetrisch")
    assert (erg.c, erg.d) == (708, 1463)


def test_main_alle_varianten(capsys):
    assert main(["--m", "2803", "--schluessel", "113", "--w", "715", "--alle-varianten"]) == 0
    aus = capsys.readouterr().out
    for v in VARIANTEN:
        assert f"# Variante: {v}" in aus
    assert "c = 708" in aus and "c = 265" in aus


def test_main_klausur(capsys):
    assert main(["--p", "43", "--q", "57", "--e", "221", "--w", "1511"]) == 0
    aus = capsys.readouterr().out
    assert "ACHTUNG: 57 ist keine Primzahl" in aus
    assert aus.rstrip().endswith("d = 821")


def test_main_fehlende_eingabe(capsys):
    assert main(["--e", "3"]) == 1
    assert "HINWEIS" in capsys.readouterr().out
