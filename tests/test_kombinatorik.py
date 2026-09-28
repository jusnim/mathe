"""Tests für skripte/kombinatorik.py mit Aufgaben aus VL2, Ü3 und Ü4."""

import math

import pytest

from skripte import kombinatorik as k


# --- Ü3 A4: 4×4-Gitter, 70 Wege; C(8,4) = 7·2·5 = 70 -----------------------

def test_ue3_a4_gitter_4x4():
    d = k.berechne_gitter(4, 4)
    assert d.wert == 70
    assert (d.oben, d.rechts) == (4, 4)
    assert d.binom.gekuerzt == [7, 2, 5]  # wie in der Musterlösung


def test_ue3_a4_beispielweg_und_liste():
    d = k.berechne_gitter(4, 4, liste=True, beispiel="rorooorr")
    assert d.beispiel == "r o r o o o r r"
    assert len(d.wege) == 70 and len(set(d.wege)) == 70


def test_ue3_a4_variante_punkte():
    # Variante: 4×4 Gitterpunkte = 3×3 Kästchen → C(6,3) = 20
    assert k.berechne_gitter(4, 4, "punkte").wert == 20
    assert k.berechne_gitter(5, 5, "punkte").wert == 70


# --- Ü3 A5: m×n-Gitter, C(m+n, n) = C(m+n, m) ------------------------------

@pytest.mark.parametrize("m,n", [(1, 1), (2, 3), (3, 5), (6, 2), (0, 4)])
def test_ue3_a5_gitter_allgemein(m, n):
    assert k.berechne_gitter(m, n).wert == math.comb(m + n, n) == math.comb(m + n, m)
    assert k.berechne_gitter(m + 1, n + 1, "punkte").wert == math.comb(m + n, n)


# --- Ü4 A2: d₄ = 9, P(4) = 3/8, Liste --------------------------------------

def test_ue4_a2_derangements_s4():
    d = k.berechne_derangement(4, liste=True)
    assert d.wert == 9
    assert d.wahrscheinlichkeit == k.Fraction(3, 8)
    loesung = {"2341", "2413", "2143", "3142", "3421", "3412", "4123", "4312", "4321"}
    assert {"".join(map(str, p)) for p in d.liste} == loesung
    text = k.text_derangement(d)
    assert "4!(1/2 − 1/6 + 1/24) = 4!·9/24 = 9" in text
    assert "Ergebnis: d_4 = 9" in text


@pytest.mark.parametrize("n,dn", [(0, 1), (1, 0), (2, 1), (3, 2), (5, 44), (6, 265)])
def test_derangement_werte(n, dn):
    assert k.derangement_zahl(n) == dn


# --- VL2 2.6: p(5) = 7, ungerade = verschieden = 3, Ferrers ----------------

def test_vl2_partitionen_von_5():
    d = k.berechne_partition(5)
    assert d.wert == 7
    assert [k.partition_text(p) for p in d.liste] == [
        "5", "4+1", "3+2", "3+1+1", "2+2+1", "2+1+1+1", "1+1+1+1+1"]


def test_vl2_partitionen_mit_einschraenkung():
    assert k.berechne_partition(5, "ungerade").wert == 3
    assert k.berechne_partition(5, "verschieden").wert == 3
    d = k.berechne_partition(12, "max-teil", 4)
    assert d.probe == ("max-anzahl", d.wert)


def test_vl2_ferrers_27():
    d = k.berechne_ferrers([8, 6, 6, 3, 2, 1, 1])
    assert d.konjugiert == [7, 5, 4, 3, 3, 3, 1, 1]


# --- VL2 2.2: Kombinationen mit Wiederholung, Bijektion --------------------

@pytest.mark.parametrize("auswahl,plaetze", [
    ([1, 1, 1, 2, 4, 4, 5], [4, 6, 7, 10]),
    ([2, 4, 4, 5, 5, 5, 5], [1, 3, 4, 7]),
    ([2, 4, 4, 4, 4, 4, 4], [1, 3, 4, 11]),
])
def test_vl2_bijektion_wiederholung(auswahl, plaetze):
    d = k.berechne_kombination(5, 7, True, auswahl)
    assert d.bijektion["A2"] == plaetze
    assert d.wert == math.comb(11, 4) == 330


# --- Ü3 A6: alternierende Summe = 1; Folgerungen a–c; Ü3 A7 Vandermonde -----

@pytest.mark.parametrize("kk", range(1, 12))
def test_ue3_a6_alternierende_summe(kk):
    assert k.berechne_summe("alternierend", kk).wert == 1


@pytest.mark.parametrize("n", range(1, 10))
def test_vl2_folgerungen(n):
    assert k.berechne_summe("alle", n).wert == 2 ** n
    assert k.berechne_summe("wechselnd", n).wert == 0
    d = k.berechne_summe("quadrate", n)
    assert d.wert == d.formelwert == math.comb(2 * n, n)


@pytest.mark.parametrize("m,n,r", [(3, 4, 2), (5, 5, 5), (6, 3, 3), (7, 9, 4)])
def test_ue3_a7_vandermonde(m, n, r):
    d = k.berechne_summe("vandermonde", n, m, r)
    assert d.wert == d.formelwert == math.comb(m + n, r)


# --- Grundformeln und Eingabeprüfung --------------------------------------

def test_grundformeln():
    assert k.berechne_geordnet(5, 3).wert == 60
    assert k.berechne_geordnet(5, 3, True).wert == 125
    assert k.berechne_permutation(4).wert == 24
    assert k.berechne_kombination(8, 4).wert == 70
    assert k.berechne_binom(10, 7).k_rechnung == 3
    assert k.berechne_lehrsatz(3, 2, -1).wert == 1
    assert k.berechne_pascal(4).zeilen[-1] == [1, 4, 6, 4, 1]


def test_eingabepruefung():
    assert k.berechne_binom(3, 5).wert == 0 and k.berechne_binom(3, 5).hinweise
    assert k.berechne_geordnet(2, 3).hinweise
    assert k.berechne_gitter(-1, 2).hinweise
    assert k.berechne_derangement(-3).hinweise
    assert k.berechne_kombination(5, 3, True, [1, 6, 2]).hinweise


def test_cli_laeuft(capsys):
    assert k.main(["gitter", "4", "4", "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert "Ergebnis: 70 Wege" in out and "Ergebnis: 20 Wege" in out
    k.main(["partition", "5"])
    assert "Ergebnis: p(5) = 7" in capsys.readouterr().out
