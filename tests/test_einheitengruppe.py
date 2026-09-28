"""Tests für skripte/einheitengruppe.py."""

import pytest

from skripte.einheitengruppe import (
    VARIANTEN,
    einheiten,
    gruppentafel,
    inversen,
    loesungsweg,
    main,
    ordne,
    ordnungen,
)
from skripte.euklid import KeinInversesFehler


def _einheiten(m, variante):
    return ordne(einheiten(m).einheiten, m, variante, "paare")


# --- ℤ_12* und ℤ_18* ---------------------------------------------------------

@pytest.mark.parametrize("m, symm, std, null", [
    (12, [1, -1, 5, -5], [1, 5, 7, 11], [2, 3, 4, 6, 8, 9, 10]),
    (18, [1, -1, 5, -5, 7, -7], [1, 5, 7, 11, 13, 17], [2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 16]),
    (6, [1, -1], [1, 5], [2, 3, 4]),          # ℤ_6*: nur 1, 5
])
def test_einheiten_und_nullteiler(m, symm, std, null):
    assert _einheiten(m, "symmetrisch") == symm
    assert _einheiten(m, "standard") == std
    e = einheiten(m)
    assert [a for a, _ in e.nullteiler] == null
    assert all(a * b % m == 0 and b % m != 0 for a, b in e.nullteiler)
    assert len(e.einheiten) == e.phi_m


def test_inverse_z12_z18():
    # Symmetrisch: alle selbstinvers in ℤ_12*
    assert {(z.a, z.inv) for z in inversen(12)} == {(1, 1), (-1, -1), (5, 5), (-5, -5)}
    # ℤ_18*: 5·(−7) ≡ 1, (−5)·7 ≡ 1
    inv18 = {z.a: z.inv for z in inversen(18)}
    assert inv18 == {1: 1, -1: -1, 5: -7, -5: 7, 7: -5, -7: 5}
    # Standard: 7·13 = 91 ≡ 1 (mod 18)
    inv18s = {z.a: z.inv for z in inversen(18, "standard")}
    assert inv18s[7] == 13 and inv18s[13] == 7


def test_gruppentafel_z12():
    t = gruppentafel(12)
    assert t.elemente == [1, -1, 5, -5]
    assert t.werte == [[1, -1, 5, -5], [-1, 1, -5, 5], [5, -5, 1, -1], [-5, 5, -1, 1]]


def test_gruppentafel_z18():
    t = gruppentafel(18)
    assert t.elemente == [1, -1, 5, -5, 7, -7]
    assert t.werte == [
        [1, -1, 5, -5, 7, -7],
        [-1, 1, -5, 5, -7, 7],
        [5, -5, 7, -7, -1, 1],
        [-5, 5, -7, 7, 1, -1],
        [7, -7, -1, 1, -5, 5],
        [-7, 7, 1, -1, 5, -5],
    ]
    assert "Gruppentafel von (ℤ_18*, ·)" in loesungsweg(18, "tafel")


def test_tafel_z6_mit_nullteilern():
    t = gruppentafel(6, "standard", mit_nullteilern=True)
    assert t.elemente == [1, 2, 3, 4, 5]
    assert t.werte == [[1, 2, 3, 4, 5], [2, 4, 0, 2, 4], [3, 0, 3, 0, 3],
                       [4, 2, 0, 4, 2], [5, 4, 3, 2, 1]]


# --- Ordnungen in ℤ_17* -----------------------------------------------------

@pytest.mark.parametrize("variante, inv", [
    ("symmetrisch", [-8, 6, -4, 7]),
    ("standard", [9, 6, 13, 7]),            # Probe: 2·9, 3·6, 4·13, 5·7
])
def test_ordnungen_z17(variante, inv):
    o = ordnungen(17, variante, elemente=[2, 3, 4, 5])
    assert [z.ordnung for z in o.zeilen] == [8, 16, 4, 16]
    assert [z.inv for z in o.zeilen] == inv
    assert [z.a for z in o.zeilen if z.primitiv] == [3, 5]
    assert len(o.primitive) == 8


def test_potenzlisten_z17():
    o = ordnungen(17, elemente=[2, 3, 4, 5])
    assert o.zeilen[0].potenzen == [2, 4, 8, -1, -2, -4, -8, 1]
    assert o.zeilen[1].potenzen[:8] == [3, -8, -7, -4, 5, -2, -6, -1]
    assert o.zeilen[2].potenzen == [4, -1, -4, 1]
    assert o.zeilen[3].potenzen[:8] == [5, 8, 6, -4, -3, 2, -7, -1]


# --- ℤ_7* -------------------------------------------------------------------

@pytest.mark.parametrize("variante", VARIANTEN)
def test_ordnungen_z7(variante):
    o = ordnungen(7, variante, elemente=[2, 3])
    assert [z.ordnung for z in o.zeilen] == [3, 6]
    assert o.zeilen[1].primitiv and not o.zeilen[0].primitiv
    assert o.primitive == [3, 5]
    if variante == "symmetrisch":
        assert o.zeilen[1].potenzen == [3, 2, -1, -3, -2, 1]
    else:
        assert o.zeilen[0].potenzen == [2, 4, 1]


# --- Inverse modulo 11 ------------------------------------------------------

def test_inverse_mod11():
    zeilen = inversen(11, "symmetrisch", "zahlen")
    assert [z.a for z in zeilen] == [1, 2, 3, 4, 5, -5, -4, -3, -2, -1]
    assert [z.inv for z in zeilen] == [1, -5, 4, 3, -2, 2, -3, -4, 5, -1]
    std = inversen(11, "standard", "zahlen")
    assert [z.inv for z in std] == [1, 6, 4, 3, 9, 2, 8, 7, 5, 10]


# --- Eingabeprüfung und Kommandozeile ---------------------------------------

def test_keine_einheit():
    with pytest.raises(KeinInversesFehler):
        inversen(12, elemente=[4])
    o = ordnungen(12, elemente=[4])
    assert not o.zeilen[0].einheit
    assert "keine Einheit" in loesungsweg(12, "ordnungen", elemente=[4, 5])


def test_modul_zu_klein(capsys):
    assert main(["--m", "1"]) == 1
    assert "HINWEIS" in capsys.readouterr().out
    with pytest.raises(ValueError):
        einheiten(1)


def test_primzahl_ist_koerper():
    assert "Körper" in loesungsweg(11, "einheiten")


@pytest.mark.parametrize("variante", VARIANTEN)
def test_main_alle_modi(capsys, variante):
    assert main(["--m", "18", "--variante", variante, "--inverse-weg", "7"]) == 0
    out = capsys.readouterr().out
    assert f"Variante: {variante}" in out and "Ergebnis:" in out and "Annahme:" in out


def test_main_alle_varianten(capsys):
    assert main(["--m", "12", "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert "Variante: symmetrisch" in out and "Variante: standard" in out
