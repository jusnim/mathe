"""Tests für skripte/linearer_code.py mit Aufgaben aus VL4, VL6, Ü5 und Ü6."""

import pytest

from skripte.linearer_code import (
    Auftrag,
    EingabeFehler,
    baue_code,
    codiere,
    decodiere,
    hamming_H,
    klassenfuehrer,
    lies_matrix,
    lies_zeile,
    loesungsweg,
    main,
    mindestabstand,
    paritaet_G,
    standardarray,
    syndrom,
    wiederholung_G,
)


def M(text, q=2):
    return lies_matrix(text, q)[0]


def W(text, q=2):
    return lies_zeile(text, q)[0]


def dec(code, r):
    tab = klassenfuehrer(code)
    return decodiere(code, W(r, code.q), mindestabstand(code), tab)


# --- Ü5 5_1: Wiederholungscode in F_3^5 (Variante B) -----------------------

def test_ue5_1_wiederholung():
    code = baue_code(3, G=wiederholung_G(3, 5), variante="B")
    assert code.H == M("21000;20100;20010;20001", 3)
    assert codiere(code, [2]).c == [2, 2, 2, 2, 2]
    assert syndrom(code, W("22120", 3)).s == [0, 2, 0, 1]
    assert syndrom(code, W("11022", 3)).s == [0, 2, 1, 1]
    assert syndrom(code, W("11111", 3)).s == [0, 0, 0, 0]
    assert mindestabstand(code).d == 5


# --- Ü5 5_2: Paritätscode in F_3^5 (Variante B) -----------------------------

def test_ue5_2_paritaet():
    code = baue_code(3, G=paritaet_G(3, 5, "B"), variante="B")
    assert code.G == M("10001;01001;00101;00011", 3)
    assert code.H == [[2, 2, 2, 2, 1]]
    assert codiere(code, W("2120", 3)).c == W("21202", 3)
    assert syndrom(code, W("22120", 3)).s == [2]
    assert syndrom(code, W("11022", 3)).s == [1]
    assert syndrom(code, W("11111", 3)).s == [0]
    assert mindestabstand(code).d == 2


# --- Ü5 5_3: H_2(2) und H_2(4) ----------------------------------------------

def test_ue5_3_hamming_2_2():
    code = baue_code(2, H=hamming_H(2, 2), variante="A")
    assert code.H == M("101;011")
    assert code.G == [[1, 1, 1]]


def test_ue5_3_hamming_2_4():
    quelle = M("100000001111111;"
               "010001110001111;"
               "001010110110011;"
               "000111011010101")
    assert hamming_H(2, 4) == quelle


# --- Ü5 5_4: Decodieren in H_2(3) mit dem H der Aufgabe --------------------

def test_ue5_4_decodieren():
    code = baue_code(2, H=M("1001101;0101011;0010111"))
    d1, d2, d3 = dec(code, "1101100"), dec(code, "1111111"), dec(code, "1111000")
    assert d1.syn.s == [1, 0, 1] and d1.c == W("1101000")
    assert d1.e == W("0000100")               # x = 5
    assert d2.ist_codewort and d2.c == W("1111111")
    assert d3.syn.s == [0, 0, 1] and d3.c == W("1101000")
    assert d3.e == W("0010000")               # x = 3


# --- Ü5 5_5: H_3(3) ---------------------------------------------------------

def test_ue5_5_hamming_3_3():
    quelle = M("1000011111111;0101100111222;0011212012012", 3)
    assert hamming_H(3, 3) == quelle
    code = baue_code(3, H=quelle)
    assert (code.n, code.k, mindestabstand(code).d) == (13, 10, 3)
    d1 = dec(code, "0000000000112")
    assert d1.syn.s == [1, 2, 2]
    assert d1.e == [0] * 12 + [1]             # x = 13, y = 1
    assert d1.c == W("0000000000111", 3)
    d2 = dec(code, "0000000000122")
    assert d2.syn.s == [2, 1, 0]
    # Die Quelle schreibt x = 13. Richtig ist x = 11, denn (1 2 0)ᵀ ist die
    # 11. Spalte von H (PRIORISIERUNG.md, Abschnitt 5). c₂ stimmt mit der Quelle.
    assert d2.e == [0] * 10 + [2, 0, 0]       # x = 11, y = 2
    assert d2.c == W("0000000000222", 3)
    text = loesungsweg(Auftrag(q=3, H=quelle, woerter=[W("0000000000122", 3)]))
    assert "x = 11" in text and "x = 13" in text   # Hinweis auf den Quellenfehler


# --- Ü6 6_4 und VL4: H_3(2) -------------------------------------------------

def test_ue6_4_hamming_3_2():
    code = baue_code(3, H=hamming_H(3, 2), variante="A")
    assert code.H == M("1011;0112", 3)
    assert code.G == M("2210;2101", 3)
    assert (code.n, code.k, mindestabstand(code).d) == (4, 2, 3)


def test_vl4_hamming_3_2_decodieren():
    code = baue_code(3, H=M("1011;0112", 3))
    d1 = dec(code, "1111")
    assert d1.syn.s == [0, 1] and d1.e == W("0100", 3) and d1.c == W("1011", 3)
    d2 = dec(code, "0011")
    assert d2.syn.s == [2, 0] and d2.e == W("2000", 3) and d2.c == W("1011", 3)


def test_vl4_hamming_2_3_G():
    code = baue_code(2, H=M("1001011;0101110;0010111"), variante="A")
    assert code.G == M("1101000;0110100;1110010;1010001")
    abst = mindestabstand(code)
    assert abst.d == 3 and abst.min_zeilengewicht == 3


# --- Ü6 6_8: Klassenführer und Decodieren -----------------------------------

def test_ue6_8_klassenfuehrer():
    code = baue_code(3, H=M("1011;0112", 3))
    tab = klassenfuehrer(code)
    erwartet = [("0000", "00"), ("1000", "10"), ("0100", "01"), ("0010", "11"),
                ("0001", "12"), ("2000", "20"), ("0200", "02"), ("0020", "22"),
                ("0002", "21")]
    assert [(f.e, f.s) for f in tab] == [(W(e, 3), W(s, 3)) for e, s in erwartet]
    d1 = dec(code, "2200")
    assert d1.syn.s == [2, 2] and d1.e == W("0020", 3) and d1.c == W("2210", 3)
    d2 = dec(code, "0121")
    assert d2.syn.roh == [3, 5] and d2.syn.s == [0, 2]
    assert d2.e == W("0200", 3) and d2.c == W("0221", 3)


# --- VL6: Nebenklassen des [4,2,2]_2-Codes ---------------------------------

def test_vl6_nebenklassen():
    code = baue_code(2, G=M("1011;0101"))
    assert mindestabstand(code).d == 2
    tab = klassenfuehrer(code, "erste")
    assert [f.e for f in tab] == [W("0000"), W("1000"), W("0100"), W("0010")]
    arr = standardarray(code, tab)
    assert arr[0] == [W("0000"), W("1011"), W("0101"), W("1110")]
    assert arr[1] == [W("1000"), W("0011"), W("1101"), W("0110")]
    assert arr[2] == [W("0100"), W("1111"), W("0001"), W("1010")]
    assert arr[3] == [W("0010"), W("1001"), W("0111"), W("1100")]
    # Decodier-Tabelle VL6: y = 0001 → 0101 (Nachricht 01), y = 0100 → 0000 (00)
    assert dec(code, "0001").c == W("0101") and dec(code, "0001").a == [0, 1]
    assert dec(code, "0100").c == W("0000") and dec(code, "0100").a == [0, 0]


def test_vl6_fuehrer_letzte():
    code = baue_code(2, G=M("1011;0101"))
    tab = klassenfuehrer(code, "letzte")
    assert W("0001") in [f.e for f in tab]      # statt 0100


# --- Varianten A, B, C (mod-11-Code aus VL6, Klausur A5) --------------------

H11 = "1111111111;1 2 3 4 5 6 7 8 9 10"


def test_variante_C_mod11():
    code = baue_code(11, H=M(H11, 11), variante="C")
    assert code.umrechnung.matrix_sys == M("9 8 7 6 5 4 3 2 1 0;3 4 5 6 7 8 9 10 0 1", 11)
    assert code.G[0] == W("1 0 0 0 0 0 0 0 2 8", 11)
    assert code.G[7] == W("0 0 0 0 0 0 0 1 9 1", 11)
    assert codiere(code, W("11500005", 11)).c == W("1 1 5 0 0 0 0 5 4 6", 11)


@pytest.mark.parametrize("variante", ["A", "B", "C"])
def test_klausur_a5_alle_varianten(variante):
    code = baue_code(11, H=M(H11, 11), variante=variante)
    assert all(v == 0 for g in code.G for v in syndrom(code, g).s)
    r1 = W("1 1 5 0 0 0 0 5 1 1", 11)
    r2 = W("1 1 5 0 0 0 0 7 3 3", 11)
    abst = mindestabstand(code)
    assert abst.d == 3
    tab = klassenfuehrer(code)
    d1 = decodiere(code, r1, abst, tab)
    assert d1.syn.s == [3, 0] and not d1.sicher      # nicht decodierbar
    d2 = decodiere(code, r2, abst, tab)
    assert d2.syn.s == [9, 10] and d2.sicher
    assert d2.c == W("1 1 5 0 0 2 0 7 3 3", 11)        # x = 6, y = 9


@pytest.mark.parametrize("variante", ["A", "B", "C"])
def test_varianten_hamming_3_2(variante):
    code = baue_code(3, H=hamming_H(3, 2, variante), variante=variante)
    assert mindestabstand(code).d == 3
    for g in code.G:
        assert syndrom(code, g).s == [0, 0]


def test_variante_A_wiederholung():
    code = baue_code(3, G=wiederholung_G(3, 5), variante="A")
    assert code.H == M("10002;01002;00102;00012", 3)


# --- Eingabeprüfung und Kommandozeile ---------------------------------------

def test_q_keine_primzahl():
    with pytest.raises(EingabeFehler):
        baue_code(4, H=[[1, 0, 1]])


def test_cli(capsys):
    assert main(["--q", "3", "--typ", "hamming", "--m", "2", "--tabelle",
                 "--r", "2200", "--r", "0121"]) == 0
    out = capsys.readouterr().out
    assert "Ergebnis: r1 = 2200 wird zu c = 2210 decodiert." in out
    assert "Ergebnis: r2 = 0121 wird zu c = 0221 decodiert." in out
    assert main(["--q", "4", "--H", "101"]) == 1
    assert "keine Primzahl" in capsys.readouterr().out


def test_cli_alle_varianten(capsys):
    assert main(["--q", "3", "--typ", "paritaet", "--n", "5", "--a", "2120",
                 "--alle-varianten"]) == 0
    out = capsys.readouterr().out
    assert out.count("Annahme (Variante") == 3
    assert "c = 21202" in out
