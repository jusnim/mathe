"""Tests für skripte/turnierplan.py."""

import pytest

from skripte.turnierplan import (
    VARIANTEN, gegner_von_2m, loesung, main, paare, probe, turnierplan,
)

# Plan für 8 Mannschaften (Spalten = Runden 1 … 7).
PLAN_8 = {
    1: [(1, 7), (2, 6), (3, 5), (4, 8)],
    2: [(1, 8), (2, 7), (3, 6), (4, 5)],
    3: [(1, 2), (3, 7), (4, 6), (5, 8)],
    4: [(1, 3), (2, 8), (4, 7), (5, 6)],
    5: [(1, 4), (2, 3), (5, 7), (6, 8)],
    6: [(1, 5), (2, 4), (3, 8), (6, 7)],
    7: [(1, 6), (2, 5), (3, 4), (7, 8)],
}

# Beispiel 2m = 4
PLAN_4 = {1: [(1, 3), (2, 4)], 2: [(1, 4), (2, 3)], 3: [(1, 2), (3, 4)]}

# Beispiel 2m = 6
PLAN_6 = {
    1: [(1, 5), (2, 4), (3, 6)],
    2: [(1, 6), (2, 5), (3, 4)],
    3: [(1, 2), (3, 5), (4, 6)],
    4: [(1, 3), (2, 6), (4, 5)],
    5: [(1, 4), (2, 3), (5, 6)],
}


@pytest.mark.parametrize("n, erwartet", [(8, PLAN_8), (4, PLAN_4), (6, PLAN_6)])
def test_plan_standard(n, erwartet):
    assert paare(turnierplan(n)) == erwartet


def test_plan_summe_2m():
    # Die Regel „x + y = 2m“ ist ungenau. Der Plan ist trotzdem der Standard-Plan.
    assert paare(turnierplan(8, "summe-2m")) == PLAN_8
    assert "nur Runde r = 1" in loesung(8, "summe-2m", False, None)


def test_plan_rest_null():
    # Gleiche Runden, nur heißt Runde 7 jetzt Runde 0.
    p = paare(turnierplan(8, "rest-null"))
    assert p[0] == PLAN_8[7]
    assert {k: v for k, v in p.items() if k} == {k: v for k, v in PLAN_8.items() if k != 7}
    assert list(p) == [0, 1, 2, 3, 4, 5, 6]


# Gegner von Mannschaft 8 in Runde r = 1 … 7: 4, 1, 5, 2, 6, 3, 7.
GEGNER_8 = {1: 4, 2: 1, 3: 5, 4: 2, 5: 6, 6: 3, 7: 7}


@pytest.mark.parametrize("variante", VARIANTEN)
@pytest.mark.parametrize("r, x", sorted(GEGNER_8.items()))
def test_gegner_von_8(variante, r, x):
    assert gegner_von_2m(8, r, variante).mannschaft == x


def test_gegner_repraesentant():
    # 28 ≡ 7 (mod 7) bei standard, 28 ≡ 0 bei rest-null.
    assert gegner_von_2m(8, 7).x == 7
    assert gegner_von_2m(8, 7, "rest-null").x == 0
    assert gegner_von_2m(8, 0, "rest-null").mannschaft == 7


@pytest.mark.parametrize("n", [2, 4, 6, 8, 10, 12, 20])
@pytest.mark.parametrize("variante", VARIANTEN)
def test_probe_allgemein(n, variante):
    assert probe(turnierplan(n, variante)) == []


def test_ungerade_anzahl():
    plan = turnierplan(7)
    assert plan.n == 8 and plan.n_eingabe == 7
    assert "ungerade" in loesung(7, "standard", False, None)


def test_ausgabe_cli(capsys):
    main(["--mannschaften", "8", "--gegner", "--alle-varianten"])
    out = capsys.readouterr().out
    assert out.count("Ergebnis:") == 6
    assert "Runde 1: 8 gegen 4" in out


def test_ungueltige_runde(capsys):
    main(["--mannschaften", "8", "--runde", "9"])
    assert "gibt es nicht" in capsys.readouterr().out
