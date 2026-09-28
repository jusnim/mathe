# Tests für die Web-Oberfläche: Formulare lassen sich für jedes Skript bauen, Rechnen klappt.
from skripte import web


def test_jedes_skript_hat_ein_schema():
    for s in web.alle_schemata():
        assert "fehler" not in s, s["name"]
        assert s["felder"] or s["unterbefehle"], s["name"]


def test_ausfuehren_crt():
    r = web.ausfuehren("crt", ["-k", "5", "13", "-k", "4", "15"])
    assert r["code"] == 0
    assert "109" in r["stdout"]
    assert r["befehl"] == "python -m skripte.crt -k 5 13 -k 4 15"


def test_unbekanntes_skript_wird_abgelehnt():
    r = web.ausfuehren("../../etc", [])
    assert r["code"] != 0
    assert "Unbekanntes Werkzeug" in r["stderr"]


def test_jedes_werkzeug_hat_eine_gruppe():
    for s in web.alle_schemata():
        assert s["gruppe"], s["name"]
