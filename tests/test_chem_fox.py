import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from poliml.chem import formula_and_mw  # noqa: E402
from poliml.fox import fox_tg, implied_weight_fraction  # noqa: E402


def test_formulas():
    assert formula_and_mw("C=CC1=CC=CC=C1")[0] == "C8H8"
    assert formula_and_mw("FC(F)=C(F)F")[0] == "C2F4"
    assert formula_and_mw("C=C(C)C(=O)OCCC[Si](OC)(OC)OC")[0] == "C10H20O5Si"
    assert formula_and_mw("OC(=O)C1=CC(=CC=C1)C(=O)O")[0] == "C8H6O4"
    assert formula_and_mw("O=C1OC(=O)C2=CC=CC=C12")[0] == "C8H4O3"


def test_mw():
    assert math.isclose(formula_and_mw("C=CC1=CC=CC=C1")[1], 104.15, abs_tol=0.05)
    assert math.isclose(formula_and_mw("CC(C)(C1=CC=C(O)C=C1)C1=CC=C(O)C=C1")[1], 228.29, abs_tol=0.05)


def test_fox_limits_and_roundtrip():
    assert math.isclose(fox_tg(1.0, 100, -50), 100)
    assert math.isclose(fox_tg(0.0, 100, -50), -50)
    assert math.isclose(fox_tg(0.5, 80, 80), 80)
    tg = fox_tg(0.3, 105, -54)
    assert math.isclose(implied_weight_fraction(tg, 105, -54), 0.3)
    assert implied_weight_fraction(10, 50, 50) is None
    assert fox_tg(0.5, 100, -50) < (100 + -50) / 2 + 20  # sanity: between the homopolymer Tgs
