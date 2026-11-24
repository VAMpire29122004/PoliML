import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from poliml.copolymer import (azeotrope_f1, classify, mayo_lewis_F1, mole_to_weight_fraction,  # noqa: E402
                              predict_tg_from_feed)


def test_ideal_random_follows_feed():
    for f in (0.1, 0.5, 0.9):
        assert math.isclose(mayo_lewis_F1(f, 1.0, 1.0), f)


def test_perfect_alternation_is_half():
    for f in (0.1, 0.5, 0.9):
        assert math.isclose(mayo_lewis_F1(f, 0.0, 0.0), 0.5)


def test_endpoints():
    assert mayo_lewis_F1(0.0, 0.3, 0.7) == 0.0
    assert mayo_lewis_F1(1.0, 0.3, 0.7) == 1.0


def test_azeotrope_composition_is_a_fixed_point():
    r1, r2 = 0.4, 0.6
    az = azeotrope_f1(r1, r2)
    assert math.isclose(az, 0.4)
    assert math.isclose(mayo_lewis_F1(az, r1, r2), az)


def test_no_azeotrope_when_one_above_one():
    assert azeotrope_f1(2.0, 0.5) is None


def test_illustrative_textbook_like_pair():
    # approximate, illustrative values only (not repository data)
    az = azeotrope_f1(0.52, 0.46)
    assert 0.4 < az < 0.6


def test_classify():
    assert classify(1.0, 1.0)[0] == "ideal"
    assert classify(0.05, 0.05)[0] == "alternating"
    assert classify(3.0, 2.0)[0] == "blocky"
    assert classify(2.0, 0.3)[0] == "prefers_1"
    assert classify(0.3, 2.0)[0] == "prefers_2"


def test_mole_to_weight():
    assert math.isclose(mole_to_weight_fraction(0.5, 100.0, 200.0), 1 / 3)


def test_predict_tg_between_homopolymers_and_equal_tg():
    out = predict_tg_from_feed(0.5, 0.5, 0.5, 100.0, 100.0, 100.0, -50.0)
    assert -50.0 < out["tg_c"] < 100.0
    same = predict_tg_from_feed(0.3, 0.5, 2.0, 104.0, 100.0, 80.0, 80.0)
    assert math.isclose(same["tg_c"], 80.0)


def test_invalid_inputs():
    for bad in (lambda: mayo_lewis_F1(1.5, 1, 1), lambda: mayo_lewis_F1(0.5, -1, 1)):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("expected ValueError")
