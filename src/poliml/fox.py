"""Fox equation for the glass transition of random copolymers.

1/Tg = w_a / Tg_a + (1 - w_a) / Tg_b   (temperatures in kelvin, w_a a weight fraction)
"""
from __future__ import annotations

K = 273.15


def fox_tg(w_a: float, tg_a_c: float, tg_b_c: float) -> float:
    """Predicted copolymer Tg (deg C) at weight fraction w_a of component A."""
    if not 0.0 <= w_a <= 1.0:
        raise ValueError("w_a must be between 0 and 1")
    ta, tb = tg_a_c + K, tg_b_c + K
    return 1.0 / (w_a / ta + (1.0 - w_a) / tb) - K


def implied_weight_fraction(tg_c: float, tg_a_c: float, tg_b_c: float) -> float | None:
    """Weight fraction of A that would give tg_c under Fox (None if Tg_a == Tg_b).

    Values outside [0, 1] mean the reported Tg cannot come from Fox for these homopolymers.
    """
    ta, tb, t = tg_a_c + K, tg_b_c + K, tg_c + K
    denom = 1.0 / ta - 1.0 / tb
    if abs(denom) < 1e-12:
        return None
    return (1.0 / t - 1.0 / tb) / denom
