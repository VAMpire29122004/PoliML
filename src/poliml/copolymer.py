"""Terminal-model (Mayo-Lewis) copolymerization helpers.

Conventions
-----------
* Monomer 1 / monomer 2 and r1 / r2 follow the row order of data/reactivity_ratios.csv:
  r1 = k11 / k12, r2 = k22 / k21.
* f1 is the *mole* fraction of monomer 1 in the feed, F1 the mole fraction of
  monomer 1 in the polymer formed at that instant (low conversion).
"""
from __future__ import annotations

import numpy as np

from .fox import fox_tg


def _check_fraction(x: float, name: str) -> None:
    if not 0.0 <= x <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1, got {x}")


def _check_ratios(r1: float, r2: float) -> None:
    if r1 < 0 or r2 < 0:
        raise ValueError("reactivity ratios cannot be negative")


def mayo_lewis_F1(f1: float, r1: float, r2: float) -> float:
    """Instantaneous copolymer composition F1 for feed composition f1 (mole fractions)."""
    _check_fraction(f1, "f1")
    _check_ratios(r1, r2)
    if f1 in (0.0, 1.0):
        return float(f1)
    f2 = 1.0 - f1
    num = r1 * f1**2 + f1 * f2
    den = r1 * f1**2 + 2.0 * f1 * f2 + r2 * f2**2
    return num / den


def azeotrope_f1(r1: float, r2: float) -> float | None:
    """Feed composition where F1 == f1, or None if the system has no azeotrope.

    An azeotrope exists only when r1 and r2 are both below 1 or both above 1.
    """
    _check_ratios(r1, r2)
    if (r1 < 1 and r2 < 1) or (r1 > 1 and r2 > 1):
        return (1.0 - r2) / (2.0 - r1 - r2)
    return None


def classify(r1: float, r2: float, tol: float = 0.15) -> tuple[str, str]:
    """Return (code, plain-language description) for a pair of reactivity ratios.

    This describes tendencies of the terminal model; it is not a yes/no answer
    to "will these monomers copolymerize" (almost any pair does to some extent).
    """
    _check_ratios(r1, r2)
    if abs(r1 * r2 - 1.0) <= tol:
        if abs(r1 - 1.0) <= 0.25:
            return "ideal", "ideal/random: polymer composition tracks the feed (F1 ~ f1)"
        pref = 1 if r1 > 1 else 2
        return "ideal", f"ideal but unequal: monomer {pref} is consistently preferred at every feed"
    if r1 < 1 and r2 < 1:
        az = azeotrope_f1(r1, r2)
        return "alternating", f"alternating tendency; azeotrope at f1 = {az:.2f}"
    if r1 > 1 and r2 > 1:
        az = azeotrope_f1(r1, r2)
        return "blocky", f"tendency to blocks/homopolymer; unstable azeotrope at f1 = {az:.2f}"
    if r1 >= 1 >= r2:
        return "prefers_1", "monomer 1 is incorporated faster; polymer is richer in 1 than the feed, so the feed drifts toward 2"
    return "prefers_2", "monomer 2 is incorporated faster; polymer is richer in 2 than the feed, so the feed drifts toward 1"


def composition_curve(r1: float, r2: float, n: int = 101) -> tuple[np.ndarray, np.ndarray]:
    """(f1, F1) arrays for plotting the composition diagram."""
    f1 = np.linspace(0.0, 1.0, n)
    return f1, np.array([mayo_lewis_F1(x, r1, r2) for x in f1])


def mole_to_weight_fraction(F1: float, mw1: float, mw2: float) -> float:
    """Weight fraction of monomer 1 units from their mole fraction (repeat-unit MW = monomer MW)."""
    _check_fraction(F1, "F1")
    return F1 * mw1 / (F1 * mw1 + (1.0 - F1) * mw2)


def predict_tg_from_feed(f1: float, r1: float, r2: float, mw1: float, mw2: float,
                         tg1_c: float, tg2_c: float) -> dict:
    """Feed composition -> polymer composition -> Fox Tg (deg C), at low conversion."""
    F1 = mayo_lewis_F1(f1, r1, r2)
    w1 = mole_to_weight_fraction(F1, mw1, mw2)
    return {"F1_mole": F1, "w1": w1, "tg_c": fox_tg(w1, tg1_c, tg2_c)}
