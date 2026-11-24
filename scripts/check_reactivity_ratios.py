"""Validate data/reactivity_ratios.csv and summarize what it implies.

    python scripts/check_reactivity_ratios.py

Fails (exit 1) on structural errors. Missing values are reported, not errors.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from poliml.copolymer import azeotrope_f1, classify  # noqa: E402


def main() -> None:
    rr = pd.read_csv(ROOT / "data" / "reactivity_ratios.csv", keep_default_na=False)
    mono = set(pd.read_csv(ROOT / "data" / "monomers_chain_growth.csv")["monomer_id"])
    errors = []
    if not rr["pair_id"].is_unique:
        errors.append("duplicate pair_id")
    for _, r in rr.iterrows():
        for col in ("monomer_1_id", "monomer_2_id"):
            if r[col] not in mono:
                errors.append(f"{r.pair_id}: unknown {col} '{r[col]}'")
    filled = rr[(rr.r1 != "") & (rr.r2 != "")].copy()
    for _, r in filled.iterrows():
        try:
            r1, r2 = float(r.r1), float(r.r2)
        except ValueError:
            errors.append(f"{r.pair_id}: r1/r2 not numeric")
            continue
        if r1 < 0 or r2 < 0:
            errors.append(f"{r.pair_id}: negative reactivity ratio")
        if not r.source:
            errors.append(f"{r.pair_id}: r values given without a source")
    print(f"{len(filled)} of {len(rr)} pairs filled")
    for _, r in filled.iterrows():
        try:
            r1, r2 = float(r.r1), float(r.r2)
            code, text = classify(r1, r2)
            print(f"  {r.pair_id}: r1={r1:g}, r2={r2:g} -> {code}: {text}")
        except ValueError:
            pass
    if errors:
        print("\nERRORS:")
        for e in errors:
            print("  -", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
