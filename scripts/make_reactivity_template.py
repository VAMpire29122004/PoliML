"""Create data/reactivity_ratios.csv, pre-filled with the monomer pairs in your copolymer table.

    python scripts/make_reactivity_template.py            # refuses to overwrite an existing file
    python scripts/make_reactivity_template.py --force    # overwrite

You then fill r1, r2 and the source columns by hand from a cited reference.
Convention: r1 belongs to monomer_1 (r1 = k11/k12), r2 to monomer_2 (r2 = k22/k21).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reactivity_ratios.csv"
COLUMNS = ["pair_id", "monomer_1_id", "monomer_2_id", "r1", "r2", "temperature_c", "medium",
           "method", "source", "verified", "notes"]


def main() -> None:
    if OUT.exists() and "--force" not in sys.argv:
        sys.exit(f"{OUT.relative_to(ROOT)} already exists; use --force to overwrite")
    co = pd.read_csv(ROOT / "data" / "copolymers_chain_growth.csv", keep_default_na=False)
    pairs = co[(co.monomer_a_id != "") & (co.monomer_b_id != "")][["monomer_a_id", "monomer_b_id", "polymer_id"]]
    rows, seen = [], set()
    for _, r in pairs.iterrows():
        key = frozenset((r.monomer_a_id, r.monomer_b_id))
        if key in seen:
            continue
        seen.add(key)
        rows.append({"pair_id": f"{r.monomer_a_id}__{r.monomer_b_id}", "monomer_1_id": r.monomer_a_id,
                     "monomer_2_id": r.monomer_b_id, "r1": "", "r2": "", "temperature_c": "", "medium": "",
                     "method": "", "source": "", "verified": "False",
                     "notes": f"needed for {r.polymer_id}"})
    pd.DataFrame(rows, columns=COLUMNS).to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)} with {len(rows)} pairs to fill")


if __name__ == "__main__":
    main()
