"""Check the chain-growth copolymer Tg values against the Fox equation.

    python scripts/fox_validation.py

For each copolymer whose two homopolymers both have a Tg in data/homopolymers.csv:
  * bounds check : a random copolymer's Tg must lie between the two homopolymer Tgs;
  * implied w_A  : the weight fraction of A that Fox would need to reproduce the
                   reported Tg (outside 0-1 = cannot be explained by Fox);
  * if `weight_fraction_a` is filled in data/copolymers_chain_growth.csv, the Fox
    prediction and its error are also reported.
Writes data/fox_validation.csv and reports/figures/fox_bounds.png.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from poliml.fox import fox_tg, implied_weight_fraction  # noqa: E402


def main() -> None:
    homo = pd.read_csv(ROOT / "data" / "homopolymers.csv", keep_default_na=False, na_values=[""])
    co = pd.read_csv(ROOT / "data" / "copolymers_chain_growth.csv", keep_default_na=False, na_values=[""])
    tg_by_monomer = homo.dropna(subset=["tg_c"]).set_index("monomer_id")["tg_c"].to_dict()

    rows = []
    for _, r in co.iterrows():
        a, b = r["monomer_a_id"], r["monomer_b_id"]
        rec = {"copolymer_id": r["polymer_id"], "name": r["name"], "monomer_a_id": a, "monomer_b_id": b,
               "tg_a_c": tg_by_monomer.get(a, np.nan), "tg_b_c": tg_by_monomer.get(b, np.nan),
               "tg_reported_c": r["tg_c"], "weight_fraction_a": r.get("weight_fraction_a", np.nan)}
        missing = [m for m, v in ((a, rec["tg_a_c"]), (b, rec["tg_b_c"])) if pd.isna(v)]
        if missing or pd.isna(rec["tg_reported_c"]):
            rec["status"] = "not_testable"
            rec["comment"] = "no homopolymer Tg for: " + ", ".join(m or "(unlinked)" for m in missing)
            rows.append(rec)
            continue
        lo, hi = sorted((rec["tg_a_c"], rec["tg_b_c"]))
        rec.update(tg_lower_c=lo, tg_upper_c=hi)
        inside = lo <= rec["tg_reported_c"] <= hi
        rec["within_homopolymer_bounds"] = inside
        rec["deviation_outside_bounds_c"] = 0.0 if inside else (
            rec["tg_reported_c"] - hi if rec["tg_reported_c"] > hi else rec["tg_reported_c"] - lo)
        w = implied_weight_fraction(rec["tg_reported_c"], rec["tg_a_c"], rec["tg_b_c"])
        rec["implied_weight_fraction_a"] = w
        if pd.notna(rec["weight_fraction_a"]):
            rec["fox_tg_c"] = fox_tg(rec["weight_fraction_a"], rec["tg_a_c"], rec["tg_b_c"])
            rec["fox_error_c"] = rec["fox_tg_c"] - rec["tg_reported_c"]
        rec["status"] = "consistent" if inside else "inconsistent"
        rec["comment"] = "" if inside else "reported Tg outside the range of the homopolymer Tgs"
        rows.append(rec)

    out = pd.DataFrame(rows)
    cols = ["copolymer_id", "name", "monomer_a_id", "monomer_b_id", "tg_a_c", "tg_b_c", "tg_reported_c",
            "tg_lower_c", "tg_upper_c", "within_homopolymer_bounds", "deviation_outside_bounds_c",
            "implied_weight_fraction_a", "weight_fraction_a", "fox_tg_c", "fox_error_c", "status", "comment"]
    for c in cols:
        if c not in out:
            out[c] = np.nan
    out = out[cols]
    out.to_csv(ROOT / "data" / "fox_validation.csv", index=False, float_format="%.4g")

    t = out[out.status != "not_testable"]
    print(f"{len(co)} copolymers: {len(t)} testable, {len(out) - len(t)} not testable")
    print(f"consistent with homopolymer bounds: {(t.status == 'consistent').sum()} / {len(t)}")
    print(t[["copolymer_id", "tg_a_c", "tg_b_c", "tg_reported_c", "implied_weight_fraction_a", "status"]]
          .to_string(index=False))
    plot(t)


def plot(t: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = t.sort_values("tg_reported_c").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(t) + 1.6))
    for i, r in t.iterrows():
        ok = r["status"] == "consistent"
        ax.plot([r["tg_lower_c"], r["tg_upper_c"]], [i, i], color="0.7", lw=6, solid_capstyle="round", zorder=1)
        ax.scatter(r["tg_reported_c"], i, s=45, color="tab:green" if ok else "tab:red", zorder=2)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels(t["copolymer_id"].str.upper())
    ax.set_xlabel("Tg (°C)")
    ax.set_title("Reported copolymer Tg vs. range of its homopolymer Tgs\n"
                 "(grey bar = Fox bounds; red = outside the bounds)", fontsize=10)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    path = ROOT / "reports" / "figures" / "fox_bounds.png"
    fig.savefig(path, dpi=160)
    print("wrote", path.relative_to(ROOT))


if __name__ == "__main__":
    main()
