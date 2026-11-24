"""Load the cleaned PoliML tables and answer the lookups the app needs (no Streamlit here)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

TABLES = {
    "monomers": "monomers_chain_growth.csv",
    "monomers_step": "monomers_step_growth.csv",
    "homopolymers": "homopolymers.csv",
    "copolymers": "copolymers_chain_growth.csv",
    "copolymers_step": "copolymers_step_growth.csv",
    "fox": "fox_validation.csv",
    "reactivity": "reactivity_ratios.csv",
}


def default_data_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data"


def load_tables(data_dir: Path | None = None) -> dict[str, pd.DataFrame]:
    data_dir = Path(data_dir) if data_dir else default_data_dir()
    out = {}
    for key, fname in TABLES.items():
        path = data_dir / fname
        if path.exists():
            out[key] = pd.read_csv(path, keep_default_na=False, na_values=[""])
    return out


def monomer_mw(monomers: pd.DataFrame) -> dict[str, float]:
    """Molecular weight per monomer_id: from the structure when available, else the listed value."""
    mw = monomers["mw_from_smiles"].where(monomers["mw_from_smiles"].notna(), monomers["mw_listed"])
    return {i: float(v) for i, v in zip(monomers["monomer_id"], mw) if pd.notna(v)}


def homopolymer_tg(homopolymers: pd.DataFrame) -> dict[str, float]:
    h = homopolymers[(homopolymers["monomer_id"].notna()) & (homopolymers["tg_c"].notna())]
    return dict(zip(h["monomer_id"], h["tg_c"].astype(float)))


def monomers_for_copolymer_tool(monomers: pd.DataFrame, homopolymers: pd.DataFrame) -> pd.DataFrame:
    """Monomers that have both a molecular weight and a homopolymer Tg."""
    mw, tg = monomer_mw(monomers), homopolymer_tg(homopolymers)
    ids = [i for i in monomers["monomer_id"] if i in mw and i in tg]
    return monomers[monomers["monomer_id"].isin(ids)].sort_values("name").reset_index(drop=True)


def lookup_reactivity(reactivity: pd.DataFrame | None, id1: str, id2: str):
    """Return (r1, r2, row) for id1/id2 in that order, or None if no filled row exists.

    A row stored in the opposite order is swapped (r1 <-> r2).
    """
    if reactivity is None or reactivity.empty:
        return None
    for _, row in reactivity.iterrows():
        try:
            a, b = float(row["r1"]), float(row["r2"])
        except (TypeError, ValueError):
            continue
        if pd.isna(a) or pd.isna(b):
            continue
        if row["monomer_1_id"] == id1 and row["monomer_2_id"] == id2:
            return a, b, row
        if row["monomer_1_id"] == id2 and row["monomer_2_id"] == id1:
            return b, a, row
    return None


def all_flags(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """One long table of every review flag across the cleaned tables."""
    rows = []
    spec = {"monomers": ("monomer_id", "monomers"), "monomers_step": ("monomer_id", "step-growth monomers"),
            "homopolymers": ("polymer_id", "homopolymers"), "copolymers": ("polymer_id", "copolymers"),
            "copolymers_step": ("polymer_id", "step-growth copolymers")}
    for key, (idcol, label) in spec.items():
        df = tables.get(key)
        if df is None or "review_flags" not in df:
            continue
        for _, r in df[df["review_flags"].notna()].iterrows():
            for flag in str(r["review_flags"]).split(" | "):
                if flag.strip():
                    rows.append({"table": label, "id": r[idcol], "name": r["name"], "flag": flag.strip()})
    return pd.DataFrame(rows, columns=["table", "id", "name", "flag"])
