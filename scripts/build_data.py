"""Build cleaned PoliML tables from the raw CSVs.

    python scripts/build_data.py

Reads  data/raw/*.csv and data/curation/structures.csv
Writes data/*.csv  (cleaned, numeric, linked by ids, with review_flags)

Nothing in the raw data is silently corrected: suspicious values are kept as
given and described in the `review_flags` column instead.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from poliml.chem import formula_and_mw  # noqa: E402

RAW = ROOT / "data" / "raw"
OUT = ROOT / "data"
CUR = ROOT / "data" / "curation"

# names used in one table that refer to a monomer named differently in another
ALIASES = {
    "2-Hydroxyethyl Methacrylate": "Hydroxyethyl Methacrylate",
    "Ethylhexyl Acrylate": "2-Ethylhexyl Acrylate",
    "Trifluoroethyl Methacrylate": "2,2,2-Trifluoroethyl Methacrylate",
    "Perfluoro Acrylate": "Perfluoroalkyl Acrylate",
    "Dimethyl Acrylamide": "N,N-Dimethylacrylamide",
    "Hexamethylene Diamine": "Hexamethylenediamine",
    "1,4-Butanediol": "Butanediol",
    "Para-phenylene Diamine": "Phenylenediamine",
}
ID_OVERRIDES = {"2,2,2-Trifluoroethyl Methacrylate": "trifluoroethyl_methacrylate"}

ABBR_RE = re.compile(r"\s*\(([A-Za-z0-9\-,]{2,14})\)\s*$")
NUM_RE = re.compile(r"\d+(?:,\d{3})*(?:\.\d+)?")


# ----------------------------------------------------------------- helpers
def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def split_abbr(name: str) -> tuple[str, str]:
    """'Polystyrene (PS)' -> ('Polystyrene', 'PS')."""
    m = ABBR_RE.search(name)
    if m:
        return name[: m.start()].strip(), m.group(1)
    return name.strip(), ""


def parse_range(text) -> tuple[float, float] | tuple[None, None]:
    if not isinstance(text, str):
        return (None, None)
    nums = [float(n.replace(",", "")) for n in NUM_RE.findall(text)]
    if not nums:
        return (None, None)
    return (min(nums), max(nums))


def crystallinity_fields(text) -> dict:
    lo, hi = parse_range(text)
    low = str(text).lower()
    if "amorphous" in low:
        cls = "amorphous"
    elif "high" in low:
        cls = "highly_crystalline"
    elif lo is not None:
        cls = "semicrystalline"
    else:
        cls = ""
    return {"crystallinity_text": text, "crystallinity_min_pct": lo, "crystallinity_max_pct": hi,
            "crystallinity_class": cls}


def mfi_fields(text) -> dict:
    lo, hi = parse_range(text)
    return {"mfi_min_g_10min": lo, "mfi_max_g_10min": hi, "mfi_note": "" if lo is not None else text}


def lists(text) -> str:
    return "; ".join(p.strip() for p in str(text).split(",")) if isinstance(text, str) else ""


def add_flag(flags: list[str], text: str) -> None:
    flags.append(text)


def structure_columns(df: pd.DataFrame, formula_col: str, mw_col: str) -> pd.DataFrame:
    """Add formula/MW computed from SMILES and a listed-vs-computed check."""
    fs, mws = [], []
    for smi in df["smiles"]:
        if isinstance(smi, str) and smi:
            f, mw = formula_and_mw(smi)
            fs.append(f)
            mws.append(round(mw, 2))
        else:
            fs.append(None)
            mws.append(np.nan)
    df["formula_from_smiles"] = fs
    df["mw_from_smiles"] = mws
    return df


def listed_vs_structure_flags(row, formula_col="formula_listed", mw_col="mw_listed") -> list[str]:
    out = []
    if pd.isna(row["formula_from_smiles"]):
        if not isinstance(row.get("smiles"), str) or not row["smiles"]:
            out.append("no structure assigned")
        return out
    if pd.isna(row[formula_col]) and pd.isna(row[mw_col]):
        return out  # added stub row: nothing listed to compare against
    bad_formula = row[formula_col] != row["formula_from_smiles"]
    bad_mw = pd.notna(row[mw_col]) and abs(row[mw_col] - row["mw_from_smiles"]) > 0.15
    if bad_formula or bad_mw:
        out.append(
            f"listed formula/MW ({row[formula_col]}, {row[mw_col]}) do not match structure "
            f"({row['formula_from_smiles']}, {row['mw_from_smiles']})"
        )
    return out


# ------------------------------------------------------------ chain monomers
def build_chain_monomers(structures: pd.DataFrame) -> pd.DataFrame:
    raw = pd.read_csv(RAW / "monomers.csv")
    names, abbrs = zip(*[split_abbr(n) for n in raw["Monomer"]])
    df = pd.DataFrame({
        "name": names,
        "abbreviation": abbrs,
        "formula_listed": raw["Formula"],
        "physical_state": raw["Physical State"],
        "mw_listed": raw["Molecular_Weight"],
        "density_g_cm3": raw["Density"],
        "boiling_point_c": raw["Boiling_Point"],
        "reactivity_ratio_legacy": raw["Reactivity_Ratio"],
        "image_file": raw["Structure"],
        "data_origin": "user",
    })
    stubs = structures[(structures.kind == "chain") & ~structures.name.isin(df.name)]
    stub_rows = pd.DataFrame({"name": stubs.name.values, "data_origin": "stub_from_structure"})
    df = pd.concat([df, stub_rows], ignore_index=True)

    st = structures[structures.kind == "chain"][["name", "smiles", "psmiles", "structure_note"]]
    df = df.merge(st, on="name", how="left")
    df["monomer_id"] = [ID_OVERRIDES.get(n, slug(n)) for n in df["name"]]
    df["polymerization_type"] = "chain_growth"
    df = structure_columns(df, "formula_listed", "mw_listed")

    manual = {
        "Tetrafluoroethylene": "density 4.1 looks like a gas density in g/L; other rows are g/cm3 (liquid TFE is about 1.5 g/cm3 - verify)",
        "Ethylene": "density 0.57 looks like liquid at the boiling point; condition unspecified",
        "N-Vinylpyrrolidone": "boiling point 90 C is probably quoted at reduced pressure (verify)",
        "Vinyl Ether": "name is ambiguous (ethyl vinyl ether vs divinyl ether)",
    }
    flags = []
    for _, r in df.iterrows():
        f = listed_vs_structure_flags(r)
        if pd.notna(r["boiling_point_c"]) and isinstance(r["physical_state"], str):
            if r["physical_state"] == "Gas" and r["boiling_point_c"] > 25:
                f.append(f"state is Gas but boiling point is {r['boiling_point_c']:g} C")
            if r["physical_state"] == "Liquid" and r["boiling_point_c"] < 15:
                f.append(f"state is Liquid but boiling point is {r['boiling_point_c']:g} C")
        if r["name"] in manual:
            f.append(manual[r["name"]])
        if r["physical_state"] == "Gas" and r["name"] not in manual:
            f.append("density of a gaseous monomer: measurement condition unspecified")
        if isinstance(r["structure_note"], str) and r["structure_note"] and "Added because" not in r["structure_note"]:
            f.append("structure: " + r["structure_note"])
        if r["data_origin"] == "stub_from_structure":
            f.append(r["structure_note"])
        flags.append(" | ".join(f))
    df["review_flags"] = flags

    cols = ["monomer_id", "name", "abbreviation", "polymerization_type", "smiles", "psmiles",
            "formula_listed", "formula_from_smiles", "mw_listed", "mw_from_smiles", "physical_state",
            "density_g_cm3", "boiling_point_c", "reactivity_ratio_legacy", "image_file", "data_origin",
            "review_flags"]
    return df[cols]


# ------------------------------------------------------------- step monomers
def build_step_monomers(structures: pd.DataFrame, chain: pd.DataFrame) -> pd.DataFrame:
    raw = pd.read_csv(RAW / "monomer_step_growth.csv")
    df = pd.DataFrame({
        "name": raw["Monomer"],
        "monomer_type": raw["Monomer Type"],
        "formula_listed": raw["Formula"],
        "physical_state": raw["Physical State"],
        "mw_listed": raw["Molecular Weight (g/mol)"],
        "density_g_cm3": raw["Density (g/cm³)"],
        "melting_point_c": raw["Melting Point (°C)"],
        "boiling_point_c": pd.to_numeric(raw["Boiling Point (°C)"], errors="coerce"),
        "boiling_point_note": np.where(pd.to_numeric(raw["Boiling Point (°C)"], errors="coerce").isna(),
                                       raw["Boiling Point (°C)"].astype(str).str.lower(), ""),
        "functional_groups": raw["Functional Groups"],
        "image_file": raw["Structure"],
    })
    st = structures[structures.kind == "step"][["name", "smiles", "structure_note"]]
    df = df.merge(st, on="name", how="left")
    df["monomer_id"] = [slug(n) for n in df["name"]]
    df = structure_columns(df, "formula_listed", "mw_listed")

    chain_ids = dict(zip(chain["name"], chain["monomer_id"]))
    df["also_in_chain_growth_table"] = [chain_ids.get(n, "") for n in df["name"]]

    manual = {
        "1,4-Dioxane": "a solvent, not a polymerizable monomer",
        "Methyl Methacrylate": "chain-growth (vinyl) monomer, not a step-growth monomer",
        "Vinyl Acetate": "chain-growth (vinyl) monomer, not a step-growth monomer",
        "Acrylonitrile": "chain-growth (vinyl) monomer, not a step-growth monomer",
        "Methacrylonitrile": "chain-growth (vinyl) monomer, not a step-growth monomer",
        "Ethyleneimine": "polymerizes by ring-opening (chain-type), not classic step-growth",
        "Caprolactam": "polymerizes by ring-opening (nylon 6), not classic step-growth",
        "Formaldehyde": "functional group is -CH=O (aldehyde); listed as -CHO",
    }
    flags = []
    for _, r in df.iterrows():
        f = listed_vs_structure_flags(r)
        if r["name"] in manual:
            f.append(manual[r["name"]])
        if isinstance(r["structure_note"], str) and r["structure_note"]:
            f.append("structure: " + r["structure_note"])
        flags.append(" | ".join(f))
    df["review_flags"] = flags
    cols = ["monomer_id", "name", "monomer_type", "smiles", "formula_listed", "formula_from_smiles",
            "mw_listed", "mw_from_smiles", "physical_state", "density_g_cm3", "melting_point_c",
            "boiling_point_c", "boiling_point_note", "functional_groups", "also_in_chain_growth_table",
            "image_file", "review_flags"]
    return df[cols]


# ------------------------------------------------------------- polymer tables
def common_polymer_fields(raw: pd.DataFrame, name_col: str) -> pd.DataFrame:
    names, abbrs = zip(*[split_abbr(n) for n in raw[name_col]])
    ids = [a.lower().replace("-", "_") if a else slug(n) for n, a in zip(names, abbrs)]
    assert len(set(ids)) == len(ids), f"duplicate polymer ids: {ids}"
    df = pd.DataFrame({
        "polymer_id": ids,
        "name": names,
        "abbreviation": abbrs,
        "physical_state": raw["Physical State"],
        "density_g_cm3": raw["Density (g/cm³)"],
        "tg_c": raw["Tg (°C)"],
        "tm_c": raw["Tm (°C)"],
    })
    cry = pd.DataFrame([crystallinity_fields(t) for t in raw["Crystallinity (%)"]])
    mfi = pd.DataFrame([mfi_fields(t) for t in raw["MFI (g/10 min)"]])
    df = pd.concat([df, cry, mfi], axis=1)
    df["processing_techniques"] = [lists(t) for t in raw["Processing Techniques"]]
    df["applications"] = [lists(t) for t in raw["Applications"]]
    df["image_file"] = raw["Structure"]
    return df


def generic_polymer_flags(r) -> list[str]:
    f = []
    if r["crystallinity_class"] == "amorphous" and pd.notna(r["tm_c"]):
        f.append(f"Tm {r['tm_c']:g} C listed but crystallinity is 'amorphous (0%)'")
    if pd.notna(r["tm_c"]) and pd.notna(r["tg_c"]) and r["tm_c"] < r["tg_c"]:
        f.append("Tm is below Tg")
    return f


def link(name, lookup: dict[str, str]) -> str:
    if not isinstance(name, str):
        return ""
    return lookup.get(ALIASES.get(name, name), "")


def build_homopolymers(chain: pd.DataFrame) -> pd.DataFrame:
    raw = pd.read_csv(RAW / "homopolymers.csv")
    df = common_polymer_fields(raw, "Homopolymer")
    mw = [parse_range(t) for t in raw["Molecular Weight (g/mol)"]]
    df["mw_min_g_mol"] = [a for a, _ in mw]
    df["mw_max_g_mol"] = [b for _, b in mw]
    lookup = dict(zip(chain["name"], chain["monomer_id"]))
    df["monomer_name_listed"] = raw["Monomer"]
    df["monomer_id"] = [link(n, lookup) for n in raw["Monomer"]]
    ps = dict(zip(chain["monomer_id"], chain["psmiles"]))
    mf = dict(zip(chain["monomer_id"], chain["formula_from_smiles"]))
    df["repeat_unit_psmiles"] = [ps.get(i, "") if isinstance(ps.get(i), str) else "" for i in df["monomer_id"]]
    listed = [re.sub(r"^\((.*)\)n$", r"\1", f) for f in raw["Formula"]]
    df["repeat_unit_formula_listed"] = listed
    df["formula_consistent_with_monomer_structure"] = [
        (listed_f == mf.get(i)) if isinstance(mf.get(i), str) else "" for listed_f, i in zip(listed, df["monomer_id"])
    ]

    manual = {
        "pnipaam": "Tg 32 C equals PNIPAAm's LCST (~32 C in water); the dry glass transition is reported well above 100 C - verify",
        "pdmaema": "Tg 55 C is higher than commonly reported values for PDMAEMA (about 20 C) - verify",
        "ppfa": "linked to 'Perfluoroalkyl Acrylate' but the formula (C5H3FO2)n does not match that monomer; generic name",
    }
    flags = []
    for _, r in df.iterrows():
        f = generic_polymer_flags(r)
        if r["monomer_id"] == "":
            f.append(f"monomer '{r['monomer_name_listed']}' not found in monomer table")
        if r["formula_consistent_with_monomer_structure"] is False:
            f.append("listed repeat-unit formula differs from the monomer structure")
        if r["polymer_id"] in manual:
            f.append(manual[r["polymer_id"]])
        flags.append(" | ".join(f))
    df["review_flags"] = flags
    cols = ["polymer_id", "name", "abbreviation", "monomer_name_listed", "monomer_id", "repeat_unit_psmiles",
            "repeat_unit_formula_listed", "formula_consistent_with_monomer_structure", "physical_state",
            "mw_min_g_mol", "mw_max_g_mol", "density_g_cm3", "tg_c", "tm_c", "crystallinity_text",
            "crystallinity_min_pct", "crystallinity_max_pct", "crystallinity_class", "mfi_min_g_10min",
            "mfi_max_g_10min", "mfi_note", "processing_techniques", "applications", "image_file",
            "review_flags"]
    return df[cols]


def build_copolymers(raw_file: str, lookup: dict[str, str], extra_flags: dict[str, str] | None = None,
                     add_composition_column: bool = False) -> pd.DataFrame:
    raw = pd.read_csv(RAW / raw_file)
    df = common_polymer_fields(raw, "Copolymer")
    df["monomer_a_name"] = raw["Monomer A"]
    df["monomer_b_name"] = raw["Monomer B"]
    df["monomer_a_id"] = [link(n, lookup) for n in raw["Monomer A"]]
    df["monomer_b_id"] = [link(n, lookup) for n in raw["Monomer B"]]
    if add_composition_column:
        df["weight_fraction_a"] = np.nan  # fill with sourced values to enable Fox testing
    extra_flags = extra_flags or {}
    flags = []
    for _, r in df.iterrows():
        f = generic_polymer_flags(r)
        for side in ("a", "b"):
            if r[f"monomer_{side}_id"] == "":
                f.append(f"monomer '{r[f'monomer_{side}_name']}' not in monomer table")
        if r["polymer_id"] in extra_flags:
            f.append(extra_flags[r["polymer_id"]])
        flags.append(" | ".join(f))
    df["review_flags"] = flags
    cols = ["polymer_id", "name", "abbreviation", "monomer_a_name", "monomer_b_name", "monomer_a_id",
            "monomer_b_id"] + (["weight_fraction_a"] if add_composition_column else []) + [
            "physical_state", "density_g_cm3", "tg_c", "tm_c", "crystallinity_text", "crystallinity_min_pct",
            "crystallinity_max_pct", "crystallinity_class", "mfi_min_g_10min", "mfi_max_g_10min", "mfi_note",
            "processing_techniques", "applications", "image_file", "review_flags"]
    return df[cols]


# --------------------------------------------------------------------- main
def write(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, index=False, float_format="%.6g")
    print(f"wrote data/{name}: {df.shape[0]} rows x {df.shape[1]} columns")


def main() -> None:
    structures = pd.read_csv(CUR / "structures.csv").fillna("")
    chain = build_chain_monomers(structures)
    step = build_step_monomers(structures, chain)
    homo = build_homopolymers(chain)

    chain_lookup = dict(zip(chain["name"], chain["monomer_id"]))
    co_chain = build_copolymers(
        "copolymers.csv", chain_lookup, add_composition_column=True,
        extra_flags={
            "abs": "ABS is a terpolymer / grafted blend, not a binary random copolymer; Tg of the rubber and SAN phases differ",
            "sbr": "composition not given (commercial SBR is typically ~23% styrene by weight)",
        },
    )
    step_lookup = dict(zip(step["name"], step["monomer_id"]))
    co_step = build_copolymers(
        "copolymers_step_growth.csv", step_lookup,
        extra_flags={
            "pps": "listed monomers (p-phenylenediamine + sulfur dichloride) are not the usual PPS route "
                   "(p-dichlorobenzene + sodium sulfide) - verify",
        },
    )
    write(chain, "monomers_chain_growth.csv")
    write(step, "monomers_step_growth.csv")
    write(homo, "homopolymers.csv")
    write(co_chain, "copolymers_chain_growth.csv")
    write(co_step, "copolymers_step_growth.csv")

    # integrity checks
    for df, key in [(chain, "monomer_id"), (step, "monomer_id"), (homo, "polymer_id"),
                    (co_chain, "polymer_id"), (co_step, "polymer_id")]:
        assert df[key].is_unique, f"duplicate {key}"
    assert set(homo["monomer_id"]) - {""} <= set(chain["monomer_id"])
    n_flag = {n: int((d["review_flags"] != "").sum()) for n, d in
              [("chain monomers", chain), ("step monomers", step), ("homopolymers", homo),
               ("chain copolymers", co_chain), ("step copolymers", co_step)]}
    print("rows with review flags:", n_flag)


if __name__ == "__main__":
    main()
