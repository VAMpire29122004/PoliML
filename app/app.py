"""PoliML showcase app.

    streamlit run app/app.py

Only uses data and code that exist in this repository. Property prediction from
structure (Kaggle-trained models) and curve analysis are on the roadmap, not here.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from poliml.chem import formula_and_mw  # noqa: E402
from poliml.copolymer import (azeotrope_f1, classify, composition_curve, mayo_lewis_F1,  # noqa: E402
                              predict_tg_from_feed)
from poliml.data_access import (all_flags, homopolymer_tg, load_tables, lookup_reactivity,  # noqa: E402
                                monomer_mw, monomers_for_copolymer_tool)

st.set_page_config(page_title="PoliML", page_icon="🧪", layout="wide")


@st.cache_data
def get_tables() -> dict[str, pd.DataFrame]:
    return load_tables()


def draw_structure(smiles: str, size=(320, 240)):
    """RDKit drawing if RDKit is installed, else None."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=size) if mol is not None else None


def show_structure(smiles: str) -> None:
    img = draw_structure(smiles)
    if img is not None:
        st.image(img)
    else:
        st.caption("Install RDKit to draw structures.")


tables = get_tables()
monomers, homo = tables["monomers"], tables["homopolymers"]
fox = tables.get("fox")
flags = all_flags(tables)

st.title("🧪 PoliML")
st.caption("Curated polymer reference data, consistency checks and copolymer tools")

tab_over, tab_lib, tab_co, tab_fox, tab_dq, tab_struct = st.tabs(
    ["Overview", "Reference library", "Copolymer explorer", "Fox check", "Data quality", "Structure tool"])

# ---------------------------------------------------------------- Overview
with tab_over:
    testable = fox[fox["status"] != "not_testable"] if fox is not None else pd.DataFrame()
    c = st.columns(5)
    c[0].metric("Chain-growth monomers", int((monomers["data_origin"] == "user").sum()))
    c[1].metric("Step-growth monomers", len(tables["monomers_step"]))
    c[2].metric("Homopolymers", len(homo))
    c[3].metric("Copolymers", len(tables["copolymers"]) + len(tables["copolymers_step"]))
    c[4].metric("Flagged data issues", len(flags))

    st.markdown(
        """
**What this is.** A hand-curated polymer reference dataset (monomers, homopolymers, copolymers) linked by
structure (SMILES / PSMILES), with automated consistency checks and copolymer-composition tools.

**What it does today**
- Browse the linked tables and see monomer structures (*Reference library*)
- Predict copolymer composition and Tg from a monomer feed using Mayo-Lewis + the Fox equation (*Copolymer explorer*)
- Test every copolymer Tg against the Fox bounds of its homopolymers (*Fox check*)
- See every place the source data disagrees with itself (*Data quality*)
- Check a SMILES string for formula and molecular weight (*Structure tool*)

**On the roadmap (not in this app yet)**
- Structure → property models (Tg, density, ...) trained on a public benchmark dataset
- Tensile, DSC and rheology curve analysis
- Sourced pairwise reactivity ratios for every copolymer pair
"""
    )
    if len(testable):
        n_ok = int((testable["status"] == "consistent").sum())
        st.info(f"Fox check: {n_ok} of {len(testable)} testable copolymers have a Tg between the Tgs of their "
                f"two homopolymers. The rest point to data issues or to structure the model doesn't capture.")

# ------------------------------------------------------- Reference library
with tab_lib:
    choice = st.radio("Table", ["Chain-growth monomers", "Step-growth monomers", "Homopolymers",
                                "Copolymers (chain-growth)", "Copolymers (step-growth)"], horizontal=True)
    key = {"Chain-growth monomers": "monomers", "Step-growth monomers": "monomers_step",
           "Homopolymers": "homopolymers", "Copolymers (chain-growth)": "copolymers",
           "Copolymers (step-growth)": "copolymers_step"}[choice]
    df = tables[key]
    q = st.text_input("Search", "")
    if q:
        mask = df.astype(str).apply(lambda col: col.str.contains(q, case=False, regex=False)).any(axis=1)
        df = df[mask]
    st.dataframe(df.drop(columns=["image_file"], errors="ignore"), hide_index=True)
    st.caption(f"{len(df)} rows")

    if key in ("monomers", "monomers_step"):
        st.subheader("Monomer detail")
        names = dict(zip(tables[key]["monomer_id"], tables[key]["name"]))
        sel = st.selectbox("Monomer", list(names), format_func=lambda i: names[i], key=f"sel_{key}")
        row = tables[key][tables[key]["monomer_id"] == sel].iloc[0]
        left, right = st.columns(2)
        with left:
            st.markdown(f"**{row['name']}**")
            if pd.notna(row.get("smiles")):
                st.code(row["smiles"], language=None)
                show_structure(row["smiles"])
            else:
                st.caption("No structure assigned.")
            if key == "monomers" and pd.notna(row.get("psmiles")):
                st.markdown("Repeat unit (PSMILES)")
                st.code(row["psmiles"], language=None)
        with right:
            st.markdown(
                f"Formula (listed / from structure): `{row['formula_listed']}` / `{row['formula_from_smiles']}`  \n"
                f"MW (listed / from structure): `{row['mw_listed']}` / `{row['mw_from_smiles']}`")
            if pd.notna(row.get("review_flags")):
                for f in str(row["review_flags"]).split(" | "):
                    st.warning(f)

# ------------------------------------------------------ Copolymer explorer
with tab_co:
    st.subheader("Feed → copolymer composition → Tg")
    usable = monomers_for_copolymer_tool(monomers, homo)
    names = dict(zip(usable["monomer_id"], usable["name"]))
    ids = list(names)
    mw, tg = monomer_mw(monomers), homopolymer_tg(homo)

    col_a, col_b = st.columns(2)
    with col_a:
        id1 = st.selectbox("Monomer 1", ids, index=ids.index("styrene") if "styrene" in ids else 0,
                           format_func=lambda i: names[i])
    with col_b:
        default2 = ids.index("methyl_methacrylate") if "methyl_methacrylate" in ids else 1
        id2 = st.selectbox("Monomer 2", ids, index=default2, format_func=lambda i: names[i])

    if id1 == id2:
        st.warning("Choose two different monomers.")
    else:
        found = lookup_reactivity(tables.get("reactivity"), id1, id2)
        if found:
            r1_default, r2_default, rr_row = found
            src = rr_row["source"] if pd.notna(rr_row["source"]) else "no source given"
            st.success(f"Reactivity ratios for this pair come from data/reactivity_ratios.csv ({src}).")
        else:
            r1_default, r2_default = 1.0, 1.0
            st.info("No sourced reactivity ratios for this pair in the repository yet. Set r1 and r2 manually; "
                    "the defaults (1.0 / 1.0) are placeholders, not literature values.")

        ctrl, plots = st.columns([1, 2])
        with ctrl:
            r1 = st.number_input("r1", min_value=0.0, value=float(r1_default), step=0.05, format="%.3f",
                                 key=f"r1_{id1}_{id2}")
            r2 = st.number_input("r2", min_value=0.0, value=float(r2_default), step=0.05, format="%.3f",
                                 key=f"r2_{id1}_{id2}")
            f1 = st.slider("Mole fraction of monomer 1 in the feed", min_value=0.0, max_value=1.0, value=0.5,
                           step=0.01)
            code, text = classify(r1, r2)
            st.markdown(f"**{code.replace('_', ' ')}**: {text}")
            res = predict_tg_from_feed(f1, r1, r2, mw[id1], mw[id2], tg[id1], tg[id2])
            st.metric("Polymer composition, mole fraction of 1", f"{res['F1_mole']:.3f}")
            st.metric("Weight fraction of 1", f"{res['w1']:.3f}")
            st.metric("Fox Tg (°C)", f"{res['tg_c']:.1f}")
            for pid in (id1, id2):
                hrow = homo[homo["monomer_id"] == pid].iloc[0]
                if pd.notna(hrow["review_flags"]) and "Tg" in str(hrow["review_flags"]):
                    st.warning(f"{hrow['name']}: {hrow['review_flags']}")

        with plots:
            f_grid, F_grid = composition_curve(r1, r2)
            fig = Figure(figsize=(5, 4))
            ax = fig.subplots()
            ax.plot([0, 1], [0, 1], "--", color="0.6", label="F1 = f1")
            ax.plot(f_grid, F_grid, color="tab:blue", label="Mayo-Lewis")
            ax.scatter([f1], [res["F1_mole"]], color="tab:red", zorder=3, label="selected feed")
            az = azeotrope_f1(r1, r2)
            if az is not None:
                ax.scatter([az], [az], marker="x", color="k", zorder=3, label=f"azeotrope f1={az:.2f}")
            ax.set_xlabel("f1 (feed, mole fraction)")
            ax.set_ylabel("F1 (polymer, mole fraction)")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.legend(fontsize=8)
            fig.tight_layout()

            tg_curve = [predict_tg_from_feed(x, r1, r2, mw[id1], mw[id2], tg[id1], tg[id2])["tg_c"]
                        for x in np.linspace(0, 1, 101)]
            fig2 = Figure(figsize=(5, 4))
            ax2 = fig2.subplots()
            ax2.plot(np.linspace(0, 1, 101), tg_curve, color="tab:green")
            ax2.scatter([f1], [res["tg_c"]], color="tab:red", zorder=3)
            ax2.set_xlabel("f1 (feed, mole fraction)")
            ax2.set_ylabel("Tg (°C)")
            ax2.set_title(f"{names[id2]} ← → {names[id1]}", fontsize=9)
            fig2.tight_layout()
            p1, p2 = st.columns(2)
            with p1:
                st.pyplot(fig)
            with p2:
                st.pyplot(fig2)
        st.caption("Terminal-model (Mayo-Lewis) composition at low conversion; Tg from the Fox equation using the "
                   "homopolymer Tgs in this repository. Composition drifts as conversion increases, and Fox ignores "
                   "sequence effects, so treat results as estimates.")

# --------------------------------------------------------------- Fox check
with tab_fox:
    st.subheader("Do reported copolymer Tgs fall between their homopolymers' Tgs?")
    if fox is None:
        st.warning("data/fox_validation.csv not found. Run scripts/fox_validation.py.")
    else:
        t = fox[fox["status"] != "not_testable"].sort_values("tg_reported_c").reset_index(drop=True)
        fig = Figure(figsize=(8, 0.4 * len(t) + 1.5))
        ax = fig.subplots()
        for i, r in t.iterrows():
            ok = r["status"] == "consistent"
            ax.plot([r["tg_lower_c"], r["tg_upper_c"]], [i, i], color="0.7", lw=6, solid_capstyle="round")
            ax.scatter(r["tg_reported_c"], i, s=40, color="tab:green" if ok else "tab:red", zorder=3)
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels(t["copolymer_id"].str.upper(), fontsize=8)
        ax.set_xlabel("Tg (°C)")
        ax.set_title("Grey bar = range of homopolymer Tgs; red dot = reported Tg outside it", fontsize=10)
        ax.grid(axis="x", alpha=0.3)
        fig.tight_layout()
        st.pyplot(fig)
        bad = t[t["status"] == "inconsistent"][["copolymer_id", "tg_a_c", "tg_b_c", "tg_reported_c",
                                                "deviation_outside_bounds_c", "implied_weight_fraction_a"]]
        st.markdown("**Outside the bounds**")
        st.dataframe(bad, hide_index=True)
        st.caption("A random copolymer's Tg should lie between its homopolymers' Tgs. Misses can come from "
                   "non-random structure (e.g. ABS), small rounding differences, or errors in the source values. "
                   "Copolymer compositions are not in the data, so a full Fox prediction is not possible yet.")

# ------------------------------------------------------------ Data quality
with tab_dq:
    st.subheader("Where the source data disagrees with itself")
    st.caption("Raw values are never edited; each issue is described instead.")
    options = sorted(flags["table"].unique())
    chosen = st.multiselect("Tables", options, default=options)
    shown = flags[flags["table"].isin(chosen)]
    st.dataframe(shown, hide_index=True)
    st.caption(f"{len(shown)} flags")

# ---------------------------------------------------------- Structure tool
with tab_struct:
    st.subheader("Formula and molecular weight from a SMILES string")
    smi = st.text_input("SMILES", "C=CC1=CC=CC=C1")
    if smi:
        try:
            formula, mol_wt = formula_and_mw(smi)
            st.markdown(f"Formula: `{formula}`  \nMolecular weight: `{mol_wt:.2f}` g/mol")
            show_structure(smi)
        except Exception as exc:  # parser errors, unsupported syntax
            st.error(f"Could not read that SMILES: {exc}. (Aromatic lowercase SMILES need RDKit installed.)")
