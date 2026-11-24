import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from poliml.data_access import (all_flags, homopolymer_tg, load_tables, lookup_reactivity,  # noqa: E402
                                monomer_mw, monomers_for_copolymer_tool)


def test_tables_load_and_link():
    t = load_tables()
    assert len(t["monomers"]) == 46 and len(t["homopolymers"]) == 36
    assert set(t["homopolymers"]["monomer_id"].dropna()) - {""} <= set(t["monomers"]["monomer_id"])


def test_copolymer_tool_monomers_have_tg_and_mw():
    t = load_tables()
    m = monomers_for_copolymer_tool(t["monomers"], t["homopolymers"])
    mw, tg = monomer_mw(t["monomers"]), homopolymer_tg(t["homopolymers"])
    assert len(m) > 10
    assert all(i in mw and i in tg for i in m["monomer_id"])


def test_reactivity_lookup_swaps_order():
    rr = pd.DataFrame([{"monomer_1_id": "a", "monomer_2_id": "b", "r1": "0.5", "r2": "2.0"}])
    assert lookup_reactivity(rr, "a", "b")[:2] == (0.5, 2.0)
    assert lookup_reactivity(rr, "b", "a")[:2] == (2.0, 0.5)
    assert lookup_reactivity(rr, "a", "c") is None
    blank = pd.DataFrame([{"monomer_1_id": "a", "monomer_2_id": "b", "r1": "", "r2": ""}])
    assert lookup_reactivity(blank, "a", "b") is None
    assert lookup_reactivity(None, "a", "b") is None


def test_flags_table():
    f = all_flags(load_tables())
    assert len(f) > 20 and {"table", "id", "name", "flag"} <= set(f.columns)
