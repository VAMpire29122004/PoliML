# PoliML data

All tables are produced from `data/raw/` by `python scripts/build_data.py`.
Raw files are never edited. Suspicious values are kept as given and described in `review_flags`.

| File | Rows | Key | Contents |
|---|---|---|---|
| `monomers_chain_growth.csv` | 46 | `monomer_id` | 44 original monomers + 2 stubs (vinylidene fluoride, N,N-dimethylacrylamide) added because other tables reference them |
| `monomers_step_growth.csv` | 25 | `monomer_id` | step-growth monomers |
| `homopolymers.csv` | 36 | `polymer_id` | linked to monomers by `monomer_id`; carries the repeat-unit PSMILES |
| `copolymers_chain_growth.csv` | 25 | `polymer_id` | linked to two monomers; empty `weight_fraction_a` column to fill with sourced compositions |
| `copolymers_step_growth.csv` | 15 | `polymer_id` | linked to monomers where they exist in the step-growth monomer table |
| `fox_validation.csv` | 25 | `copolymer_id` | output of `scripts/fox_validation.py` |

## Conventions
- Temperatures in deg C (`tg_c`, `tm_c`, `boiling_point_c`); density in g/cm3; molecular weight in g/mol.
- Range strings were split into numeric columns: `mw_min/max_g_mol`, `crystallinity_min/max_pct`, `mfi_min/max_g_10min`.
  Non-numeric MFI entries ("Limited data", "Not applicable", ...) are kept in `mfi_note`.
- `crystallinity_class` is one of `amorphous`, `semicrystalline`, `highly_crystalline`.
- List columns (`processing_techniques`, `applications`) are separated by `; `.
- `smiles`, `psmiles`, and the formula/MW computed from them were drafted for this project and
  should be spot-checked against a trusted source (PubChem). PSMILES use `*` for chain connection points.
- `formula_listed` / `mw_listed` are the original values; `formula_from_smiles` / `mw_from_smiles` are computed.
- `reactivity_ratio_legacy` is the original single value per monomer. Reactivity ratios (r1, r2) belong to a
  monomer *pair* and are not well defined per monomer; replace with a pairwise table before using Mayo-Lewis.
- Gaseous-monomer densities have unspecified measurement conditions.

## Known issues surfaced by validation (see `review_flags`)
- Listed formula/MW disagree with the structure: trifluoroethylene, lauryl acrylate, diallyl maleate,
  N-vinylcaprolactam (formula and MW), dodecyl methacrylate (MW 270.41 vs 254.41 from C16H30O2).
- PNIPAAm Tg of 32 C matches its LCST, not its glass transition. PDMAEMA Tg may be high. Verify both.
- Polymers listed as amorphous but given a Tm (PMMA, PC, PU, PI, PAI, PBI, PSU, PEC, polyarylates).
- Several "step-growth monomers" are chain-growth or ring-opening monomers, or a solvent (1,4-dioxane).
