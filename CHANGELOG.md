# Changelog

## 0.2.0 - unreleased
- Add Streamlit showcase app (`app/app.py`) and `poliml.data_access` helpers.
- Add `poliml.copolymer`: Mayo-Lewis composition, azeotrope, tendency classification, mole-to-weight
  conversion, and feed -> polymer composition -> Fox Tg.
- Add `data/reactivity_ratios.csv` template (pairwise r1/r2, pre-filled with the copolymer pairs; values to be added with sources).
- Add `scripts/make_reactivity_template.py` and `scripts/check_reactivity_ratios.py`.
- Add packaging (`pyproject.toml`), CI, license and `.gitignore`.

## 0.1.0
- Cleaned and linked reference tables for chain- and step-growth monomers, homopolymers and copolymers.
- SMILES / PSMILES curation table with formula and MW cross-checks and `review_flags`.
- Fox-equation consistency check of copolymer Tg (`scripts/fox_validation.py`).
