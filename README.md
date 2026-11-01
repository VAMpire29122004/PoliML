# PoliML

Polymer reference data, validation tools, and (in progress) property prediction from
polymer structure (SMILES / PSMILES) with a Streamlit app.

## Status
- Curated, cleaned reference tables for chain-growth and step-growth monomers and polymers (`data/`)
- Fox-equation consistency check of copolymer Tg values (`scripts/fox_validation.py`)
- Planned: SMILES -> property models, pairwise reactivity-ratio (Mayo-Lewis) tool, tensile / DSC / rheology curve analysis

## Quick start
```bash
pip install -r requirements.txt
python scripts/build_data.py
python scripts/fox_validation.py
python -m pytest        # or run the tests in tests/
```

## Data
See `data/DATA_DICTIONARY.md`. Values are hand-curated and flagged where unverified (`review_flags`).
Large third-party datasets (e.g. Kaggle) are not committed; see `scripts/download_data.md` (to be added).
