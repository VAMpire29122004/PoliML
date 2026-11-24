# PoliML

Polymer reference data, validation tools, and (in progress) property prediction from polymer
structure (SMILES / PSMILES) with a Streamlit app.

## Status
- Curated, cleaned reference tables for chain-growth and step-growth monomers and polymers (`data/`)
- Fox-equation consistency check of copolymer Tg values (`scripts/fox_validation.py`)
- Mayo-Lewis copolymer composition tools (`poliml.copolymer`); the pairwise reactivity-ratio table
  (`data/reactivity_ratios.csv`) is a template awaiting sourced values
- Streamlit showcase app (`app/app.py`) over the data and tools above
- Planned: SMILES -> property models, tensile / DSC / rheology curve analysis

## Run the app
```bash
pip install -r requirements.txt
streamlit run app/app.py
```
Tabs: overview, reference library (with structure drawings), copolymer explorer (Mayo-Lewis + Fox Tg),
Fox consistency check, data-quality flags, SMILES formula/MW tool.
Deploy on Streamlit Community Cloud with main file path `app/app.py`
(if the build fails on `rdkit`, remove that line from `requirements.txt`; structure drawing is then skipped).

## Quick start
```bash
pip install -e ".[dev]"        # add ",chem" to install RDKit
python scripts/build_data.py
python scripts/fox_validation.py
python scripts/check_reactivity_ratios.py
pytest
```

```python
from poliml.copolymer import predict_tg_from_feed
# illustrative numbers only
predict_tg_from_feed(f1=0.5, r1=0.5, r2=0.5, mw1=104.15, mw2=100.12, tg1_c=100, tg2_c=105)
```

## Data
See `data/DATA_DICTIONARY.md`. Values are hand-curated and flagged where unverified (`review_flags`).
The Kaggle dataset is not committed; see `scripts/download_data.md`.
SMILES/PSMILES were drafted during the project and should be spot-checked against PubChem.
