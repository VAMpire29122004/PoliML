# External data

The Open Polymer Prediction (Kaggle, 2025) training data is **not** stored in this repository.
Competition data usually restricts redistribution, so read the competition's data terms first.

1. Download `train.csv` from the competition's Data tab on Kaggle.
2. Place it at `data/external/train.csv` (this folder is gitignored).

Columns: `id`, `SMILES` (PSMILES with `*` connection points) and the targets `Tg`, `FFV`, `Tc`,
`Density`, `Rg`. Each polymer is typically labelled for only one or two targets.
