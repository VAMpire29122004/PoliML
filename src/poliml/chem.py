"""Formula / molecular-weight helpers.

Uses RDKit when it is installed. Otherwise falls back to a tiny parser that
supports Kekule SMILES over the organic subset (C N O F Cl Br I S P Si), simple
bracket atoms such as [Si], and '*' dummy atoms. Aromatic (lower-case) SMILES
need RDKit.
"""
from __future__ import annotations

import re
from collections import Counter

ATOMIC_WEIGHTS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904,
}
DEFAULT_VALENCE = {"C": 4, "N": 3, "O": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1, "S": 2, "P": 3, "Si": 4}
_TOKEN = re.compile(r"\[([A-Z][a-z]?)[^\]]*\]|Cl|Br|Si|[BCNOSPFI]|\*|=|#|-|\(|\)|\d")
_BOND_ORDER = {"-": 1, "=": 2, "#": 3}


def _parse(smiles: str) -> Counter:
    atoms: list[str] = []
    bracket: list[bool] = []
    bonds: list[int] = []
    stack: list[int] = []
    rings: dict[str, tuple[int, int]] = {}
    prev = None
    order = 1
    pos = 0
    while pos < len(smiles):
        m = _TOKEN.match(smiles, pos)
        if not m:
            raise ValueError(f"Unsupported SMILES syntax at position {pos}: {smiles!r}")
        tok = m.group(0)
        pos = m.end()
        if tok in _BOND_ORDER:
            order = _BOND_ORDER[tok]
        elif tok == "(":
            stack.append(prev)
        elif tok == ")":
            prev = stack.pop()
        elif tok.isdigit():
            if tok in rings:
                j, o = rings.pop(tok)
                o = max(o, order)
                bonds[prev] += o
                bonds[j] += o
            else:
                rings[tok] = (prev, order)
            order = 1
        else:  # atom
            element = m.group(1) if tok.startswith("[") else tok
            idx = len(atoms)
            atoms.append(element)
            bracket.append(tok.startswith("["))
            bonds.append(0)
            if prev is not None:
                bonds[prev] += order
                bonds[idx] += order
            order = 1
            prev = idx
    if rings or stack:
        raise ValueError(f"Unbalanced rings or branches in {smiles!r}")
    counts: Counter = Counter()
    for el, is_br, b in zip(atoms, bracket, bonds):
        if el == "*":
            continue
        counts[el] += 1
        if not is_br:
            h = max(DEFAULT_VALENCE[el] - b, 0)
            if h:
                counts["H"] += h
    return counts


def _hill(counts: Counter) -> str:
    keys = sorted(counts)
    if "C" in counts:
        keys = ["C"] + (["H"] if "H" in counts else []) + [k for k in keys if k not in ("C", "H")]
    return "".join(f"{k}{counts[k] if counts[k] > 1 else ''}" for k in keys)


def formula_and_mw(smiles: str) -> tuple[str, float]:
    """Return (Hill formula, molecular weight in g/mol) for a SMILES string."""
    try:  # pragma: no cover - exercised only when RDKit is installed
        from rdkit import Chem
        from rdkit.Chem.Descriptors import MolWt
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit could not parse {smiles!r}")
        return CalcMolFormula(mol), float(MolWt(mol))
    except ImportError:
        counts = _parse(smiles)
        mw = sum(ATOMIC_WEIGHTS[k] * n for k, n in counts.items())
        return _hill(counts), round(mw, 3)
