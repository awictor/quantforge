"""Categorical feature encoding.

Turn categorical labels into numbers a model can use:

  * label encoding -- map each distinct category to an integer index (ordinal
    position, arbitrary but consistent).
  * one-hot encoding -- expand each category to a 0/1 indicator column, so no
    false ordinal relationship is implied.

Fit the encoder on the training categories, then apply the same mapping to test
data; unseen categories map to ``-1`` (label) or an all-zero row (one-hot). Pure
standard library.
"""


def fit_label_encoder(values):
    """Fit a label encoder: sorted distinct categories -> integer indices.

    Returns ``{"categories": [...]}`` where the list position is the code.
    """
    if len(values) == 0:
        raise ValueError("need at least one value")
    return {"categories": sorted(set(values))}


def label_encode(encoder, values):
    """Encode ``values`` to integer codes; unseen categories map to -1."""
    index = {c: i for i, c in enumerate(encoder["categories"])}
    return [index.get(v, -1) for v in values]


def label_decode(encoder, codes):
    """Invert label codes back to categories; -1 (unseen) maps to ``None``."""
    cats = encoder["categories"]
    return [cats[c] if 0 <= c < len(cats) else None for c in codes]


def one_hot_encode(encoder, values):
    """One-hot encode ``values`` against a fitted label encoder.

    Each row is a length-``len(categories)`` list of 0/1; a known category sets one
    entry to 1 (rows sum to 1), an unseen category yields an all-zero row.
    """
    cats = encoder["categories"]
    index = {c: i for i, c in enumerate(cats)}
    rows = []
    for v in values:
        row = [0.0] * len(cats)
        if v in index:
            row[index[v]] = 1.0
        rows.append(row)
    return rows
