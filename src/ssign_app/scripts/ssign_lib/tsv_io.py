"""Shared TSV loaders used by the pipeline scripts.

Replaces a near-duplicate "open a TSV, index its rows by some column"
pattern that grew up in cross_validate_predictions, enrichment_testing,
and the t5_passenger helper. The three variants differed only in:
- which column to key on (always ``locus_tag`` in the strict variants;
  cross_validate_predictions falls back to ``protein_id`` / ``seq_id``);
- whether a missing input file is an error or silently returns ``{}``.

Both axes are now parameters here.
"""

from __future__ import annotations

import csv
import os
from typing import Sequence


def load_tsv_by_key(
    path: str | os.PathLike[str],
    key_columns: Sequence[str] = ("locus_tag",),
    missing_ok: bool = True,
) -> dict[str, dict[str, str]]:
    """Read a TSV; return ``{key: row_dict}``.

    ``key_columns`` is a sequence; the FIRST one present in the TSV's
    header is used as the key. Lets callers express tolerant fallbacks
    (e.g. ``("locus_tag", "protein_id", "seq_id")`` — common for tool
    outputs whose ID-column name varies). Rows whose key value is empty
    or whitespace-only are skipped silently.

    When ``missing_ok=True`` (default) a non-existent or empty ``path``
    returns ``{}``. When ``False``, a missing file raises
    ``FileNotFoundError`` and an empty path raises ``ValueError``.

    Returns an empty dict if the TSV has no header or none of
    ``key_columns`` appears in the header.
    """
    if not path:
        if missing_ok:
            return {}
        raise ValueError("load_tsv_by_key: path must be non-empty")
    path = os.fspath(path)
    if not os.path.exists(path):
        if missing_ok:
            return {}
        raise FileNotFoundError(path)

    out: dict[str, dict[str, str]] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        key_col = next((c for c in key_columns if c in fieldnames), None)
        if key_col is None:
            return {}
        for row in reader:
            tag = (row.get(key_col) or "").strip()
            if tag:
                out[tag] = row
    return out
