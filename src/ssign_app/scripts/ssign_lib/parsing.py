"""Tolerant primitive parsers shared across pipeline scripts.

These exist because tool outputs occasionally encode primitives in surprising
ways (float strings for integers, range strings like ``"22-23"`` for SignalP
cleavage-site positions, empty cells for "no value"). Centralising the
tolerance rules in one place keeps the scripts that consume those outputs
from each defining their own slightly-different parser.
"""

from __future__ import annotations


def parse_int_or_none(value: object, allow_range: bool = False) -> int | None:
    """Parse an integer from a tool-output cell; return ``None`` on failure.

    Behaviour:
    - Empty / ``None`` / whitespace-only → ``None``.
    - Plain integer string (``"22"``) → ``22``.
    - Float string (``"22.0"``) → ``22`` (truncating via ``int(float(...))``).
    - With ``allow_range=True``, range strings (``"22-23"``) are split on
      the first ``-`` and the leading token parsed. Used for SignalP
      cleavage-site positions, which the wrapper sometimes emits as
      ``"22-23"`` and sometimes as ``"22"``.
    - Any other malformed value → ``None``.

    ``value`` is typed ``object`` because tool-output rows from ``csv.DictReader``
    surface as ``str`` but tests sometimes pass ``None`` directly.
    """
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    token = s.split("-", 1)[0].strip() if allow_range else s
    try:
        return int(float(token))
    except (ValueError, TypeError):
        return None
