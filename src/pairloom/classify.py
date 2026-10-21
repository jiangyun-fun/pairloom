"""Read-pair classification and configurable BAM tagging.

Classifications are mutually exclusive and, with ``total``, exhaustive::

    total = mergeable + left_aligned_only + right_aligned_only
            + both_aligned + neither_aligned

``left`` and ``right`` refer to R1 and R2 (read number), never to physical
reference-coordinate order.
"""

from __future__ import annotations

from typing import Any

from pairloom.config import PairLoomConfig
from pairloom.exceptions import ConfigurationError

TOTAL = "total"
MERGEABLE = "mergeable"
LEFT_ALIGNED_ONLY = "left_aligned_only"
RIGHT_ALIGNED_ONLY = "right_aligned_only"
BOTH_ALIGNED = "both_aligned"
NEITHER_ALIGNED = "neither_aligned"

CLASSIFICATION_ORDER: tuple[str, ...] = (
    TOTAL,
    MERGEABLE,
    LEFT_ALIGNED_ONLY,
    RIGHT_ALIGNED_ONLY,
    BOTH_ALIGNED,
    NEITHER_ALIGNED,
)

_KNOWN = frozenset(CLASSIFICATION_ORDER)


class StatsCounter:
    """Integer counters for the six classification categories."""

    def __init__(self) -> None:
        self._counts = {name: 0 for name in CLASSIFICATION_ORDER}

    def add(self, classification: str) -> None:
        if classification not in self._counts:
            raise ValueError(f"unknown classification {classification!r}")
        self._counts[classification] += 1
        self._counts[TOTAL] += 1

    def counts(self) -> dict[str, int]:
        """Ordered counts including the total (which must equal the sum)."""
        result = {name: self._counts[name] for name in CLASSIFICATION_ORDER}
        parts = self._counts[TOTAL]
        summed = sum(result[name] for name in CLASSIFICATION_ORDER if name != TOTAL)
        if parts != summed:  # pragma: no cover - guards internal invariant
            raise AssertionError(
                f"classification invariant violated: total {parts} != sum {summed}"
            )
        return result


def assign_status_tag(record: Any, cfg: PairLoomConfig, classification: str) -> None:
    """Attach the configured status tag to one record.

    A pre-existing tag whose value is not a PairLoom classification string is
    treated as an unrelated value: PairLoom refuses to overwrite it and fails
    loudly, suggesting a different ``--status-tag``.
    """
    tag = cfg.status_tag
    if record.has_tag(tag):
        existing = record.get_tag(tag)
        if isinstance(existing, str):
            if existing == classification:
                record.set_tag(tag, classification, value_type="Z")
                return
            if existing not in _KNOWN:
                raise ConfigurationError(
                    f"read {record.query_name!r} already carries tag {tag} with "
                    f"unrelated value {existing!r}; PairLoom does not silently "
                    f"overwrite it. Choose another --status-tag or remove the tag."
                )
        else:
            raise ConfigurationError(
                f"read {record.query_name!r} carries tag {tag} with non-string "
                f"value {existing!r}; refusing to overwrite. Choose another "
                "--status-tag."
            )
    record.set_tag(tag, classification, value_type="Z")
