"""Conservative MAPQ policy after alignment adjustment.

Documented rules (README "MAPQ policy"), applied per read:

============ =====================================================
Condition    Adjusted MAPQ
============ =====================================================
ambiguous    (>= 2 equivalent-best placements remain) min(orig, 3)
truncated    (candidate discovery hit its cap)         min(orig, 10)
unchanged    placement identical to the BAM record     orig
exact +      changed placement, verified exact match,  60
corroborated original MAPQ >= 30
otherwise    changed placement without corroboration   min(orig, 25)
============ =====================================================

MAPQ 60 is therefore never assigned merely because an alignment is an exact
match: the read must also be unambiguous in the searched window, discovery
must not have been truncated, and the original aligner must have considered
the read high quality. The ``--force-mapq60-exact`` compatibility flag
overrides this for exact matches and is reported with an explicit warning
about overstating mapping confidence.
"""

from __future__ import annotations

MAPQ_AMBIGUOUS_CAP = 3
MAPQ_TRUNCATED_CAP = 10
MAPQ_CHANGED_CAP = 25
CORROBORATING_ORIG_MAPQ = 30


def adjusted_mapq(
    orig_mapq: int,
    n_equiv: int,
    truncated: bool,
    exact: bool,
    placement_changed: bool,
    force_mapq60_exact: bool = False,
) -> int:
    if force_mapq60_exact and exact:
        return 60
    if n_equiv > 1:
        return min(orig_mapq, MAPQ_AMBIGUOUS_CAP)
    if truncated:
        return min(orig_mapq, MAPQ_TRUNCATED_CAP)
    if not placement_changed:
        return orig_mapq
    if exact and orig_mapq >= CORROBORATING_ORIG_MAPQ:
        return 60
    return min(orig_mapq, MAPQ_CHANGED_CAP)
