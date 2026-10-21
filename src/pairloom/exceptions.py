"""Exception hierarchy for PairLoom.

Every failure path in the package raises one of these exceptions.  Nothing is
silently swallowed: unexpected conditions terminate processing with a nonzero
exit code and an actionable message (see the fail-loud policy in the README).
"""

from __future__ import annotations


class PairLoomError(Exception):
    """Base class for all PairLoom errors."""


class ConfigurationError(PairLoomError):
    """Invalid CLI arguments or inconsistent run configuration."""


class InputError(PairLoomError):
    """Unreadable, missing, or corrupt input data (BAM, FASTA)."""


class DataInvariantError(PairLoomError):
    """A record or intermediate structure violates a documented invariant.

    Examples: malformed CIGARs, templates with more than two primary records,
    reference-coordinate arithmetic that falls outside the reference.
    """


class OutputError(PairLoomError):
    """A required output could not be written, validated, or published."""


class WorkerError(PairLoomError):
    """A multiprocessing worker failed; the whole run is aborted."""
