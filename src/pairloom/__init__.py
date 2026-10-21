"""PairLoom: pair-aware BAM realignment, overlap analysis, and read merging."""

from pairloom.config import PairLoomConfig
from pairloom.exceptions import (
    ConfigurationError,
    DataInvariantError,
    InputError,
    OutputError,
    PairLoomError,
    WorkerError,
)

__version__ = "0.1.0"

__all__ = [
    "PairLoomConfig",
    "PairLoomError",
    "ConfigurationError",
    "InputError",
    "DataInvariantError",
    "OutputError",
    "WorkerError",
    "__version__",
]
