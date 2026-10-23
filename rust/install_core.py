#!/usr/bin/env python
"""Build the Rust core and install it into the active Python environment.

Run via `pixi run core` from the repository root. Requires cargo on PATH.
"""

from __future__ import annotations

import shutil
import subprocess
import sysconfig
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> None:
    subprocess.run(
        ["cargo", "build", "--release", "--manifest-path", str(ROOT / "Cargo.toml")],
        check=True,
    )
    purelib = Path(sysconfig.get_paths()["purelib"])
    dest = purelib / "pairloom_core"
    dest.mkdir(exist_ok=True)
    suffix = sysconfig.get_config_var("EXT_SUFFIX")
    shutil.copy(
        ROOT / "target/release/libpairloom_core.so",
        dest / f"pairloom_core{suffix}",
    )
    (dest / "__init__.py").write_text(
        "from .pairloom_core import *\nfrom .pairloom_core import __doc__\n"
    )
    print(f"installed pairloom_core -> {dest}")


if __name__ == "__main__":
    main()
