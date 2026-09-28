#!/usr/bin/env python3
"""Keep the public repository free of large data and model files.

Fails if a tracked file is larger than MAX_BYTES (apart from the few large
inputs that were already published), or if a tracked file looks like a
potential model or an uncompressed trajectory.

    python3 tests/regress/check_hygiene.py
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MAX_BYTES = 2 * 1024 * 1024

# Inputs larger than MAX_BYTES that were already in v2.0.0.
ALLOWED_LARGE = {
    "LAMMPS_input/Beam_ArPt/ArPt_slab.data",
    "initialisation/SurfGenerator/SpectralCharacterisation/Figures/fig1.eps",
}

# Potential models and raw trajectories never belong in the repository.
FORBIDDEN_SUFFIXES = (".pb", ".pth", ".pt", ".lammpstrj", ".restart")


def tracked_files():
    out = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "-z"],
        check=True, capture_output=True,
    ).stdout.decode()
    return [p for p in out.split("\0") if p]


def main():
    problems = []
    for path in tracked_files():
        target = ROOT / path
        if not target.exists():
            continue
        if path.endswith(FORBIDDEN_SUFFIXES):
            problems.append(f"forbidden file type: {path}")
        size = target.stat().st_size
        if size > MAX_BYTES and path not in ALLOWED_LARGE:
            problems.append(f"too large ({size / 1e6:.1f} MB): {path}")
    for msg in problems:
        print(msg, file=sys.stderr)
    if problems:
        return 1
    print("no large or forbidden files tracked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
