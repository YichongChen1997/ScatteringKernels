#!/usr/bin/env python3
"""Require an entry in tests/reference/CHANGES.md for any change to the references.

Compares the current commit with BASE (a branch, tag or commit). If any file
under tests/reference/ changed, including FROZEN.sha256, then
tests/reference/CHANGES.md must have changed too.

    python3 tests/regress/check_reference_changes.py origin/main
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF_DIR = "tests/reference/"
CHANGES = "tests/reference/CHANGES.md"


def changed_files(base):
    out = subprocess.run(
        ["git", "-C", str(ROOT), "diff", "--name-only", f"{base}...HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout
    return [p for p in out.splitlines() if p]


def main():
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    files = changed_files(sys.argv[1])
    touched = [p for p in files if p.startswith(REF_DIR) and p != CHANGES]
    if touched and CHANGES not in files:
        print(f"{len(touched)} reference files changed but {CHANGES} did not:", file=sys.stderr)
        for p in touched:
            print(f"  {p}", file=sys.stderr)
        return 1
    print(f"{len(touched)} reference files changed; CHANGES.md "
          f"{'updated' if CHANGES in files else 'not needed'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
