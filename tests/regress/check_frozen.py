#!/usr/bin/env python3
"""Check that the files used by the published papers are unchanged.

tests/reference/FROZEN.sha256 lists one "<sha256>  <path>" pair per line,
the same format as `shasum -a 256`. Every listed file must exist and match.

    python3 tests/regress/check_frozen.py            # check
    python3 tests/regress/check_frozen.py --write    # rebuild the list (see below)

Rebuilding the list is only allowed together with an entry in
tests/reference/CHANGES.md that says which file changed and why.
"""

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIST = ROOT / "tests" / "reference" / "FROZEN.sha256"

# Directories whose tracked files are frozen: the inputs, sources and tables
# behind the published papers.
FROZEN_DIRS = (
    "LAMMPS_input",
    "LAMMPS_src",
    "SPARTA_input",
    "SPARTA_src",
    "tools",
    "examples",
    "initialisation",
)
FROZEN_FILES = ("FlowChart.png", "LICENSE.md")

# LaTeX build files and editor settings are clutter, not paper inputs.
NOT_FROZEN_SUFFIXES = (
    ".aux", ".bbl", ".blg", ".fdb_latexmk", ".fls", ".log", ".synctex.gz",
    "-eps-converted-to.pdf",
)
NOT_FROZEN_PARTS = (".vscode",)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tracked_files():
    out = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "-z"],
        check=True, capture_output=True,
    ).stdout.decode()
    return [p for p in out.split("\0") if p]


def is_frozen(path):
    if any(part in NOT_FROZEN_PARTS for part in Path(path).parts):
        return False
    if path.endswith(NOT_FROZEN_SUFFIXES):
        return False
    return path in FROZEN_FILES or path.split("/", 1)[0] in FROZEN_DIRS


def write_list():
    paths = sorted(p for p in tracked_files() if is_frozen(p))
    lines = [f"{sha256(ROOT / p)}  {p}\n" for p in paths]
    LIST.parent.mkdir(parents=True, exist_ok=True)
    LIST.write_text("".join(lines))
    print(f"wrote {len(lines)} entries to {LIST.relative_to(ROOT)}")
    return 0


def check_list():
    if not LIST.exists():
        print(f"missing {LIST.relative_to(ROOT)}", file=sys.stderr)
        return 1
    bad = []
    n = 0
    for line in LIST.read_text().splitlines():
        if not line.strip():
            continue
        digest, path = line.split("  ", 1)
        n += 1
        target = ROOT / path
        if not target.exists():
            bad.append(f"missing  {path}")
        elif sha256(target) != digest:
            bad.append(f"changed  {path}")
    listed = {line.split("  ", 1)[1] for line in LIST.read_text().splitlines() if line.strip()}
    added = [p for p in tracked_files() if is_frozen(p) and p not in listed]
    bad += [f"new file in a frozen directory, not in the list: {p}" for p in added]
    for msg in bad:
        print(msg, file=sys.stderr)
    if bad:
        print(f"{len(bad)} problems with the {n} frozen files", file=sys.stderr)
        return 1
    print(f"all {n} frozen files unchanged, none added")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true",
                        help="rebuild FROZEN.sha256 from the tracked files")
    args = parser.parse_args()
    return write_list() if args.write else check_list()


if __name__ == "__main__":
    sys.exit(main())
