#!/usr/bin/env python3
"""Compare the outputs of the legacy tools with a stored reference set.

    compare_outputs.py REF_OUTPUTS NEW_OUTPUTS [--exact] [--rtol 1e-6]
                       [--atol 0] [--skip tool,tool,...]

REF_OUTPUTS is a stored outputs/ directory, for example
tests/reference/R3_channel_mini/mini/outputs, with one subdirectory per tool.
NEW_OUTPUTS has the same layout; tests/regress/run_legacy_tools.sh writes it
as WORKDIR/tools.

Every stored reference file must exist in NEW_OUTPUTS and match:

  --exact   byte for byte.
  default   line by line and word by word. Words that are not numbers must be
            equal; numbers a and b must satisfy
            |a - b| <= atol + rtol * max(|a|, |b|).
            nan and -nan count as equal; inf must have the same sign.

Files listed in outputs/NOT_STORED.txt (too large to store) are checked by
their sha256 when they were produced, and reported as skipped otherwise.
Without --exact a different sha256 of such a file is only a warning, since
it cannot be compared with a tolerance.

Exit codes are compared with the "exit_code" column of the tool_status.tsv
next to each directory (override with --ref-status / --new-status). Tools
that the runner did not run (exit code NA, note "skipped ...") are left out,
like those named with --skip.

A few outputs depend on undefined behaviour in the legacy code and change
with the compiler (UNDEFINED_COLUMNS and UNDEFINED_EXIT below; known issues
KI-4 and KI-7 in docs/legacy_tools.md). Without --exact a difference there
is only a warning.

Files in NEW_OUTPUTS that the reference does not have are listed as extra;
with --exact they count as a mismatch. Symbolic links (the input links the
runner puts in each tool directory) are ignored.

Exit status: 0 if everything matches, 1 on any mismatch, 2 on bad usage.
"""

import argparse
import hashlib
import math
import os
import re
import sys

# A number: nan/inf (not inside a word) or a decimal with optional exponent.
NUMBER = re.compile(
    r"((?<![A-Za-z_])[-+]?(?:nan|inf(?:inity)?)(?![A-Za-z_])"
    r"|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)",
    re.IGNORECASE,
)
MAX_DETAILS = 3  # differences shown per file

# Outputs that depend on undefined behaviour in the legacy tools, so they
# can differ between compilers. Only checked strictly with --exact.
UNDEFINED_COLUMNS = {
    # word numbers (from 1) on every line of the file
    "pp_meas_GasOverTime/gasOverTime.txt":
        ({6}, "KI-4: pressure from the uninitialised stressTemp"),
}
UNDEFINED_EXIT = {
    # crash codes seen: 139 = segmentation fault, 134 = abort from the
    # C++ library's own bounds check (GCC 15 and later at -O0)
    "pp_meas_Traj": ({"134", "139"}, "KI-7: reads past the end of a vector"),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def regular_files(top):
    """Relative paths of the regular files under top (links left out)."""
    found = []
    for root, dirs, files in os.walk(top):
        dirs.sort()
        for name in sorted(files):
            path = os.path.join(root, name)
            if os.path.islink(path) or not os.path.isfile(path):
                continue
            found.append(os.path.relpath(path, top))
    return found


def read_not_stored(ref_outputs):
    """{path relative to outputs/: (sha256, size)} from NOT_STORED.txt."""
    path = os.path.join(ref_outputs, "NOT_STORED.txt")
    entries = {}
    if not os.path.isfile(path):
        return entries
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            digest, size, rel = line.split(None, 2)
            entries[rel.strip()] = (digest, int(size))
    return entries


def read_status(path):
    """{tool: (exit code, note)} from a tool_status.tsv, or None."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    header = rows[0]
    for name in ("exit_code", "exit"):
        if name in header:
            col = header.index(name)
            break
    else:
        raise SystemExit(f"{path}: no exit_code column in {header}")
    note = header.index("note") if "note" in header else None
    status = {}
    for row in rows[1:]:
        if len(row) > col:
            text = row[note] if note is not None and len(row) > note else ""
            status[row[0]] = (row[col], text)
    return status


def to_float(text):
    try:
        return float(text)
    except ValueError:
        return None


class Diff:
    """Largest difference seen between two numbers that were not equal."""

    def __init__(self):
        self.abs = 0.0
        self.rel = 0.0

    def add(self, a, b):
        if math.isfinite(a) and math.isfinite(b):
            d = abs(a - b)
            self.abs = max(self.abs, d)
            scale = max(abs(a), abs(b))
            if scale > 0:
                self.rel = max(self.rel, d / scale)


def numbers_close(a, b, rtol, atol):
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def compare_word(w_ref, w_new, rtol, atol, diff):
    """None if the two words agree, else a short reason."""
    if w_ref == w_new:
        return None
    p_ref = NUMBER.split(w_ref)
    p_new = NUMBER.split(w_new)
    if len(p_ref) != len(p_new):
        return f"'{w_ref}' vs '{w_new}'"
    for i, (a, b) in enumerate(zip(p_ref, p_new)):
        if i % 2 == 0:  # text between numbers
            if a != b:
                return f"'{w_ref}' vs '{w_new}'"
            continue
        x, y = to_float(a), to_float(b)
        if x is None or y is None:
            if a != b:
                return f"'{w_ref}' vs '{w_new}'"
            continue
        diff.add(x, y)
        if not numbers_close(x, y, rtol, atol):
            return f"{a} vs {b}"
    return None


def compare_text(ref_path, new_path, rtol, atol, loose_words=()):
    """Returns (differences, differences in loose_words only, Diff)."""
    with open(ref_path, "rb") as f:
        ref_lines = f.read().decode("latin-1").splitlines()
    with open(new_path, "rb") as f:
        new_lines = f.read().decode("latin-1").splitlines()
    problems = []
    loose = []
    diff = Diff()
    if len(ref_lines) != len(new_lines):
        problems.append(f"{len(ref_lines)} lines in the reference, {len(new_lines)} here")
    for n, (l_ref, l_new) in enumerate(zip(ref_lines, new_lines), start=1):
        if l_ref == l_new:
            continue
        w_ref, w_new = l_ref.split(), l_new.split()
        if len(w_ref) != len(w_new):
            problems.append(f"line {n}: {len(w_ref)} words vs {len(w_new)}")
            continue
        for k, (a, b) in enumerate(zip(w_ref, w_new), start=1):
            if k in loose_words:
                if compare_word(a, b, rtol, atol, Diff()):
                    loose.append(n)
                continue
            reason = compare_word(a, b, rtol, atol, diff)
            if reason:
                problems.append(f"line {n} word {k}: {reason}")
                break
    return problems, loose, diff


class Report:
    def __init__(self, exact):
        self.exact = exact
        self.failures = []
        self.warnings = []
        self.close = []
        self.skipped = []
        self.n_identical = 0

    def fail(self, msg):
        self.failures.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)


def compare_tool(tool, args, not_stored, report):
    ref_dir = os.path.join(args.ref_outputs, tool)
    new_dir = os.path.join(args.new_outputs, tool)
    ref_files = regular_files(ref_dir)
    big = {rel.split("/", 1)[1]: v for rel, v in not_stored.items()
           if rel.split("/", 1)[0] == tool}
    if not os.path.isdir(new_dir):
        report.fail(f"{tool}: no output directory {new_dir}")
        return
    for rel in ref_files:
        name = f"{tool}/{rel}"
        ref_path = os.path.join(ref_dir, rel)
        new_path = os.path.join(new_dir, rel)
        if not os.path.isfile(new_path):
            report.fail(f"missing   {name}")
            continue
        with open(ref_path, "rb") as f1, open(new_path, "rb") as f2:
            same = f1.read() == f2.read()
        if same:
            report.n_identical += 1
            continue
        if args.exact:
            report.fail(f"differs   {name} (bytes)")
            continue
        words, why = UNDEFINED_COLUMNS.get(name, (set(), ""))
        problems, loose, diff = compare_text(ref_path, new_path, args.rtol, args.atol, words)
        if loose:
            report.warn(f"{name}: word {sorted(words)} differs on {len(loose)} lines ({why})")
        if problems:
            shown = "; ".join(problems[:MAX_DETAILS])
            more = f" (+{len(problems) - MAX_DETAILS} more)" if len(problems) > MAX_DETAILS else ""
            report.fail(f"differs   {name}: {shown}{more}")
        else:
            apart = f"apart from word {sorted(words)}, " if loose else ""
            report.close.append(
                f"close     {name} ({apart}max abs diff {diff.abs:.3g}, "
                f"max rel diff {diff.rel:.3g})")
    for rel, (digest, size) in sorted(big.items()):
        name = f"{tool}/{rel}"
        new_path = os.path.join(new_dir, rel)
        if not os.path.isfile(new_path):
            report.skipped.append(f"skipped   {name} (hash only, not produced)")
            continue
        new_size = os.path.getsize(new_path)
        if new_size == size and sha256(new_path) == digest:
            report.n_identical += 1
            continue
        msg = f"{name} (hash only): sha256/size differ ({new_size} vs {size} bytes)"
        if args.exact:
            report.fail(f"differs   {msg}")
        else:
            report.warn(f"not compared {msg}; cannot use a tolerance on it")
    known = set(ref_files) | set(big)
    for rel in regular_files(new_dir):
        if rel not in known:
            msg = f"extra     {tool}/{rel}"
            (report.fail if args.exact else report.warn)(msg)


def compare_status(tools, ref, new, report, exact):
    """Compares exit codes; ref and new come from read_status."""
    if ref is None:
        report.warn("no reference tool_status.tsv; exit codes not compared")
        return
    if new is None:
        report.fail("no tool_status.tsv next to the new outputs")
        return
    for tool in tools:
        if tool not in ref:
            report.warn(f"{tool}: not in the reference tool_status.tsv")
            continue
        if tool not in new:
            report.fail(f"exit code {tool}: not run")
            continue
        ref_code, new_code = ref[tool][0], new[tool][0]
        if ref_code == new_code:
            continue
        codes, why = UNDEFINED_EXIT.get(tool, (set(), ""))
        msg = f"exit code {tool}: {new_code}, reference {ref_code}"
        if not exact and ref_code in codes and new_code in codes:
            report.warn(f"{msg} ({why})")
        else:
            report.fail(msg)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        epilog="See the top of this file for details.")
    parser.add_argument("ref_outputs")
    parser.add_argument("new_outputs")
    parser.add_argument("--exact", action="store_true", help="compare bytes")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=0.0,
                        help="absolute tolerance; keep 0 for SI outputs of order 1e-13 or less")
    parser.add_argument("--skip", default="", help="comma-separated tools to leave out")
    parser.add_argument("--ref-status", help="reference tool_status.tsv")
    parser.add_argument("--new-status", help="tool_status.tsv of the new run")
    args = parser.parse_args()

    for d in (args.ref_outputs, args.new_outputs):
        if not os.path.isdir(d):
            parser.error(f"not a directory: {d}")
    ref_status = read_status(args.ref_status or os.path.join(
        os.path.dirname(os.path.abspath(args.ref_outputs)), "tool_status.tsv"))
    new_status = read_status(args.new_status or os.path.join(
        os.path.dirname(os.path.abspath(args.new_outputs)), "tool_status.tsv"))

    # Tools left out: by --skip, or not run by the runner (exit code NA,
    # note "skipped...").
    skip = {t.strip(): "--skip" for t in args.skip.split(",") if t.strip()}
    for tool, (code, note) in (new_status or {}).items():
        if code == "NA" and note.startswith("skipped") and tool not in skip:
            skip[tool] = note[len("skipped"):].lstrip(": ") or "not run"
    tools = sorted(t for t in os.listdir(args.ref_outputs)
                   if os.path.isdir(os.path.join(args.ref_outputs, t)))
    not_stored = read_not_stored(args.ref_outputs)
    report = Report(args.exact)
    compared = [t for t in tools if t not in skip]
    # Every tool in the reference tool_status.tsv must have stored outputs,
    # so that a deleted reference directory cannot pass unnoticed.
    for tool in sorted(set(ref_status or {}) - set(tools) - set(skip)):
        report.fail(f"reference outputs of {tool} missing from {args.ref_outputs}")
    for tool in sorted(set(skip) & set(tools)):
        report.skipped.append(f"skipped   {tool} ({skip[tool]})")
    for tool in compared:
        compare_tool(tool, args, not_stored, report)
    new_tools = {t for t in os.listdir(args.new_outputs)
                 if os.path.isdir(os.path.join(args.new_outputs, t))}
    for tool in sorted(new_tools - set(tools)):
        report.warn(f"extra tool directory {tool} (not in the reference)")
    compare_status(compared, ref_status, new_status, report, args.exact)

    mode = "byte for byte" if args.exact else f"rtol {args.rtol:g}, atol {args.atol:g}"
    print(f"compare {args.new_outputs} with {args.ref_outputs} ({mode})")
    for line in report.close + report.skipped:
        print("  " + line)
    for line in report.warnings:
        print("  warning: " + line)
    for line in report.failures:
        print("  FAIL " + line)
    print(f"{len(compared)} tools: {report.n_identical} files identical, "
          f"{len(report.close)} within tolerance, {len(report.failures)} problems, "
          f"{len(report.skipped)} skipped")
    return 1 if report.failures else 0


if __name__ == "__main__":
    sys.exit(main())
