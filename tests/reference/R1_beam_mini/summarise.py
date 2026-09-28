#!/usr/bin/env python3
"""Small helpers for make.sh. Only reads its inputs and writes the named output.

  summarise.py dumpinfo DUMP
      print "sha256 size_bytes frames" of a LAMMPS text dump
  summarise.py gzip SRC DST MAX_BYTES
      gzip SRC to DST (level 9, mtime 0, no file name). If the result is not
      below MAX_BYTES, keep only the first N frames, with N the largest count
      that fits, and print how many frames were kept.
  summarise.py expected OUT.json [meta.KEY=VALUE ...] NAME=DIR [NAME=DIR ...]
      collect the dict printed by extract_alpha.py (last line of
      DIR/extract_alpha_stdout.txt) into OUT.json
"""
import ast
import gzip
import hashlib
import io
import json
import os
import sys

FRAME_TAG = b"ITEM: TIMESTEP"
EXPECTED_KEYS = ("alpha_t", "alpha_n", "eac", "err_t", "err_n", "events", "eps", "theta")
RUN_KEYS = ("V_INC", "THETA", "PHI", "T_WALL", "n_insert", "n_loops", "TAILSTEPS")


def frames_of(data):
    """Byte offsets where each frame starts."""
    starts, pos = [], 0
    while True:
        pos = data.find(FRAME_TAG, pos)
        if pos < 0:
            return starts
        if pos == 0 or data[pos - 1:pos] == b"\n":
            starts.append(pos)
        pos += len(FRAME_TAG)


def gz_bytes(data):
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", compresslevel=9, fileobj=buf, mtime=0) as g:
        g.write(data)
    return buf.getvalue()


def cmd_dumpinfo(path):
    with open(path, "rb") as f:
        data = f.read()
    print(hashlib.sha256(data).hexdigest(), len(data), len(frames_of(data)))


def cmd_gzip(src, dst, max_bytes):
    with open(src, "rb") as f:
        data = f.read()
    starts = frames_of(data)
    total = len(starts)
    packed = gz_bytes(data)
    kept = total
    if len(packed) >= max_bytes:
        lo, hi = 1, total - 1          # largest frame count whose gzip fits
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if len(gz_bytes(data[:starts[mid]])) < max_bytes:
                lo = mid
            else:
                hi = mid - 1
        kept = lo
        data = data[:starts[kept]]
        packed = gz_bytes(data)
    with open(dst, "wb") as f:
        f.write(packed)
    print(f"frames_kept {kept} of {total}; uncompressed_sha256 "
          f"{hashlib.sha256(data).hexdigest()} uncompressed_bytes {len(data)} gz_bytes {len(packed)}")


def run_vars(command_file):
    """The -var NAME VALUE pairs of the stored LAMMPS command line."""
    with open(command_file) as f:
        words = f.read().split()
    return {words[i + 1]: words[i + 2] for i, w in enumerate(words[:-2]) if w == "-var"}


def cmd_expected(out, pairs):
    """NAME=DIR adds a condition; meta.KEY=VALUE adds a string to "meta"."""
    meta, conditions = {}, {}
    for pair in pairs:
        name, value = pair.split("=", 1)
        if name.startswith("meta."):
            meta[name[len("meta."):]] = value
            continue
        with open(os.path.join(value, "extract_alpha_stdout.txt")) as f:
            last = f.read().strip().splitlines()[-1]
        row = ast.literal_eval(last)
        entry = {k: row[k] for k in EXPECTED_KEYS}
        lvars = run_vars(os.path.join(value, "command.txt"))
        entry.update({k: int(lvars[k]) for k in RUN_KEYS})
        conditions[name] = entry
    with open(out, "w") as f:
        json.dump({"meta": meta, "conditions": conditions}, f, indent=2)
        f.write("\n")


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    cmd = sys.argv[1]
    if cmd == "dumpinfo":
        cmd_dumpinfo(sys.argv[2])
    elif cmd == "gzip" and len(sys.argv) == 5:
        cmd_gzip(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    elif cmd == "expected" and len(sys.argv) >= 4:
        cmd_expected(sys.argv[2], sys.argv[3:])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
