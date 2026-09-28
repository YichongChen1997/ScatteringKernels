#!/usr/bin/env python3
"""Crop the Ar-Pt beam slab laterally to NCELL x NCELL fcc unit cells.

usage: python3 crop_slab.py SRC DST [NCELL] [A]

SRC is LAMMPS_input/Beam_ArPt/ArPt_slab.data from v2.0.0 (103 x 103 cells,
169744 Pt atoms, atom_style full). Every atom with 0 <= x-xlo < NCELL*A and
0 <= y-ylo < NCELL*A is kept, in all 8 layers, in source order. Ids are
renumbered 1..N and the molecule id is set equal to the atom id, as in the
source. Atom types, charges, coordinates and velocities are copied as text,
so no number is reformatted. Masses and the z bounds come from the source
header. The lattice constant A = 3.92 Angstrom keeps the fcc periodicity
across the new periodic x and y boundaries.

The output depends only on the bytes of SRC and on the arguments.
"""
import sys

DEFAULT_NCELL = 10
DEFAULT_A = 3.92
X_COL, Y_COL = 4, 5   # atom_style full: id mol type q x y z


def section(lines, name):
    """Index of the first line whose first word is `name`."""
    for i, line in enumerate(lines):
        words = line.split()
        if words and words[0] == name:
            return i
    raise SystemExit(f"crop_slab: no '{name}' section in the source")


def body(lines, start):
    """Non-blank lines of the section whose title is on line `start`."""
    out = []
    for line in lines[start + 2:]:
        if not line.strip():
            break
        out.append(line)
    return out


def header_value(lines, key):
    for line in lines:
        if line.rstrip().endswith(key):
            return line.split()
    raise SystemExit(f"crop_slab: no '{key}' line in the source header")


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    ncell = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_NCELL
    a = float(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_A
    length = ncell * a

    with open(src) as f:
        lines = f.read().splitlines()
    ia, iv, im = section(lines, "Atoms"), section(lines, "Velocities"), section(lines, "Masses")
    atoms = [l.split() for l in body(lines, ia)]
    vels = {l.split()[0]: l.split()[1:] for l in body(lines, iv)}
    masses = body(lines, im)
    ntypes = header_value(lines, "atom types")[0]
    xlo = float(header_value(lines, "xlo xhi")[0])
    ylo = float(header_value(lines, "ylo yhi")[0])
    zlo, zhi = header_value(lines, "zlo zhi")[:2]

    keep = [t for t in atoms
            if float(t[X_COL]) - xlo < length and float(t[Y_COL]) - ylo < length]
    out = [f"LAMMPS data file: ArPt_slab.data (v2.0.0) cropped to {ncell}x{ncell} cells by crop_slab.py",
           "", f"{len(keep)}\tatoms", "", f"{ntypes}\tatom types", "",
           f"{xlo:g}\t{xlo + length:.2f}\txlo xhi",
           f"{ylo:g}\t{ylo + length:.2f}\tylo yhi",
           f"{zlo}\t{zhi}\tzlo zhi", "", "Masses", ""]
    out += masses
    out += ["", "Atoms", ""]
    for k, t in enumerate(keep, start=1):
        out.append("\t".join([str(k), str(k)] + t[2:]))
    out += ["", "Velocities", ""]
    for k, t in enumerate(keep, start=1):
        out.append("\t".join([str(k)] + vels[t[0]]))
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"{len(keep)} atoms, L = {length:.2f} Angstrom, written to {dst}")


if __name__ == "__main__":
    main()
