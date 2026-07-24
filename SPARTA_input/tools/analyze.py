#!/usr/bin/env python3
"""Post-process runs/*/surf_force.dump.
  cylinder: C_D and the stagnation heat flux q_stag
  plate:    C_L, C_D and C_m about the leading edge, from per-element forces

In the compute surf output fx and fy are forces in N per metre of span in 2D,
while etot is a net energy flux per unit area in W/m^2. The heat flux is etot
itself and the integrated power is sum(etot * L_elem).

usage:
  python3 tools/analyze.py runs/cyl_kn10_diffuse
  python3 tools/analyze.py --all           # scan runs/, write summary.tsv
"""
import argparse
import glob
import math
import os
import re
import sys

M_AR = 6.6335e-26
D_REF = 0.1          # reference length: cylinder diameter = plate chord
U_INF = 7000.0
LE = (-0.0, 0.0)     # plate leading edge stays at the origin after rotation


def read_dump(path):
    """parse a 2D SPARTA surf dump: [(id,x1,y1,x2,y2,fx,fy,etot)]"""
    rows, insurf = [], False
    with open(path) as f:
        for line in f:
            if line.startswith("ITEM: SURFS"):
                insurf = True
                continue
            if line.startswith("ITEM:"):
                insurf = False
                continue
            if insurf and line.strip():
                v = line.split()
                rows.append(tuple(float(x) for x in v[:8]))
    return rows


def nrho_from_cmd(rdir):
    with open(os.path.join(rdir, "cmd.sh")) as f:
        m = re.search(r"-var NRHO ([0-9.eE+-]+)", f.read())
    return float(m.group(1))


def analyze(rdir, verbose=True):
    dump = os.path.join(rdir, "surf_force.dump")
    if not os.path.exists(dump):
        return None
    rows = read_dump(dump)
    if not rows:
        return None
    nrho = nrho_from_cmd(rdir)
    rho = nrho * M_AR
    qref = 0.5 * rho * U_INF**2
    Fx = sum(r[5] for r in rows)
    Fy = sum(r[6] for r in rows)
    Q = sum(r[7] * math.hypot(r[3] - r[1], r[4] - r[2]) for r in rows)
    cd = Fx / (qref * D_REF)
    cl = Fy / (qref * D_REF)
    # moment about the LE, or the centre for the cylinder: M = sum(xc*fy - yc*fx)
    M = sum((0.5 * (r[1] + r[3]) - LE[0]) * r[6] - (0.5 * (r[2] + r[4]) - LE[1]) * r[5]
            for r in rows)
    cm = M / (qref * D_REF**2)
    # stagnation element = smallest x
    stag = min(rows, key=lambda r: 0.5 * (r[1] + r[3]))
    q_stag = stag[7]
    res = dict(case=os.path.basename(rdir), CD=cd, CL=cl, Cm=cm,
               Qtot_W_per_m=Q, q_stag_W_m2=q_stag, n_elem=len(rows))
    if verbose:
        print(f"{res['case']:28s} CD={cd:8.4f} CL={cl:+8.4f} Cm={cm:+8.4f} "
              f"q_stag={q_stag:11.4g} W/m^2 Qtot={Q:10.4g} W/m")
        print(f"  sum(fx)={Fx:.6g} N/m")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", nargs="?", help="a single run directory")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if args.all:
        out = ["case\tCD\tCL\tCm\tq_stag_W_m2\tQtot_W_per_m"]
        for rdir in sorted(glob.glob(os.path.join(here, "runs", "*"))):
            if not os.path.isdir(rdir):
                continue
            r = analyze(rdir, verbose=True)
            if r:
                out.append(f"{r['case']}\t{r['CD']:.5g}\t{r['CL']:.5g}\t"
                           f"{r['Cm']:.5g}\t{r['q_stag_W_m2']:.5g}\t{r['Qtot_W_per_m']:.5g}")
        path = os.path.join(here, "runs", "summary.tsv")
        with open(path, "w") as f:
            f.write("\n".join(out) + "\n")
        print("\nwrote", path)
    elif args.rundir:
        analyze(args.rundir)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
