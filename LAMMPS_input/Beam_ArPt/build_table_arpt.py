#!/usr/bin/env python3
"""Build the alpha_t, alpha_n lookup tables read by surf_collide cll/md.

Input is calibration_curve_ArPt.csv from extract_alpha.py --all, with columns
eps,theta,alpha_t,alpha_n,eac,err_t,err_n,events. The full 5x5 measured grid is
required: eps in {0.4,1,2,5,10} eV and theta in {15,30,45,60,75} deg. The 25
measured points enter the table as they are, so no separable form
alpha(eps,theta) = g(eps)h(theta) is assumed.

Padding of the measured grid:
  eps = 0.05 eV holds the 0.4 eV row, for re-emitted thermal molecules
  eps above 10 eV is held at the 10 eV values by the lookup
  theta = 0 holds the 15 deg column, alpha_t being undefined at normal incidence
  theta = 90 is either a linear extrapolation of the 60 and 75 deg values clipped
    to [0,1], or a hold of the 75 deg value. The difference between the two
    variants measures the uncertainty of that extrapolation.
  everything is clipped to [0,1], which the CL kernel requires. alpha_n is
  negative at some grazing conditions and is clipped to zero there, biasing the
  normal accommodation upward. Clipped points are reported one by one.

Writes alpha_arpt_table.txt (extrap) and alpha_arpt_table_hold.txt here, and
with --install copies them into SPARTA_input/tables/.

usage: python3 build_table_arpt.py [--csv F] [--outdir D] [--install]
"""
import argparse
import csv
import os
import shutil
from datetime import date

EPS_MEAS = [0.4, 1.0, 2.0, 5.0, 10.0]
TH_MEAS = [15.0, 30.0, 45.0, 60.0, 75.0]
EPS_PAD_LO = 0.05
TABLES_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "SPARTA_input", "tables"))
CANON_NAMES = {"extrap": "alpha_arpt_table.txt", "hold": "alpha_arpt_table_hold.txt"}


def read_calibration(path):
    """read the CSV into {(eps,theta): row} and check the 5x5 grid is complete"""
    data = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = (float(row["eps"]), float(row["theta"]))
            data[key] = {k: float(v) for k, v in row.items()}
    missing = [(e, t) for e in EPS_MEAS for t in TH_MEAS if (e, t) not in data]
    if missing:
        raise SystemExit(f"[FATAL] {len(missing)} grid points missing: {missing}\n"
                         f"        check that all 25 beam runs finished, then rerun "
                         f"extract_alpha.py --all")
    return data


def clip01(x, label, clips):
    if x < 0.0 or x > 1.0:
        clips.append(f"{label}: {x:+.4f} -> {min(max(x, 0.0), 1.0):.1f}")
    return min(max(x, 0.0), 1.0)


def build_grids(data):
    """the 25 measured points as 5x5 nested lists [ie][it], plus the clip report"""
    clips = []
    at = [[clip01(data[(e, t)]["alpha_t"], f"a_t({e}eV,{t:.0f}deg)", clips)
           for t in TH_MEAS] for e in EPS_MEAS]
    an = [[clip01(data[(e, t)]["alpha_n"], f"a_n({e}eV,{t:.0f}deg)", clips)
           for t in TH_MEAS] for e in EPS_MEAS]
    return at, an, clips


def pad_table(rows, grazing):
    """5x5 -> 6x7 with the 0.05 eV row, the 0 deg column and the 90 deg column"""
    out = []
    for r in [rows[0]] + rows:          # 0.05 eV row copies the 0.4 eV one
        if grazing == "extrap":
            g90 = r[-1] + (r[-1] - r[-2])  # 75 deg plus the 60 to 75 deg slope
            g90 = min(max(g90, 0.0), 1.0)
        else:
            g90 = r[-1]
        out.append([r[0]] + list(r) + [g90])   # 0 deg holds 15 deg
    return out


def write_table(path, eps_grid, th_grid, at, an, grazing, src):
    n_e, n_t = len(eps_grid), len(th_grid)
    with open(path, "w") as f:
        f.write(f"# {os.path.basename(path)}: Ar-Pt accommodation coefficients "
                f"alpha_t(eps,theta), alpha_n(eps,theta) [grazing={grazing}]\n")
        f.write(f"# from {os.path.basename(src)} by build_table_arpt.py, "
                f"{date.today().isoformat()}. LAMMPS monoenergetic beam, 2500 atoms per\n"
                f"# condition, coefficients flux-weighted over the returned atoms only\n")
        f.write("# 5 eps x 5 theta measured points, no separable form assumed. Padded with "
                "eps=0.05\n# (holds 0.4), theta=0 (holds 15) and theta=90 "
                "(extrapolated and clipped, or holds 75)\n")
        f.write("# format: n_eps n_theta / eps_grid[eV] / theta_grid[deg] / "
                "alpha_t[ie][ith] / alpha_n[ie][ith]\n")
        f.write(f"{n_e} {n_t}\n")
        f.write(" ".join(f"{x:.4f}" for x in eps_grid) + "\n")
        f.write(" ".join(f"{x:.1f}" for x in th_grid) + "\n")
        f.write("# alpha_t\n")
        for r in at:
            f.write(" ".join(f"{v:.5f}" for v in r) + "\n")
        f.write("# alpha_n\n")
        for r in an:
            f.write(" ".join(f"{v:.5f}" for v in r) + "\n")


def verify_tokens(path, n_e, n_t):
    """read the table the way cll/md does: 2 + n_e + n_t + 2*n_e*n_t tokens"""
    tokens = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tokens.extend(s.split())
    expect = 2 + n_e + n_t + 2 * n_e * n_t
    assert len(tokens) == expect, f"{path}: {len(tokens)} tokens, expected {expect}"
    assert int(tokens[0]) == n_e and int(tokens[1]) == n_t


def bilinear(eps, th, eps_grid, th_grid, tbl):
    """bilinear lookup with hold outside the grid, as in cll/md"""
    def frac(x, grid):
        if x <= grid[0]:
            return 0, 0.0
        if x >= grid[-1]:
            return len(grid) - 2, 1.0
        for i in range(len(grid) - 1):
            if x <= grid[i + 1]:
                return i, (x - grid[i]) / (grid[i + 1] - grid[i])
        return len(grid) - 2, 1.0
    i, fe = frac(eps, eps_grid)
    j, ft = frac(th, th_grid)
    return ((1 - fe) * (1 - ft) * tbl[i][j] + fe * (1 - ft) * tbl[i + 1][j]
            + (1 - fe) * ft * tbl[i][j + 1] + fe * ft * tbl[i + 1][j + 1])


def install(outdir, made):
    os.makedirs(TABLES_DIR, exist_ok=True)
    for grazing, src_path in made.items():
        dst = os.path.join(TABLES_DIR, CANON_NAMES[grazing])
        shutil.copy2(src_path, dst)
        print(f"[INSTALL] {dst}")


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--csv", default=os.path.join(here, "calibration_curve_ArPt.csv"))
    ap.add_argument("--outdir", default=here)
    ap.add_argument("--install", action="store_true",
                    help="copy the tables into SPARTA_input/tables/")
    a = ap.parse_args()

    data = read_calibration(a.csv)
    at5, an5, clips = build_grids(data)
    if clips:
        print(f"[CLIP] {len(clips)} values outside [0,1] were clipped:")
        for c in clips:
            print(f"       {c}")

    eps_grid = [EPS_PAD_LO] + EPS_MEAS
    th_grid = [0.0] + TH_MEAS + [90.0]
    made = {}
    for grazing in ("extrap", "hold"):
        at = pad_table(at5, grazing)
        an = pad_table(an5, grazing)
        path = os.path.join(a.outdir, CANON_NAMES[grazing])
        write_table(path, eps_grid, th_grid, at, an, grazing, a.csv)
        verify_tokens(path, len(eps_grid), len(th_grid))
        made[grazing] = path
        print(f"[OK] {path}  ({len(eps_grid)}x{len(th_grid)}, grazing={grazing}, tokens checked)")
        if grazing == "extrap":
            print("spot checks (bilinear lookup, hold outside the grid):")
            for e, t, lbl in [(5.0, 45.0, "5 eV, 45 deg (measured point)"),
                              (10.0, 45.0, "10 eV, 45 deg (measured point)"),
                              (13.0, 45.0, "13 eV, 45 deg (held at 10 eV)"),
                              (10.0, 83.7, "10 eV, 83.7 deg (grazing extrapolation)"),
                              (0.05, 45.0, "0.05 eV, 45 deg (held at 0.4 eV)")]:
                print(f"  {lbl:<34} alpha_t={bilinear(e, t, eps_grid, th_grid, at):.3f}"
                      f"  alpha_n={bilinear(e, t, eps_grid, th_grid, an):.3f}")

    if a.install:
        install(a.outdir, made)
    else:
        print("[NOTE] not installed, use --install to write into SPARTA_input/tables/")


if __name__ == "__main__":
    main()
