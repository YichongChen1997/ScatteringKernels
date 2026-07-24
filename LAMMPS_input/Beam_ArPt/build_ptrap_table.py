#!/usr/bin/env python3
"""Build the trapping probability table p_trap(eps,theta) read by cll/md/trap.

Input is calibration_curve_ArPt.csv from extract_alpha.py --all, where the events
column holds N_return out of the N_INSERT = 2500 atoms of each condition, so
p_trap = 1 - N_return/2500 as in build_event_library.py.

Where every atom is trapped, N_return = 0 and p_trap = 1. The alpha table carries
the complete accommodation limit at those points, but cll/md/trap never uses it
there because the explicit trapping branch takes over first. The returned-only
coefficients act on the untrapped part of the flux alone.

Padding follows build_table_arpt.py: eps = 0.05 eV holds the 0.4 eV row,
theta = 0 holds the 15 deg column, theta = 90 is either a linear extrapolation of
the 60 and 75 deg values or a hold of the 75 deg value, and everything is clipped
to [0,1].

Writes ptrap_arpt_table.txt (extrap) and ptrap_arpt_table_hold.txt here, and with
--install copies them into SPARTA_input/tables/. Table format:
  n_eps n_theta / eps_grid[eV] / theta_grid[deg] / ptrap[ie][ith]

usage: python3 build_ptrap_table.py [--csv F] [--outdir D] [--install]
"""
import argparse
import csv
import os
import shutil
from datetime import date

EPS_MEAS = [0.4, 1.0, 2.0, 5.0, 10.0]
TH_MEAS = [15.0, 30.0, 45.0, 60.0, 75.0]
EPS_PAD_LO = 0.05
N_INSERT = 2500
TABLES_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "SPARTA_input", "tables"))
NAMES = {"extrap": "ptrap_arpt_table.txt", "hold": "ptrap_arpt_table_hold.txt"}


def read_calibration(path: str) -> dict:
    """read the CSV into {(eps,theta): N_return} and check the 5x5 grid is complete"""
    data = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = (float(row["eps"]), float(row["theta"]))
            data[key] = int(float(row["events"]))
    missing = [(e, t) for e in EPS_MEAS for t in TH_MEAS if (e, t) not in data]
    if missing:
        raise SystemExit(f"[FATAL] {len(missing)} grid points missing: {missing}")
    return data


def build_grid(data: dict) -> list:
    """p_trap on the 25 measured points as a 5x5 nested list [ie][it]"""
    return [[1.0 - data[(e, t)] / N_INSERT for t in TH_MEAS] for e in EPS_MEAS]


def pad_table(rows: list, grazing: str) -> list:
    """5x5 -> 6x7 with the 0.05 eV row, the 0 deg column and the 90 deg column"""
    out = []
    for r in [rows[0]] + rows:
        if grazing == "extrap":
            g90 = r[-1] + (r[-1] - r[-2])
            g90 = min(max(g90, 0.0), 1.0)
        else:
            g90 = r[-1]
        out.append([r[0]] + list(r) + [g90])
    return out


def write_table(path: str, eps_grid: list, th_grid: list, pt: list,
                grazing: str, src: str) -> None:
    n_e, n_t = len(eps_grid), len(th_grid)
    with open(path, "w") as f:
        f.write(f"# {os.path.basename(path)}: Ar-Pt trapping probability "
                f"p_trap(eps,theta) [grazing={grazing}]\n")
        f.write(f"# from {os.path.basename(src)} by build_ptrap_table.py, "
                f"{date.today().isoformat()}, as p_trap = 1 - N_return/{N_INSERT}\n")
        f.write("# fifth argument of surf_collide cll/md/trap, padded as the alpha table\n")
        f.write("# format: n_eps n_theta / eps_grid[eV] / theta_grid[deg] / "
                "ptrap[ie][ith]\n")
        f.write(f"{n_e} {n_t}\n")
        f.write(" ".join(f"{x:.4f}" for x in eps_grid) + "\n")
        f.write(" ".join(f"{x:.1f}" for x in th_grid) + "\n")
        f.write("# ptrap\n")
        for r in pt:
            f.write(" ".join(f"{v:.6f}" for v in r) + "\n")


def verify_tokens(path: str, n_e: int, n_t: int) -> None:
    """read the table the way cll/md/trap does: 2 + n_e + n_t + n_e*n_t tokens in [0,1]"""
    tokens = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tokens.extend(s.split())
    expect = 2 + n_e + n_t + n_e * n_t
    assert len(tokens) == expect, f"{path}: {len(tokens)} tokens, expected {expect}"
    assert int(tokens[0]) == n_e and int(tokens[1]) == n_t
    vals = [float(x) for x in tokens[2 + n_e + n_t:]]
    assert all(0.0 <= v <= 1.0 for v in vals), f"{path}: p_trap outside [0,1]"


def main() -> None:
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--csv", default=os.path.join(here, "calibration_curve_ArPt.csv"))
    ap.add_argument("--outdir", default=here)
    ap.add_argument("--install", action="store_true",
                    help="copy the tables into SPARTA_input/tables/")
    a = ap.parse_args()

    data = read_calibration(a.csv)
    pt5 = build_grid(data)
    print("p_trap on the measured grid (eps by row, theta by column):")
    print("        " + " ".join(f"{t:6.0f}deg" for t in TH_MEAS))
    for e, r in zip(EPS_MEAS, pt5):
        print(f"{e:5.1f}eV " + " ".join(f"{v:8.4f}" for v in r))

    eps_grid = [EPS_PAD_LO] + EPS_MEAS
    th_grid = [0.0] + TH_MEAS + [90.0]
    made = {}
    for grazing in ("extrap", "hold"):
        pt = pad_table(pt5, grazing)
        path = os.path.join(a.outdir, NAMES[grazing])
        write_table(path, eps_grid, th_grid, pt, grazing, a.csv)
        verify_tokens(path, len(eps_grid), len(th_grid))
        made[grazing] = path
        print(f"[OK] {path}  ({len(eps_grid)}x{len(th_grid)}, grazing={grazing}, "
              f"tokens checked)")

    if a.install:
        os.makedirs(TABLES_DIR, exist_ok=True)
        for grazing, src_path in made.items():
            dst = os.path.join(TABLES_DIR, NAMES[grazing])
            shutil.copy2(src_path, dst)
            print(f"[INSTALL] {dst}")
    else:
        print("[NOTE] not installed, use --install to write into SPARTA_input/tables/")


if __name__ == "__main__":
    main()
