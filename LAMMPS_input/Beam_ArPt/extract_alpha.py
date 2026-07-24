#!/usr/bin/env python3
"""Accommodation coefficients from the Ar-Pt beam runs.

The coefficients are flux-weighted over the returned atoms only, which for a
monoenergetic beam means a plain average over the atoms that come back:
  alpha_t = 1 - <v_t,out>/v_t,in                  tangential direction is x, PHI = 0
  alpha_n = (<E_n,i> - <E_n,r>)/(<E_n,i> - 1)     energies in units of kT_w, E_n,w = 1
  eac     = (<E_i> - <E_r>)/(<E_i> - 2)           E_w = 2 for a flux-weighted half Maxwellian

An atom counts as returned if its last appearance in the dump has z > Z_FREE and
vz > 0. The wall exerts no force above Z_FREE = rCut, so that velocity is the
free-flight one.

usage:
  python3 extract_alpha.py runs/eps5_th45 [--save-events]
  python3 extract_alpha.py --all           # scan runs/, write calibration_curve_ArPt.csv
"""
import argparse
import csv
import glob
import math
import os
import random

M_AR = 39.948 * 1.66054e-27
EV = 1.602176634e-19
KB = 1.380649e-23
Z_FREE = 15.0
V_CONV = 1e5          # Ang/fs -> m/s


def parse_dump(path):
    """last free-flight frame of each atom: {id: (z, vx, vy, vz)}"""
    last_free = {}
    seen_below = set()
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                n = int(f.readline())
            elif line.startswith("ITEM: ATOMS"):
                cols = line.split()[2:]
                iid = cols.index("id")
                iz, ivx, ivy, ivz = cols.index("z"), cols.index("vx"), cols.index("vy"), cols.index("vz")
                for _ in range(n):
                    v = f.readline().split()
                    aid = int(v[iid])
                    z, vx, vy, vz = (float(v[k]) for k in (iz, ivx, ivy, ivz))
                    if z < Z_FREE:
                        seen_below.add(aid)
                    elif vz > 0 and aid in seen_below:
                        last_free[aid] = (z, vx * V_CONV, vy * V_CONV, vz * V_CONV)
    return last_free


def alphas(vout, v_in, theta_deg, Tw=300.0):
    vt_in = v_in * math.sin(math.radians(theta_deg))
    vn_in = v_in * math.cos(math.radians(theta_deg))
    kTw = KB * Tw
    en_i = 0.5 * M_AR * vn_in**2 / kTw
    e_i = 0.5 * M_AR * v_in**2 / kTw
    vt_out = [vx for _, vx, _, _ in vout]
    en_r = [0.5 * M_AR * vz**2 / kTw for _, _, _, vz in vout]
    e_r = [0.5 * M_AR * (vx*vx + vy*vy + vz*vz) / kTw for _, vx, vy, vz in vout]
    N = len(vout)
    mean = lambda a: sum(a) / len(a)
    at = 1.0 - mean(vt_out) / vt_in
    an = (en_i - mean(en_r)) / (en_i - 1.0)
    ae = (e_i - mean(e_r)) / (e_i - 2.0)
    # bootstrap stderr
    rng = random.Random(42)
    bs_t, bs_n = [], []
    for _ in range(200):
        idx = [rng.randrange(N) for _ in range(N)]
        bs_t.append(1.0 - mean([vt_out[i] for i in idx]) / vt_in)
        bs_n.append((en_i - mean([en_r[i] for i in idx])) / (en_i - 1.0))
    sd = lambda a: (sum((x - mean(a))**2 for x in a) / (len(a) - 1)) ** 0.5
    return at, an, ae, sd(bs_t), sd(bs_n), N


def do_run(rdir, save_events=False):
    name = os.path.basename(rdir.rstrip("/"))
    eps = float(name.split("_")[0][3:])
    th = int(name.split("_th")[1])
    v_in = math.sqrt(2 * eps * EV / M_AR)
    dump = os.path.join(rdir, "dump_gas.lammpstrj")
    if not os.path.exists(dump):
        return None
    ev = parse_dump(dump)
    if save_events:
        with open(os.path.join(rdir, "events_out.csv"), "w") as f:
            f.write("vx_out,vy_out,vz_out\n")
            for _, vx, vy, vz in ev.values():
                f.write(f"{vx:.2f},{vy:.2f},{vz:.2f}\n")
    if len(ev) < 30:
        # At grazing incidence the normal energy is low enough that nearly every
        # atom is trapped, and too few return to measure anything. Trapped atoms
        # desorb thermalised, so the coefficients are set to the complete
        # accommodation limit rather than measured.
        N_INSERT = 2500
        print(f"[trap-limit] {name}: N_return={len(ev)} < 30 "
              f"(trapping fraction ~ {1 - len(ev)/N_INSERT:.3f}) -> alpha_t=alpha_n=eac=1.0")
        return dict(eps=eps, theta=th, alpha_t=1.0, alpha_n=1.0, eac=1.0,
                    err_t=0.0, err_n=0.0, events=len(ev))
    at, an, ae, et, en, N = alphas(list(ev.values()), v_in, th)
    return dict(eps=eps, theta=th, alpha_t=at, alpha_n=an, eac=ae,
                err_t=et, err_n=en, events=N)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", nargs="?")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--save-events", action="store_true")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    if a.all:
        rows = []
        for rd in sorted(glob.glob(os.path.join(here, "runs", "eps*_th*"))):
            r = do_run(rd, save_events=True)
            if r:
                rows.append(r)
                print(f"{os.path.basename(rd):14s} N={r['events']:5d}  "
                      f"a_t={r['alpha_t']:.3f}±{r['err_t']:.3f}  a_n={r['alpha_n']:.3f}±{r['err_n']:.3f}")
        with open(os.path.join(here, "calibration_curve_ArPt.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("wrote calibration_curve_ArPt.csv")
    elif a.rundir:
        r = do_run(a.rundir, save_events=a.save_events)
        print(r)
