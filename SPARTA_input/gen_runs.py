#!/usr/bin/env python3
"""Write one run directory per case: 2 geometries x 3 Knudsen numbers x 7 wall models.

Kn = lambda_VHS(T_inf)/D with D = 0.1 m, the cylinder diameter or the plate chord, and
lambda_VHS = [sqrt(2) pi d_ref^2 n (T_ref/T_inf)^(omega-1/2)]^-1. The argon VHS
parameters below are those of ar.vss.

With --seeds N the same cases are repeated under runs_ver/ with different random seeds,
which is what tools/analyze_ci.py and tools/analyze_dvsm.py read for their confidence
intervals.
"""
import argparse
import math
import os
import stat

HERE = os.path.dirname(os.path.abspath(__file__))

# Kn -> free-stream number density
D = 0.1
T_INF, T_REF, DREF, OMEGA = 200.0, 273.15, 4.11e-10, 0.81
SIGMA_FACTOR = math.sqrt(2.0) * math.pi * DREF**2 * (T_REF / T_INF) ** (OMEGA - 0.5)

def nrho_for_kn(kn):
    lam = kn * D
    return 1.0 / (SIGMA_FACTOR * lam)

KN_CASES = [0.1, 1.0, 10.0]
GEOMS = {"cyl": "in.cyl", "plate": "in.plate"}
# diffuse is the fully accommodating baseline, cll_lo and cll_hi bracket the constant
# CLL coefficients. The two cll/md cases share one table and differ only at grazing
# angles: one extrapolates to 90 deg and clips to [0,1], the other holds the 75 deg
# value, so the pair brackets that extrapolation. cllmdtrap adds the trapping branch
# and dvsm replays the MD events directly.
BCS = {
    "diffuse":      "-var WALLMODEL diffuse -var SIGMA 1.0",
    "cll_lo":       "-var WALLMODEL cll -var ALPHAT 0.20 -var ALPHAN 0.05",
    "cll_hi":       "-var WALLMODEL cll -var ALPHAT 0.87 -var ALPHAN 0.70",
    "cllmd_extrap": "-var WALLMODEL cllmd -var ALPHATABLE ../../tables/alpha_arpt_table.txt",
    "cllmd_hold":   "-var WALLMODEL cllmd -var ALPHATABLE ../../tables/alpha_arpt_table_hold.txt",
    "cllmdtrap":    ("-var WALLMODEL cllmdtrap -var ALPHATABLE ../../tables/alpha_arpt_table.txt"
                     " -var PTRAPTABLE ../../tables/ptrap_arpt_table.txt"),
    "dvsm":         "-var WALLMODEL dvsm -var ALPHATABLE ../../events/dvsm_events.txt",
}
NP_MPI = 4
SPA = "${SPARTA_BIN:-spa_mpi}"
BASE_SEED = 12345


def write_case(rdir, deck, nrho, bcargs, seed=None):
    os.makedirs(rdir, exist_ok=True)
    seed_arg = f"-var SEED {seed} " if seed is not None else ""
    cmd = (f"#!/bin/sh\ncd \"$(dirname \"$0\")\"\n"
           f"mpirun -np {NP_MPI} {SPA} -in ../../{deck} "
           f"-var NRHO {nrho:.6g} {bcargs} {seed_arg}"
           f"> stdout.log 2>&1\n")
    sh = os.path.join(rdir, "cmd.sh")
    with open(sh, "w") as f:
        f.write(cmd)
    os.chmod(sh, os.stat(sh).st_mode | stat.S_IEXEC)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=1,
                    help="replicates per case; >1 also fills runs_ver/")
    a = ap.parse_args()

    manifest = ["case\tgeom\tkn\tnrho_m3\tbc\tdir"]
    for geom, deck in GEOMS.items():
        for kn in KN_CASES:
            nrho = nrho_for_kn(kn)
            for bc, bcargs in BCS.items():
                name = f"{geom}_kn{kn:g}_{bc}"
                write_case(os.path.join(HERE, "runs", name), deck, nrho, bcargs)
                manifest.append(f"{name}\t{geom}\t{kn:g}\t{nrho:.6g}\t{bc}\truns/{name}")
                for i in range(2, a.seeds + 1):
                    rep = f"rep_{name}_s{i}"
                    write_case(os.path.join(HERE, "runs_ver", rep), deck, nrho, bcargs,
                               seed=BASE_SEED + i)
                    manifest.append(f"{rep}\t{geom}\t{kn:g}\t{nrho:.6g}\t{bc}\truns_ver/{rep}")

    with open(os.path.join(HERE, "runs", "manifest.tsv"), "w") as f:
        f.write("\n".join(manifest) + "\n")

    print(f"generated {len(manifest) - 1} runs")
    for kn in KN_CASES:
        print(f"  Kn={kn:g}  ->  n_inf = {nrho_for_kn(kn):.4g} m^-3")


if __name__ == "__main__":
    main()
