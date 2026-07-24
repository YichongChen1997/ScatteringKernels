#!/usr/bin/env python3
"""Separate the energy-resolved coefficients from the shape of the kernel by
comparing the separable CLL fit against direct resampling of the MD events.

Three boundary conditions, seed-paired over 4 seeds (12345 + s2/s3/s4):
  diffuse       fully accommodated baseline
  cllmd_extrap  energy-resolved coefficients with the CLL shape (separable)
  dvsm          resampled MD events (non-separable, kernel-free reference)

Decomposition of each observable O, per geometry and Kn:
  Δ_full  = O_dvsm  − O_diffuse    full effect of the MD information
  Δ_coeff = O_cllmd − O_diffuse    part captured by the separable layer
  Δ_shape = O_dvsm  − O_cllmd      shape residual
  X = Δ_coeff/Δ_full (captured fraction),  Y = Δ_shape/O_cllmd
  Δ_shape_pts = (O_dvsm − O_cllmd)/|O_diffuse|, on the common base, so that
                Δ_full = Δ_coeff + Δ_shape_pts
Everything is paired by seed first, then given as a 95% CI over the 4
replicates (t(3)=3.182).

Observables: cylinder CD and q_stag, plate CD and Cm.
Writes runs_ver/dvsm_decomposition.tsv and a report on stdout.

usage:
  python3 tools/analyze_dvsm.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_ci import metrics, nrho_of, stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

OBS = {"cyl": ["CD", "q"], "plate": ["CD", "Cm"]}
KNS = ["0.1", "1", "10"]


def reps(case: str) -> list:
    dirs = [f"{HERE}/runs/{case}"] + \
           [f"{HERE}/runs_ver/rep_{case}_s{i}" for i in (2, 3, 4)]
    nrho = nrho_of(case)
    out = []
    for d in dirs:
        has = os.path.exists(os.path.join(d, "surf_force.dump"))
        out.append(metrics(d, nrho) if has else None)
    return out


def main() -> None:
    lines = ["geom\tkn\tobs\tdiffuse\tcllmd\tdvsm"
             "\tDfull_pct\tDfull_ci\tDcoeff_pct\tDcoeff_ci"
             "\tDshape_pct\tDshape_ci\tDshapepts_pct\tDshapepts_ci"
             "\tcapture_X_pct"]
    for geom, obs_list in OBS.items():
        for kn in KNS:
            r_dif = reps(f"{geom}_kn{kn}_diffuse")
            r_cll = reps(f"{geom}_kn{kn}_cllmd_extrap")
            r_dvs = reps(f"{geom}_kn{kn}_dvsm")
            for obs in obs_list:
                trip = [(a[obs], b[obs], c[obs])
                        for a, b, c in zip(r_dif, r_cll, r_dvs)
                        if a and b and c]
                if len(trip) < 2:
                    print(f"[skip] {geom} kn{kn} {obs}: "
                          f"{len(trip)} paired replicates")
                    continue
                dfull = [100 * (v - d) / abs(d) for d, _, v in trip]
                dcoef = [100 * (c - d) / abs(d) for d, c, _ in trip]
                dshap = [100 * (v - c) / abs(c) for _, c, v in trip]
                dspts = [100 * (v - c) / abs(d) for d, c, v in trip]
                mu_f, ci_f = stats(dfull)
                mu_c, ci_c = stats(dcoef)
                mu_s, ci_s = stats(dshap)
                mu_p, ci_p = stats(dspts)
                x_pct = 100 * mu_c / mu_f if abs(mu_f) > 1e-12 else float("nan")
                d0 = sum(t[0] for t in trip) / len(trip)
                c0 = sum(t[1] for t in trip) / len(trip)
                v0 = sum(t[2] for t in trip) / len(trip)
                lines.append(
                    f"{geom}\t{kn}\t{obs}\t{d0:.5g}\t{c0:.5g}\t{v0:.5g}"
                    f"\t{mu_f:+.2f}\t{ci_f:.2f}\t{mu_c:+.2f}\t{ci_c:.2f}"
                    f"\t{mu_s:+.2f}\t{ci_s:.2f}\t{mu_p:+.2f}\t{ci_p:.2f}"
                    f"\t{x_pct:.1f}")
                print(f"{geom:5s} kn{kn:4s} {obs:2s}: "
                      f"Δfull={mu_f:+7.2f}%±{ci_f:.2f}  "
                      f"Δcoeff={mu_c:+7.2f}%±{ci_c:.2f}  "
                      f"Δshape_pts={mu_p:+6.2f}±{ci_p:.2f}  "
                      f"(vs cllmd {mu_s:+6.2f}%±{ci_s:.2f})  "
                      f"capture X={x_pct:5.1f}%")
    out = f"{HERE}/runs_ver/dvsm_decomposition.tsv"
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
