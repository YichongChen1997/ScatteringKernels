#!/usr/bin/env python3
"""Replicate statistics for the verification runs.

- CI: four independent replicates per case, base (seed 12345) plus
  rep_*_s{2,3,4}, with t(3, 0.975) = 3.182. Deviations against diffuse are
  paired by seed before the CI is taken.
- convergence: conv_{geom}_kn1_cllmd_{dt2,npc2} against the matching base run.
- degeneracy: degen_{geom}_kn1_cllmdconst against {geom}_kn1_cll_hi. A constant
  table in cll/md should reproduce the native cll statistics.

Writes runs_ver/ci_summary.tsv and a report on stdout.

usage:
  python3 tools/analyze_ci.py
"""
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import read_dump, M_AR, D_REF, U_INF, LE

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T975_3 = 3.182


def metrics(rdir, nrho):
    rows = read_dump(os.path.join(rdir, "surf_force.dump"))
    if not rows:
        return None
    qref = 0.5 * nrho * M_AR * U_INF**2
    Fx = sum(r[5] for r in rows)
    Fy = sum(r[6] for r in rows)
    M = sum((0.5 * (r[1] + r[3])) * r[6] - (0.5 * (r[2] + r[4])) * r[5] for r in rows)
    stag = min(rows, key=lambda r: 0.5 * (r[1] + r[3]))
    return dict(CD=Fx / (qref * D_REF), CL=Fy / (qref * D_REF),
                Cm=M / (qref * D_REF**2), q=stag[7])


def nrho_of(case):
    kn = float(case.split("_")[1][2:])
    import math as m
    lam = kn * 0.1
    sf = m.sqrt(2) * m.pi * (4.11e-10)**2 * (273.15 / 200.0)**(0.81 - 0.5)
    return 1.0 / (sf * lam)


def stats(vals):
    n = len(vals)
    mu = sum(vals) / n
    if n < 2:
        return mu, float("nan")
    sd = (sum((v - mu)**2 for v in vals) / (n - 1)) ** 0.5
    return mu, T975_3 * sd / math.sqrt(n)


def replicates(case):
    dirs = [f"{HERE}/runs/{case}"] + [f"{HERE}/runs_ver/rep_{case}_s{i}" for i in (2, 3, 4)]
    nrho = nrho_of(case)
    out = []
    for d in dirs:
        m_ = metrics(d, nrho) if os.path.exists(os.path.join(d, "surf_force.dump")) else None
        out.append(m_)
    return out


def main():
    cases = sorted(os.path.basename(d) for d in glob.glob(f"{HERE}/runs/*") if os.path.isdir(d))
    lines = ["case\tCD_mean\tCD_ci95\tCm_mean\tCm_ci95\tq_mean\tq_ci95\tn_rep"]
    R = {}
    for c in cases:
        reps = [r for r in replicates(c) if r]
        if not reps:
            continue
        R[c] = reps
        cd, cdc = stats([r["CD"] for r in reps])
        cm, cmc = stats([r["Cm"] for r in reps])
        q, qc = stats([r["q"] for r in reps])
        lines.append(f"{c}\t{cd:.5g}\t{cdc:.2g}\t{cm:.5g}\t{cmc:.2g}\t{q:.5g}\t{qc:.2g}\t{len(reps)}")
    with open(f"{HERE}/runs_ver/ci_summary.tsv", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote runs_ver/ci_summary.tsv ({len(R)} cases, up to 4 replicates)\n")

    print("moment deviation, plate cllmd_extrap vs diffuse, seed-paired 95% CI")
    for kn in ["kn0.1", "kn1", "kn10"]:
        a = replicates(f"plate_{kn}_cllmd_extrap")
        b = replicates(f"plate_{kn}_diffuse")
        devs = [100 * (x["Cm"] - y["Cm"]) / abs(y["Cm"]) for x, y in zip(a, b) if x and y]
        mu, ci = stats(devs)
        sig = "significant" if abs(mu) > ci else "not significant!"
        print(f"  {kn:6s}: ΔCm = {mu:+7.2f}% ± {ci:.2f}%  (n={len(devs)}) {sig}")

    print("\nconvergence, kn1 cllmd_extrap against base")
    for geom in ["cyl", "plate"]:
        base = metrics(f"{HERE}/runs/{geom}_kn1_cllmd_extrap", nrho_of(f"{geom}_kn1_x"))
        for tag in ["dt2", "npc2"]:
            d = f"{HERE}/runs_ver/conv_{geom}_kn1_cllmd_{tag}"
            m_ = metrics(d, nrho_of(f"{geom}_kn1_x")) if os.path.exists(os.path.join(d, "surf_force.dump")) else None
            if m_ and base:
                print(f"  {geom}/{tag}: ΔCD={100*(m_['CD']-base['CD'])/abs(base['CD']):+.2f}%  "
                      f"Δq={100*(m_['q']-base['q'])/abs(base['q']):+.2f}%  "
                      f"ΔCm={100*(m_['Cm']-base['Cm'])/abs(base['Cm']):+.2f}%")

    print("\ndegeneracy, cll/md with a constant table vs native cll_hi, kn1")
    for geom in ["cyl", "plate"]:
        d = f"{HERE}/runs_ver/degen_{geom}_kn1_cllmdconst"
        m_ = metrics(d, nrho_of(f"{geom}_kn1_x")) if os.path.exists(os.path.join(d, "surf_force.dump")) else None
        reps = replicates(f"{geom}_kn1_cll_hi")
        if m_ and reps:
            cd, cdc = stats([r["CD"] for r in reps if r])
            q, qc = stats([r["q"] for r in reps if r])
            print(f"  {geom}: CD {m_['CD']:.4f} vs {cd:.4f}±{cdc:.4f} "
                  f"({'ok' if abs(m_['CD']-cd) < max(cdc, 0.02*abs(cd)) else 'differs!'});  "
                  f"q {m_['q']:.4g} vs {q:.4g}±{qc:.2g} "
                  f"({'ok' if abs(m_['q']-q) < max(qc, 0.03*abs(q)) else 'differs!'})")


if __name__ == "__main__":
    main()
