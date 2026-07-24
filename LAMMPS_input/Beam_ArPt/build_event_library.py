#!/usr/bin/env python3
"""Pack the outgoing velocities of the 25 beam runs into the dvsm/md event library.

Events come from runs/eps*_th*/events_out.csv, written by extract_alpha.py
--save-events, in the beam frame with +z along the outward normal and +x along
the incident tangential direction (PHI = 0), in m/s. The trapping probability is
p_trap = 1 - N_return/N_INSERT with N_INSERT = 2500 atoms per condition.

Output is a token stream with # comment lines ignored, as surf_collide dvsm/md
reads it:
  n_eps n_theta
  eps_grid (eV)
  theta_grid (deg)
  per bin, eps outer and theta inner:
    eps theta N_return p_trap
    N_return lines of vx vy vz

usage: python3 build_event_library.py [--install]
  --install writes into SPARTA_input/events/dvsm_events.txt
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EPS_GRID = [0.4, 1.0, 2.0, 5.0, 10.0]
THETA_GRID = [15, 30, 45, 60, 75]
N_INSERT = 2500
INSTALL_PATH = os.path.join(
    HERE, "..", "..", "SPARTA_input", "events", "dvsm_events.txt")


def dirname_for(eps: float, th: int) -> str:
    eps_s = str(int(eps)) if eps == int(eps) else str(eps)
    return os.path.join(HERE, "runs", f"eps{eps_s}_th{th}")


def read_events(rdir: str) -> list[tuple[float, float, float]]:
    path = os.path.join(rdir, "events_out.csv")
    if not os.path.exists(path):
        sys.exit(f"FATAL: missing {path}")
    rows = []
    with open(path) as f:
        header = f.readline()
        if not header.startswith("vx_out"):
            sys.exit(f"FATAL: unexpected header in {path}: {header!r}")
        for line in f:
            line = line.strip()
            if not line:
                continue
            vx, vy, vz = (float(x) for x in line.split(","))
            if vz <= 0:
                sys.exit(f"FATAL: non-outgoing event vz={vz} in {path}")
            rows.append((vx, vy, vz))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--install", action="store_true")
    args = ap.parse_args()

    out = os.path.abspath(INSTALL_PATH) if args.install \
        else os.path.join(HERE, "dvsm_events.txt")

    lines = [
        "# DVSM event library: Ar-Pt beam MD outgoing velocities (m/s)",
        "# frame: +z outward normal, +x incident tangential (PHI=0)",
        "# per bin: eps(eV) theta(deg) N_return p_trap, then N_return x (vx vy vz)",
        f"{len(EPS_GRID)} {len(THETA_GRID)}",
        " ".join(str(e) for e in EPS_GRID),
        " ".join(str(t) for t in THETA_GRID),
    ]
    total = 0
    for eps in EPS_GRID:
        for th in THETA_GRID:
            ev = read_events(dirname_for(eps, th))
            n = len(ev)
            if n > N_INSERT:
                sys.exit(f"FATAL: eps{eps}_th{th} has {n} > N_INSERT events")
            p_trap = 1.0 - n / N_INSERT
            lines.append(f"{eps} {th} {n} {p_trap:.6f}")
            lines.extend(f"{vx:.2f} {vy:.2f} {vz:.2f}" for vx, vy, vz in ev)
            total += n
            print(f"eps{eps}_th{th}: N_return={n:4d} p_trap={p_trap:.4f}")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}: {len(EPS_GRID)}x{len(THETA_GRID)} bins, "
          f"{total} events total")


if __name__ == "__main__":
    main()
