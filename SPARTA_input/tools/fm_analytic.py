#!/usr/bin/env python3
"""Free-molecular limit, used as an analytic check on the SPARTA runs.

Panel formulae from Schaaf & Chambre (1961); the cylinder coefficients follow
from integrating them around the circumference (standard result, as in the
appendix of Bird 1994).

Speed ratio S = U/√(2RT∞), s_n = S sinθ with θ the angle between the free
stream and the surface; σ_n, σ_t are the normal and tangential momentum
accommodation coefficients and T_w the wall temperature.

  C_p = (1/S²)·[ (2-σ_n)/√π · s_n · e^(−s_n²)
        + { (2−σ_n)(s_n² + 1/2) + σ_n/2 · √(π T_w/T∞) · s_n }·(1 + erf s_n) ]
  C_τ = σ_t cosθ /(S√π) · [ e^(−s_n²) + √π s_n (1 + erf s_n) ]

Leeward side: s_n → −s_n. σ_n=σ_t=0 gives specular reflection, σ_n=σ_t=1 fully
diffuse reflection.

usage:
  python3 tools/fm_analytic.py      # default case U=7000, T=200, Tw=300
"""
import math

R_AR = 208.13          # J/(kg K)
M_AR = 6.6335e-26      # kg


def _elem(S, theta, sn_sign, sig_n, sig_t, TwTi):
    """(Cp, Ct) for one panel. theta = angle between the free stream and the
    surface (rad); sn_sign = +1 windward, -1 leeward."""
    sn = sn_sign * S * math.sin(theta)
    e = math.exp(-sn * sn)
    F = 1.0 + math.erf(sn)
    cp = (1.0 / S**2) * ((2.0 - sig_n) / math.sqrt(math.pi) * sn * e
                         + ((2.0 - sig_n) * (sn * sn + 0.5)
                            + 0.5 * sig_n * math.sqrt(math.pi * TwTi) * sn) * F)
    ct = sig_t * math.cos(theta) / (S * math.sqrt(math.pi)) * (e + math.sqrt(math.pi) * sn * F)
    return cp, ct


def plate_coeffs(S, aoa_deg, sig_n=1.0, sig_t=1.0, TwTi=1.5):
    """Two-sided flat plate: (C_L, C_D, C_m_LE), referenced to 0.5ρU²c and
    0.5ρU²c² for the moment. Uniform loading puts the centre of pressure at
    mid-chord, so C_m_LE = −0.5·C_N with C_N the body-axis normal force."""
    a = math.radians(aoa_deg)
    cpw, ctw = _elem(S, a, +1.0, sig_n, sig_t, TwTi)   # windward
    cpl, ctl = _elem(S, a, -1.0, sig_n, sig_t, TwTi)   # leeward
    cn = cpw - cpl            # body-axis normal, positive towards the lee side
    ca = ctw + ctl            # body-axis axial, along the chord
    cl = cn * math.cos(a) - ca * math.sin(a)
    cd = cn * math.sin(a) + ca * math.cos(a)
    # same convention as analyze.py: z-moment about the LE, anticlockwise positive
    # centre of pressure at mid-chord, r=(c/2)(cosa,-sina), F=(Fd,Fl) → Mz = x·Fy − y·Fx
    xm, ym = 0.5 * math.cos(a), -0.5 * math.sin(a)
    cm_le = xm * cl - ym * cd
    return cl, cd, cm_le


def cylinder_cd(S, sig_n=1.0, sig_t=1.0, TwTi=1.5, N=2000):
    """Cylinder drag coefficient, referenced to 0.5ρU²D. The panel formulae are
    integrated over the windward half and doubled by symmetry; the leeward half
    enters through the s_n<0 branch."""
    D_acc = 0.0
    for k in range(N):
        phi = math.pi * (k + 0.5) / N          # angle between the panel normal and −x, 0..π
        theta = 0.5 * math.pi - phi            # angle between the free stream and the surface
        sign = 1.0 if theta >= 0 else -1.0
        cp, ct = _elem(S, abs(theta), sign, sig_n, sig_t, TwTi)
        # upper and lower halves combined
        D_acc += (cp * math.cos(phi) + ct * math.sin(phi)) * (math.pi / N)
    return D_acc  # already normalised by the diameter


if __name__ == "__main__":
    U, T, Tw = 7000.0, 200.0, 300.0
    S = U / math.sqrt(2.0 * R_AR * T)
    print(f"S = {S:.2f}  (U={U}, T_inf={T}, Ar)")
    for tag, sn, st in [("diffuse σ=1", 1.0, 1.0), ("specular σ=0", 0.0, 0.0),
                        ("cll_lo proxy (0.05,0.20)", 0.05, 0.20),
                        ("cll_hi proxy (0.70,0.87)", 0.70, 0.87)]:
        cl, cd, cm = plate_coeffs(S, 20.0, sn, st, Tw / T)
        cdc = cylinder_cd(S, sn, st, Tw / T)
        print(f"{tag:26s}  plate: CL={cl:+.4f} CD={cd:.4f} Cm_LE={cm:+.4f} | cyl: CD={cdc:.4f}")
    print("\nNote: σ_n/σ_t are momentum accommodation coefficients, not the CLL α_n")
    print("(normal energy accommodation), so the cll_lo/hi rows are indicative only.")
    print("For a strict check compare the diffuse and specular limits with Kn=10 (< 5%).")
