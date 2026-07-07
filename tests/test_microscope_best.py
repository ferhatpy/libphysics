#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_microscope.py
==================
Compound Microscope System via ABCD Ray-Transfer Matrix Method.

Based on:
    Kloos, Gerhard. "Matrix Methods for Optical Layout." SPIE, 2007.
    Chapter 1: An Introduction to Tools and Concepts.

Physics summary
---------------
A compound microscope is modelled as a *lens doublet* (Sec. 1.7):
  - Objective lens  : focal length f1 > 0
  - Eyepiece  lens  : focal length f2 > 0
  - Optical interval: E > 0   with d = f1 + E + f2  (Eq. 1.57)

  System matrix of the doublet (Eq. 1.55):
      S = L2 · T(d) · L1

  Effective focal length (Eq. 1.58):
      f_eff = -f1 * f2 / E          (divergent system, f < 0)

  Full imaging chain (Eq. 1.47/1.49):
      S_total = T(b) · L2 · T(d) · L1 · T(g)
      Imaging condition: s12 = 0

  Total magnification (classical formula):
      M_total = -(E · d_n) / (f1 · f2)
      d_n ≈ 250 mm  (near-point / standard viewing distance)

References
----------
    Kloos2007  Eq. 1.7   Translation matrix T
    Kloos2007  Eq. 1.17  Refraction matrix R
    Kloos2007  Eq. 1.34  Thin lens matrix L
    Kloos2007  Eq. 1.47  S = B · A · G  (image · lens · object)
    Kloos2007  Eq. 1.50  Imaging condition: 1/g + 1/b = 1/f
    Kloos2007  Eq. 1.55  Lens-doublet system matrix
    Kloos2007  Eq. 1.57  d = f1 + E + f2
    Kloos2007  Eq. 1.58  f = -f1*f2/E  (compound microscope / divergent doublet)
"""

import sys
import os

lstPaths = ["../src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)

from libsympy import *
from sympy.abc import*
from optics import *

print(sys.version)
print(sys.path)

# =============================================================================
# Settings
# =============================================================================

class sets:
    """
    Execution settings.

    flow keys
    ---------
    100  get_formulary           – print the optics formulary
    200  doublet_system_matrix   – symbolic S for compound microscope doublet
    300  imaging_condition       – derive s12=0, solve for image distance
    400  magnification           – total lateral magnification
    500  numerical_example       – plug in real microscope numbers
    600  ray_trace               – trace two characteristic rays
    """
    global dictflow, test_all

    def __init__(self):
        pass

    input_dir  = "input/optics"
    output_dir = "output/optics"

    test_all = {0: False, 1: True}[0]

    dictflow = dict(
        ch1={
            100: "get_formulary",
            200: "doublet_system_matrix",
            300: "imaging_condition",
            400: "magnification",
            500: "numerical_example",
            600: "ray_trace",
        }
    )

    # Choose which tests to run:
    flow = [dictflow["ch1"][i] for i in [200, 300, 400, 500, 600]]
    if test_all:
        flow = flatten([list(dictflow[i].values()) for i in dictflow.keys()])

print("Running: {0}".format(sets.flow))

# =============================================================================
# 100 – get_formulary
# =============================================================================

if "get_formulary" in sets.flow:
    print("\n" + "="*60)
    print("ABCD Formulary")
    print("="*60)
    oopti.__init__()
    oopti.get_formulary()

# =============================================================================
# 200 – Doublet system matrix  (Kloos2007 Eq. 1.55)
# =============================================================================

if "doublet_system_matrix" in sets.flow:
    print("\n" + "="*60)
    print("200 · Compound Microscope – Doublet System Matrix")
    print("="*60)

    oopti.__init__()
    oopti.verbose = False

    # --- symbols -------------------------------------------------------
    f1, f2 = symbols('f_1 f_2', positive=True)
    E      = symbols('E',       positive=True)   # optical interval (tube length)
    d      = symbols('d',       positive=True)   # physical lens separation

    abcd = oopti.ABCD

    # Thin-lens matrix  L(f) = [[1, 0], [-1/f, 1]]  (Kloos2007 Eq. 1.34)
    L1 = abcd.thin_lens.abcd(f1).rhs.doit()   # objective
    L2 = abcd.thin_lens.abcd(f2).rhs.doit()   # eyepiece

    # Translation matrix  T(d) = [[1, d], [0, 1]]  (Kloos2007 Eq. 1.7)
    Td = abcd.T(d).rhs.doit()

    # Doublet system matrix:  S = L2 · T(d) · L1  (Kloos2007 Eq. 1.52 / 1.55)
    S_doublet = MatMul(L2, Td, L1).doit()

    S_sym = symbols('S_M')
    print("\nS = L2 · T(d) · L1  (Kloos2007 Eq. 1.55):")
    display(Eq(S_sym, UnevaluatedExpr(S_doublet)))

    s11 = S_doublet[0, 0]
    s12 = S_doublet[0, 1]
    s21 = S_doublet[1, 0]
    s22 = S_doublet[1, 1]

    print("\nMatrix entries (simplified):")
    print("  s11 =", simplify(s11))
    print("  s12 =", simplify(s12))
    print("  s21 =", simplify(s21))
    print("  s22 =", simplify(s22))

    # Verify det S = 1  (both media in air, n1=n2=1)
    det_S = simplify(s11*s22 - s12*s21)
    print("\n  det(S) =", det_S, "  (should be 1)")
    assert det_S == 1, "Determinant check failed!"

    # --- Apply d = f1 + E + f2  (Kloos2007 Eq. 1.57) -------------------
    d_expr = f1 + E + f2
    S_micro = S_doublet.subs(d, d_expr)
    print("\nAfter d = f1 + E + f2  (Kloos2007 Eq. 1.57):")
    display(Eq(S_sym, UnevaluatedExpr(simplify(S_micro))))

    # Effective focal length from s21  (Kloos2007 Eq. 1.56 / 1.58)
    s21_micro = simplify(S_micro[1, 0])
    f_eff_sym = simplify(-1 / s21_micro)
    f_eff_lhs = symbols('f_eff')
    print("\n  f_eff = -f1·f2/E  (Kloos2007 Eq. 1.58):")
    display(Eq(f_eff_lhs, f_eff_sym))

    assert simplify(f_eff_sym - (-f1*f2/E)) == 0, "Focal length formula mismatch!"
    print("  ✓ Confirmed: f_eff = -f1·f2/E  (divergent system)")

# =============================================================================
# 300 – Imaging condition  s12 = 0  (Kloos2007 Eq. 1.47–1.50)
# =============================================================================

if "imaging_condition" in sets.flow:
    print("\n" + "="*60)
    print("300 · Imaging Condition  s12 = 0")
    print("="*60)
    """
    Full imaging chain:
        S_total = T(b) · L2 · T(d) · L1 · T(g)
    
    g = object distance from objective
    b = image  distance to   image plane (behind eyepiece)
    Imaging condition:  S_total[0,1] = 0  (Kloos2007 Eq. 1.48)
    """
    oopti.__init__()
    oopti.verbose = False

    f1, f2  = symbols('f_1 f_2', positive=True)
    E       = symbols('E',       positive=True)
    d, g, b = symbols('d g b',   positive=True)

    abcd = oopti.ABCD

    L1 = abcd.thin_lens.abcd(f1).rhs.doit()
    L2 = abcd.thin_lens.abcd(f2).rhs.doit()
    Td = abcd.T(d).rhs.doit()
    Tg = abcd.T(g).rhs.doit()
    Tb = abcd.T(b).rhs.doit()

    # Full system matrix  S = T(b) · L2 · T(d) · L1 · T(g)
    S_total   = MatMul(Tb, L2, Td, L1, Tg).doit()
    s12_total = simplify(S_total[0, 1])

    s12_lhs = symbols('s_{12}')
    print("\nImaging condition  s12 = 0:")
    display(Eq(s12_lhs, s12_total))

    # With d = f1 + E + f2
    d_expr   = f1 + E + f2
    s12_sub  = simplify(s12_total.subs(d, d_expr))
    print("\nAfter substituting d = f1 + E + f2:")
    display(Eq(s12_lhs, s12_sub))

    # Solve for b given g
    b_sol = solve(s12_sub, b)
    print("\nImage distance b(g, f1, f2, E)  from s12=0:")
    for sol in b_sol:
        display(Eq(b, simplify(sol)))

# =============================================================================
# 400 – Total lateral magnification
# =============================================================================

if "magnification" in sets.flow:
    print("\n" + "="*60)
    print("400 · Total Lateral Magnification")
    print("="*60)
    """
    For the imaging system S_total = T(b)·L2·T(d)·L1·T(g) with s12=0:
        lateral magnification M = s11_total

    Classical microscope formula:
        M_obj   = -v_obj / u_obj  ≈ -E / f1
        M_eye   = d_n / f2               (angular magnification of eyepiece)
        M_total = -(E · d_n) / (f1 · f2)
    """
    oopti.__init__()
    oopti.verbose = False

    f1, f2  = symbols('f_1 f_2', positive=True)
    E       = symbols('E',       positive=True)
    d, g, b = symbols('d g b',   positive=True)
    dn      = symbols('d_n',     positive=True)   # near-point ≈ 250 mm

    abcd = oopti.ABCD
    L1 = abcd.thin_lens.abcd(f1).rhs.doit()
    L2 = abcd.thin_lens.abcd(f2).rhs.doit()
    Td = abcd.T(d).rhs.doit()
    Tg = abcd.T(g).rhs.doit()
    Tb = abcd.T(b).rhs.doit()

    S_total   = MatMul(Tb, L2, Td, L1, Tg).doit()
    s11_total = simplify(S_total[0, 0])

    M_lhs = symbols('M')
    print("\nLateral magnification M = s11  (on the imaging condition):")
    display(Eq(M_lhs, s11_total))

    # Substitute d = f1+E+f2, g = f1 (object at front focal of objective)
    d_expr    = f1 + E + f2
    s11_sub   = simplify(s11_total.subs(d, d_expr).subs(g, f1))
    print("\nM with g=f1 and d=f1+E+f2 (paraxial approximation):")
    display(Eq(M_lhs, s11_sub))

    # Classical formula
    M_classical_lhs = symbols('M_total')
    M_classical_rhs = -(E * dn) / (f1 * f2)
    print("\nClassical formula  M_total = -(E·d_n)/(f1·f2):")
    display(Eq(M_classical_lhs, M_classical_rhs))
    print("  where d_n ≈ 250 mm is the standard near-point viewing distance.")

# =============================================================================
# 500 – Numerical example with real microscope parameters
# =============================================================================

if "numerical_example" in sets.flow:
    print("\n" + "="*60)
    print("500 · Numerical Example – Real Microscope Parameters")
    print("="*60)
    """
    Typical 40× objective / 10× eyepiece compound microscope (DIN standard):
        f1  =   4 mm   objective focal length
        f2  =  25 mm   eyepiece  focal length
        E   = 160 mm   optical interval (DIN tube length)
        d_n = 250 mm   near-point distance

    Expected total magnification ≈ 400×
    """
    import numpy as np

    # Parameters [mm]
    nf1 = 4.0
    nf2 = 25.0
    nE  = 160.0
    ndn = 250.0
    nd  = nf1 + nE + nf2

    print(f"\n  Microscope parameters:")
    print(f"    f1  = {nf1} mm   (objective focal length)")
    print(f"    f2  = {nf2} mm   (eyepiece  focal length)")
    print(f"    E   = {nE}  mm   (optical interval / tube length)")
    print(f"    d   = f1+E+f2 = {nd} mm  (physical lens separation)")
    print(f"    d_n = {ndn} mm   (near-point distance)")

    # Effective focal length
    nf_eff = -nf1 * nf2 / nE
    print(f"\n  Effective focal length:")
    print(f"    f_eff = -f1·f2/E = {nf_eff:.3f} mm  (negative → divergent ✓)")

    # Thin-lens matrices (numpy)
    def T_mat(t):
        return np.array([[1.0, t], [0.0, 1.0]])

    def L_mat(f):
        return np.array([[1.0, 0.0], [-1.0/f, 1.0]])

    # Doublet matrix
    S_num = L_mat(nf2) @ T_mat(nd) @ L_mat(nf1)
    print(f"\n  Doublet system matrix S = L2·T(d)·L1:")
    print(f"    [{S_num[0,0]:+.6f}  {S_num[0,1]:+.6f}]")
    print(f"    [{S_num[1,0]:+.6f}  {S_num[1,1]:+.6f}]")
    print(f"  det(S) = {np.linalg.det(S_num):.8f}  (should be 1.0)")

    # Object / image distances
    # Object just beyond f1 so objective forms real intermediate image at v = f1+E
    ng  = 1.0 / (1.0/nf1 - 1.0/(nf1 + nE))   # thin-lens: 1/g + 1/v = 1/f1, v=f1+E
    nv  = nf1 + nE                             # intermediate image at v from objective
    # Eyepiece used as magnifier: object at its front focal plane → virtual image at ∞
    # or image at near-point: 1/b + 1/f2 = 1/0 ... use standard angular magnification

    M_obj   = -nv / ng
    M_eye   = ndn / nf2
    M_total = M_obj * M_eye
    M_class = -(nE * ndn) / (nf1 * nf2)

    print(f"\n  Object distance from objective:  g = {ng:.4f} mm")
    print(f"  Intermediate image distance:     v = {nv:.4f} mm  (from objective)")

    print(f"\n  Magnifications:")
    print(f"    M_obj   = -v/g         = {M_obj:.2f}×")
    print(f"    M_eye   = d_n/f2       = {M_eye:.1f}×")
    print(f"    M_total = M_obj × M_eye = {M_total:.1f}×")
    print(f"    M_class = -(E·d_n)/(f1·f2) = {M_class:.1f}×  (classical formula)")

    # Principal planes of the doublet (Kloos2007 Eq. 1.38 / 1.39)
    s11_n = S_num[0, 0]
    s21_n = S_num[1, 0]
    s22_n = S_num[1, 1]
    h1_n  = (1.0 - s11_n) / s21_n
    h2_n  = (1.0 - s22_n) / s21_n
    print(f"\n  Principal plane positions (from 2nd ref. plane of L1):")
    print(f"    h1 = {h1_n:.4f} mm   (first  principal plane)")
    print(f"    h2 = {h2_n:.4f} mm   (second principal plane)")

    # Rayleigh resolution limit
    nlambda = 550e-6   # 550 nm green light [mm]
    nNA     = 0.65     # typical 40× objective numerical aperture
    delta_r = 0.61 * nlambda / nNA
    print(f"\n  Resolution limit (Rayleigh criterion):")
    print(f"    λ  = {nlambda*1e6:.0f} nm")
    print(f"    NA = {nNA}")
    print(f"    δr = 0.61·λ/NA = {delta_r*1e3:.4f} µm  = {delta_r*1e6:.1f} nm")

# =============================================================================
# 600 – Paraxial ray tracing through the microscope
# =============================================================================

if "ray_trace" in sets.flow:
    print("\n" + "="*60)
    print("600 · Paraxial Ray Tracing (Central Theorem, Kloos2007 Sec. 1.9)")
    print("="*60)
    """
    Trace two characteristic rays (Kloos2007 Sec. 1.9):
      Ray a (marginal):  enters on-axis at small angle β
      Ray b (chief):     enters at maximum height, parallel to axis

    Optical path:  Object → T(g) → L1(f1) → T(d) → L2(f2) → T(b) → Image
    """
    import numpy as np
    import matplotlib
    # matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # Parameters [mm]
    nf1 = 4.0
    nf2 = 25.0
    nE  = 160.0
    nd  = nf1 + nE + nf2
    ndn = 250.0

    # Object distance (thin-lens formula for objective)
    ng = 1.0 / (1.0/nf1 - 1.0/(nf1 + nE))
    nv = nf1 + nE       # intermediate image from objective

    # Image distance for eyepiece: object at front focal → virtual image at ∞
    # Use finite image at near-point for plotting
    nb_eye = 1.0 / (1.0/nf2 + 1.0/ndn)   # real image through eyepiece

    def T_mat(t):
        return np.array([[1.0, t], [0.0, 1.0]])

    def L_mat(f):
        return np.array([[1.0, 0.0], [-1.0/f, 1.0]])

    def trace(ray0, steps):
        """Trace ray through ordered (label, matrix, dz) steps.
        Returns (z_list, y_list) of ray height at each plane."""
        z_pts = [0.0]
        y_pts = [ray0[0]]
        z = 0.0
        r = ray0.copy()
        for label, M, dz in steps:
            r  = M @ r
            z += dz
            z_pts.append(z)
            y_pts.append(r[0])
        return np.array(z_pts), np.array(y_pts)

    # Steps: (name, matrix, z-advance)
    steps = [
        ("T_g",  T_mat(ng),   ng),
        ("L1",   L_mat(nf1),  0.0),
        ("T_d",  T_mat(nd),   nd),
        ("L2",   L_mat(nf2),  0.0),
        ("T_b",  T_mat(nb_eye), nb_eye),
    ]

    # Initial ray vectors [y, beta]
    ray_a = np.array([0.0,  0.04])   # marginal ray: on-axis, angle
    ray_b = np.array([0.08, 0.0 ])   # chief   ray: off-axis, parallel

    z_a, h_a = trace(ray_a, steps)
    z_b, h_b = trace(ray_b, steps)

    print(f"\n  Marginal ray (a): y0={ray_a[0]}, β0={ray_a[1]}")
    for zi, hi in zip(z_a, h_a):
        print(f"    z={zi:8.3f} mm  →  y={hi:+.5f} mm")

    print(f"\n  Chief ray (b): y0={ray_b[0]}, β0={ray_b[1]}")
    for zi, hi in zip(z_b, h_b):
        print(f"    z={zi:8.3f} mm  →  y={hi:+.5f} mm")

    # -------- Plot --------
    os.makedirs(sets.output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 5))
    fig.patch.set_facecolor('#0d1117')
    ax.set_facecolor('#0d1117')

    ax.plot(z_a, h_a, color='#58a6ff', linewidth=2.0,
            marker='o', markersize=5, label='Marginal ray (a)')
    ax.plot(z_b, h_b, color='#f97583', linewidth=2.0,
            marker='s', markersize=5, label='Chief ray (b)')

    # Optical axis
    ax.axhline(0, color='#8b949e', linewidth=0.8, linestyle='--', alpha=0.7)

    # Lens planes
    lens_z = [ng, ng + nd]
    lens_labels = [f"Objective L₁\n(f₁={nf1} mm)", f"Eyepiece L₂\n(f₂={nf2} mm)"]
    ymin, ymax = ax.get_ylim()
    for lz, ll in zip(lens_z, lens_labels):
        ax.axvline(lz, color='#3fb950', linewidth=1.5, linestyle=':', alpha=0.8)
        ax.text(lz, max(h_a.max(), h_b.max())*1.05,
                ll, ha='center', va='bottom', fontsize=8.5, color='#3fb950')

    # Intermediate image marker
    z_inter = ng
    ax.axvline(z_inter, color='#e3b341', linewidth=1.0, linestyle='--', alpha=0.5)

    # Labels and style
    ax.set_xlabel("z  [mm]", fontsize=12, color='#c9d1d9')
    ax.set_ylabel("Ray height y  [mm]", fontsize=12, color='#c9d1d9')
    M_tot = -(nE * ndn) / (nf1 * nf2)
    ax.set_title(
        f"Compound Microscope – Paraxial Ray Tracing\n"
        f"f₁={nf1} mm  ·  f₂={nf2} mm  ·  E={nE} mm  ·  "
        f"M ≈ {M_tot:.0f}×  (Kloos2007)",
        fontsize=12, color='#c9d1d9', pad=10
    )
    ax.tick_params(colors='#8b949e')
    for spine in ax.spines.values():
        spine.set_edgecolor('#30363d')

    legend = ax.legend(fontsize=11, facecolor='#161b22',
                       edgecolor='#30363d', labelcolor='#c9d1d9')
    ax.grid(True, alpha=0.2, color='#30363d')
    plt.tight_layout()

    fig_path = os.path.join(sets.output_dir, "microscope_ray_trace.png")
    plt.savefig(fig_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.show()
    print(f"\n  Ray-trace plot saved → {fig_path}")

print("\n" + "="*60)
print("test_microscope.py  DONE")
print("="*60)
