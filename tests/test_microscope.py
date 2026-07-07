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
    600  ray_trace               – interactive ray tracing with sliders (replaces static plot)
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
# 600 – Interactive Ray Tracing with Sliders (replaces static plot)
# =============================================================================

if "ray_trace" in sets.flow:
    print("\n" + "="*60)
    print("600 · Interactive Ray Tracing with Sliders")
    print("="*60)
    print("A window will open with sliders for f1, f2, E, and object height.")
    print("Adjust the sliders to see the rays and image height update in real time.")
    print("Close the figure window to continue.\n")

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    # ------------------------------------------------------------------
    # Helper functions (numerical)
    # ------------------------------------------------------------------

    def T_mat(t):
        return np.array([[1.0, t], [0.0, 1.0]])

    def L_mat(f):
        return np.array([[1.0, 0.0], [-1.0/f, 1.0]])

    def compute_system(f1, f2, E, d_n=250.0):
        d = f1 + E + f2
        g = 1.0 / (1.0/f1 - 1.0/(f1 + E))
        M_obj = -(f1 + E) / g
        M_eye = d_n / f2
        M_total = M_obj * M_eye
        L1 = L_mat(f1)
        L2 = L_mat(f2)
        Td = T_mat(d)
        S_dbl = L2 @ Td @ L1
        return g, d, M_obj, M_eye, M_total, S_dbl

    def trace_ray(ray0, steps):
        z_pts = [0.0]
        y_pts = [ray0[0]]
        r = ray0.copy()
        z = 0.0
        for M, dz in steps:
            r = M @ r
            z += dz
            z_pts.append(z)
            y_pts.append(r[0])
        return np.array(z_pts), np.array(y_pts)

    def trace_microscope(f1, f2, E, h_obj, d_n=250.0):
        g, d, M_obj, M_eye, M_total, S_dbl = compute_system(f1, f2, E, d_n)
        steps = [
            (T_mat(g),      g),
            (L_mat(f1),     0.0),
            (T_mat(d),      d),
            (L_mat(f2),     0.0),
            (T_mat(d_n),    d_n),   # screen at near‑point after eyepiece
        ]
        ray_marg = np.array([0.0, 0.04])
        z_marg, y_marg = trace_ray(ray_marg, steps)
        ray_chief = np.array([h_obj, 0.0])
        z_chief, y_chief = trace_ray(ray_chief, steps)
        return z_marg, y_marg, z_chief, y_chief, M_total, g, d

    # ------------------------------------------------------------------
    # Initial parameters
    # ------------------------------------------------------------------
    INIT_f1 = 4.0
    INIT_f2 = 25.0
    INIT_E  = 160.0
    INIT_h  = 0.1
    d_n = 250.0

    # ------------------------------------------------------------------
    # Create figure and axes
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 6))
    plt.subplots_adjust(bottom=0.35)          # space for sliders

    # Light theme for readability
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Slider positions: [left, bottom, width, height]
    # Stack from bottom up: h_obj, E, f2, f1
    slider_config = [
        ('f1', 'f₁ (mm)', 0.28, 1.0, 100.0, 0.5, INIT_f1),
        ('f2', 'f₂ (mm)', 0.22, 10.0, 100.0, 1.0, INIT_f2),
        ('E',  'E (mm)',  0.16, 100.0, 300.0, 5.0, INIT_E),
        ('h',  'h_obj (mm)', 0.10, 0.01, 1.0, 0.01, INIT_h),
    ]

    sliders = {}
    for key, label, bottom, vmin, vmax, step, init in slider_config:
        ax_slider = plt.axes([0.15, bottom, 0.65, 0.03])
        slider = Slider(
            ax_slider,
            label,
            vmin, vmax,
            valinit=init,
            valstep=step,
            valfmt='%1.2f'
        )
        # Make labels and ticks clearly visible
        ax_slider.xaxis.label.set_fontsize(11)
        ax_slider.xaxis.label.set_color('black')
        ax_slider.tick_params(colors='black', labelsize=9)
        sliders[key] = slider

    # ------------------------------------------------------------------
    # Update function
    # ------------------------------------------------------------------
    def update_plot(val):
        f1 = sliders['f1'].val
        f2 = sliders['f2'].val
        E  = sliders['E'].val
        h  = sliders['h'].val

        ax.clear()
        ax.set_facecolor('white')

        z_marg, y_marg, z_chief, y_chief, M_total, g, d = trace_microscope(f1, f2, E, h)

        image_height = M_total * h

        # Plot rays
        ax.plot(z_marg, y_marg, color='blue', linewidth=2.0,
                marker='o', markersize=4, label='Marginal ray')
        ax.plot(z_chief, y_chief, color='red', linewidth=2.0,
                marker='s', markersize=4, label='Chief ray')

        # Optical axis
        ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', alpha=0.7)

        # Lens planes
        lens_z = [g, g + d]
        lens_labels = [f"Objective L₁\n(f₁={f1:.1f} mm)", f"Eyepiece L₂\n(f₂={f2:.1f} mm)"]
        ymax = max(abs(y_marg.max()), abs(y_marg.min()), abs(y_chief.max()), abs(y_chief.min()))
        ymax = max(ymax, 0.5) * 1.2
        for lz, ll in zip(lens_z, lens_labels):
            ax.axvline(lz, color='green', linewidth=1.5, linestyle=':', alpha=0.8)
            ax.text(lz, ymax*0.95, ll, ha='center', va='bottom', fontsize=8.5, color='green')

        # Intermediate image plane
        z_inter = g + f1 + E
        ax.axvline(z_inter, color='orange', linewidth=1.0, linestyle='--', alpha=0.5)

        # Info box
        M_classical = -(E * d_n) / (f1 * f2)
        info_text = (f"f₁ = {f1:.2f} mm  |  f₂ = {f2:.2f} mm  |  E = {E:.1f} mm\n"
                     f"g = {g:.2f} mm  |  d = {d:.2f} mm  |  M_total = {M_total:.2f}×\n"
                     f"Object height = {h:.3f} mm  |  Image height = {image_height:.3f} mm")
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                fontsize=10, color='black', verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray'))

        # Labels and style
        ax.set_xlabel("z  [mm]", fontsize=12, color='black')
        ax.set_ylabel("Ray height y  [mm]", fontsize=12, color='black')
        ax.set_title(f"Compound Microscope – Interactive Ray Tracing\n"
                     f"M ≈ {M_classical:.0f}×  (classical formula)",
                     fontsize=12, color='black', pad=10)
        ax.tick_params(colors='black')
        for spine in ax.spines.values():
            spine.set_edgecolor('black')
        ax.grid(True, alpha=0.3, color='gray')
        ax.set_ylim(-ymax, ymax)

        ax.legend(fontsize=11, facecolor='white', edgecolor='gray')

        fig.canvas.draw_idle()

    # Connect sliders
    for slider in sliders.values():
        slider.on_changed(update_plot)

    # Initial draw
    update_plot(None)

    # Show the interactive window
    plt.show()

    print("Interactive session ended.")

print("\n" + "="*60)
print("test_microscope.py  DONE")
print("="*60)