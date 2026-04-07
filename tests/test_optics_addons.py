# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Tue Dec  3 15:21:15 2024

@author: yubuntu
"""
import numpy as np
from scipy.constants import pi
from scipy.integrate import dblquad
from numba import jit
import matplotlib.pyplot as plt

# Constants
I = 1j  # Imaginary unit
k = 2 * pi / 632e-6  # Wave number (example wavelength 632e-6 mm)
l = 632e-6  # Wavelength in mm
z = 800.0   # Propagation distance in mm
Uapr, Eli = 1.0, 1.0  # Amplitudes (example values)
screen_factor = 2.5  # Scaling factor for visualization

# Define the Fresnel integrand


@jit(nopython=True)
def fresnel_integrand(x, y, z, x0, y0):
    """
    Uapr*Eli*exp( I*k/(2*z)*((x-x0)**2+(y-y0)**2) )
    """
    return (
        Uapr
        * Eli
        * np.exp(I*k/(2*z)*((x-x0)**2+(y-y0)**2))
    ).real  # Use `.real` for compatibility with Numba

# Double integral for the Fresnel integral


def compute_fresnel_integral(x, y, z, x0min, x0max, y0min, y0max):
    integral, _ = dblquad(
        lambda x0, y0: fresnel_integrand(x, y, z, x0, y0),
        y0min, y0max,  # y0 limits
        lambda _: x0min, lambda _: x0max  # x0 limits
    )
    prefactor = (np.exp(I * k * z) / (I * l * z))
    return np.abs(prefactor * integral)**2  # Intensity (absolute square)


# Observation plane (grid)
n = 10  # Number of grid points per axis
nLx, nLy = 2.0, 2.0  # Half-width of the observation plane
x_vals = np.linspace(-nLx*screen_factor, nLx*screen_factor, n)
y_vals = np.linspace(-nLy*screen_factor, nLy*screen_factor, n)
X, Y = np.meshgrid(x_vals, y_vals)

# Compute the Fresnel integral across the observation plane
Z = np.zeros_like(X)
x0min, x0max = -nLx*screen_factor, nLx*screen_factor  # Integration limits for x0
y0min, y0max = -nLy*screen_factor, nLy*screen_factor  # Integration limits for y0

for i in range(n):
    for j in range(n):
        Z[i, j] = compute_fresnel_integral(
            X[i, j], Y[i, j], z, x0min, x0max, y0min, y0max)

# Normalize intensity
brightness = 1.0
#Z = Z / np.max(Z)

# Plotting 2D Diffraction Intensity
fig = plt.figure(figsize=(5, 5))
ax1 = fig.add_subplot(111)
ax1.imshow(Z, cmap=plt.cm.gray, interpolation='bilinear',
           origin='lower', vmin=np.min(Z), vmax=brightness * np.max(Z))
ax1.set_xticks(np.linspace(0, n, 5))
ax1.set_xticklabels([-nLx * screen_factor, -nLx * screen_factor *
                    0.5, 0, nLx * screen_factor * 0.5, nLx * screen_factor])
ax1.set_yticks(np.linspace(0, n, 5))
ax1.set_yticklabels([-nLy * screen_factor, -nLy * screen_factor *
                    0.5, 0, nLy * screen_factor * 0.5, nLy * screen_factor])
plt.savefig("fresnel_diffraction.png", format="png",
            dpi=600, bbox_inches='tight')
plt.show()
