#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 07:53:39 2025

@author: g
"""

import numpy as np
from scipy.special import eval_legendre
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt

def legendre_smooth_conserve_total_volume(h_original, theta, R, sigma, l_max=None):
    """
    Smooths an axisymmetric height profile h(theta) on a sphere using
    Legendre polynomial filtering, conserving the total volume defined by
    Integral[ (R+h)^3 * sin(theta) ].

    Args:
        h_original (np.ndarray): 1D array of height values at theta points.
        theta (np.ndarray): 1D array of polar angles (radians) from 0 to pi,
                             corresponding to h_original. Must be sorted.
        R (float): The base radius of the sphere.
        sigma (float): Controls the smoothing strength applied in the Legendre domain
                       to the function f=(R+h)^3. Smaller sigma means more smoothing.
                       Note: Smoothing f=(R+h)^3 might require different sigma/l_max
                       than smoothing h directly for similar visual result on h.
        l_max (int, optional): Maximum Legendre mode (degree) to use for f=(R+h)^3.
                               Defaults to N-1, where N is the number of points.

    Returns:
        np.ndarray: The smoothed height profile h_smoothed(theta).
    """
    N = len(theta)
    if len(h_original) != N:
        raise ValueError("h_original and theta must have the same length.")
    if R <= 0:
        raise ValueError("Base radius R must be positive.")
    if theta[0] != 0 or not np.isclose(theta[-1], np.pi):
         print("Warning: theta should ideally span exactly [0, pi] for best results.")

    if np.any(R + h_original < 0):
        print("Warning: R + h_original is negative at some points. "
              "Cube root might produce unexpected results if f_smoothed goes negative.")

    if l_max is None:
        l_max = N - 1
    if l_max >= N:
         print(f"Warning: l_max ({l_max}) >= N ({N}). Reducing l_max to {N-1}.")
         l_max = N-1

    # --- 1. Define the function f = (R+h)^3 ---
    f_original = (R + h_original)**3

    x = np.cos(theta) # Argument for Legendre polynomials
    sin_theta = np.sin(theta)

    # --- 2. Calculate Legendre Coefficients b_l for f(theta) ---
    b_l = np.zeros(l_max + 1)
    print(f"Calculating Legendre coefficients for f=(R+h)^3 up to l={l_max}...")
    for l in range(l_max + 1):
        P_l = eval_legendre(l, x)
        # Integrand for coefficient calculation: f(theta) * P_l(cos(theta)) * sin(theta)
        integrand = f_original * P_l * sin_theta
        # Integrate over theta using trapezoidal rule
        integral_val = trapezoid(integrand, x=theta)
        # Formula for coefficients: b_l = ( (2l+1)/2 ) * integral
        b_l[l] = (2 * l + 1) / 2.0 * integral_val
        # if l % 10 == 0:
        #      print(f"  Calculated b_{l}")
    print("Finished calculating coefficients for f.")


    # --- 3. Define and Apply Spectral Filter W(l) ---
    l_values = np.arange(l_max + 1)
    # Gaussian spectral filter W(l) = exp(-l(l+1) * sigma^2 / 2)
    if sigma <= 0:
         print("Warning: sigma must be positive for smoothing. No smoothing applied.")
         W_l = np.ones_like(l_values, dtype=float)
    else:
         W_l = np.exp(-l_values * (l_values + 1) * sigma**2 / 2.0)

    # Ensure W(0) = 1 for conservation of the integral of f*sin(theta)
    if not np.isclose(W_l[0], 1.0):
         print("Warning: W(0) is not 1! Volume conservation might be affected. Setting W(0)=1.")
         W_l[0] = 1.0

    # Apply the filter
    b_l_smoothed = b_l * W_l

    # --- 4. Reconstruct Smoothed Profile f_smoothed ---
    f_smoothed = np.zeros_like(f_original)
    print("Reconstructing smoothed profile f_smoothed...")
    # Sum contributions from each Legendre mode
    for l in range(l_max + 1):
        P_l = eval_legendre(l, x)
        f_smoothed += b_l_smoothed[l] * P_l
        # if l % 10 == 0:
        #      print(f"  Added contribution from l={l} for f")

    # --- 5. Recover h_smoothed from f_smoothed ---
    # Use np.cbrt for real cube root, handles potential negative f_smoothed gracefully
    # although f_smoothed should ideally remain positive.
    h_smoothed = np.cbrt(f_smoothed) - R

    # Optional: Check if f_smoothed went negative (undesirable physically)
    if np.any(f_smoothed < 0):
        neg_indices = np.where(f_smoothed < 0)[0]
        print(f"Warning: f_smoothed=(R+h_smoothed)^3 became negative at {len(neg_indices)} points."
              " Check if smoothing is too aggressive or data has issues.")
        # Depending on the application, you might want to clip h_smoothed
        # h_smoothed = np.maximum(h_smoothed, -R + epsilon) # Ensure R+h > epsilon


    print("Finished reconstruction.")
    return h_smoothed

def calculate_total_volume(h, theta, R):
    """Calculates the total volume V = (2pi/3) * Integral[ (R+h)^3 * sin(theta) ]."""
    if np.any(R + h < 0):
         print("Warning during volume calculation: R + h < 0 found.")

    f = (R + h)**3
    integrand = f * np.sin(theta)
    # Integrate f(theta) * sin(theta) d(theta) from 0 to pi
    integral_val = trapezoid(integrand, x=theta)
    # Full volume is (2 * pi / 3) * integral
    volume = (2 * np.pi / 3.0) * integral_val
    return volume

# --- Example Usage ---
if __name__ == "__main__":
    # --- Setup Grid ---
    N_points = 181 # Number of points from North Pole (0) to South Pole (pi)
    theta_grid = np.linspace(0, np.pi, N_points)

    # --- Define Base Radius ---
    R_sphere = 10.0 # Example: Base radius is significant

    # --- Create Sample Data (e.g., a cosine bell bump + smaller bumps) ---
    center_lat_deg = 0  # Equator
    width_deg = 40
    center_theta = np.deg2rad(90 - center_lat_deg)
    width_rad = np.deg2rad(width_deg)

    h_profile_original = np.zeros_like(theta_grid)
    mask = np.abs(theta_grid - center_theta) < (width_rad / 2.0)
    h_profile_original[mask] = 1.5 * (1 + np.cos(2 * np.pi * (theta_grid[mask] - center_theta) / width_rad)) # Height is not << R

    # Add another smaller bump
    center_theta_2 = np.deg2rad(30)
    width_rad_2 = np.deg2rad(20)
    mask2 = np.abs(theta_grid - center_theta_2) < (width_rad_2 / 2.0)
    h_profile_original[mask2] += 0.8 * (1 + np.cos(2 * np.pi * (theta_grid[mask2] - center_theta_2) / width_rad_2))

    # Add some noise
    noise_level = 0.1
    h_profile_original += noise_level * (np.random.rand(N_points) - 0.5) * (R_sphere * 0.05) # Noise relative to R
    # Ensure R+h > 0 , maybe clip h if it goes too negative
    h_profile_original = np.maximum(h_profile_original, -R_sphere * 0.95)


    # --- Define Smoothing Parameter ---
    # Sigma controls smoothing of f=(R+h)^3. You might need to adjust this
    # experimentally compared to smoothing h directly.
    # Let's try smoothing features around 15 degrees.
    smoothing_angle_deg = 15.0
    sigma_smooth = np.deg2rad(smoothing_angle_deg)

    # --- Perform Smoothing ---
    # Reduce l_max slightly more aggressively maybe, as f=(R+h)^3 has more power in high modes
    h_profile_smoothed = legendre_smooth_conserve_total_volume(
        h_profile_original,
        theta_grid,
        R = R_sphere,
        sigma=sigma_smooth,
        l_max = N_points // 2 # Example max L
    )

    # --- Check Volume Conservation ---
    volume_original = calculate_total_volume(h_profile_original, theta_grid, R=R_sphere)
    volume_smoothed = calculate_total_volume(h_profile_smoothed, theta_grid, R=R_sphere)

    print("-" * 30)
    print(f"Base Sphere Radius R: {R_sphere}")
    print(f"Original Total Volume: {volume_original:.8f}")
    print(f"Smoothed Total Volume: {volume_smoothed:.8f}")
    print(f"Relative Volume Change: {(volume_smoothed - volume_original) / volume_original:.6e}")
    print("-" * 30)

    # --- Plot Results ---
    plt.figure(figsize=(12, 7))

    plt.subplot(2, 1, 1)
    plt.plot(np.rad2deg(theta_grid), R_sphere + h_profile_original, '.-', label='Original R+h', alpha=0.7)
    plt.plot(np.rad2deg(theta_grid), R_sphere + h_profile_smoothed, 'r-', label=f'Smoothed R+h (sigma={smoothing_angle_deg:.1f} deg)', linewidth=2)
    plt.ylabel('Total Radius r(theta)')
    plt.title(f'Axisymmetric Profile Smoothing Conserving V = (2pi/3) Integral[(R+h)^3 sin(theta)] (R={R_sphere})')
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(np.rad2deg(theta_grid), h_profile_original, '.-', label='Original h', alpha=0.7)
    plt.plot(np.rad2deg(theta_grid), h_profile_smoothed, 'r-', label=f'Smoothed h', linewidth=2)
    plt.xlabel('Polar Angle theta (degrees)')
    plt.ylabel('Height Profile h(theta)')
    plt.legend()
    plt.grid(True)
    # plt.ylim(bottom=min(0, plt.ylim()[0])) # Ensure y-axis starts at or below 0 if h can be negative

    plt.tight_layout()
    plt.show()