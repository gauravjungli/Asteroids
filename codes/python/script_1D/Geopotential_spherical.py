#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 16:45:55 2025

@author: g
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.constants import G

def get_radius_at_angle(theta, C_potential, M_total, omega):
    """
    Solves for the radius r at a given angle (colatitude) theta for a
    specific geopotential constant C.

    The equation is:
    (0.5 * omega^2 * sin(theta)^2) * r^3 + C * r + G * M_total = 0
    This is a cubic equation in the form ar^3 + br^2 + cr + d = 0
    """
    
    sin_theta = np.sin(theta)
    
    # Coefficients of the cubic equation for r
    a = 0.5 * omega**2 * sin_theta**2
    b = 0.0
    c = C_potential
    d = G * M_total
    
    # At the pole (theta=0), sin(theta)=0, and the equation simplifies
    # if np.isclose(a, 0):
    #     if C_potential == 0:
    #         # Should not happen in a real scenario
    #         return np.nan 
    #     r = -d / c
    #     return r

    coeffs = [a, b, c, d]
    
    # Find the roots of the cubic equation
    roots = np.roots(coeffs)
    
    # We need the single, positive, real root.
    # Physical solutions must be real and positive.
    real_positive_roots = [r.real for r in roots if np.isreal(r) and r.real > 0]
    
    if len(real_positive_roots) == 1:
        print("one positive root",real_positive_roots)
        return real_positive_roots[0]
    elif len(real_positive_roots) > 1:
    #    print("more than one positive root",theta)
        # This case is physically unlikely for this problem setup,
        # but we'll return the smallest positive root (closest to the sphere)
        return min(real_positive_roots)
    else:
        # No physical solution found
        print("No root found")
        return np.nan

def calculate_fluid_volume(C_potential, M_total, omega, R_sphere):
    """
    Calculates the total volume of the fluid layer for a given
    geopotential C. It does this by integrating the volume
    between the surface r(theta) and the solid sphere R_sphere.
    
    Volume V = integral(phi=0 to 2pi) integral(theta=0 to pi) integral(r=R_sphere to r(theta)) [ r^2 * sin(theta) dr d(theta) d(phi) ]
    V = (2*pi/3) * integral(theta=0 to pi) [ (r(theta)^3 - R_sphere^3) * sin(theta) d(theta) ]
    """
    
    # Define the integrand function for the integration
    def integrand(theta):
        r_at_theta = get_radius_at_angle(theta, C_potential, M_total, omega)
        
        if np.isnan(r_at_theta) or r_at_theta < R_sphere:
            # If the calculated radius is inside the solid sphere,
            # the fluid volume at this angle is 0.
           # print("small value for r",r_at_theta)
            return 0.0
            
        return (r_at_theta**3 - R_sphere**3) * np.sin(theta)

    # Perform the numerical integration from pole (0) to pole (pi)
    try:
        volume, error = quad(integrand, 0, np.pi)
        return (2 * np.pi/3) * volume
    except Exception as e:
        print(f"Integration failed for C={C_potential}: {e}")
        return np.nan

def volume_error_function(C_potential, M_total, omega, R_sphere, target_volume):
    """
    This is the error function we want to minimize (find the root of).
    It returns the difference between the calculated volume and the target volume.
    """
    if isinstance(C_potential, (list, np.ndarray)):
        C_scalar = C_potential[0]
    else:
        C_scalar = C_potential
    
    calculated_volume = calculate_fluid_volume(C_scalar, M_total, omega, R_sphere)
    
    if np.isnan(calculated_volume):
        # Return a large error if integration fails
        return 1e30 
        
    error = calculated_volume - target_volume
    print(f"Trying C = {C_potential[0]:.5e}, Calculated Vol = {calculated_volume:.5e}, Target Vol = {target_volume:.5e}, Error = {error:.5e}")
    return error

def plot_shape(C_final, M_total, omega, R_sphere):
    """
    Plots a 2D cross-section of the solid sphere and the
    calculated equipotential fluid surface.
    """
    print("\nGenerating plot...")
    
    # --- Plot the Fluid Surface ---
    # Generate angles from 0 to 2*pi for a full circle
    thetas_plot = np.linspace(1.570796326794896653e-03,3.140021857262998317e+00, 1000)
    
    # Calculate the radius at each angle
    radii_plot = [max(get_radius_at_angle(t, C_final, M_total, omega) -R_sphere,(0.0075)**2*R_sphere) for t in thetas_plot]
    
    # Convert from polar (r, theta) to Cartesian (x, z)
    # Note: We use theta for colatitude in physics, but for plotting (x,y)
    # it's easier to think of it as the standard angle.
    # Let's be careful. theta=0 is north pole (z-axis).
    # x = r * sin(theta)
    # z = r * cos(theta)
    
    x_surface = [r * np.sin(t) for r, t in zip(radii_plot, thetas_plot)]
    z_surface = [r * np.cos(t) for r, t in zip(radii_plot, thetas_plot)]

    # --- Plot the Solid Sphere ---
    x_sphere = [R_sphere * np.sin(t) for t in thetas_plot]
    z_sphere = [R_sphere * np.cos(t) for t in thetas_plot]
    
    # --- Create the Plot ---
    #plt.figure(figsize=(10, 10))
    #plt.plot(x_surface, z_surface, 'b-', label=f'Fluid Surface (Equipotential)')
    #plt.plot(x_sphere, z_sphere, 'k--', label=f'Solid Sphere (Radius = {R_sphere:,.0f} m)')

    my_array = np.array(radii_plot)/R_sphere/0.0075
    plt.plot(thetas_plot,my_array,linewidth=2)
    np.save('/home/g/Asteroids/output/check_5/run1/height.npy',my_array)
    # --- Formatting ---
   # plt.title('Geopotential Surface of a Rotating Sphere', fontsize=16)
   # plt.xlabel('X (distance from rotation axis, m)', fontsize=12)
   # plt.ylabel('Z (distance along rotation axis, m)', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.7)
    
    # Use equal axis scaling to see the shape without distortion
    #plt.axis('equal') 
    
    # Calculate and display the equatorial bulge
    r_equator = get_radius_at_angle(np.pi/2, C_final, M_total, omega)
    r_pole = get_radius_at_angle(0, C_final, M_total, omega)
    bulge = r_equator - r_pole
    print(f"Equatorial Radius (r_eq): {r_equator:,.2f} m")
    print(f"Polar Radius (r_pole):     {r_pole:,.2f} m")
    print(f"Equatorial Bulge (r_eq - r_pole): {bulge:,.2f} m")
    
    # plt.text(0.05, 0.05,
    #          f'Equatorial Radius: {r_equator:,.1f} m\n'
    #          f'Polar Radius: {r_pole:,.1f} m\n'
    #          f'Equatorial Bulge: {bulge:,.1f} m',
    #          transform=plt.gca().transAxes,
    #          bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.8))
             
    plt.show()


def main():
    # --- 1. DEFINE PHYSICAL CONSTANTS ---
    # (Using Earth and its oceans as a default example)
    

    
    # Density of the fluid (e.g., seawater)
    rho_fluid = 1250  # kg/m^3
    
    # Radius of the solid sphere (e.g., average radius of solid Earth)
    R_sphere = 250*(1-0.0075)  # m (6368 km)
    
    # Mass of the solid sphere (e.g., solid Earth)
    M_sphere = 4/3*np.pi*rho_fluid*R_sphere**3  # kg 
    
    # Angular velocity of the sphere (e.g., Earth's rotation)
    # omega = 2 * pi / (Period in seconds)
    omega = 2 * np.pi / (10* 3600)  # rad/s
    
    # Mass of the fluid (e.g., Earth's oceans)
    # This is the 'M' from your request.
    M_fluid =  4*np.pi*rho_fluid*R_sphere**3*0.0075 # kg

    # --- 2. CALCULATE DERIVED VALUES ---
    M_total = M_sphere + M_fluid
    target_volume = M_fluid / rho_fluid

    print("--- Geopotential Surface Calculator ---")
    print(f"Solid Sphere Mass (M_sphere): {M_sphere:.3e} kg")
    print(f"Fluid Mass (M_fluid):         {M_fluid:.3e} kg")
    print(f"Total Mass (M_total):         {M_total:.3e} kg")
    print(f"Solid Sphere Radius (R_sphere): {R_sphere:,.0f} m")
    print(f"Fluid Density (rho_fluid):    {rho_fluid:.1f} kg/m^3")
    print(f"Target Fluid Volume (V_f):    {target_volume:.3e} m^3")
    print(f"Angular Velocity (omega):     {omega:.3e} rad/s")

    # --- 3. FIND THE GEOPOTENTIAL CONSTANT C ---
    
    # We need a good initial guess for the potential C
    # Let's calculate the radius of a non-rotating sphere with the same volume
    R_approx_no_rotation = (M_total* 3 / (4 * np.pi*rho_fluid))**(1/3)
    
    # The potential C will be close to the potential of this non-rotating sphere
    # C = U_g + U_c. Let's guess U_c is small.
    C_guess = -G*M_sphere/(1.0075*R_sphere) -0.5*omega**2*(1.0075*R_sphere)**2
    
    print(f"\nInitial guess for C: {C_guess:.3e}")
    print("Solving for the correct geopotential constant C... (this may take a moment)")
    
    # Use fsolve to find the root of the volume_error_function
    # This finds the value of C where (calculated_volume - target_volume) = 0
    C_final =  fsolve(volume_error_function,
                     C_guess,
                     args=(M_sphere, omega, R_sphere, target_volume))
    #
    if isinstance(C_final, np.ndarray):
        C_final = C_final[0] # fsolve returns an array

    print(f"Solution found!")
    print(f"Final Geopotential Constant (C): {C_final:.5e} J/kg")
    
    # --- 4. PLOT THE FINAL SHAPE ---
    plot_shape(C_final, M_sphere, omega, R_sphere)

if __name__ == "__main__":
    main()