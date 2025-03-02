#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 12:06:35 2024

@author: Kumar Gaurav, Gemini and Marco Fenucci
"""

import numpy as np
from gaurav import Output_File, Height, Landslides, shape_gen
from collisions import G, wobblecalcf
import math

"""Main YORP function that updates the spin state."""

def YORP(target, parameters, impacttime, oldtime, myomega):
    print("Simulating YORP")
    if  not target.YORP:
        return
    # which omega to use because it is being changed by the yorp
    wobblecalcf(target, impacttime, oldtime)

    t_yorp = oldtime
    file = Output_File(parameters, "output", ["output.yorp"])
    

    h_yorp = float(parameters["h_yorp"])
    while impacttime > t_yorp+1:
        if impacttime-t_yorp < h_yorp:
            h_yorp = impacttime-t_yorp
        t_yorp += h_yorp

    # Check for collisional reorientation if P >= 1000 h
        if (2*np.pi / target.omega[2]) / 3600 >= 1000.0:
            reorient_flag = reor_probability(parameters, target)

            if reorient_flag:
                print("Too slow spinning causing collisonal reorientation")
                spin_state_new(target)
                target.coeff_f, target.coeff_g = shape_gen(target.K,parameters['Nature'])

        rk4(parameters, target)

        with open(file, "a") as output_file:
            output_file.write(f"{t_yorp:12.6e} { target.omega[2]:12.8e} {target.obliq:12.8e}\n")

        omegaLimit = (G*4/3 * math.pi*target.dens)**0.5
        
       
        if target.omega[2] > 0.9*omegaLimit:
            print("Too fast spinning causing landslides")
            myomega.append([t_yorp, target.omega[2]])
            # parameters["uni_h"]=min(max((target.omega[2]-0.9*omegaLimit)/(omegaLimit)*(0.2/float(parameters["epsilon"])),1),10)
            Height(parameters, target)
            Landslides(target, parameters,t_yorp,myomega)
    myomega.append([t_yorp, target.omega[2]])
    print("Omega after the yorp effect:", target.omega[2])

# %%


def rk4(parameters, target):
    """
    Fourth-order Runge-Kutta integration for a system of two equations.

    Args:
        y (np.ndarray): Initial values of the two state variables.
        h (float): Step size.
        rpar (np.ndarray): Array of 5 parameters for the 'yorp_vf' function.

    Returns:
        np.ndarray: Estimated values of the state variables after the step.
    """

    # Initialize arrays for the intermediate k values (f1, f2, f3, f4 in Fortran)
    h = float(parameters["h_yorp"])
    f = np.zeros((4, 2))
    y = np.zeros(2)
    y[0], y[1] = target.omega[2], target.obliq
    # Compute the k values using the 'yorp_vf' function (replace with your actual function)
    f[0] = yorp_vf(y, parameters, target)
    f[1] = yorp_vf(y + h * f[0] / 2.0, parameters, target)
    f[2] = yorp_vf(y + h * f[1] / 2.0, parameters, target)
    f[3] = yorp_vf(y + h * f[2], parameters, target)

    # Update the state variables using the weighted average of the k values
    ytph = y + h * (f[0] + 2.0 * f[1] + 2.0 * f[2] + f[3]) / 6.0

    #  NOTE: this is done to handle the cases that goes asymptotically too fast!
    #       Indeed, in the step of 50 years is too long and they are decelerating very fast,
    #       we simply set a rotation period of 1000 h and put \gamma to the closest asymptotic
    #       value.
    # NOTE: the failure of the step is recognized because \omega becomes negative, which does not
    #       make any physical sense!

    if abs(ytph[1]) < 1e-10:
        ytph[1] = 0.0

    if 2*np.pi / ytph[0]/3600 > 1000.0:
        ytph[0] = 2*np.pi / (1000.0 * 3600)

    if ytph[0] < 0.0:
        ytph[0] = 2*np.pi / (1000.0 * 3600)
        d = np.abs([y[1], y[1] - np.pi / 2.0, y[1] - np.pi])
        i = np.argmin(d)
        ytph[1] = [0.0, np.pi/2.0, np.pi][i]

    target.omega[2] = ytph[0]
    target.obliq = ytph[1]

# %%


def yorp_vf(x, parameters, target):
    """
    Calculates the YORP effect vector field (vf) based on the state variables (x) and parameters (rpar).

    Args:
        x (np.ndarray): Array of two state variables:
            - x[0]: Angular velocity (omega) in 1/d
            - x[1]: Obliquity (gamma) in radians

    Returns:
        np.ndarray: Array of two components of the vector field:
            - vf[0]: Time derivative of omega
            - vf[1]: Time derivative of gamma
    """

    # Normalize gamma to [0, pi]
    x[1] = x[1] + 180 if x[1] < 0 else x[1]-180 if x[1] > 180 else x[1]

    # Compute f and g (using the comp_f_g function, which needs to be implemented separately)
    f = target.f_spline(x[1])
    g = target.g_spline(x[1])

    # Compute rescaling factor
    fact = float(parameters["c_YORP"]) * (2500.0 / target.dens) * \
        (2.5 / target.sma)**2.0 * (2000.0 / target.d)**2.0

    vf = np.array([
        f / (1e+6)  * fact * target.coeff_f,
        g*180/np.pi / x[0] / (1e+6)  * fact * target.coeff_g
    ])
    

    # Special cases: Freeze evolution of omega for extreme rotation periods
    if 2*np.pi / (x[0])/3600 > 1e3:
        vf = np.zeros(2)
        
    if parameters['Nature']=='Increasing':
        vf[0]=abs(vf[0])
    elif parameters['Nature']=='Decreasing':
        vf[0]=-abs(vf[0])
        
    return vf

# %%


def reor_probability(parameters, target):
    """
    Determines the probability of reorientation based on parameters and random number generation.

    Args:
        omega (float): Angular velocity (assumed in 1/d).
        D (float): Diameter (assumed in meters).
        h_yorp (float): Timestep for the YORP effect integration (assumed in days).

    Returns:
        bool: True if reorientation occurs, False otherwise.
    """

    # Parameters
    B = 84.5e3
    beta1 = 5.0 / 6.0
    beta2 = 4.0 / 3.0
    omega0 = 2*np.pi / (5.0 * 3600)
    D0 = 2.0

    # Calculate characteristic timescale for collisional reorientation
    t_reor = float(parameters["c_reor"]) * B * \
        (target.omega[2] / omega0)**beta1 * (target.d / D0)**beta2

    # Generate a random number between 0 and 1
    rand_num = np.random.rand()

    # Decide if reorientation happens based on probability
    reorient_flag = rand_num < 1.0 - \
        np.exp(-float(parameters["h_yorp"]) / t_reor)

    return reorient_flag

# %%


def spin_state_new(target):
    """
    generate a new initial condition for the spin state.
   !          The obliquity is such that cos\gamma is equidistributed,
   !          while the rotation period is Gaussian with peak at 10 h and
   !          standard deviation of 4 h. Values smaller than 4 h and larger than 24 h 
   !          are discarded

    Returns:
        tuple: A tuple containing:
            - gamma_new (float): New obliquity in degree.
            - omega_new (float): New rotation rate in 1/s.
    """

    # Parameter
    P_peak = 10.0  # hours

    # Generate a new value of gamma (obliquity) with a cosine distribution
    rand_num = np.random.rand()
    target.obliq = np.arccos(2.0 * rand_num - 1.0)*180/np.pi

    # Generate rotation period from a Maxwellian distribution (approximated using Gaussians)
    # Three standard normal random numbers
    rand_gauss = np.random.normal(size=3)
    # Maxwellian distribution
    P = P_peak * np.sqrt(rand_gauss[0]**2 +
                         rand_gauss[1]**2 + rand_gauss[2]**2)

    # Convert period (P) to rotation rate (omega_new) in 1/day
    target.omega[2] = 1.0 / (P * 3600)
