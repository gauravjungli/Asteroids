#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 12:06:35 2024

@author: Kumar Gaurav, Gemini and Marco Fenucci
"""

import numpy as np
from landslides import  Landslides
from Height import Height
from collisions import G, wobblecalcf
import math
import subprocess
from IO import Output_File
from scipy.interpolate import  CubicSpline
import sys

"""Main YORP function that updates the spin state."""

def YORP( parameters, target, impacttime, oldtime, myomega):
 #   print("Simulating YORP")
    check_time = math.ceil(oldtime / 1000) * 1000
    if  not target.YORP:
        return
    # which omega to use because it is being changed by the yorp
    wobblecalcf(target, impacttime, oldtime)

    t_yorp = oldtime
    print_time = oldtime
    file = Output_File(target, "output", ["output.yorp"])
    
    h_yorp = float(parameters["h_yorp"])
    
    omega_crit, omega_shed =  Critical_omega(target)
    
    if (omega_shed<omega_crit):
        
        print("Shedding starts before the landslide.")
    
    while impacttime > t_yorp+1:
        
        if (t_yorp>print_time):
            print_time +=1000
            with open(file, "a") as output_file:
                output_file.write(f"{t_yorp:12.6e} { target.omega[2]:12.8e} {target.obliq:12.8e} {omega_crit:12.8e} {omega_shed:12.8e}\n")
            
        
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
        
       
        if not target.collision:
            print(min(99,round(t_yorp/impacttime*100)),file=sys.__stdout__,flush=True)
            sys.stdout.flush()
       
        if t_yorp > check_time and target.omega[2] > 1.0*omega_crit:
            
            target.fast_rotation_flag = True
            
            myomega.append([t_yorp, target.omega[2]])
            Height(target)
            Landslides( parameters,target,t_yorp,myomega)
            
            target.fast_rotation_flag = False
            check_time = t_yorp + 1000
            
            omega_crit, omega_shed =  Critical_omega(target)
            
            
            if (omega_shed<omega_crit):
                
                print("Shedding starts before the landslide.")
                
    
            
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

    # Compute f and g 
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

#%%

def shape_gen(K):
    """
    generate random coefficients to determine the functions f,g. To
    this purpose, we use the statistics presented in Capek & Vokrouhlicky 2004

    Args:
        K (float): Input value used for determining probabilities.

    Returns:
        tuple: A tuple containing the generated values of coeff1 and coeff2.
    """

    K_t = 0.005
    max_g = 1.8 / 1.1
    min_g = 0.4 / 1.1
    std_g = abs(1.1 - 1.8 / 1.1) / 3.0
    max_f = 3.0 / 2.0
    min_f = -3.0 / 2.0
    std_f = 0.5**2.0

    if K <= K_t:
        # 80% probability to reach 0/180 (g, coeff2)
        # 40% probability to accelerate (f, coeff1)
        coeff2 = np.random.normal(1.0, std_g)
        coeff2 = np.clip(coeff2, min_g, max_g)  # Ensure coeff2 is within bounds

        if np.random.rand() > 0.8:
            coeff2 = -coeff2  # Switch to reaching 90 degrees

        # ! If the asymptotic state is 90 (i.e. coeff2 < 0), then we always decelerate
        # ! the rotation rate! See Capek & Vokrouhlicky 2004, Fig 7.
        if coeff2 < 0.0:
            # Always decelerate if asymptotic state is 90
            sgn = -1.0
        else:
            # 60/40 probability to decelerate/accelerate for 0/180
            sgn = 1.0 if np.random.rand() <= 0.6 else -1.0 

        coeff1 = sgn * (np.random.normal(1.0, std_f))
        coeff1 = np.clip(coeff1, min_f, max_f)

    else:  # K > K_t
        # 100% probability to reach 0/180 (g, coeff2)
        # 50% probability to accelerate (f, coeff1)

        sgn = 1.0 if np.random.rand() > 0.5 else -1.0  # Choose sign for coeff1
        coeff1 = sgn * (np.random.normal(1.0, std_f))
        coeff1 = np.clip(coeff1, min_f, max_f)
       
        coeff2 = np.random.normal(1.0, std_g)
        coeff2 = np.clip(coeff2, min_g, max_g)
    print(f"The value of the coefficients are {coeff1} and {coeff2}")
    return coeff1, coeff2 

#%% 

"""old version of YORP which uses orbit9, a fortran code, Please use the newer version of the YORP implemented in YORP.py """

def Yorp(target,parameters,impacttime,oldtime,myomega):
    
    wobblecalcf(target,impacttime,oldtime)  #which omega to use because it is being changed by the yorp
    while impacttime > oldtime + 10:
        subprocess.run(["make"], cwd="../OrbFit/tests/gaurav",stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)
        try:
            yark = []
            with open("../OrbFit/tests/gaurav/yarkovsky.in", "r") as file:
                yark = [list( line.strip().split()) for line in file]
        except FileNotFoundError:
            print("Failed in importing the data from yarkovsky.in")
            return
        yark[0][5] = target.obliq
        yark[0][6] = 2 * math.pi / (target.omega[2] * 3600)
        try:
            with open("../OrbFit/tests/gaurav/yarkovsky.in", "w") as file:
                for row in yark:
                    file.write("\t".join(map(str, row)) + "\n")
        except IOError:
            print("Failed in exporting the yarkovsky.in")
            return
        exitcode = subprocess.run(["./orbit9.x"], cwd="../OrbFit/tests/gaurav").returncode
        if exitcode != 0:
            print("Failed in running orbit9")
            return
        try:
            omegOrb = []
            with open("../OrbFit/tests/gaurav/clo0.yorp", "r") as file:
                omegOrb = [list( line.strip().split()) for line in file]
        except FileNotFoundError:
            print("Failed in importing the data from orbit9")
            return
        omegaLimit = (G*4/3* math.pi*target.dens)**0.5
        index = min(round((impacttime - oldtime) / 50), 2000)
        target.omega[2] = 2 * math.pi / (float(omegOrb[index][1]) * 3600)
        target.obliq = float(omegOrb[index][2])
        oldtime = oldtime + 10**5
        myomega.append([min(impacttime,oldtime), target.omega[2]])
        if target.omega[2] > 0.9*omegaLimit:
            print("Too fast spinning causing landslides")
            Height(parameters,target)
            Landslides(target,parameters,min(impacttime,oldtime),myomega)
        
        
    print("Omega after the yorp effect:", target.omega[2])


#%%

def read_f_g_spline(target):
    """
    Reads data from files, creates spline interpolations for f and g functions.

    This function assumes the following file structure:

    - input/yorp_f.txt: Contains gamma (in degrees) and corresponding f values.
    - input/yorp_g.txt: Contains gamma (in degrees) and corresponding g values.

    Returns:
        tuple: Two CubicSpline objects representing the spline interpolations for f and g.
    """

    # --- Read data for the f function ---
    f = []
    mydir=Output_File(target,"input",["yorp_f.txt"])
    with open(mydir, 'r') as file:
        for line in file:
            x, y = map(float, line.split(","))
            f.append([x,y])
          
    f=np.array(f)

    # Create cubic spline interpolation for f
    f_spline = CubicSpline(f[:,0], f[:,1])

    g= []
    mydir=Output_File(target,"input",["yorp_g.txt"])
    with open(mydir, 'r') as file:
        for line in file:
            x, y = map(float, line.split(","))
            g.append([x,y])
    
    g=np.array(g)

    # Create cubic spline interpolation for g
    g_spline= CubicSpline(g[:,0], g[:,1])

    return f_spline, g_spline


#%%


def Critical_omega(target):
    """
    Computes valid rotation rates (omega) for the condition -mu * f_n = |f_theta|.
    Assumes u_phi = 0 and n = 0.
    
    Parameters:
    - theta: 1D numpy array of polar angles in radians.
    - R: 1D numpy array of radial distances corresponding to theta.
    - mu: Coefficient of friction (scalar).
    - bn0: Static normal force component (scalar or matching array).
    - btheta0: Static tangential force component (scalar or matching array).
    - R_phi: Radius of curvature of the phi-curve (scalar or matching array).
    
    Returns:
    - A dictionary containing arrays for the computed curvatures and the valid 
      omega roots for both the positive and negative f_theta cases.
    """
    mydir = Output_File (target,"output",["base.txt"])

    base = np.loadtxt(mydir,delimiter=",",dtype=float)
    # 1. Numerically differentiate R with respect to theta
    # np.gradient uses second-order accurate central differences
    Res= target.res
    theta = base[2:Res-2,0]
    R = base[2:Res-2,1]*target.d/2
    R_prime = base[2:Res-2,3]*target.d/2#np.gradient(R, theta)
    mu = np.tan(math.radians(target.delta))
    b_n = target.rgrav[2:Res-2]
    b_theta = target.tgrav[2:Res-2]
    # 2. Compute the normalization factor N
    N = np.sqrt(R**2 + R_prime**2)
    rho = target.dens
    # 3. Compute base curvatures (kappa_phi_n and kappa_phi_g)
    # Using errstate to suppress warnings at theta = 0 or pi where sin(theta) = 0
    with np.errstate(divide='ignore', invalid='ignore'):
        # Corrected denominator including R
        denom = R * np.sin(theta) * N
        
        R_cos_prime = R_prime * np.cos(theta) - R * np.sin(theta)
        R_sin_prime = R_prime * np.sin(theta) + R * np.cos(theta)
        
        # Since n=0, the full curvatures equal their base components
        kappa_phi_n = -R_cos_prime / denom
        kappa_phi_g = R_sin_prime / denom
        R_phi = R*np.sin(theta)
        
        c_0 =target.cohesion_cons
        c_1 = target.cohesion_linear
    # Initialize output dictionary with NaNs for the physical roots
    results = {
        'theta': theta,
        'R': R,
        'R_prime': R_prime,
        'kappa_phi_n': kappa_phi_n,
        'kappa_phi_g': kappa_phi_g,
        'omega_case_A': np.full_like(theta, 1, dtype=float),
        'omega_case_B': np.full_like(theta, 1, dtype=float)
    }
    
    # ---------------------------------------------------------
    # Case A: Assuming f_theta >= 0 (-mu * f_n = f_theta)
    # ---------------------------------------------------------
    
    
    num_A = c_0*c_1/rho -(mu * b_n + b_theta)
    den_A = kappa_phi_g + mu * kappa_phi_n
    
    with np.errstate(divide='ignore', invalid='ignore'):
        omega_sq_A = num_A / ((R_phi**2) * den_A)
        f_theta_A = b_theta + (R_phi**2) * kappa_phi_g * omega_sq_A
        
        # Mask for valid physical solutions
        valid_mask_A = (omega_sq_A >= 0) & (f_theta_A >= -1e-9) & ~np.isnan(omega_sq_A)
        results['omega_case_A'][valid_mask_A] = np.sqrt(omega_sq_A[valid_mask_A])

    # ---------------------------------------------------------
    # Case B: Assuming f_theta < 0 (-mu * f_n = -f_theta)
    # ---------------------------------------------------------
    num_B = c_0*c_1/rho - mu * b_n + b_theta
    den_B = - kappa_phi_g + mu * kappa_phi_n
    
    with np.errstate(divide='ignore', invalid='ignore'):
        omega_sq_B = num_B / ((R_phi**2) * den_B)
        f_theta_B = b_theta + (R_phi**2) * kappa_phi_g * omega_sq_B
        
        # Mask for valid physical solutions
        valid_mask_B = (omega_sq_B >= 0) & (f_theta_B < 1e-9) & ~np.isnan(omega_sq_B)
        results['omega_case_B'][valid_mask_B] = np.sqrt(omega_sq_B[valid_mask_B])
        
    crit_omega_A = np.min(results['omega_case_A'])
    crit_omega_B = np.min(results['omega_case_B'])
    
    shed_omega = np.sqrt(np.maximum(0,-b_n/kappa_phi_n))/R_phi
    
    
    
    return min(crit_omega_A,crit_omega_B),min(shed_omega)

