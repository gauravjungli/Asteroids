#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 22:04:49 2024

@author: g
"""

import numpy as np
from scipy.special import jv, jvp,lpn, hyp2f1  # Hypergeometric function  # Bessel function of first kind and its derivative and Legendre polynomial
from scipy.optimize import root_scalar, brentq
from numba import jit
import time
import math
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from scipy.interpolate import make_interp_spline, CubicSpline

def find_threshold_time(E0, t0,theta,target):
    """Find the time at which energy falls below 10% of its initial value."""

    threshold = 0.1 * E0  # 10% of initial energy
    
    # Define search range
    left, right = t0, 100*t0  # Adjust upper bound as needed
   
   # Perform binary search
    while right - left > 1e-6:
       mid = (left + right) / 2
       if energy(t=mid,theta=theta,target=target) > threshold:
           left = mid  # Move right if energy is still above threshold
       else:
           right = mid  # Move left if energy is below threshold
   
    return (left + right) / 2 -t0  # Return the mid-point as the estimated time
    


""" 
The function gives the peak seismic energy due to impact at any location theta. 
Time t corresoponds to the peak  time."""

def  energy(t,theta, target,E=1.0):
    n, m = target.roots.shape

    k_d=2*np.pi*target.f/target.Q
    G = 3/2 
    legendre = lpn(n - 1, theta)[0]  # Compute all needed Legendre polynomials at once

    # Vectorized computation for B
    i_vals = np.arange(n).reshape(n, 1)  # Shape (n, 1) to match (n, m)
    B = (2 * i_vals + 1) / (1 - i_vals * (i_vals + 1) / (target.d / 2 * target.roots) ** 2)
    
    # Compute exponential term only once per root
    exp_term = np.exp(-target.roots**2 * target.k_s * t)

    # Vectorized sum over i and j
    G += np.sum(B * legendre[:, None] * exp_term)  # Using broadcasting

    # Final scaling
    G *= E * np.exp(-k_d * t) / (2 * np.pi * (target.d / 2) ** 3)

    return G
       

"""Next three functions are used for computing the peak time. The function is broken into three subfunctions to speed up
 the execution using Numba. """

# This function will not be compiled by Numba since it uses lpn from scipy
def compute_legendre(n, theta):
    return lpn(n, theta)[0]

# JIT-compiled function for performance
@jit(nopython = True)
def max_time_optimized(t, theta, legendre, roots, k_d, k_s, d):
    n = roots.shape[0]
    m = roots.shape[1]
    value = 3 / 2 * k_d
    k_s_t = k_s * t
    
    for i in range(n):
        legendre_i = legendre[i]
        for j in range(m):
            root_ij = roots[i, j]
            root_ij_squared = root_ij**2
            B = (2*i+1) / (1 - i * (i + 1) / (d / 2 * root_ij)**2)
            exp_term = np.exp(-root_ij_squared * k_s_t)
            value += B * legendre_i * exp_term * (k_d + root_ij_squared * k_s)
    
    return value

# This function computes all necessary values and passes them to the optimized function
def max_time(t, theta, target):
    # Precompute legendre outside of Numba-compiled function
    legendre = compute_legendre(target.roots.shape[0], theta)
    
    # Extract necessary values from target
    roots = target.roots
    k_d = np.pi * 2 * target.f / target.Q
    k_s = target.k_s
    d = target.d
    
    # Call the optimized function
    return max_time_optimized(t, theta, legendre, roots, k_d, k_s, d)



def compute_energy(target):
    N =target.N 
    theta = np.cos(target.theta) 
    t = np.zeros(N)
    E = np.zeros(N)
    t_lan =  np.zeros(N)
    for i in range(N):
        lower = t[i-1] if i>0 else 1e-1
        t[i]= root_scalar(max_time,bracket=[lower,100],args=(theta[i],target),method='brentq').root
        E[i] = energy(theta=theta[i],t=t[i], target=target)
        t_lan[i] =  find_threshold_time(E0=E[i], t0=t[i],theta=theta[i],target=target)
        print ("Time               Energy             cos(theta)")
        print (t[i],E[i],theta[i],t_lan[i])
    avg_lan = np.average(t_lan)   
    return E,avg_lan
    


def compute_height(target=None,impactor=None):
    
    
    # Define constants and variables (replace these with actual values)
  
    pi = math.pi
    R = target.d/2      
    theta = target.theta      
    rho = target.dens            
    #rgrav=make_interp_spline(x, y)             
    E = 1/2*impactor.M*(impactor.vel**2)*target.efficiency*target.energy
    v_p = target.wave_speed        
    
    peak_p = v_p*np.sqrt(2*rho*E)
    normal_p = rho*np.abs(target.grav + target.omega[2]**2*R*np.sin(theta)**2)
    tangential_p = rho*np.abs(target.omega[2]**2*R*np.sin(theta)*np.cos(theta))
    height = (peak_p*np.tan(target.delta*pi/180)-target.cohesion_cons)/(np.tan(target.delta*pi/180)*normal_p
                   -tangential_p +target.cohesion_linear)
    avg_height = np.trapz(height*np.sin(theta),theta)/2
    
    return avg_height/R

 
""" Currently not in use. Need to change implementation. May be wrong result. """

def frequency(target):

    # Define constants and variables (replace these with actual values)
    G = 6.67430e-11    # gravitational constant
    pi = math.pi
    R = target.d/2            # example radius value
    rho = target.dens          # example density value
    phi = 0.6          # example phi value
    Z = 6.0            # example Z value
    lambda_ = 1.0e+9      # example lambda value
    mu_0 = 1.0e+9         # example mu_0 value
    Co = 2.0           # example Co value
    
    # Compute terms in the numerator
    numerator_part1 = math.sqrt(5) * math.sqrt(3 * lambda_ + 4 * mu_0) * math.sqrt(rho)
    numerator_part2 = 2**(5/6) * pi**(1/3) * (lambda_ + 2 * mu_0)**(1/3) * R
    
    # Argument for the hypergeometric function
    hypergeom_arg = (2 * G * pi * rho**2 * R**2) / (2 * G * pi * R**2 * rho**2 + 3 * Co)
    hypergeom_val = hyp2f1(1/6, 1/2, 3/2, hypergeom_arg)
    
    # Complete the numerator
    numerator = numerator_part1 * numerator_part2 * hypergeom_val
    
    # Compute terms in the denominator
    denominator_part1 = 2 * phi**(1/3) * Z**(1/3) * mu_0**(1/3) * (lambda_ + mu_0)**(1/3)
    denominator_part2 = math.sqrt(13 * lambda_ + 20 * mu_0) * (2 * G * pi * R**2 * rho**2 + 3 * Co)**(1/6)
    
    # Complete the denominator
    denominator = denominator_part1 * denominator_part2
    
    # Final result
    result = numerator / denominator
    
    return 1/result
    
    print("Result:", result)


def bessel_eq(x, n):
    return 2 * x  * jvp(n + 0.5, x ) - jv(n + 0.5, x )

""" Searching bounds for the root of Bessel's equation to be used in the root_scalar function"""

def dynamic_bounds_search(n, lower, step=1.0, max_iter=100):

    upper = lower + step
    iter_count = 0
    # Increment upper bound until sign changes
    while bessel_eq(lower, n) * bessel_eq(upper, n) > 0:
        lower = upper
        upper += step
        iter_count += 1
        if iter_count > max_iter:
            raise RuntimeError("Exceeded maximum iterations while searching for bounds.")
    return lower, upper


""" Serial version of finding roots. New parallel version is implemented. Not in use currently"""

def find_roots(n, num_roots=5, step=1.0):
    roots = []
    lower = max(0.1, np.sqrt(n * (n + 1)))  # Start from a small value
    for _ in range(num_roots):
        # Dynamically find the next upper bound where the sign changes
        lower, upper = dynamic_bounds_search(n, lower, step)
        # Use the brentq method to find the root in the interval
        sol = root_scalar(bessel_eq, args=(n,), bracket=[lower, upper], method='brentq')
        roots.append(sol.root)
        # Set new lower bound for the next root search
        lower = upper
    return np.array(roots)



def find_intervals(n, num_roots, step):
    """Find unique intervals where each root is located."""
    intervals = []
    lower = max(0.1, np.sqrt(n * (n + 1)))  # Start from a small value
    for _ in range(num_roots):
        lower, upper = dynamic_bounds_search(n, lower, step)
        intervals.append((lower, upper))
        lower = upper  # Update lower for the next search
    return intervals


def root_worker(n, interval):
    """Find the root in a given interval."""
    lower, upper = interval
    sol = root_scalar(bessel_eq, args=(n,), bracket=[lower, upper], method='brentq')
    return sol.root


def find_roots_parallel(n, num_roots=5, step=1.0):
    # Step 1: Generate unique intervals for each root
    intervals = find_intervals(n, num_roots, step)
    
    # Step 2: Use parallel processing to find roots in these intervals
    roots = []
    with ThreadPoolExecutor(max_workers=24) as executor:
        futures = [executor.submit(root_worker, n, interval) for interval in intervals]
        for future in futures:
            roots.append(future.result())
    
    return np.array(roots)


def compute_roots(i, m, d):
    roots = find_roots_parallel(i, num_roots=m) / (d / 2)
    return i, roots

""" Parallel version for computing roots of the Bessel equation. """
def parallel_root_computation(n, m, d):
    roots = np.zeros((n + 1, m))  # Initialize the roots array
    start =time.time()
    with ProcessPoolExecutor(max_workers=24) as executor:
        # Submit all tasks in parallel
        futures = [executor.submit(compute_roots, i, m, d) for i in range(0, n + 1)]
        
        # Collect results and place them in the correct row of `roots`
        for future in futures:
            i, result = future.result()
            roots[i, :] = result
    end =time.time()
    print(f"Time taken in finding roots:{end-start}")
    return roots