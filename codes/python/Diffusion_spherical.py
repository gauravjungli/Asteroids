#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 22:04:49 2024

@author: g
"""

import numpy as np
from scipy.special import jv, jvp,legendre_p_all  # Hypergeometric function  # Bessel function of first kind and its derivative and Legendre polynomial
from scipy.optimize import root_scalar
from numba import jit
import time
import math
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from scipy.interpolate import make_interp_spline
import traceback



def find_threshold_time(E0,t0,theta,target,rgrav):
    """Find the time at which energy falls below 10% of its initial value."""
 
    Gamma = 0.25
    f = target.f

    threshold = target.dens*(Gamma*(rgrav-(target.omega[2]*np.sin(theta))**2*target.d/2)/2/np.pi/f)**2/2
    
    # Define search range
    left, right = t0, 10000*t0  # Adjust upper bound as needed
   
   # Perform binary search
    while right - left > 1e-1:
       mid = (left + right) / 2
       if E0*energy(t=mid,theta=np.cos(theta),target=target) > threshold:
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
    legendre = legendre_p_all(n - 1, theta)[0]  # Compute all needed Legendre polynomials at once

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
    return legendre_p_all(n, theta)[0]

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
    
   # print(value,t,theta,k_s,k_d)
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

    for i in range(N):
        upper =  (i+1)*100*(target.d/500)
        if i==0:
            lower = upper
            max_iters = 100
            count = 0
            while max_time(lower, theta[i], target) > 0 and count < max_iters:
                upper = lower
                lower = lower/2
                count+=1
        else:
            
            lower = t[i-1]/2 
        
        try:
            t[i]= root_scalar(max_time,bracket=[lower,upper],args=(theta[i],target),method='brentq').root
            E[i] = energy(theta=theta[i],t=t[i], target=target)
        except ValueError as e:
            t[i]=1
            print(t[i],E[i])
            print("No root found",i,target.f,target.wave_speed_P,target.d,lower)
            print(e)
            print("--- Full Traceback ---")
            traceback.print_exc()
            print("----------------------")

    return E,t

    

def compute_time(target,impactor):
    N =target.N 
    theta = target.theta
    t = target.t_max
    t_lan =  np.zeros(N)
    E = 1/2*impactor.M*(impactor.vel**2)*target.efficiency
    for i in range(N):
        t_lan[i] =  find_threshold_time( E0=E, t0=t[i],theta=theta[i],target=target,rgrav=target.grav) #needs improvement for non-spherical case

    avg_lan = np.average(t_lan)  
    print (f"Seismic shaking time is {avg_lan}")
    return avg_lan


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
   # print(f"Time taken in finding roots:{end-start}")
    return roots

#%%


