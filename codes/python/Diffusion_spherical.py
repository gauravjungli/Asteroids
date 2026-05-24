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


def compute_height(parameters,target=None,impactor=None,grid=None,dimension='1D',crater_depth=None):
    
    
    # Define constants and variables (replace these with actual values)
    
    pi = math.pi
    R = target.d/2      
    theta = target.theta      
    rho = target.dens   
    Nx = int(parameters['X Resolution']) 
    Ny = int(parameters['Y Resolution']) 

    E = 1/2*impactor.M*(impactor.vel**2)*target.efficiency*target.energy
    v_p = target.wave_speed_P  
    v_s = target.wave_speed_S 
    mu_t = v_s**2*rho
    lambda_t = (v_p**2-2*v_s**2)*rho
    
    
    if  dimension != '1D':
        
        grid2 =np.zeros((Nx,4))
        grid2[:,0] = grid[::Ny,0]
        grid2[:,1] = grid[::Ny,4]
        grid2[:,2] =grid[::Ny,3]
        grid2[:,3] = grid[::Ny,5]
       
    else:
        grid2=grid
        
    r_grav = make_interp_spline(grid2[2:-2,0], target.rgrav[2:-2])
    t_grav = make_interp_spline(grid2[2:-2,0], target.tgrav[2:-2]) 
    f_base = make_interp_spline(grid2[2:-2,0], grid2[2:-2,1])
    f_dbase = make_interp_spline(grid2[2:-2,0], grid2[2:-2,3])


    base = f_base(theta)*R
    dbase = f_dbase(theta)*R
    metric = np.sqrt(base**2 + dbase**2)
    J =  metric*base*np.sin(theta)
    kappa_phi_g = 1/J*(dbase*np.sin(theta)+base*np.cos(theta))
    kappa_phi_n = -1/J*(dbase*np.cos(theta)-base*np.sin(theta))
    rgrav = r_grav(theta)
    tgrav = t_grav(theta)
    normal_p = rho*np.abs(rgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_n)
    tangential_p = rho*np.abs(tgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_g)
    
    if target.failure_mode == "P-wave":
        gamma_max = np.sqrt(2*E/(lambda_t+2*mu_t))*(lambda_t*np.tan(target.delta*pi/180)+mu_t*np.tan(np.pi/4+target.delta*pi/360))
    else:
        gamma_max = v_s*np.sqrt(2*rho*E)/np.cos(target.delta*pi/180)
        
    alpha = (np.tan(target.delta*pi/180)*normal_p -tangential_p +target.cohesion_linear)

    height = (gamma_max-target.cohesion_cons)/alpha
    Gamma = np.abs(2*np.pi*target.f/rgrav*np.sqrt(2*E/target.dens))
    
    if height[-1]<0:
        
        print("No global sesmic shaking")
        
        return 0, 0
        
    
    avg_Gamma = np.trapezoid(Gamma*J*height,theta)/np.trapezoid(J*height,theta)
    if dimension != '1D':
        theta = np.append(0,theta)
        max_height = np.max(height)
        height =np.append(max_height,height)
        height = np.clip(height, a_min=0, a_max=crater_depth)
        
        theta_new, phi_new = transform_spherical_coords(grid[:,0], grid[:,1], impactor.Theta, impactor.Phi)
        h_interp=make_interp_spline(theta,height)
        height_new = h_interp(theta_new)
        height_new = np.clip(height_new, a_min=0, a_max=max_height)

        return height_new/R, avg_Gamma  #Contains the flag to check whether the grid lies inside the crater
    else:
      #  f_h = make_interp_spline(target.theta, height)
        avg_height = np.trapezoid(height*J,theta)/(np.trapezoid(J,theta))
        H = avg_height/R
        
        return H, avg_Gamma



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


def spherical_S2_to_S1(theta2, phi2, theta_pole, phi_pole):
    """
    Converts spherical coordinates from a rotated system (S2) to a primary system (S1).

    The S2 system's North Pole (theta2=0) is located at (theta_pole, phi_pole)
    in the S1 system.

    Args:
        theta2 (float or np.ndarray): Polar angle(s) in S2 (radians, 0 to pi).
        phi2 (float or np.ndarray): Azimuthal angle(s) in S2 (radians, 0 to 2*pi).
        theta_pole (float): Polar angle of the S2 pole in S1 (radians, 0 to pi).
        phi_pole (float): Azimuthal angle of the S2 pole in S1 (radians, 0 to 2*pi).

    Returns:
        tuple: (theta1, phi1)
            theta1 (float or np.ndarray): Polar angle(s) in S1 (radians, 0 to pi).
            phi1 (float or np.ndarray): Azimuthal angle(s) in S1 (radians, 0 to 2*pi).
    """
    theta2 = np.asarray(theta2)
    phi2 = np.asarray(phi2)

    # --- Handle Edge Cases ---
    # S2 pole is at S1 North Pole (systems aligned)
    
    if np.isclose(theta_pole, 0.0):
        # print("Edge Case: theta_pole is near 0. Systems aligned.")
        # We might need to consider a potential rotation around Z if phi_pole matters
        # Assuming standard alignment where S2's prime meridian aligns with S1's
        # along the great circle connecting poles (which is undefined here).
        # Safest assumption is identity if theta_pole is 0.
        phi1 = (phi2 + phi_pole) % (2 * np.pi) # Align prime meridians based on phi_pole? Or just phi1=phi2?
                                              # Let's assume the rotation M below handles this correctly even near 0.
                                              # Reverting to calculation below, as M should become Rz(phi_pole)
        # return theta2, (phi2 + phi_pole) % (2 * np.pi) # Simple rotation around Z
        pass # Let the main calculation handle it, check matrix M below

    # S2 pole is at S1 South Pole
    elif np.isclose(theta_pole, np.pi):
        # print("Edge Case: theta_pole is near pi. Systems anti-aligned.")
        theta1 = np.pi - theta2
        phi1 = (phi_pole + np.pi + phi2) % (2 * np.pi) # Need phi_pole for orientation, +pi flip, +phi2 offset
                                                     # Let's re-derive: cartesian x1=-x2, y1=-y2, z1=-z2 requires rotation matrix diag(-1,-1,-1) IF phi_pole=0
                                                     # If phi_pole != 0, the S2 x/y axes are rotated first.
                                                     # x1 = -sin(theta2)cos(phi2+phi_pole?) No this seems wrong.
                                                     # Let's recalculate phi1 = atan2(y1,x1)
                                                     # v1 = Rz(phi_pole) @ Ry(pi) @ v2 ?? No.
                                                     # Let's use the formula derived: theta1=pi-theta2, phi1 = (phi2+pi) % (2*pi).
                                                     # Where does phi_pole come in? The definition of phi2=0 depends on it.
                                                     # Let's trust the matrix M calculation below for consistency.
                                                     # If theta_pole=pi, matrix M should handle it.
        # return np.pi - theta2, (phi2 + np.pi) % (2 * np.pi) # Ignoring phi_pole influence here, might be wrong?
        pass # Let the main calculation handle it


    # --- General Calculation ---

    # 1. Convert S2 spherical to S2 Cartesian (assuming r=1)
    x2 = np.sin(theta2) * np.cos(phi2)
    y2 = np.sin(theta2) * np.sin(phi2)
    z2 = np.cos(theta2)
    v2 = np.stack([x2, y2, z2], axis=0)
    # Ensure v2 has shape (3,) or (3, N)
    if v2.ndim == 1:
       v2 = v2[:, np.newaxis] # Make it (3, 1) for matmul

    # 2. Construct the transformation matrix M (S2 basis in S1 coords)
    ct = np.cos(theta_pole)
    st = np.sin(theta_pole)
    cp = np.cos(phi_pole)
    sp = np.sin(phi_pole)

    # S2 basis vectors (i2, j2, k2) expressed in S1 coordinates
    # Note: Need to handle st=0 case if we used the j2 = normalize(k1 x k2) definition directly.
    # Using the derived components avoids division by zero implicitly.
    i2 = np.array([ct * cp, ct * sp, -st])
    j2 = np.array([-sp, cp, 0])
    k2 = np.array([st * cp, st * sp, ct])

    # Matrix M has i2, j2, k2 as columns
    M = np.stack([i2, j2, k2], axis=1) # Shape (3, 3)

    # 3. Apply Rotation: v1 = M @ v2
    # M is (3, 3), v2 is (3, N). Result v1 is (3, N)
    v1 = M @ v2

    # 4. Convert S1 Cartesian to S1 Spherical
    x1 = v1[0, :]
    y1 = v1[1, :]
    z1 = v1[2, :]

    # Calculate r1 for robustness, although it should be 1 (or original r)
    r1 = np.sqrt(x1**2 + y1**2 + z1**2)
    # Avoid division by zero or NaNs if r1 is very small (input was origin)
    # Also clip z1/r1 for numerical stability with arccos
    safe_r1 = np.where(np.isclose(r1, 0.0), 1.0, r1) # Use 1 if r=0
    z1_over_r1 = np.clip(z1 / safe_r1, -1.0, 1.0)
    theta1 = np.arccos(z1_over_r1)

    # Phi1: atan2 handles quadrants correctly. Result is in (-pi, pi]
    phi1_raw = np.arctan2(y1, x1)
    # Adjust to [0, 2*pi)
    phi1 = phi1_raw % (2 * np.pi)

    # If input was scalar, return scalar
    if theta2.ndim == 0:
        return theta1.item(), phi1.item()
    else:
        return theta1, phi1
    
    
def spherical_to_cartesian(theta, phi, r=1):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.stack([x, y, z], axis=-1)

def cartesian_to_spherical(cartesian_coords):
    x, y, z = cartesian_coords[..., 0], cartesian_coords[..., 1], cartesian_coords[..., 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.clip(z / r, -1.0, 1.0))  # Ensuring theta is in [0, pi]
    phi = np.arctan2(y, x)
    return np.stack([theta, phi], axis=-1)

def rotation_matrix(theta_0, phi_0):
    Rz = np.array([[np.cos(phi_0), -np.sin(phi_0), 0],
                   [np.sin(phi_0), np.cos(phi_0), 0],
                   [0, 0, 1]])
    
    Ry = np.array([[np.cos(theta_0), 0, np.sin(theta_0)],
                   [0, 1, 0],
                   [-np.sin(theta_0), 0, np.cos(theta_0)]])
    
    return Ry @ Rz  # Apply rotation about z-axis first, then y-axis

def rotate_spherical_coordinates(theta, phi, theta_0, phi_0, r=1):
    cartesian_coords = spherical_to_cartesian(theta, phi, r)
    R = rotation_matrix(theta_0, phi_0)
    rotated_cartesian = np.einsum('ij,...j->...i', R, cartesian_coords)
    return cartesian_to_spherical(rotated_cartesian)


def transform_spherical_coords(theta, phi, theta_0, phi_0):
    """
    Transforms spherical coordinates (theta, phi) to a new system (theta', phi')
    where the direction (theta_0, phi_0) becomes the new North Pole (Z' axis).

    Assumes standard physics spherical coordinates:
    - theta: polar angle from Z+ (radians, 0 to pi)
    - phi: azimuthal angle from X+ (radians, 0 to 2*pi)

    Args:
        theta (float or np.ndarray): Original polar angle(s) in radians.
        phi (float or np.ndarray): Original azimuthal angle(s) in radians.
        theta_0 (float): Polar angle of the new North Pole in the old system (radians).
        phi_0 (float): Azimuthal angle of the new North Pole in the old system (radians).

    Returns:
        tuple[float or np.ndarray, float or np.ndarray]:
            - theta_prime: Polar angle(s) in the new coordinate system (radians).
            - phi_prime: Azimuthal angle(s) in the new coordinate system (radians, 0 to 2*pi).
    """
    # Ensure inputs are numpy arrays for broadcasting
    theta = np.asarray(theta)
    phi = np.asarray(phi)

    # === Step 1: Define the New Z' axis ===
    sin_t0 = np.sin(theta_0)
    cos_t0 = np.cos(theta_0)
    sin_p0 = np.sin(phi_0)
    cos_p0 = np.cos(phi_0)

    Z_prime_vec = np.array([sin_t0 * cos_p0, sin_t0 * sin_p0, cos_t0])

    # === Step 2: Define the New X' and Y' axes ===
    # Choose X' to point "downhill" from Z' towards the original Z axis, projected.
    # Or equivalently, Z = cos(theta_0)*Z' + sin(theta_0)*X' (if phi_0=0)
    # So X' is proportional to Z - cos(theta_0)*Z'
    # Unit vector for original Z
    Z_vec = np.array([0.0, 0.0, 1.0])
    # Component of Z along Z'
    Z_proj_Zprime = cos_t0 * Z_prime_vec
    # Vector pointing from Z' towards Z (perpendicular to Z')
    X_prime_dir_temp = Z_vec - Z_proj_Zprime

    # Handle edge case where new pole is old pole (theta_0 = 0 or pi)
    # Use small epsilon to avoid division by zero if theta_0 is exactly 0 or pi
    epsilon = 1e-12
    norm_X_prime_temp = np.linalg.norm(X_prime_dir_temp)

    if norm_X_prime_temp < epsilon:
         # If Z' aligns with Z, the standard X,Y axes (rotated by phi_0) can be used.
         # X' should correspond to theta'=pi/2, phi'=0 in new system.
         # A simple choice is rotation of original X by phi_0: [cos(phi_0), sin(phi_0), 0]
         # Let's refine this using Y' first for this case.
         # Y' = Z' x X'. If Z' = Z, Y' should be related to original Y rotated by phi_0.
         # Let's use the standard cross product derived Y' vector below, it should simplify.
         Y_prime_vec = np.array([sin_p0, -cos_p0, 0.0])
         # Then X' = Y' x Z' (note order for right-handed system: X=YxZ)
         X_prime_vec = np.cross(Y_prime_vec, Z_prime_vec)
    else:
        # Normalize the derived X' direction
        X_prime_vec = X_prime_dir_temp / norm_X_prime_temp
        # Calculate Y' = Z' x X' to complete the right-handed system
        Y_prime_vec = np.cross(Z_prime_vec, X_prime_vec)


    # === Step 3: Convert Original Point(s) to Cartesian ===
    # Assume unit sphere (rho=1) as only direction matters for angles
    sin_t = np.sin(theta)
    cos_t = np.cos(theta)
    sin_p = np.sin(phi)
    cos_p = np.cos(phi)

    # Need to handle theta possibly being an array
    if theta.shape: # If theta is an array
        P_vec_x = sin_t * cos_p
        P_vec_y = sin_t * sin_p
        P_vec_z = cos_t
        # Stack into (N, 3) array if theta/phi are 1D, or handle higher dims if needed
        # For dot product using np.dot, easier if P_vec is (3, N) or similar
        P_vec = np.stack((P_vec_x, P_vec_y, P_vec_z), axis=0) # Shape (3, ...)
    else: # If theta is a scalar
         P_vec = np.array([sin_t * cos_p, sin_t * sin_p, cos_t]) # Shape (3,)


    # === Step 4: Project onto New Axes ===
    # Use np.dot which handles matrix/vector multiplication correctly
    # If P_vec is (3,) -> result is scalar
    # If P_vec is (3, N) -> result is (N,)
    x_prime = np.dot(X_prime_vec, P_vec)
    y_prime = np.dot(Y_prime_vec, P_vec)
    z_prime = np.dot(Z_prime_vec, P_vec)

    # === Step 5: Convert New Cartesian to New Spherical ===
    # Calculate new polar angle theta'
    # Clip argument to arccos to avoid domain errors due to floating point inaccuracies
    rho_prime = np.sqrt(x_prime**2 + y_prime**2 + z_prime**2) # Should be ~1.0
    cos_theta_prime_arg = np.clip(z_prime / rho_prime, -1.0, 1.0)
    theta_prime = np.arccos(cos_theta_prime_arg)

    # Calculate new azimuthal angle phi'
    phi_prime = np.arctan2(y_prime, x_prime)

    # Adjust phi_prime range from [-pi, pi] to [0, 2*pi)
    phi_prime = np.mod(phi_prime, 2 * np.pi)

    # Handle poles explicitly: phi' is degenerate at the poles (theta'=0 or pi)
    # arctan2(0, 0) returns 0, which is a fine convention.
    # Check if any input points were exactly poles where phi is degenerate
    # if theta.shape:
    #     near_pole = (theta < epsilon) | (theta > np.pi - epsilon)
    #     # Or check if output points are near new poles
    #     near_new_pole = (theta_prime < epsilon) | (theta_prime > np.pi - epsilon)
    #     # Can set phi_prime to 0 for these cases if needed, but arctan2(0,0) handles it.
    # else:
    #     if (theta < epsilon) or (theta > np.pi - epsilon):
    #         pass # phi_prime already 0 from arctan2(0,0) if input is pole

    return theta_prime, phi_prime
