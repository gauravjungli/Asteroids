#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:09:28 2026

@author: g
"""

import math
from scipy.interpolate import make_interp_spline
import numpy as np
from collisions import G
from Diffusion_spherical import  compute_energy
from IO import  Output_File
from Fit import Fit
from script_2D.Crater import Crater
from transformations import transform_spherical_coords
from Scheeres import Scheeres
from scipy.optimize import root_scalar
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import os

def add_gaussian_fillet(y, A=None, sigma=10, index=None,flag=1):
    """
    Adds a Gaussian bump near the steepest drop to ease slope discontinuity.
    A: Amplitude of the Gaussian patch (defaults to 20% of peak height)
    sigma: Width of the Gaussian patch in array indices
    offset: Shifts the Gaussian center relative to the max slope point
    """
    steep_idx = index
    
    if A is None:
        
        A = 0.1 * np.max(y[index-1:index+int(flag*sigma):flag])  # Default amplitude heuristic
        
    x = np.arange(len(y))
    
    gaussian_patch = A * np.exp(-0.5 * ((x - steep_idx) / sigma) ** 2) 
    
    y =  np.where(gaussian_patch>y,gaussian_patch,y)
    return y 

def stress_pressure(target,normal_p,tangential_p):
    
    phi = target.delta*np.pi/180
    
    avg_stress = normal_p 
    
    diff_stress = 0
    
    theta_prime = -phi/2
    
    sigma_x_prime = avg_stress + diff_stress * np.cos(2 *theta_prime ) + np.abs(tangential_p) * np.sin(2 * theta_prime)
    
    tau_xy_prime = -diff_stress * np.sin(2 * theta_prime) + tangential_p * np.cos(2 * theta_prime)
    
    return sigma_x_prime, tau_xy_prime

def stress_principal(target,normal_p,tangential_p):
    
     sin_psi = tangential_p/np.sqrt(normal_p**2 + tangential_p**2) 
    
     psi = np.arcsin(sin_psi)
    
     phi = target.delta*np.pi/180
    # # Step 1: Calculate sigma_yy
    
     sigma_yy = np.ones(target.res)*normal_p
    
     sigma_yy =  normal_p - (2 * tangential_p) / np.tan(2 * psi)
    
    # # Step 2: Calculate stresses at angle theta
     avg_stress = (normal_p + sigma_yy) / 2
     diff_stress = (normal_p - sigma_yy) / 2
    
     theta_prime = math.pi/4 + psi + phi/2
    
     sigma_x_prime = avg_stress + diff_stress * np.cos(2 *theta_prime ) + tangential_p * np.sin(2 * theta_prime)
    
     tau_xy_prime = -diff_stress * np.sin(2 * theta_prime) + tangential_p * np.cos(2 * theta_prime)
     
     return sigma_x_prime, tau_xy_prime
    


def compute_height(target=None,impactor=None,grid=None,dimension='1D',crater_depth=None):
    
    
    # Define constants and variables
    pi = math.pi
    R = target.d/2      
    theta = target.theta      
    rho = target.dens   
    
    res = target.res
    Nx = target.x_res
    Ny = target.y_res
    
    if  dimension != '1D':
        
        grid2 = np.zeros((Nx,4))
        grid2[:,0] = grid[::Ny,0]
        grid2[:,1] = grid[::Ny,4]
        grid2[:,2] =grid[::Ny,3]
        grid2[:,3] = grid[::Ny,5]
       
    else:
        
        grid2=grid
    
    if target.failure == "Scheeres":
        
        height = Scheeres(parameters,target,grid)/R
        
        epsilon = np.average(height)
        
        return  height/epsilon, epsilon, 0
        
    if impactor:
        
        E = 1/2*impactor.M*(impactor.vel**2)*target.efficiency*target.energy
        v_p = target.wave_speed_P  
        v_s = target.wave_speed_S 
        mu_t = v_s**2*rho
        lambda_t = (v_p**2-2*v_s**2)*rho
            
        r_grav = make_interp_spline(grid2[2:-2,0], target.rgrav[2:-2])
        t_grav = make_interp_spline(grid2[2:-2,0], target.tgrav[2:-2]) 
        f_base = make_interp_spline(grid2[2:-2,0], grid2[2:-2,1])
        f_dbase = make_interp_spline(grid2[2:-2,0], grid2[2:-2,3])
        base = f_base(theta)*R
        dbase = f_dbase(theta)*R
        rgrav = r_grav(theta)
        tgrav = t_grav(theta)

        metric = np.sqrt(base**2 + dbase**2)
        J =  metric*base*np.sin(theta)
        kappa_phi_g = 1/J*(dbase*np.sin(theta)+base*np.cos(theta))
        kappa_phi_n = -1/J*(dbase*np.cos(theta)-base*np.sin(theta))
    
        normal_p = rho*np.abs(rgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_n)
        tangential_p = rho*np.abs(tgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_g)

        gamma_max = 0
        
        if target.failure_mode == "P-wave":
            gamma_max = np.sqrt(2*E/(lambda_t+2*mu_t))*(lambda_t*np.tan(target.delta*pi/180)+mu_t*np.tan(np.pi/4+target.delta*pi/360))
        else:
            gamma_max = v_s*np.sqrt(2*rho*E)/np.cos(target.delta*pi/180)
        
        alpha = (np.tan(target.delta*pi/180)*normal_p -tangential_p +target.cohesion_linear)
        
        height = (gamma_max-target.cohesion_cons)/alpha
        Gamma = np.abs(2*np.pi*target.f/rgrav*np.sqrt(2*E/target.dens))
        
    else:
        
        base = grid[:,1]*R
        theta = grid[:,0]
        dbase = grid[:,3]*R
        rgrav = target.rgrav
        tgrav = target.tgrav

        metric = np.sqrt(base**2 + dbase**2)
        J =  metric*base*np.sin(theta)
        kappa_phi_g = 1/J*(dbase*np.sin(theta)+base*np.cos(theta))
        kappa_phi_n = -1/J*(dbase*np.cos(theta)-base*np.sin(theta))
    
        normal_p = -rho*(rgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_n)
        tangential_p = rho*(tgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_g)
        
        # sin_psi = tangential_p/np.sqrt(normal_p**2 + tangential_p**2) 
        
        # psi = np.arcsin(sin_psi)
        
        phi = target.delta*pi/180
        # # Step 1: Calculate sigma_yy
        
        # sigma_yy = np.ones(res)*normal_p
        
        # sigma_yy =  normal_p - (2 * tangential_p) / np.tan(2 * psi)
        
        # # Step 2: Calculate stresses at angle theta
        # avg_stress = (normal_p + sigma_yy) / 2
        # diff_stress = (normal_p - sigma_yy) / 2
        
        # theta_prime = math.pi/4 + psi + phi/2
        
        # sigma_x_prime = avg_stress + diff_stress * np.cos(2 *theta_prime ) + tangential_p * np.sin(2 * theta_prime)
        
        # tau_xy_prime = -diff_stress * np.sin(2 * theta_prime) + tangential_p * np.cos(2 * theta_prime)
        
        #normal_p = sigma_x_prime
        
        #tangential_p = tau_xy_prime
        
        #normal_p , tangential_p = stress_pressure(target,normal_p,tangential_p)
        
        c_0 = target.cohesion_cons
        c_1 = target.cohesion_linear

        
        mu = np.tan(phi)
        
        sign =np.sign(tangential_p)
       
        k = sign*tangential_p - mu*normal_p
        
        
        def mohr_coulomb(h,k):
            if k<0:
                return 0
            return k*h - c_0*(np.exp(c_1*h)-1)
        
        def mass_shedding(h,normal_p):

            return normal_p*h + c_0*(np.exp(c_1*h)-1)
        
        height = np.ones(res)*1e-4
        
        for i in range(len(k)):
            
            if -normal_p[i]>c_0*c_1:
            
                print("Mass shedding due to rotation")
                try:
                    sol = root_scalar(mass_shedding,bracket=[1e-6,100],args=(normal_p[i]))
                except:
                    print(f"Could not find shed mass for i= {i}")
                height[i] = sol.root
                
            
            if k[i]<c_0*c_1:
                
                continue
            try:
                sol = root_scalar(mohr_coulomb,bracket=[1e-4,50],args=(k[i]))
                height[i] = max( sol.root, height[i])
            except: 
                print(f"Could not find height for i = {i}")
        
        fig, ax = plt.subplots()
        ax.plot(theta,height)   
        
        height = gaussian_filter1d(height, sigma=int(res/100))
        # flag = True
        # for i in range(len(k)):
             
        #     if k[i]> c_0*c_1 and flag:
        #         height = add_gaussian_fillet(height,sigma=30,index = i,flag=1)#
        #         flag =False
        #     if k[i] < c_0*c_1 and not flag:
        #         height = add_gaussian_fillet(height,sigma=30,index = i,flag=-1)#gaussian_filter1d(height, sigma=int(res/100))
        #         flag = True
                
        ax.plot(theta,height)
        
        save_folder =  Output_File(target,"output",["height_failed"])
        os.makedirs(save_folder, exist_ok=True) 
        fig.savefig(save_folder + f"/failure_height_{target.slides}.png")
        plt.close(fig)
        J = np.ones_like(theta)
        
        Gamma = 0
        
    if impactor and height[-1] < 0:
            
        print("No global sesmic shaking")
            
        return height, 0, 0
            
    avg_Gamma = np.trapezoid(Gamma*J*height,theta)/np.trapezoid(J*height,theta)
    
    if dimension != '1D':
        theta = np.append(0,theta)
        max_height = np.max(height)
        height =np.append(max_height,height)
        height = np.clip(height, a_min=0, a_max=crater_depth)
        
        theta_new, phi_new = transform_spherical_coords(grid[:,0], grid[:,1], impactor.Theta, impactor.Phi)
        h_interp = make_interp_spline(theta,height)
        height_new = h_interp(theta_new)
        height_new = np.clip(height_new, a_min=0, a_max=max_height)

        return height_new/R, avg_Gamma  #Contains the flag to check whether the grid lies inside the crater
   
    else:
      #  f_h = make_interp_spline(target.theta, height)
        avg_height = np.trapezoid(height*J,theta)/(np.trapezoid(J,theta))
        epsilon = avg_height/R
        height =height/R/epsilon
        if impactor:
            height = np.ones(res)
        
        return  height, epsilon, avg_Gamma
        
#%%       ############################################################################################################################################################

def Height(target,impactor=None):
    
    """
        This calculates the failure height of the landslide. First check whether landslide model is included in the 
        simulation or not. If yes then check whether it is initiated by an impact or it is a rotational failure. 
        For impacts, update epsilon for energy dependent simulations and select a failure profile from Gaussian and
        uniform. Update these in the base array and save it in a text file that will be later used by the C++ module.
    """ 

    if not target.landslide:
        return
    
    mydir=Output_File (target,"output",["base.txt"])

    epsilon = target.epsilon
    base=np.loadtxt(mydir,delimiter=",",dtype=float)
    
    min_epsilon = target.min_epsilon
    max_epsilon = target.max_epsilon
    
    if impactor:
        K_0 = target.K_0
        mu_0 = target.mu_0
        dens =target.dens
        radius = target.d/2
        omega = target.omega[2]
        beta = target.beta
        target.wave_speed_P = np.sqrt(K_0+4/3*mu_0)*((radius**2/15/dens)*(4*np.pi*G*dens-2*omega**2))**(1/4)
        target.wave_speed_S = np.sqrt(mu_0)*((radius**2/15/dens)*(4*np.pi*G*dens-2*omega**2))**(1/4)  
        target.f = ((2*impactor.dens*target.efficiency)/(np.pi*beta**2*dens))**(1/3)*2*target.wave_speed_P/impactor.d
        target.k_d = 2*np.pi*target.f/target.Q
        target.energy, target.t_max = compute_energy(target)
 
        if target.failure_mode == "P-wave":
            target.wave_speed = target.wave_speed_P
        else:
            target.wave_speed = target.wave_speed_S 
    
    if target.dim == '1D':
        
        
        height, epsilon_c, target.Gamma = compute_height(target=target,impactor=impactor,grid = base,dimension='1D') 
        
        print("The average failure height and maximum acceleration are", epsilon_c, target.Gamma )
        
        if epsilon_c < min_epsilon:
            target.epsilon = epsilon_c
            return
        
        if impactor:
                
            if (target.Gamma*np.exp(-target.k_d/2*1/(G*4/3*3.14*target.dens)**(0.5)*0.3)<1):
                return
            
            print( f'the destablization height is {epsilon_c} and the impactor diameter is {impactor.d}\
                  and impactor velocity is {impactor.vel}' )
        
        elif target.fast_rotation_flag:
       
            print("This is the case of rotational failure")
    
        epsilon = np.clip(epsilon_c,min_epsilon,max_epsilon)    
        target.epsilon = epsilon 
    
        base[:,2] = -height[:]
        base = Fit(target,base_old=base)
        base[:,2] = height[:] 
    
    else:
        
        if impactor:
            crater_data, crater_depth = Crater(target,impactor) #Fourth column contains the flag to check whether the point lies inside the crater
        
            H,target.Gamma = compute_height(target=target,impactor=impactor,grid=base,dimension='2D',crater_depth = crater_depth)
            H = H#*(1-crater_data[:,3])
        
            if (target.Gamma*np.exp(-target.k_d/2*1/(G*4/3*3.14*1250)**(0.5)*0.3)<0.5):
                return
       # epsilon = (np.min(H) + np.max(H))/2
       # epsilon = np.clip(epsilon,min_epsilon,max_epsilon)

            base[:,2] = (crater_data[:,2]- H[:])/epsilon 
            base[:,3] = H[:]/epsilon  
              
      #  target.epsilon = epsilon   
    
        else:
            target.Gamma = 0
            
    np.savetxt(mydir,base,delimiter=",")
