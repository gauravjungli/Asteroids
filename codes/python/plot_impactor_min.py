#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 28 19:59:23 2025

@author: g
"""




from Target import Target
import pdb
import numpy as np
import math
import matplotlib.pyplot as plt
from Diffusion_spherical import   energy, compute_energy,max_time
from scipy.constants import gravitational_constant
from debug import debug_main, debug_plot
from scipy.optimize import root_scalar
G = gravitational_constant
#%%



def compute_energy_2(target):
    
    
    theta = np.cos(np.pi) 
    upper =  1000*(target.d/500)**2
    lower = upper
    max_iters = 100
    count = 0
    
    while max_time(lower, theta, target) > 0 and count < max_iters:
        upper = lower
        lower = lower/2
        count+=1
    
    try:
        t= root_scalar(max_time,bracket=[lower,upper],args=(theta,target),method='brentq').root
        E= energy(theta=theta,t=t, target=target)
    except ValueError as e:

        print("No root found",target.f,target.wave_speed_P,target.d,lower)
        print(e)

    return E,t

#%%
def compute_d(imp_d,target):
    
    K_0 = target.K_0
    mu_0 = target.mu_0
    dens =target.dens
    radius = target.d/2
    omega = target.omega[2]
    beta = target.beta

    imp_R = imp_d/2
    imp_dens =1500
    imp_M = 4/3*np.pi* imp_R**3*imp_dens
    imp_v = 5500
    
    target.wave_speed_P = np.sqrt(K_0+4/3*mu_0)*((radius**2/15/dens)*(4*np.pi*G*dens-2*omega**2))**(1/4)
    target.wave_speed_S = np.sqrt(mu_0)*((radius**2/15/dens)*(4*np.pi*G*dens-2*omega**2))**(1/4)  
    target.f = ((2*imp_dens*target.efficiency)/(np.pi*beta**2*dens))**(1/3)*2*target.wave_speed_P/(imp_d)
    target.k_d = 2*np.pi*target.f/target.Q
    target.energy, target.t_max = compute_energy_2(target)
    
    v_p = target.wave_speed_P  
    v_s = target.wave_speed_S 
    mu_t = v_s**2*dens
    lambda_t = (v_p**2-2*v_s**2)*dens
    
    if target.failure_mode == "P-wave":
        energy_min = (lambda_t+2*mu_t)*target.cohesion_cons**2/2/(lambda_t*np.tan(target.delta*math.pi/180)+mu_t*np.tan(np.pi/4+target.delta*math.pi/360))**2
    else:
        energy_min = target.cohesion_cons**2*np.cos(target.delta*math.pi/180)**2/(2*target.wave_speed**2*target.dens)

    energy_cons =  1/2*imp_M*(imp_v**2)*target.efficiency*target.energy
    
   # print(imp_d,energy_min - energy_cons)
    return energy_min - energy_cons


def plot_impactor(parameters):
  #  pdb.set_trace()
    parameters['run'] =1
    plt.rcParams.update({'font.size' : 16})
   # compute_deriv(target)
    #roots =parallel_root_computation(300, 50, 2)
    
    i=0
    j=0
    dexplicit =np.zeros((2,50))
    D = np.linspace(100,5000,50)
    for failure_mode in ['P-wave','S-wave']: 
        
        parameters['Failure wave'] = failure_mode
        j=0
        for dia in D:
 
            parameters['Diameter'] = dia
            target = Target(parameters)
    
            dexplicit[i,j] =root_scalar(compute_d,bracket=[0.01,1000],args=(target),method='brentq').root
            print(failure_mode,target.d,dexplicit[i,j])
            j+=1
            
        i+=1

    # Plot data
    plt.plot(D/1000,dexplicit[0,:],'-r', markersize=4, linewidth=2,label='P-waves')
    plt.plot(D/1000,dexplicit[1,:],'-b', marker='o', markersize=4, linewidth=2,label='S-waves')
    plt.xlabel("Target diameter (km)", fontsize=16)
    plt.ylabel("dmin", fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.yscale('log', base=2)
    plt.xscale('log', base=2)
    xticks =[0.1,0.2,0.5,1.0,2.0,3.0,5.0]
    labels = ['0.1','0.2','0.5','1.0','2.0','3.0','5.0']
    yticks =[0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0]
    ylabels = ['0.01','0.02','0.05','0.1','0.2','0.5','1.0','2.0']
    plt.xticks(xticks,labels)
    plt.yticks(yticks,ylabels)
    #%%
    
    
def compute_height(parameters,target):
    
    K_0 = target.K_0
    mu_0 = target.mu_0
    dens =target.dens
    radius = target.d/2
    omega = target.omega[2]
    beta = target.beta
   
    
    pi = math.pi
    R = target.d/2  
    theta = target.theta      
    rho = target.dens  
    imp_R = 0.1
    imp_dens =1500
    imp_M = 4/3*np.pi* imp_R**3*imp_dens
    imp_v = 5500
    
    target.wave_speed_P = np.sqrt(K_0+4/3*mu_0)*((radius**2/15/dens)*(4*np.pi*G*dens-2*omega**2))**(1/4)
    target.wave_speed_S = np.sqrt(mu_0)*((radius**2/15/dens)*(4*np.pi*G*dens-2*omega**2))**(1/4)  
    target.f = ((2*imp_dens*target.efficiency)/(np.pi*beta**2*dens))**(1/3)*2*target.wave_speed_P/(2*imp_R)
    target.k_d = 2*np.pi*target.f/target.Q
    target.energy, target.t_max = compute_energy(target)


    E = 1/2*imp_M*(imp_v**2)*target.efficiency*target.energy
    v_p = target.wave_speed_P  
    v_s = target.wave_speed_S 
    mu_t = v_s**2*rho
    lambda_t = (v_p**2-2*v_s**2)*rho
    
    base = 1
    dbase = 0
    metric = np.sqrt(base**2 + dbase**2)
    J =  metric*base*np.sin(theta)
    kappa_phi_g = 1/J*(dbase*np.sin(theta)+base*np.cos(theta))
    kappa_phi_n = -1/J*(dbase*np.cos(theta)-base*np.sin(theta))
    rgrav = -G*4/3*np.pi*R*rho
    tgrav = 0
    normal_p = rho*np.abs(rgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_n)
    tangential_p = rho*np.abs(tgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_g)
    
    if parameters['Failure wave']  == "P-wave":
        gamma_max = np.sqrt(2*E/(lambda_t+2*mu_t))*(lambda_t*np.tan(target.delta*pi/180)+mu_t*np.tan(np.pi/4+target.delta*pi/360))
    else:
        gamma_max = v_s*np.sqrt(2*rho*E)/np.cos(target.delta*pi/180)
        
    alpha = (np.tan(target.delta*pi/180)*normal_p -tangential_p +target.cohesion_linear)

    height = (gamma_max-target.cohesion_cons)/alpha

    return height


    
def fail_height(parameters):
    fig = plt.figure(figsize=(12, 5), dpi=300) 
    parameters['run'] =1
    plt.rcParams.update({'font.size' : 9})

    height =np.zeros((2,1000))
    i=0
    for failure_mode in ['P-wave','S-wave']: 
        
        parameters['Failure wave'] = failure_mode
        
        target = Target(parameters)

        height[i,:] = compute_height(parameters, target)
        i+=1

        # Plot data
    plt.plot(target.theta,height[0,:],'-r',  linewidth=1.5,label='P-waves')
    plt.plot(target.theta,height[1,:],'-b',  linewidth=1.5,label='S-waves')
    plt.xlabel("Colatitude")
    plt.ylabel("Failure height (m)")
    plt.xscale('log', base=2)
    xticks = [np.pi/12, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi]
    xticks2 = [np.pi/12,np.pi/8, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3, np.pi]
    labels = [r'$\dfrac{\pi}{12}$', r'$\dfrac{\pi}{6}$', r'$\dfrac{\pi}{4}$', r'$\dfrac{\pi}{3}$', r'$\dfrac{\pi}{2}$', r'$\dfrac{2\pi}{3}$', r'$\dfrac{3\pi}{4}$', r'$\dfrac{5\pi}{6}$', r'$\pi$']
    labels2 = [r'$\dfrac{\pi}{12}$',r'$\dfrac{\pi}{8}$', r'$\dfrac{\pi}{6}$', r'$\dfrac{\pi}{4}$', r'$\dfrac{\pi}{3}$', r'$\dfrac{\pi}{2}$', r'$\dfrac{2\pi}{3}$', r'$\pi$']

    plt.xticks(xticks,labels)
    plt.grid(True)
    plt.tight_layout()
  #  plt.xticks(fontsize=10)
   # plt.yticks(fontsize=10)

#%%

if __name__=="__main__":
    parameters = debug_main()
    plot_impactor(parameters)