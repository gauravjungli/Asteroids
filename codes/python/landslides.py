#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:13:39 2023

@author: kumargaurav
"""

import math
import subprocess
from circle_fit import taubinSVD
import numpy as np
import os
from collisions import G
import time
from Diffusion_spherical import  compute_height,compute_time, compute_energy
from IO import Parameter, Exparameter, Output_File
from scipy.optimize import least_squares
from gravity import Gravitycalc
from script_2D.Crater import Crater
import pdb
from scipy.interpolate import UnivariateSpline, CubicSpline  , interp1d  
from Fit import Fit, Fit_2D

  
 #%%  


def Height(parameters,target,impactor=None):
    
    """
        This calculates the failure height of the landslide. First check whether landslide model is included in the 
        simulation or not. If yes then check whether it is initiated by an impact or it is a rotational failure. 
        For impacts, update epsilon for energy dependent simulations and select a failure profile from Gaussian and
        uniform. Update these in the base array and save it in a text file that will be later used by the C++ module.
    """ 
    
    if  not target.landslide:
        return
    
    mydir=Output_File (parameters,"output",["base.txt"])

    epsilon=float(parameters["epsilon"])
    base=np.loadtxt(mydir,delimiter=",",dtype=float)
    
    min_epsilon = float(parameters['Minimum epsilon'])
    max_epsilon = float(parameters['Maximum epsilon'])
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
        target.wave_speed =target.wave_speed_P
    else:
        target.wave_speed =target.wave_speed_S 
    
    if parameters['Dimension'] == '1D':
        
        res=int(parameters["Resolution"])
        height =np.ones(res)
        
        if impactor:
            #compute landslide duration
            #target.t_lan = compute_time(target,impactor)
            H, target.Gamma = compute_height(parameters,target=target,impactor=impactor,grid = base,dimension='1D') 
            print("The average failure height, maximum acceleration and diameter are ", H, target.Gamma, impactor.d )
            if H<min_epsilon:
                print("Not enough energetic imactor. No landslide simulated")
                parameters['epsilon']=H
                return
            if (target.Gamma*np.exp(-target.k_d/2*1/(G*4/3*3.14*1250)**(0.5)*0.3)<1):
                return
            epsilon = np.clip(H,min_epsilon,max_epsilon) #change_P make it constant if you don't want impactor dependent failure height
            
            print(f'the destablization height is {epsilon} and the impactor diameter is {impactor.d}\
                  and impactor velocity is {impactor.vel}' )
        
        else:
       
            print("No impactor found to inititate landslide. Probably its a rotational failure")
            epsilon=10*min_epsilon
            target.t_lan = np.sqrt(target.d/2/target.grav)
        
        parameters['epsilon'] = epsilon
        base[:,2] = -height[:]
        base = Fit(parameters,base_old=base)
        
        base[:,2]=height[:] #np.load('/home/g/Asteroids/output/check_5/height.npy')#height[:]        //change
    
    else:
        
        crater_data, crater_depth = Crater(parameters,target,impactor) #Fourth column contains the flag to check whether the point lies inside the crater
        #pdb.set_trace()
        
        H,target.Gamma = compute_height(parameters,target=target,impactor=impactor,grid=base,dimension='2D',crater_depth = crater_depth)
        H = H#*(1-crater_data[:,3])
        if (target.Gamma*np.exp(-target.k_d/2*1/(G*4/3*3.14*1250)**(0.5)*0.3)<0.5):
            return
       # epsilon = (np.min(H) + np.max(H))/2
       # epsilon = np.clip(epsilon,min_epsilon,max_epsilon)

        base[:,2] = (crater_data[:,2]- H[:])/epsilon 
        base[:,3] = H[:]/epsilon  
        
        
      #  parameters['epsilon'] = epsilon   
    print(f"The epsilon: {parameters['epsilon']}")

    

    np.savetxt(mydir,base,delimiter=",")



#%%
"""
This is the main function which calls the Landslide module:
    1. Update YORP if sochastic YORP is simulated. 
    2. Nondimensionalize all the relevant values such as omega, and inertia.
    3. Update slide count, impact time and dump directory of the landslide module.
    4. Export updated parameters
    5. Check whether the executable exist, the directories exist and create relevant directories
    6. Run the landslide executables
    7. Import the update values in the parameters
    8. Create the best fit circle and update the parameters dict
    9. Update omega,diameter and inertia of the target asteroid.
    10. Update gravity if necessary.
    
"""

def Landslides(target,parameters,impacttime,myomega):
    
    
    if  not target.landslide:
        return
    if float(parameters['epsilon'])<float(parameters['Minimum epsilon']):
        print(f"Value of epsilon {float(parameters['epsilon'])} too small to simulate landslide")
        return
    if (target.Gamma*np.exp(-target.k_d/2*1/(G*4/3*3.14*target.dens)**(0.5)*0.3)<1):
        print("Too small time for global landslides. Not simulating global landslide")
        return

    parameters["omega"] = target.omega[2] /(G * (4/3) * math.pi * target.dens)**0.5
    parameters['Maximum acceleration'] = target.Gamma
    parameters["slides"] = int(parameters["slides"])+1
    parameters["Current diameter"] = target.d
    parameters["jinertia"]=target.jinertia[2]/(target.d/2)**5/target.dens
    parameters["jinertia1"]=target.jinertia[0]/(target.d/2)**5/target.dens
    parameters["time"]=impacttime
    parameters["k_d"] = target.k_d/(G * (4/3) * math.pi * target.dens)**0.5
    
    
    print("The value of k_d is and the frequency is", parameters["k_d"], target.f)
    if bool(parameters['verbose']):
        parameters['verbose_dir']=Output_File(parameters,"output",['data',f"landslides_{parameters['slides']}"])
       
    try:
        Exparameter(parameters)
    except IOError:
        print("Failed in exporting the data")
    
   # exitcode = subprocess.run(["/home/g/Asteroids/build/asteroid"], cwd=".", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


    # 1. Ensure Executable Path & Permissions
    executable_path=Output_File(parameters,"output",[parameters["executable"]])
    if not os.path.isfile(executable_path) or not os.access(executable_path, os.X_OK):
        raise FileNotFoundError(f"Executable not found or not executable: {executable_path}")
    
    # 2. Adjust Working Directory if Necessary
    working_directory = Output_File(parameters,"output")  # Update if the executable is elsewhere
    if not os.path.isdir(working_directory):
        raise NotADirectoryError(f"Working directory not found: {working_directory}")
    
    try:
        # 3. Capture Output for Debugging (initially)
        os.mkdir(parameters['verbose_dir'])
        result = subprocess.run([executable_path], cwd=working_directory, capture_output=True, text=True) 
    
        if result.returncode != 0:
            print("Aborting due to error in running the executable landslide")
            raise subprocess.CalledProcessError(result.returncode, executable_path, output=result.stderr)
            
        else:
            print("Ran successfully slide", parameters["slides"])
    
        # If successful, you can later redirect output to DEVNULL 
        # result = subprocess.run([executable_path], cwd=working_directory, 
        #                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"Error output:\n{e.stderr}")  
    

    Parameter(parameters,"output")
    
    if parameters['Dimension'] == '1D': #expand later for 2D
        Fit(parameters)
        Parameter(parameters,"output")
    
        
        target.epsilon+=float(parameters["epsilon"])
        if target.epsilon>=0.001 or float(parameters["omega"])>0.95: 
            print("Updating gravity")
            target.epsilon = 0
            target.rgrav, target.tgrav = Gravitycalc(parameters,target)
    
    else:
        Fit_2D(parameters)
        
    target.omega[2] = float(parameters["omega"]) * (G * (4/3) * math.pi * target.dens )**0.5
    target.roots = target.roots*target.d/float(parameters["Current diameter"])
    target.d = float( parameters["Current diameter"])
    
    r = target.d / 2
    target.jinertia[2] = float(parameters["jinertia"]) * r**5 * target.dens
    target.jinertia[0] = target.jinertia[1] = float(parameters["jinertia1"]) * r**5 * target.dens
    
    myomega.append([impacttime,target.omega[2]])
    print("Omega after the Landslides", target.omega[2])





    
