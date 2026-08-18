#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:13:39 2023

@author: kumargaurav
"""

import math
import subprocess
import numpy as np
import os
from collisions import G
import time
from IO import Parameter, Exparameter, Output_File
from gravity import Gravitycalc
import pdb
from Fit import Fit, Fit_2D

  
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

def Landslides(parameters,target,impacttime,myomega):
    
    
    if  not target.landslide:
        return
    if target.epsilon<target.min_epsilon:
        print(f"Value of epsilon {target.epsilon} too small to simulate landslide")
        return
    if not (target.fast_rotation_flag or parameters["Single run"].lower() == 'yes') and (target.Gamma*np.exp(-target.k_d/2*1/(G*4/3*3.14*target.dens)**(0.5)*0.3)<1):
        print("Too small time for global landslides. Not simulating global landslide")
        return
    
    parameters['epsilon'] =   target.epsilon 
    parameters["omega"] =  target.omega[2] /(G * (4/3) * math.pi * target.dens)**0.5
    parameters['Maximum acceleration'] = target.Gamma
    target.slides += 1 
    parameters["slides"] = target.slides
    parameters["Current diameter"] = target.d
    parameters["jinertia"] = target.jinertia[2]/(target.d/2)**5/target.dens
    parameters["jinertia1"] = target.jinertia[0]/(target.d/2)**5/target.dens
    parameters["time"] = impacttime
    parameters["k_d"] = target.k_d/(G * (4/3) * math.pi * target.dens)**0.5
    
    if (target.fast_rotation_flag):
    
        target.rgrav, target.tgrav = Gravitycalc(target)
    
    print("The value of k_d is and the frequency is", parameters["k_d"], target.f)
    if bool(parameters['verbose']):
        parameters['verbose_dir'] = Output_File(target,"output",['data',f"landslides_{parameters['slides']}"])
       
    try:
        Exparameter(parameters)
    except IOError:
        print("Failed in exporting the data")
    
   # exitcode = subprocess.run(["/home/g/Asteroids/build/asteroid"], cwd=".", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


    # 1. Ensure Executable Path & Permissions
    executable_path=Output_File(target,"output",[parameters["executable"]])
    if not os.path.isfile(executable_path) or not os.access(executable_path, os.X_OK):
        raise FileNotFoundError(f"Executable not found or not executable: {executable_path}")
    
    # 2. Adjust Working Directory if Necessary
    working_directory = Output_File(target,"output")  # Update if the executable is elsewhere
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
    target.shed_mass = float(parameters["Mass shed"])
    target.omega[2] = float(parameters["omega"]) * (G * (4/3) * math.pi * target.dens )**0.5
    
    if parameters['Dimension'] == '1D': #expand later for 2D
        Fit(target)
        Parameter(parameters,"output")
    
        
        if target.epsilon>=0.001 or float(parameters["omega"])>0.95: 
            print("Updating gravity")
            target.rgrav, target.tgrav = Gravitycalc(target)
    
    else:
        Fit_2D(parameters)
        
    target.roots = target.roots*target.d/float(parameters["Current diameter"])
    
   
    myomega.append([impacttime,target.omega[2]])
        
    print("Omega after the Landslides", target.omega[2])



#%%

