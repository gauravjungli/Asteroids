#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:36:25 2023

@author: kumargaurav and chatgpt 
This is the main script which starts one simulation for each collision history. It needs some arguments: Output folder and run.
The log of this script is saved in python_log.txt and it gives output in form of % run to the progressbar. It creates an object target
and the collisional history which is then stored in the istuff object of the impactor class. It then calls the YORP, collisions
and landslide modules to simulate history and stores the value of omega in a file. 
"""
import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'script_1D')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'script_2D')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'Failure_profile')))

from Height import Height
from collisions import   Collision, Istuff
from landslides import Landslides
from Target import Target
from IO import ExportOmega, Output_File_old, Output_File, Parameter, Cumdistr
from Initialize import Initialize
import multiprocessing 
import time
from YORP import YORP, shape_gen

#%%

if __name__=="__main__":
    
    multiprocessing.freeze_support()
    parameters={}
    
    for i  in range(1,len(sys.argv),2):
        if i<1:
            continue
        else:
            parameters[sys.argv[i]]=sys.argv[i+1]

    Parameter(parameters, "output") 

# =============================================================================
    myfile=Output_File_old(parameters,"output",["python_log.txt"])
    log_file=open(myfile, "w") 
    sys.stdout = log_file 
    sys.stderr = log_file
# =============================================================================
    start_time = time.time()
    target = Target(parameters)
    cumdistr = Cumdistr(parameters)
    Initialize(parameters,target)
    myomega=[[0,target.omega[2]]]
    
    tmaxby=float(parameters['Simulation period'])*1e+6
    oldtime = 0
    
    if parameters["Single run"].lower() == 'yes':
        
        Height(target,)
        Landslides(parameters,target,0,myomega)
        
    elif not target.collision:
        YORP(parameters,target,tmaxby,oldtime,myomega)
        
    else:       
    
        
        istuff = Istuff(parameters,target,tmaxby,cumdistr)
        
       
        print("Starting real simulation")
        sys.stdout.flush() 
        
        for i in range(len(istuff)):
    
            YORP(parameters,target,istuff[i].impacttime,oldtime,myomega)
            
            Collision(target,istuff[i],myomega)
            
            if target.collision and istuff[i].explicit:
                #print("Calling landslide")
                Height(target,istuff[i])
                Landslides(parameters,target,istuff[i].impacttime,myomega)
                if parameters["stoc_yorp"].lower()=='yes' and int(parameters["slides"])%int(parameters['N_stoch'])==0:
                    target.coeff_f, target.coeff_g = shape_gen(target.K)
       
            oldtime = istuff[i].impacttime
            #sys.__stdout__ always points to the original standard output stream. 99 is set to max because 
            #near completion time it will be very close to 100 that will stop the simulation
            
            print(min(99,round(oldtime/tmaxby*100)),file=sys.__stdout__,flush=True)
            sys.stdout.flush() 
            ExportOmega(parameters,myomega) 
            
        print("All collisions simulated. Final YORP simulations")
        YORP(target,parameters,tmaxby,oldtime,myomega)
            
    ExportOmega(parameters,myomega)    
    
    executable_file=Output_File(target,"output",[parameters["executable"]])
    
    try:
        os.remove(executable_file)
        print(f"File '{executable_file}' has been deleted successfully.")
    except FileNotFoundError:
        print(f"File '{executable_file}' not found.")
    except PermissionError:
        print(f"Permission denied: '{executable_file}'.")
    except Exception as e:
        print(f"Error occurred while trying to delete the file: {e}")


    print(100,file=sys.__stdout__,flush=True)
    end_time=time.time()
    elapsed_time=end_time-start_time
    
    print("Time taken:", elapsed_time, "seconds")
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    log_file.close()  
    
    
