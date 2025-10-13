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

from collisions import   Collision, Istuff
from landslides import Landslides, Height
from Target import Target
from IO import ExportOmega, Output_File, Parameter, Cumdistr
from Initialize import Initialize
import multiprocessing 
import time
from YORP import YORP, shape_gen


#%%



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
    myfile=Output_File(parameters,"output",["python_log.txt"])
    log_file=open(myfile, "w") 
    sys.stdout = log_file 
    sys.stderr = log_file
# =============================================================================
    start_time = time.time()
    target = Target(parameters)
    cumdistr = Cumdistr(parameters)
    Initialize(parameters,target)
    tmaxby=float(parameters['Simulation period'])*1e+6
#    fig = plt.figure(figsize=(10,6))

    istuff = Istuff(parameters,target,tmaxby,cumdistr)
    oldtime = 0
    myomega=[[0,target.omega[2]]]
    
    print("Starting real simulation")
    sys.stdout.flush() 
    for i in range(len(istuff)):

        YORP(target,parameters,istuff[i].impacttime,oldtime,myomega)
        
        
        Collision(target,istuff[i],myomega)
        
        if istuff[i].explicit:
            print("Calling landslide")
            Height(parameters,target,istuff[i])
            Landslides(target,parameters,istuff[i].impacttime,myomega)
            if parameters["stoc_yorp"].lower()=='yes' and int(parameters["slides"])%10==0:
                target.coeff_f, target.coeff_g = shape_gen(target.K)

               
             
        #print(istuff[i].d)
        oldtime = istuff[i].impacttime
        #sys.__stdout__ always points to the original standard output stream. 99 is set to max because 
        #near completion time it will be very close to 100 that will stop the simulation
        
        print(min(99,round(oldtime/tmaxby*100)),file=sys.__stdout__,flush=True)
        sys.stdout.flush() 
        ExportOmega(myomega,parameters) 
        #plt.clf()
        #plt.plot([data[0] for data in myomega],[(2*math.pi/data[1]/3600) for data in myomega])
        #plt.show()
        #plt.pause(0.5)
    print("All collisions simulated. Final YORP simulations")
    YORP(target,parameters,tmaxby,oldtime,myomega)
    ExportOmega(myomega,parameters)    
    
    executable_file=Output_File(parameters,"output",[parameters["executable"]])
    
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
    
    
