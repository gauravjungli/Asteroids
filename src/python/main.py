#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:36:25 2023

@author: kumargaurav, chatgpt and Holsapple
"""

from collisions import   Collision
from gaurav import Parameter, Initialize, Target, ExportOmega, Landslides, Cumdistr, Istuff, Height, shape_gen,Output_File, Exparameter
import multiprocessing 
import matplotlib.pyplot as plt
import time
from YORP import YORP
import sys
import os
import subprocess
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
    myfile=Output_File(parameters,"output",["python_log.txt"])
    log_file=open(myfile, "w") 
    sys.stdout = log_file
    sys.stderr = log_file
    start_time=time.time()
    target=Target(parameters)
    cumdistr=Cumdistr(parameters)
    Initialize(parameters,target)
    tmaxby=float(parameters['Simulation period'])
#    fig = plt.figure(figsize=(10,6))

    istuff = Istuff(parameters,target,tmaxby,cumdistr)
    oldtime = 0
    myomega=[[0,target.omega[2]]]
    
    for i in range(len(istuff)):

  #      Yorp(target,target1,target2,parameters,istuff[i].impacttime,oldtime,myomega)

        YORP(target,parameters,istuff[i].impacttime,oldtime,myomega)
        
        Collision(target,istuff[i],myomega)
        
        
        if istuff[i].explicit:
            print("Calling landslide")
            Height(parameters,target,istuff[i])
            Landslides(target,parameters,istuff[i].impacttime,myomega)

               
             
        #print(istuff[i].d)
        oldtime = istuff[i].impacttime
        print(round(oldtime/tmaxby*100),file=sys.__stdout__,flush=True)
        sys.stdout.flush() 
        ExportOmega(myomega,parameters) 
        #plt.clf()
        #plt.plot([data[0] for data in myomega],[(2*math.pi/data[1]/3600) for data in myomega])
        #plt.show()
        #plt.pause(0.5)
    print("All collisions simulated. Final YORP simulations")
    YORP(target,parameters,tmaxby,oldtime,myomega)
    ExportOmega(myomega,parameters)    
    print(100,file=sys.__stdout__,flush=True)
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




    end_time=time.time()
    elapsed_time=end_time-start_time
    
    print("Time taken:", elapsed_time, "seconds")
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    log_file.close()  
    
    
