#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:36:25 2023

@author: kumargaurav, chatgpt and Holsapple
"""

from functions import   Collision
from gaurav import Parameter, Initialize, Target, ExportOmega, Landslides, Cumdistr, Istuff, Height, shape_gen,Output_File
import multiprocessing 
import matplotlib.pyplot as plt
import time
from YORP import YORP
import sys
import os
import subprocess
#%%



#%%

def simulate(run):
    
    start_time=time.time()
    #parameters={}
    #Parameter(parameters,"input")
    return
    target=Target(parameters)
    cumdistr=Cumdistr(parameters)

    Initialize(parameters,target,run)
    myfile=Output_File(parameters,"output",["python_log.txt"])
    sys.stdout = open(myfile, "w")  
    tmaxby=float(parameters["Simulation period"])
    fig = plt.figure(figsize=(10,6))

    istuff = Istuff(target,tmaxby,cumdistr)
    oldtime = 0
    myomega=[[0,target.omega[2]]]
    
    for i in range(len(istuff)):
      
        

  #      Yorp(target,target1,target2,parameters,istuff[i].impacttime,oldtime,myomega)

        YORP(target,parameters,istuff[i].impacttime,oldtime,myomega)
        
        Collision(target,istuff[i],myomega)
        
        
        if istuff[i].explicit:
            Height(parameters,target,istuff[i])
            Landslides(target,parameters,istuff[i].impacttime,myomega)
            if parameters["stoc_yorp"].lower()=='yes':
                target.coeff_f, target.coeff_g = abs(shape_gen(target.K))
               
               
        #print(istuff[i].d)
        oldtime = istuff[i].impacttime
        #print("Total time till now",oldtime)
        #plt.clf()
        #plt.plot([data[0] for data in myomega],[(2*math.pi/data[1]/3600) for data in myomega])
        #plt.show()
        #plt.pause(0.5)
    
    YORP(target,parameters,tmaxby,oldtime,myomega)
    myomega.append([tmaxby,target.omega[2],target.omega[2]])
    ExportOmega(myomega,parameters,target)    
        
    elapsed_time=end_time-start_time
    print("Time taken:", elapsed_time, "seconds")
    sys.stdout.close()  



if __name__=="__main__":
    
    multiprocessing.freeze_support()
    run = list(range(1, 10))  # List of input values
    parameters={}
    Parameter(parameters,"input")
    parameters["run"]=0
    mydir=Output_File (parameters,"output")
    if os.path.exists(mydir):
        subprocess.run(["rm", "-r", mydir])
    else:
        print("No directory exists")
    start_time=time.time()
    
    simulate(1)
   # with multiprocessing.Pool() as pool:
   #    pool.map(simulate, run) 
       
    end_time=time.time()
    elapsed_time=end_time-start_time
    
    
