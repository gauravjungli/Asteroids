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
    sys.stdout = open(myfile, "w") 
    start_time=time.time()
    target=Target(parameters)
    cumdistr=Cumdistr(parameters)
    Initialize(parameters,target)
    tmaxby=float(parameters['Simulation period'])
#    fig = plt.figure(figsize=(10,6))

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
        print(round(oldtime/tmaxby*100),file=sys.__stdout__,flush=True)
        sys.stdout.flush() 
        ExportOmega(myomega,parameters,target) 
        #plt.clf()
        #plt.plot([data[0] for data in myomega],[(2*math.pi/data[1]/3600) for data in myomega])
        #plt.show()
        #plt.pause(0.5)
    
    YORP(target,parameters,tmaxby,oldtime,myomega)
    myomega.append([tmaxby,target.omega[2]])
    ExportOmega(myomega,parameters,target)    
    print(100,file=sys.__stdout__,flush=True)

    end_time=time.time()
    elapsed_time=end_time-start_time
    
    print("Time taken:", elapsed_time, "seconds")
    sys.stdout.close()  
    
    
