#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:36:25 2023

@author: kumargaurav, chatgpt and Holsapple
"""

from functions import   Collision
from gaurav import Parameter, Initialize, Target, ExportOmega, Yorp, Landslides, Cumdistr, Istuff,velave,Height
import numpy as np
import multiprocessing 
import matplotlib.pyplot as plt
import time
from functions import qstarf
import math

#%%



#%%

if __name__=="__main__":
    
    
    multiprocessing.freeze_support()
    parameters={}
    start_time=time.time()
    Parameter(parameters)
    target=Target(parameters)
    target1=Target(parameters)
    target2=Target(parameters)
    cumdistr=Cumdistr(parameters)
    Initialize(parameters,target)
    tmaxby=float(parameters["tmaxby"])
    fig = plt.figure(figsize=(10,6))
    while(True):
        istuff = Istuff(target,tmaxby,cumdistr)
        if max(node.d for node in istuff)>1:
            break
    oldtime = 0
    myomega=[[0,target.omega[2],target1.omega[2],target2.omega[2]]]
    
    
    for i in range(len(istuff)):
      
        Yorp(target,target1,target2,parameters,istuff[i].impacttime,oldtime,myomega)
        
        Collision(target1,target,target2,istuff[i],myomega)
        
        
        if istuff[i].explicit:
            Height(parameters,target2,istuff[i])
            Landslides(target2,target,target1,parameters,istuff[i].impacttime,myomega)
               
        print(istuff[i].d)
        oldtime = istuff[i].impacttime
        print("Total time till now",oldtime)
        plt.clf()
        plt.plot([data[0] for data in myomega],[(2*math.pi/data[2]/3600) for data in myomega])
        plt.plot([data[0] for data in myomega],[(2*math.pi/data[3]/3600) for data in myomega])
        plt.show()
        plt.pause(0.2)
        ExportOmega(myomega,parameters)
    Yorp(target,target1,target2,parameters,tmaxby,oldtime,myomega)
   # myomega.append([tmaxby,target.omega[2],target1.omega[2]])
    ExportOmega(myomega,parameters)
    
        
    end_time=time.time()
    elapsed_time=end_time-start_time
    print("Time taken:", elapsed_time, "seconds")

