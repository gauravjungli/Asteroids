#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:36:25 2023

@author: kumargaurav, chatgpt and Holsapple
"""

from functions import   Collision
from gaurav import Parameter, Initialize, Target, ExportOmega, Yorp, Landslides, Cumdistr, Istuff, velave, Height, shape_gen
import numpy as np
import multiprocessing 
import matplotlib.pyplot as plt
import time
from functions import qstarf
import math
from YORP import YORP
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
    # For selecting cases where there is an impactor of size greater than 1 meter atleast

    istuff = Istuff(target,tmaxby,cumdistr)
    oldtime = 0
    myomega=[[0,target.omega[2],target1.omega[2],target2.omega[2]]]
    
    for i in range(len(istuff)):
      
        
      #three targets are given for the ease of plotting. target: YORP, target1: YC and target2: YCL
  #      Yorp(target,target1,target2,parameters,istuff[i].impacttime,oldtime,myomega)
        YORP(target,target1,target2,parameters,istuff[i].impacttime,oldtime,myomega)
        
        Collision(target1,target,target2,istuff[i],myomega)
        
        
        if istuff[i].explicit:
            Height(parameters,target2,istuff[i],"uniform")
            Landslides(target2,target,target1,parameters,istuff[i].impacttime,myomega)
            if bool(int(parameters["stoc_yorp"])):
                target.coeff_f, target.coeff_g = abs(shape_gen(target.K))
               
               
        print(istuff[i].d)
        oldtime = istuff[i].impacttime
        print("Total time till now",oldtime)
        plt.clf()
        plt.plot([data[0] for data in myomega],[(2*math.pi/data[1]/3600) for data in myomega])
        plt.plot([data[0] for data in myomega],[(2*math.pi/data[3]/3600) for data in myomega])
        plt.show()
        plt.pause(0.2)

    YORP(target,target1,target2,parameters,tmaxby,oldtime,myomega)
    myomega.append([tmaxby,target.omega[2],target1.omega[2]])
    ExportOmega(myomega,parameters,"Omega_gaurav")    
        
    end_time=time.time()
    elapsed_time=end_time-start_time
    print("Time taken:", elapsed_time, "seconds")

