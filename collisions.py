#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:36:25 2023

@author: kumargaurav, chatgpt and Holsapple
"""

from functions import   Collision
from gaurav import Parameter, Initialize, Target, ExportOmega, Yorp, Landslides, Cumdistr, Istuff
import numpy as np
import multiprocessing 
import matplotlib.pyplot as plt
import time



#%%


if __name__=="__main__":
    
    start_time=time.time()
    parameters=Parameter()
    target=Target(parameters)
    cumdistr=Cumdistr(parameters)
    Initialize(parameters,target)
    tmaxby=float(parameters["tmaxby"])
    istuff = Istuff(target,tmaxby,cumdistr)
    oldtime = 0
    myomega=[]
    multiprocessing.freeze_support()
    
    for i in range(len(istuff)):
        
        Yorp(target,parameters,istuff[i].impacttime,oldtime,myomega)
        
        Collision(target,istuff[i])
        
        if istuff[i].explicit:
            Landslides(target,parameters,istuff[i].impacttime,myomega)
               
        #print(istuff[i].d)
        oldtime = istuff[i].impacttime
        #print(oldtime)
    
    Yorp(target,parameters,tmaxby,oldtime,myomega)
    ExportOmega(myomega,parameters)
    
        
    end_time=time.time()
    elapsed_time=end_time-start_time
    print("Time taken:", elapsed_time, "seconds")

