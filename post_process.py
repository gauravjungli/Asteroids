#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from pathlib import Path
import glob
import os
import re
from gaurav import Parameter

#from mpl_toolkits.mplot3d import Axes3D
parameters={}
Parameter(parameters)
omega=float(parameters["omega_in"])
delta=float(parameters["delta"])
slides=int(parameters["slides"])
epsilon=float(parameters["epsilon"])
gamma=float(parameters["gamma"])
res=int(parameters["res"])
offset=float(parameters["offset"])
dx=(math.pi-2*offset)/res

file1="output/files_"+str(format(delta,".6f"))+"_"+str(format(omega,".6f"))
# omega=np.loadtxt(file1+"/omega.txt",delimiter=" ")
# plt.clf()
# plt.plot(omega[:,0],omega[:,1])

#plt.close()
#%%    
fig = plt.figure(figsize=(6,6))
for count in range(1,slides):
    file=glob.glob(file1+"/field_"+str(count+1)+".csv",recursive=True)
    w=np.loadtxt(file[0],delimiter=",",dtype=float)
    print(file)
   # if count!=slides-1 and count!=0:
    #    continue
    #plt.clf()
    x=np.sin(w[:,0])*(1+(epsilon*w[:,1]+gamma*w[:,2]))
    y=np.cos(w[:,0])*(1+(epsilon*w[:,1]+gamma*w[:,2]))
    plt.clf()
    plt.axis('equal')
    plt.plot(x,y,'-b')
    x=-np.sin(w[:,0])*(1+(epsilon*w[:,1]+gamma*w[:,2]))
    plt.plot(x,y,'-b')
    plt.title("lanslide number="+str(count+1))
    plt.pause(0.2)
    plt.savefig(file1+"/img_"+str(count+1))
#plt.close()
 #%%   
fig = plt.figure(figsize=(6,6))  
dirFiles = os.listdir(file1+"/data") #list of directory files
dirFiles.sort(key=lambda f: int(re.sub('\D', '', f)))
os.chdir(file1+"/data")
ang_mom=[]
lin_mom=[]
count=0
for file in dirFiles:

    
    w=np.loadtxt(file,delimiter=",",dtype=float)
  #  print(file)
    
    plt.clf()
    x=(w[:,0])
    y=(w[:,1]+gamma/epsilon*w[:,2])
    # y=(w[:,3])
    ang_mom.append([count,sum(w[:,4])*dx])
    lin_mom.append([count,sum(w[:,3])*dx])
    count +=1
    if count%10!=0:
          continue
    plt.plot(x,y)
    plt.title("Time="+str(count))
    plt.pause(0.1)  
    # print(sum(w[:,1]))  
# ang_mom=np.array(ang_mom)
# plt.plot(ang_mom[:,0],ang_mom[:,1])
# lin_mom=np.array(lin_mom)
# plt.plot(lin_mom[:,0],lin_mom[:,1])
os.chdir("..")
#%%
x=w[:,0]
grav=np.loadtxt(file1+"/grav.txt")
plt.plot(x,grav[2:-2,0])
plt.plot(x,grav[2:-2,1])
plt.grid()
# for count in range(len(omega)):

#     file=glob.glob(file1+"/data/field_"+str(count)+".csv",recursive=True)
#     w=np.loadtxt(file[0],delimiter=",",dtype=float)
#     print(file)
    
#     plt.clf()
#     x=(w[:,0])
#     y=(w[:,1]*w[:,3])
#   #  if count%50!=0:
#    #     continue
#     plt.plot(x,y)
#     plt.title("Time="+str(count))
#     plt.pause(0.1)  
#plt.close()  

#%%

omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
#plt.clf()
#plt.plot(omega[:,0],(2*math.pi/omega[:,1]/3600),linewidth=2)
plt.plot(omega[:,0],(2*math.pi/omega[:,2]/3600),linewidth=2)
plt.plot(omega[:,0],(2*math.pi/omega[:,3]/3600),linewidth=2)

#plt.plot(omega[:,0],(omega[:,1]),linewidth=2)
#plt.plot(omega[:,0],(omega[:,2]),linewidth=2)