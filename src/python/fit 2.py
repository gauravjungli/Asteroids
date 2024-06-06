# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
from circle_fit import taubinSVD
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import re
import glob

parameters={}
dir = os.getcwd()
inputfile = os.path.join(dir, "parameters")
with open(inputfile, "r") as file:
    for line in file:
        line = line.strip()
        if '--' in line:
            continue
        if line:
            values = re.split(r"\s+",line)
            key = values[0]
            value = values[1]
            parameters[key] = value
      
        

#%%


omega=float(parameters["omega_in"])
delta=float(parameters["delta"])
slides=int(parameters["slides"])
res=int(parameters["res"])-4
epsilon=parameters["epsilon"]
file1="output/files_"+str(format(delta,".6f"))+"_"+str(format(omega,".6f"))



file=glob.glob(file1+"/field_"+str(slides)+".csv",recursive=True)
w=np.loadtxt(file[0],delimiter=",",dtype=float)
#print(file)

#plt.clf()
x=np.sin(w[:,0])*(1+epsilon*(w[:,1]+w[:,2]))
y=np.cos(w[:,0])*(1+epsilon*(w[:,1]+w[:,2]))

# fig = plt.figure(figsize=(6,6))
# plt.axis('equal')
# plt.plot(x,y)


point = []
for i in range(res):
    point.append([x[i],y[i]])
    point.append([-x[i],y[i]])
xc, yc, r, sigma = taubinSVD(point)
parameters["radius"]=r

with open("parameters","w") as f:
    for key in parameters.keys():
         
        f.writelines([key, "    ", str(parameters[key]), " ", "\n" ]  )
    
print ("The best fit value of r is ", r)

# z= np.linspace(0,math.pi,200)
# zx=xc+r*np.sin(z[:])
# zy=yc+r*np.cos(z[:])
# plt.plot(zx,zy)
