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
import seaborn as sns

#from mpl_toolkits.mplot3d import Axes3D
plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("rocket",7)

parameters={}
Parameter(parameters)
omega=float(parameters["omega_in"])
delta=float(parameters["Friction angle"])
slides=int(parameters["slides"])
epsilon=float(parameters["epsilon"])
Gamma=float(parameters["Gamma"])
res=int(parameters["res"])
offset=float(parameters["offset"])
dx=(math.pi-2*offset)/res
#omega=0.65
#delta=30
file1="/home/g/Asteroids/output/files_"+str(format(delta,".6f"))+"_"+str(format(omega,".6f"))
#file1="output/omega_15_0.65_0.002"
# omega=np.loadtxt(file1+"/omega.txt",delimiter=" ")
#%%

fig = plt.figure(figsize=(6,6))
for count in range(0,slides):
    file=glob.glob(file1+"/field_"+str(count+1)+".csv",recursive=True)
    w=np.loadtxt(file[0],delimiter=",",dtype=float)
    print(file)
   # if count!=slides-1 and count!=0:
    #    continue
    #plt.clf()
    x=np.sin(w[:,0])*(1+(epsilon*w[:,1]+Gamma*w[:,2]))
    y=np.cos(w[:,0])*(1+(epsilon*w[:,1]+Gamma*w[:,2]))
    plt.clf()
    plt.axis('equal')
    plt.plot(x,y,'-r',linewidth=4)
    x=-np.sin(w[:,0])*(1+(epsilon*w[:,1]+Gamma*w[:,2]))
    plt.plot(x,y,'-r',linewidth=4)
    plt.title("lanslide number="+str(count+1))
    plt.pause(0.5)
    plt.savefig(file1+"/img_"+str(count+1)+".svg",dpi=300,bbox_inches="tight")
    
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
    y=(w[:,1]+Gamma/epsilon*w[:,2])
    # y=(w[:,3])
    ang_mom.append([count,sum(w[:,4])*dx])
    lin_mom.append([count,sum(w[:,3])*dx])
    count +=1
    if count%1!=0:
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
fig = plt.figure(figsize=(14,6)) 
x=w[:,0]
grav=np.loadtxt(file1+"/grav.txt",delimiter=" ")
#plt.clf()
plt.plot(x,grav[:,1],linewidth=2,markersize=8)

plt.grid()
plt.xlabel(r'$\theta$')
plt.ylabel('Non-dimensionalized Normal Gravity')
plt.xlim([0,3.14])
plt.ylim([-0.25,0.25])
plt.minorticks_on()
plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
plt.tick_params(labelsize=14)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)

plt.savefig('output/grav_20_static.svg', dpi=300,bbox_inches="tight")


#%%
#fig = plt.figure(figsize=(14,6)) 
omega=np.loadtxt(file1+"/omega_L.txt",delimiter="\t")


#plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,markersize=8,color=colors[2],linestyle='dotted',label='L')

plt.xlabel('Time (Myr)')
plt.ylabel('Time Period (hr)')
plt.xlim([0,1])
plt.ylim([2.5,4.5])
plt.minorticks_on()
plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
plt.tick_params(labelsize=14)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.legend(loc='best')
plt.legend(fontsize=14) 
plt.savefig('output/Omega.svg', dpi=300,bbox_inches="tight")

#%%
#script for plot 1
file1="output/omega_15_0.65_0.002"
fig = plt.figure(figsize=(6,6)) 
omega=np.loadtxt(file1+"/omega_CL.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,2]/3600)),linewidth=2,markersize=8,color=colors[0],linestyle='dashdot',label='C')

omega=np.loadtxt(file1+"/omega_CLY.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,1]/3600)),linewidth=2,markersize=8,color=colors[6],linestyle='dashed',label='Y')

omega=np.loadtxt(file1+"/omega_L.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='o',mfc='w',markersize=8,color=colors[2],linestyle='None',label='L')


omega=np.loadtxt(file1+"/omega_CLY.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,(2*math.pi/omega[:,2]/3600),linewidth=2,marker='s',mfc='w',markersize=8,color=colors[3],linestyle='None',label='CY')

omega=np.loadtxt(file1+"/omega_CL.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[4],linestyle='None',label='CL')

omega=np.loadtxt(file1+"/omega_LY.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,markersize=8,color=colors[5],linestyle='solid',label='LY')

omega=np.loadtxt(file1+"/omega_CLY.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,markersize=8,color=colors[1],linestyle='dotted',label='CLY')

yticks = np.arange(2.5,4.6,0.5)
#plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
plt.yticks(yticks)
plt.xlabel('Time (Myr)')
plt.ylabel('Time Period (hr)')
plt.xlim([0,1])
plt.ylim([2.48,4.5])
plt.minorticks_on()
plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
plt.tick_params(labelsize=14)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.legend(loc='best')
plt.legend(fontsize=14) 
plt.savefig('output/Omega.svg', dpi=300,bbox_inches="tight")

#%%

#script for plot 2

file1="output/saved_data/files_15.000000_0.200000"

fig = plt.figure(figsize=(6,6)) 
omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,2]/3600)),linewidth=2,markersize=8,color=colors[0],linestyle='dashdot',label='CY')

omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,3]/3600)),linewidth=2,marker='o',mfc='w',markersize=8,color=colors[6],linestyle='None',label='CLY')


file1="output/saved_data/files_15.000000_0.800000"

omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,2]/3600)),linewidth=2,markersize=8,color=colors[1],linestyle='solid',label='CY')

omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,3]/3600)),linewidth=2,markersize=8,color=colors[4],linestyle='dotted',label='CLY')


yticks = [2.5,3,4,6,8,10,12,14]
#plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
plt.yticks(yticks)
plt.xlabel('Time (Myr)')
plt.ylabel('Time Period (hr)')
plt.xlim([0,1])
plt.ylim([2.3,14])
plt.minorticks_on()
plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
plt.tick_params(labelsize=14)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.legend(loc='best')
plt.legend(fontsize=14) 
plt.savefig('output/Omega_1.svg', dpi=300,bbox_inches="tight")


#%%

#script for plot 3

file1="/Users/kumargaurav/Asteroid_data/saved_data/files_15.000000_0.650000_1"

fig = plt.figure(figsize=(6,6)) 
omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,2]/3600)),linewidth=2,markersize=8,color=colors[0],linestyle='dashdot',label='CY')

omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,3]/3600)),linewidth=2,markersize=8,color=colors[3],linestyle='solid',label='CLY')



yticks = [3,3.5,4,4.5,5]
#plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
plt.yticks(yticks)
plt.xlabel('Time (Myr)')
plt.ylabel('Time Period (hr)')
plt.xlim([0,5])
plt.ylim([2.9,5.5])
plt.minorticks_on()
plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
plt.tick_params(labelsize=14)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.legend(loc='best')
plt.legend(fontsize=14) 
plt.savefig('output/Omega_2.svg', dpi=300,bbox_inches="tight")


#%%

#script for plot 3

file1="/Users/kumargaurav/Documents/Asteroids/output/saved_data/files_15.000000_0.650000_C"

fig = plt.figure(figsize=(6,6)) 
omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,1]/3600)),linewidth=2,markersize=8,color=colors[0],linestyle='dashdot',label='Y')

omega=np.loadtxt(file1+"/omega_T.txt",delimiter="\t")
plt.semilogy(omega[:,0]/1e+6,((2*math.pi/omega[:,2]/3600)),linewidth=2,markersize=8,color=colors[3],linestyle='solid',label='CY')



yticks = [3,4,5,10,15,20,25]
#plt.plot(omega[:,0]/1e+6,(2*math.pi/omega[:,3]/3600),linewidth=2,marker='^',mfc='w',markersize=8,color=colors[0],linestyle='solid',label='L')
plt.yticks(yticks)
plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter())
plt.gca().yaxis.set_minor_formatter(mtick.NullFormatter())
plt.xlabel('Time (Myr)')
plt.ylabel('Time Period (hr)')
plt.xlim([0,5])
plt.ylim([2.8,25.5])
plt.minorticks_on()
plt.tick_params(direction='in',right=True, top=True, left=True, bottom=True)
plt.tick_params(labelsize=14)
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.legend(loc='best')
plt.legend(fontsize=14) 
plt.savefig('output/Omega_3.svg', dpi=300,bbox_inches="tight")