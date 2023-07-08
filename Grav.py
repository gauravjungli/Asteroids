# -*- coding: utf-8 -*-
"""
Created on Sat May 23 13:31:27 2020 
Modified in June 2023
@author: Deepayan Banik & Kumar Gaurav
"""
from scipy.special import ellipk, ellipe,elliprf,elliprj
from scipy.interpolate import make_interp_spline
from scipy.constants import gravitational_constant
import numpy as np
from mpmath import ellippi
import math
import matplotlib.pyplot as plt
import time
import multiprocessing 
from itertools import repeat
import glob


parameters={}
# reading the data from the file
with open('parameters') as f:
    for line in f:
        fields = line.split("    ")
        parameters[fields[0]]=float(fields[1])
        
Res=int(parameters["res"]-4)
omega=parameters["omega_initial"]
delta=parameters["delta"]
slides=int(parameters["slides"])
epsilon=parameters["epsilon"]
gamma=parameters["gamma"]
density=parameters["density"]
rad=parameters["dia"]*1000/2

file1="output/files_"+str(format(delta,".6f"))+"_"+str(format(omega,".6f"))
file=glob.glob(file1+"/field_"+str(slides)+".csv",recursive=True)
w=np.loadtxt(file[0],delimiter=",",dtype=float)

#Ellippi = np.loadtxt('I.txt')


R=rad*np.sin(w[:,0])*(1+gamma*(w[:,1])+epsilon*w[:,2])
Z=rad*np.cos(w[:,0])*(1+gamma*(w[:,1])+epsilon*w[:,2])

fR=make_interp_spline(w[:,0],R)
fZ=make_interp_spline(w[:,0],Z)


res=10000

theta=np.linspace(0,math.pi,res)
r = fR(theta)
z = fZ(theta)

def Gravity(R, Z): # z = vertical lower limit of current disc, R, Z -> radial and vertical location of current point of evaluation
    
    r_grav=0
    z_grav=0
    #print(res)
    for i in range(1,res-1):
        a = r[i] # radius of disc being integrated
        zeta = Z - z[i] # vertical disctance of disc from point of evaluation
        delta = np.sqrt((a + R)**2 + (zeta)**2) # parameter for elliptic integrals
        k = 2 * np.sqrt(a * R) / delta  # 
        m = 2 * np.sqrt(a * R) / (a + R)
        if (R < a):
            eps = 1
        elif (R > a):
            eps = 0
        else:
            eps = 0.5
        if (k>=1 or m>=1):
            print("k = ",k," m = ",m,Z,z[i],i)
        ks = ellipk(k**2)
        es = ellipe(k**2)
        #pi= Ellippi[min(999,round(m**2*1000)),min(999,round(k**2*1000))]
        #pi=ellippi(m**2,k**2)
        pi=elliprf(0,1-k**2,1)+1/3*m**2*elliprj(0,1-k**2,1,1-m**2)

        r_grav=r_grav + np.abs((z[i]-z[i-1]))*(2 * delta * ((1 - k**2 / 2) * ks - es) / R)
        z_grav=z_grav + np.abs((z[i]-z[i-1]))*(2 * np.pi * np.sign(zeta) * eps + 2 * zeta * ((R - a)/(R + a) * pi - ks) / delta)
    return r_grav, z_grav



def main():
    
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    
  
    arguments=[(R[i]+rad*epsilon*np.sin(w[i,0]),Z[i]+ rad*epsilon*np.cos(w[i,0])) for i in range(int(Res/2))]
    
    G=gravitational_constant
    grav = pool.starmap(Gravity, arguments)
    grav=np.array(grav)
    grav1=np.array([(grav[i,0],-grav[i,1]) for i in range(round(Res/2)-1,-1,-1)])
    grav=np.vstack((grav,grav1))
    
     
   
    r_grav=-G*density*(grav[:,0]*np.sin(w[:,0])+grav[:,1]*np.cos(w[:,0]))/(4/3*np.pi*density*rad*G)
    t_grav=-G*density*(grav[:,0]*np.cos(w[:,0])-grav[:,1]*np.sin(w[:,0]))/(4/3*np.pi*density*rad*G)
    plt.plot(w[:,0],r_grav)
    plt.plot(w[:,0],t_grav)
    grav=np.hstack((r_grav.reshape(-1,1),t_grav.reshape(-1,1)))
    np.savetxt(file1+"/grav.txt",grav)
    
    # Close the pool
    pool.close()
    pool.join()
    
if __name__=="__main__":
    
    

    start_time=time.time()
    
    print(res)


    multiprocessing.freeze_support()
    main()


    end_time=time.time()
    elapsed_time=end_time-start_time
    print("Time taken:", elapsed_time, "seconds")



#%%

####### For generating tables of elliptic functions

    # k=np.linspace(0,1,100001)
    # E=pool.map(ellipe,np.array(k))
    # Ew=np.array(list(zip(np.array(k), np.array(E) )))
    # np.savetxt("E.txt",Ew)
    # content = np.loadtxt('E.txt')
    # print(content[10,1])

    # K=pool.map(ellipk,np.array(k))
    # Kw=np.array(list(zip(np.array(k), np.array(K) )))
    # np.savetxt("K.txt",Kw)
    # content = np.loadtxt('K.txt')
    # print(np.sum(content[:,1]-K))
    # print(content[10,1])
    
    
    #arguments = [(m/res, k/res) for m in range(res+1) for k in range(res+1)]
    #results = pool.starmap(calculate_ellippi, arguments)
    #I = [results[i:i+res+1] for i in range(0, len(results), res+1)]
    #np.savetxt("I.txt",I)
    
    
    #%%
    

    
