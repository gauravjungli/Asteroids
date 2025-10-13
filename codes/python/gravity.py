#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 18:16:35 2025

@author: g
"""
import time
from IO import Output_File
import numpy as np
from scipy.interpolate import make_interp_spline, CubicSpline
import multiprocessing
from collisions import G
from scipy.special import ellipk, ellipe,elliprf,elliprj
import matplotlib.pyplot as plt
#%%
"""
Calculates gravity for the axisymmetric body.
"""

def Gravitycalc(parameters):
    
    start =time.time()
    Res=int(parameters["Resolution"])
    epsilon=1e-3 #This is a different epsilon
    density=float(parameters["Density"])
    rad=float(parameters["Current diameter"])/2
    
    file=Output_File(parameters,"output",["base.txt"])

    try:
        w=np.loadtxt(file,dtype=float,delimiter=",")
    except:
            print("No file available for the fit")
            return
        

    R=rad*np.sin(w[:,0])*w[:,1]
    Z=rad*np.cos(w[:,0])*w[:,1]
    fR=make_interp_spline(w[2:Res-2,0],R[2:Res-2]) 
    fZ=make_interp_spline(w[2:Res-2,0],Z[2:Res-2])

    res=10000

    theta=np.linspace(0,np.pi,res)
    r = fR(theta)
    z = fZ(theta)
    
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    
    arguments=[(R[i]+epsilon/2*rad*np.sin(w[i,0]),Z[i]+ epsilon/2*rad*np.cos(w[i,0]),r,z) for i in range(int(Res/2))]#changed it from Res/2
    
    grav = pool.starmap(Gravity, arguments)
    grav=np.array(grav)
    grav1=np.array([(grav[i,0],-grav[i,1]) for i in range(round(Res/2)-1,-1,-1)])#change uncomment if using equator symmetry
    grav=np.vstack((grav,grav1))
    
    R_grav = -G*density*(grav[:,0]*np.sin(w[:,0])+grav[:,1]*np.cos(w[:,0]))
    T_grav = -G*density*(grav[:,0]*np.cos(w[:,0])-grav[:,1]*np.sin(w[:,0]))
    
    #Only for the axisymmetric case
    theta = w[2:Res-2,0]
    base = w[2:Res-2,1]
    height = w[2:Res-2,2]
    dbase = w[2:Res-2,3]

    metric = base**2 + dbase**2
    for i in range(2,Res-2):
        R_grav_temp = R_grav[i]
        R_grav[i] = (R_grav[i]*base[i-2]-dbase[i-2]*T_grav[i])/np.sqrt(metric[i-2])
        T_grav[i] = (T_grav[i]*base[i-2]+dbase[i-2]*R_grav_temp)/np.sqrt(metric[i-2])
    r_grav = R_grav/(4/3*np.pi*density*rad*G)
    t_grav = T_grav/(4/3*np.pi*density*rad*G)
   # plt.plot(w[:,0],r_grav)
   # plt.plot(w[:,0],t_grav)
    print("gravity updated")
    grav=np.hstack((r_grav.reshape(-1,1),t_grav.reshape(-1,1)))
    file=Output_File(parameters,"output",["grav.txt"])
    np.savetxt(file,grav)
    pool.close()
    pool.join()       
    end =time.time()
    print(f"Time taken in calculating gravity:{end-start}")
    return R_grav, T_grav



def Gravitycalc_2D(parameters,lats,lons,devs,radius,lats_comp,lons_comp,devs_comp):
    
    start =time.time()
    epsilon=0.001 #This is a different epsilon
    density=float(parameters["Density"])
    rad=radius*1000
    devs = devs*1000
    
    Res = lats_comp.shape[0]  
    R=np.sin(lats)*(rad +devs)
    Z=np.cos(lats)*(rad +devs)
    
    R_comp=np.sin(lats_comp)*(rad +devs_comp)
    Z_comp=np.cos(lats_comp)*(rad +devs_comp)
    
    fR=make_interp_spline(lats,R)
    fZ=make_interp_spline(lats,Z)

    res=10000

    theta=np.linspace(0,np.pi,res)
    r = fR(theta)
    z = fZ(theta)
    
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    
    arguments=[(R_comp[i]+epsilon/2*rad*np.sin(lats_comp[i]),Z_comp[i]+ 
                epsilon/2*rad*np.cos(lats_comp[i]),r,z) for i in range(int(Res))]
    
    grav = pool.starmap(Gravity, arguments)
    grav=np.array(grav)
    
    
    R_grav = -G*density*(grav[:,0]*np.sin(lats_comp)+grav[:,1]*np.cos(lats_comp))
    T_grav = -G*density*(grav[:,0]*np.cos(lats_comp)-grav[:,1]*np.sin(lats_comp))
    r_grav = R_grav/(4/3*np.pi*density*rad*G)
    t_grav = T_grav/(4/3*np.pi*density*rad*G)

    lats_comp_2D, lons_comp_2D = np.meshgrid(lats_comp,lons_comp,indexing='ij')
    r_grav_2D = np.tile(r_grav.reshape(-1,1),(1,len(lons_comp)))
    t_grav_2D = np.tile(t_grav.reshape(-1,1),(1,len(lons_comp)))
    gravity_data = np.stack((r_grav_2D.ravel(),t_grav_2D.ravel())).T 
    print("gravity updated")
    plt.plot(lats_comp,r_grav)
    file=Output_File(parameters,"output",["grav.txt"])
    np.savetxt(file,gravity_data)
    
    pool.close()
    pool.join()       
    end =time.time()
    print(f"Time taken in calculating gravity:{end-start}")
    return R_grav, T_grav


def Gravity(R, Z,r,z): 
    
    r_grav=0
    z_grav=0
    for i in range(1,len(r)-1):
        a = r[i] # radius of disc being integrated
        zeta = Z - z[i] # vertical disctance of disc from point of evaluation
        delta = np.sqrt((a + R)**2 + (zeta)**2) # parameter for elliptic integrals
        k = 2 * np.sqrt(a * R) / delta  
        m = 2 * np.sqrt(a * R) / (a + R)
        if (R < a):
            eps = 1
        elif (R > a):
            eps = 0
        else:
            eps = 0.5
        if (k>=1 or m>=1):
            print("k = ",k," m = ",m,Z,z[i],a,R)
            m=min(m,1-1e-6)
        ks = ellipk(k**2)
        es = ellipe(k**2)
        pi=elliprf(0,1-k**2,1)+1/3*m**2*elliprj(0,1-k**2,1,1-m**2)

        r_grav = r_grav + np.abs((z[i]-z[i-1]))*(2 * delta * ((1 - k**2 / 2) * ks - es) / R)
        z_grav = z_grav + np.abs((z[i]-z[i-1]))*(2 * np.pi * np.sign(zeta) * eps + 2 * zeta * ((R - a)/(R + a) * pi - ks) / delta)
    return r_grav, z_grav    
    
