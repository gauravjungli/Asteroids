#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:13:39 2023

@author: kumargaurav
"""
import random
import math
import subprocess
from scipy.interpolate import make_interp_spline
import numpy as np
import glob
from circle_fit import taubinSVD
import re
import os
from scipy.special import ellipk, ellipe,elliprf,elliprj
import matplotlib.pyplot as plt
import multiprocessing 
from functions import velave,G,getdiaf,qstarf,probi,astnum, wobblecalcf
import scipy.stats

#%%  Verified

class Target:
    
    def __init__(self, parameters):
        self.d=float(parameters["dia_in"])
        self.atype=parameters["atype"]
        self.delta=float(parameters["delta"])
        if self.atype == "S-Type":
            self.mu = 0.55
            self.dens = 2500
            self.Y0, self.d0strength   =  1.44e7, 0.1
            self.nsize = 3  # strength decreases with size as 1/nsize
            self.qconst1, self.qconst2 =  1e3, 1e6
            self.k1, self.k2=0.06, 1
        else:
            # otherwise - C-Type
            self.mu = 0.41
            self.dens = 1500
            self.Y0, self.d0strength   = 1e5,  0.1
            self.nsize = 3  
            self.qconst1, self.qconst2 = 2e3, 4e5 
            self.k1, self.k2=0.15, 1
            
        self.M= (math.pi / 6) * self.dens * self.d**3
        self.jinertia= [2/5 * self.M * (self.d/2)**2]*3
        self.omega=[0,0, (G * 4 / 3 * math.pi * self.dens) ** 0.5* float(parameters[ "omega_in"])]
        self.kvg=0.3
        self.grav=G*self.M/(self.d/2)**2
        self.obliq=float(parameters["obliq"])
        

class Impactor:
    
    def __init__(self,tmaxby,low,high,cumdistr,explicit=True):
        if explicit:
            self.phi = math.acos(1 - 2 * random.random()) / 2
            self.vel = scipy.stats.maxwell.ppf(random.random(), scale=3.232)*1000
            self.d=self.dia(low,high,cumdistr)
            self.theta = 2 * math.pi * random.random()
            self.Theta = 2 * math.pi * random.random()
            self.Phi = math.acos(1 - 2 * random.random())
        else:
            self.d = math.exp((math.log(low) + math.log(high)) / 2)
            self.phi = math.pi / 4
            self.vel = velave
            self.theta = math.pi
            self.Theta = 0
            self.Phi = math.pi / 2
        if random.randint(1, 4) == 1:
            self.dens = 2500
        else:
            self.dens = 1500
        self.impacttime = random.uniform(0, tmaxby)
        self.M= (math.pi / 6) * self.dens * self.d**3
        self.explicit=explicit
    def dia(self,low,high,cumdistr):
        d=np.random.randint(low=low, high=high)
        return getdiaf(d,cumdistr)

#%%

def Nondimensionalize(parameters,w,r):
    
    dia=float(parameters["dia"])
    jinertia=float(parameters["jinertia"])
    jinertia1=float(parameters["jinertia1"])
    parameters["dia"]=dia*r
    parameters["jinertia"]=jinertia/r**5
    parameters["jinertia1"]=jinertia1/r**5
    w[:,1]=w[:,1]/r
    
#%%

def Fit(parameters):
    
    
    gamma=float(parameters["gamma"])
    
    
    try:
        w,file=Shape(parameters)
    except:
        print("cannot import shape")
        return
    
    #uncomment for the special fit
    #res=int(parameters["res"])
    #x=np.sin(w[:,0])*(1+gamma*w[:,1])
    #y=np.cos(w[:,0])*(1+gamma*w[:,1])
     
    #point = []
    #for i in range(res):
     #   point.append([x[i],y[i]])
     #   point.append([-x[i],y[i]])
    #xc, yc, r, sigma = taubinSVD(point)
    
    r=1+gamma*(min(w[:,1])+max(w[:,1]))/2
    w[:,1]=((1+gamma*w[:,1])-r)/gamma
    Nondimensionalize(parameters,w,r)
    Exparameter(parameters)
    np.savetxt(file+"/base.txt",w,delimiter=",")
    print ("The best fit value of r is ", r)

#%%

def Shape(parameters):
     
     omega=float(parameters["omega_in"])
     delta=float(parameters["delta"])
     
     file=os.path.join(os.getcwd(),"output/files_"+str(format(delta,".6f"))+"_"+str(format(omega,".6f")))
     try:
         w=np.loadtxt(file+"/base.txt",dtype=float,delimiter=",")
     except:
         print("No file available for the fit")
         return None, file
     return w,file
 
#%%   Verified

def Parameter(parameters):
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
    return parameters
                
#%%    Verified

def Initialize(parameters,target):
    
    dir = os.getcwd()
    mydir="output/files_"+str(format(target.delta,".6f"))+"_"+str(format(float(parameters["omega_in"]),".6f"))
    mydir=os.path.join(dir,mydir)
    if os.path.exists(mydir):
        subprocess.run(["rm", "-r", mydir])
    else:
        print("No directory exists")
    os.mkdir(mydir)
    parameters['jinertia1'] = target.jinertia[0] / (target.d / 2)**5 / target.dens
    parameters['jinertia'] = target.jinertia[2] / (target.d / 2)**5 / target.dens
    parameters['slides'] = 0
    parameters['time'] = 0
    parameters['omega'] = target.omega[2]/(G * 4 / 3 * math.pi * target.dens) ** 0.5
    parameters['dia']=target.d
    Exparameter(parameters)
    res=int(parameters["res"])
    base=np.zeros((res,2))
    offset=float(parameters["offset"])
    dx=(math.pi-2*offset)/res
    for i in range(res):
        base[i,0] = offset + dx * (i+ 0.5)
    np.savetxt(mydir+"/base.txt",base,delimiter=",")
    Gravitycalc(parameters)

#%%  Verified

def Exparameter(parameters):
        
    with open("parameters","w") as f:
        for key in parameters.keys():
            f.writelines(["-"*50,"\n"])
            f.writelines(f"{key.ljust(25)}{parameters[key]}\n")
        f.writelines(["-"*50])
   
#%%

def Gravitycalc(parameters):
    
    Res=int(parameters["res"])
    epsilon=float(parameters["epsilon"])
    gamma=float(parameters["gamma"])
    density=float(parameters["density"])
    rad=float(parameters["dia"])/2
    
    try:
        w,file=Shape(parameters)
    except TypeError:
        print("No shape avialable to calulate gravity,calculating for the spherical body")
        return
        
    R=rad*np.sin(w[:,0])*(1+gamma*(w[:,1]))
    Z=rad*np.cos(w[:,0])*(1+gamma*(w[:,1]))
    
    fR=make_interp_spline(w[:,0],R)
    fZ=make_interp_spline(w[:,0],Z)

    res=10000

    theta=np.linspace(0,math.pi,res)
    r = fR(theta)
    z = fZ(theta)
    
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    
    arguments=[(R[i]+epsilon*np.sin(w[i,0]),Z[i]+ epsilon*np.cos(w[i,0]),r,z) for i in range(int(Res/2))]
    
    grav = pool.starmap(Gravity, arguments)
    grav=np.array(grav)
    grav1=np.array([(grav[i,0],-grav[i,1]) for i in range(round(Res/2)-1,-1,-1)])
    grav=np.vstack((grav,grav1))
    
    r_grav=-G*density*(grav[:,0]*np.sin(w[:,0])+grav[:,1]*np.cos(w[:,0]))/(4/3*np.pi*density*rad*G)
    t_grav=-G*density*(grav[:,0]*np.cos(w[:,0])-grav[:,1]*np.sin(w[:,0]))/(4/3*np.pi*density*rad*G)
   # plt.plot(w[:,0],r_grav)
   # plt.plot(w[:,0],t_grav)
    print("gravity updated")
    grav=np.hstack((r_grav.reshape(-1,1),t_grav.reshape(-1,1)))
    np.savetxt(file+"/grav.txt",grav)
    
    pool.close()
    pool.join()       
    
#%%


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

        r_grav=r_grav + np.abs((z[i]-z[i-1]))*(2 * delta * ((1 - k**2 / 2) * ks - es) / R)
        z_grav=z_grav + np.abs((z[i]-z[i-1]))*(2 * np.pi * np.sign(zeta) * eps + 2 * zeta * ((R - a)/(R + a) * pi - ks) / delta)
    return r_grav, z_grav    
    
#%% verified

def ExportOmega(myomega,parameters):
    
    omega=float(parameters["omega_in"])
    delta=float(parameters["delta"])
    mydir=os.getcwd()
    file=os.path.join(mydir,"output/files_"+str(format(delta,".6f"))+"_"+str(format(omega,".6f"))+"/omega_T.txt")
    
    resultExport = ""
    with open(file, "w") as file:
        resultExport = file.write("\n".join(["\t".join(map(str, omega)) for omega in myomega]))
    if resultExport == -1:
        print("Failed in exporting the data")

#%%

def Landslides(target,parameters,impacttime,myomega):
    
   # target.d = target.d / (1 + float(parameters["uni_h"]) * float(parameters["epsilon"]))
    parameters["omega"] = target.omega[2] / (G * (4/3) * math.pi * target.dens)**0.5
    parameters["slides"] = int(parameters["slides"])+1
    parameters["dia"] = target.d
    parameters["jinertia"]=target.jinertia[2]/(target.d/2)**5/target.dens
    parameters["jinertia1"]=target.jinertia[0]/(target.d/2)**5/target.dens
    try:
        Exparameter(parameters)
    except IOError:
        print("Failed in exporting the data")
    
    exitcode = subprocess.run(["./gaurav"], cwd=".", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
    if exitcode != 0:
        print("Aborting due to error in running the executable gaurav")
        return
    else:
        print("Ran successfully slide", parameters["slides"])
    
    Parameter(parameters)
    Fit(parameters)
    Parameter(parameters)
    
    target.omega[2] = float(parameters["omega"]) * (G * (4/3) * math.pi * target.dens )**0.5
    target.d = float( parameters["dia"])
    
    r = target.d / 2
    target.jinertia[2] = float(parameters["jinertia"]) * r**5 * target.dens
    target.jinertia[0] = target.jinertia[1] = float(parameters["jinertia1"]) * r**5 * target.dens
    
    Gravitycalc(parameters)
    
    myomega.append([impacttime, target.omega[2]])
    print("Omega after the Landslides", target.omega[2])

#%% Verified

def Yorp(target,parameters,impacttime,oldtime,myomega):
    
    wobblecalcf(target,impacttime,oldtime)  #which omega to use because it is being changed by the yorp
    while impacttime > oldtime + 10:
        subprocess.run(["make"], cwd="/Users/kumargaurav/Documents/OrbFit/tests/gaurav",stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)
        try:
            yark = []
            with open("/Users/kumargaurav/Documents/OrbFit/tests/gaurav/yarkovsky.in", "r") as file:
                yark = [list( line.strip().split()) for line in file]
        except FileNotFoundError:
            print("Failed in importing the data from yarkovsky.in")
            return
        yark[0][5] = target.obliq
        yark[0][6] = 2 * math.pi / (target.omega[2] * 3600)
        try:
            with open("/Users/kumargaurav/Documents/OrbFit/tests/gaurav/yarkovsky.in", "w") as file:
                for row in yark:
                    file.write("\t".join(map(str, row)) + "\n")
        except IOError:
            print("Failed in exporting the yarkovsky.in")
            return
        exitcode = subprocess.run(["./orbit9.x"], cwd="/Users/kumargaurav/Documents/OrbFit/tests/gaurav").returncode
        if exitcode != 0:
            print("Failed in running orbit9")
            return
        try:
            omegOrb = []
            with open("/Users/kumargaurav/Documents/OrbFit/tests/gaurav/clo0.yorp", "r") as file:
                omegOrb = [list( line.strip().split()) for line in file]
        except FileNotFoundError:
            print("Failed in importing the data from orbit9")
            return
        omegaLimit = (G*4/3* math.pi*target.dens)**0.5
        index = min(round((impacttime - oldtime) / 50), 2000)
        target.omega[2] = 2 * math.pi / (float(omegOrb[index][1]) * 3600)
        target.obliq = float(omegOrb[index][2])
        if target.omega[2] > 0.90 * omegaLimit:
            print("Too fast spinning causing landslides")
            Landslides(target,parameters,impacttime,myomega)
        oldtime = oldtime + 10**5
        myomega.append([impacttime, target.omega[2]])
    print("Omega after the yorp effect:", target.omega[2])


#%% Verified

def Cumdistr(parameters):
    mydir=os.getcwd()
    pathcum = os.path.join(mydir, "cumpopulation", f"{parameters['cumdistr']}.txt")
    cumdistr = np.loadtxt(pathcum)
    cumdistr[:,0]=cumdistr[:,0]*1000
    return cumdistr

#%% verified


def Istuff(target,tmaxby,cumdistr):
    dstarave = qstarf(target, math.pi / 4, velave)[2]
    print("dstar", dstarave)
    prob = probi * tmaxby * (target.d/2) ** 2

    numgtd = round(astnum(dstarave,cumdistr)[0])
    dialittle = cumdistr[-1][0]
    dexplicit = max(dialittle, 0.05 * dstarave)
    binexplicit = astnum(dexplicit,cumdistr)[1]
    dexplicit = cumdistr[binexplicit + 1][0]
    nexplicit = round(astnum(dexplicit,cumdistr)[0])

    nexpimpactors = round(prob * nexplicit)
    dimplicit =  0.1 * dexplicit
    binimplicit = astnum(dimplicit,cumdistr)[1]

    istuff=[]
    for j in range(nexpimpactors):
        istuff.append(Impactor(tmaxby,numgtd,nexplicit,cumdistr,True))


    for j in range(binexplicit+1, min(binimplicit,len(cumdistr)-1)):
        istuff.append(Impactor(tmaxby,cumdistr[j,0],cumdistr[j+1,0],cumdistr,False))

    istuff.sort(key=lambda x: x.impacttime)
    return istuff


