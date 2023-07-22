#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:14:37 2023

@author: kumargaurav, ChatGpt, Holsapple
"""

import math
import random
from scipy.interpolate import make_interp_spline
from scipy.constants import gravitational_constant
import numpy as np


G=gravitational_constant
velave=5.5e3
probi = 2.85e-24

#%%   Verified
    
def astnum(dia,cumdistr):
    f=make_interp_spline(np.flip(cumdistr[:,0]),np.flip(cumdistr[:,1]))
    interp=math.ceil(f(dia))
    Bin=len(cumdistr)
    Bin = next((pos for pos, val in enumerate(cumdistr) if val[0] < dia), None)
    return [interp, Bin - 1]

#%%  verified

def Collision(target,impactor):
    
    zeta = zetaf(impactor.phi, target.d, target.atype)
    delamomentum = target.d/2 * impactor.M * impactor.vel * math.sin(impactor.phi) * zeta * np.array([-math.sin(impactor.theta)
                            * math.cos(impactor.Theta) * math.cos(impactor.Phi) - math.cos(impactor.theta) * math.sin(impactor.Theta), 
                            -math.sin(impactor.theta) * math.sin(impactor.Theta) * math.cos(impactor.Phi) + math.cos(impactor.theta) *
                            math.cos(impactor.Theta),math.sin(impactor.theta) * math.sin(impactor.Phi)])
    delomega = np.divide(delamomentum , target.jinertia)  
    qstar = qstarf(target, impactor.phi, impactor.vel)[0]
    qq = (impactor.M / 2) * impactor.vel ** 2 / target.M
    qratio = qq / qstar
    delomegdrain =  omegdrainf(target,impactor) if qratio < 0.5 else [0, 0, 0]  # 3 components, parallel to omeg negative in function
    delomega = delomega + [0,0,delomegdrain]
    target.omega = target.omega + delomega
    
    print("Omega after the collision", target.omega[2])
        
#%%  Verified

def getdiaf(num,cumdistr):
    pos = next(i for i, val in enumerate(cumdistr[:, 1]) if val > num)
    low = cumdistr[pos]
    high = cumdistr[pos - 1]
    slope = math.log(high[1] / low[1]) / math.log(high[0] / low[0])
    di = high[0] * (num / high[1]) ** (1 / slope)
    return di

#%% verified

def omegdrainf(target,impactor): #need to find Y upon distribution of size

    pi2=target.grav*impactor.d/2/impactor.vel**2
    pi3=target.Y0/(target.dens*impactor.vel**2)
    vstar=target.kvg*impactor.vel*(pi2)**(1/(2+target.mu))
    rc = vstar / math.sqrt((8/3) * math.pi * target.dens * G)
    piv = target.k1 * (pi2 * (target.dens / impactor.dens) ** (-1/3) + (target.k2 * pi3) ** ((2 + target.mu) / 2)) ** (-(3 * target.mu) / (2 + target.mu))
    massc = piv * impactor.M
    masseject = 0.6 * massc
    ra = rc * (masseject / impactor.M) ** (1 / (3 * target.mu))
    delomegdrain = - (5/12) * (3 * target.mu) * target.omega[2] * (impactor.M / target.M) * (target.d/2 / ra) ** (-3 * target.mu)
    drain = delomegdrain if target.d/2 > rc else 0
    return drain

#%% Verified

def qstarf(target, phi, vel):
    slope = 3 * target.mu / (1 - 2 *target.nsize)
    qstars = min((target.qconst1 * (target.d/2 ) ** slope), target.qconst1)
    qstarg = target.qconst2 * (target.d/2 / 5e5) ** (3 * target.mu)
    qstar = ((qstars  + qstarg ) * 
             (math.cos(phi) / math.cos(math.radians(45))) ** (-3 * target.mu) * 
             (vel / velave) ** (2 - 3 * target.mu))
    
    massstar = 2 * qstar * target.M / vel ** 2
    dstar = ((6 / math.pi) * massstar / target.dens) ** (1 / 3)
    return [qstar, massstar, dstar]

#%%  verified

def wobblecalcf(target,impacttime,time):
    
    if target.omega[2]==0:
        print("No rotation")
        return
    
    wobble = np.arctan(math.sqrt(target.omega[0]**2 + target.omega[1]**2)/abs(target.omega[2]))
   
    if wobble < math.radians(1):
        target.omega[0] = math.sqrt(target.omega[0]**2 + target.omega[1]**2) * math.cos(random.uniform(0, 2 * math.pi))
        target.omega[1] = math.sqrt(target.omega[0]**2 + target.omega[1]**2) * math.sin(random.uniform(0, 2 * math.pi))
        print("No wobble")
        return
    
    trelax =math.log(wobble / math.radians(1.0)) / 2 * 1.1 * 1e-3 / (target.d/1000) ** 2 / np.linalg.norm(target.omega) ** 3
    print("relaxation time",trelax)
    print("Time between impacts",impacttime-time)
    if (trelax + time < impacttime):  
        target.omega = [0, 0, np.sign(target.omega[2]) * np.sqrt(np.sum(np.multiply(np.square(target.jinertia),np.square(target.omega))))/target.jinertia[2]]
        return

    tsince = impacttime - time
    gamma = (math.radians(1.0) / wobble) ** (tsince / trelax)

    if gamma > 1 or gamma < 0:
        print("Oops! Error in wobblecalcf. {timesince, trelax, wobble, gamma, omeg} =",tsince, trelax, wobble, gamma, target.omega)
        return

    alpha2 = math.sqrt(np.sum(np.multiply(np.square(target.jinertia),np.square(target.omega))))/ math.sqrt(
            target.omega[2]**2 * (target.jinertia[2]**2 + ((target.jinertia[0]**2 * target.omega[0]**2 + target.jinertia[1]**2 * target.omega[1]**2)
                                               * math.tan(gamma * wobble)**2) / (target.omega[0]**2 + target.omega[1]**2)))

    alpha1 = math.tan(gamma * wobble) /math.tan(wobble)* alpha2

    if alpha2.imag > 0:
        print("Wobble wrong omeg =", target.omega)
        return

    target.omega = [alpha1*target.omega[0],alpha1*target.omega[1],alpha2*target.omega[2]]


#%%  Verified

def zetaf(phi, d, name):
    
    def multiptf(x):
        
        lst=[[1e-1, 4], [10, 4], [1e3, 1], [2e3, 1], [1e6, 1.0]]
        if x>=lst[-1][0]:
            return lst[-1][1]
        elif x<=lst[0][0]:
            return lst[0][1]
        else:
            Bin=next(i for i, val in enumerate(lst) if val[0] > x)
            exp = math.log(lst[Bin][1] / lst[Bin-1][1]) / math.log(lst[Bin][0] / lst[Bin-1][0])
            return lst[Bin-1][1] * (x / lst[Bin-1][0]) ** exp

    
    if name == "S-Type":
        zeta = 0.41 * math.cos(phi)**2
    else:
        zeta = 0.8 * (1 - 2 * phi / math.pi)
    
    return multiptf(d)*zeta
 
