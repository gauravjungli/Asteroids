#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:14:37 2023

@author: Kumar Gaurav, ChatGpt, Holsapple
"""

import math

from scipy.interpolate import make_interp_spline
from scipy.constants import gravitational_constant
from IO import Output_File
import random
G=gravitational_constant
probi = 2.85e-24

import numpy as np

#%%
""" 
    Is the class of impactors. Its attributes are as follows:
        explict: IS the flag to check whether an impact is big enough to cause landslides
        impacttime: is the time of impact of each impact
        M: is the masss of each impactor
        dia: is the function that generates the random number which represents the number of asteroid greater
            than a specific dia. This dia is returned by the function getdiaf
        d: is the diameter of the impactor
        theta, Theta, Phi: are the angles of the impact.
 
"""

class Impactor:
    
    def __init__(self,tmaxby,low,high,cumdistr,file,explicit=True):
        if explicit:
            self.phi   = math.acos(1 - 2 * np.random.random()) / 2
            self.vel   = self.velocity(file) 
            self.d     = self.dia(low,high,cumdistr)  
            self.theta = 2 * math.pi * np.random.random()
            self.Phi = 2 * math.pi * np.random.random()
            self.Theta = np.clip(math.acos(1 - 2 * np.random.random()),np.pi/6,5*np.pi/6)
        else:
            self.d     = math.exp((math.log(low) + math.log(high)) / 2)
            self.phi   = math.pi / 4
            self.vel   = self.velocity(file) 
            self.theta = math.pi
            self.Phi = 0
            self.Theta   = math.pi / 2
        if np.random.randint(1, 4) == 1:
            self.dens  = 2500
        else:
            self.dens  = 1500
        self.impacttime= np.random.uniform(0, tmaxby)
        self.M         = (math.pi / 6) * self.dens * self.d**3
        self.explicit  = explicit
        
    def dia(self,low,high,cumdistr):
        d=np.random.randint(low=low, high=high)
        return getdiaf(d,cumdistr)
    
    def velocity(self,file):
        ranges = []
        probabilities = []
        
        with open(file, 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) == 3:
                    lower, upper, prob = map(float, parts)
                    ranges.append((lower, upper))
                    probabilities.append(prob)

# Step 2: Choose a range based on the distribution
        selected_range = random.choices(ranges, weights=probabilities, k=1)[0]

# Step 3: Sample a value uniformly from the selected range
        sample = random.uniform(selected_range[0], selected_range[1])
        
        return sample*1e+3





#%%   
""" It gives the number of asteroid greater than a specific dia. It also returns the bin in which this specific 
    dia belongs to. 
    1. f: is the spine fit to the data given for the population density.  
    2. interp: finds the number of asteroid using the spline f.
    3. Bin: Contains the bin in which the asteroid with the given dia resides.
    """
    
def astnum(dia,cumdistr):
    f=make_interp_spline(np.flip(cumdistr[:,0]),np.flip(cumdistr[:,1]))
    interp=math.ceil(f(dia))
    Bin = next((pos for pos, val in enumerate(cumdistr) if val[0] < dia), None)
    return [interp, Bin - 1]

#%%
""" Main function that implements spin change due to collisions also obliquity change added. 
    zeta: Is the efficiency of the of angular momentum transfer
    delamomentum: Is the angular momentum transffered to the body
    delomega: Is the change in the angular velocity
    myomega: Stores the omega values at all time instant. It is used in post processing.
"""

def Collision(target,impactor,myomega):
    
    if  not target.collision:
        return
    print(f" obliquity before impact is {target.obliq} and angular velocity is {target.omega}")
    zeta = zetaf(impactor.phi, target.d, target.atype)
    amomentum =  np.multiply(target.omega,target.jinertia)
    delamomentum = target.d/2 * impactor.M * impactor.vel * math.sin(impactor.phi) * zeta * np.array([-math.sin(impactor.theta)
                            * math.cos(impactor.Phi) * math.cos(impactor.Theta) - math.cos(impactor.theta) * math.sin(impactor.Phi), 
                            -math.sin(impactor.theta) * math.sin(impactor.Phi) * math.cos(impactor.Theta) + math.cos(impactor.theta) *
                            math.cos(impactor.Phi),math.sin(impactor.theta) * math.sin(impactor.Theta)])
    delomega = np.divide(delamomentum , target.jinertia)  
    post_amomentum = np.add(amomentum,delamomentum)
    ran_angle = 2*math.pi*np.random.rand()
    normal = [ np.sin(target.obliq/180*np.pi)*np.cos(ran_angle),
              np.sin(target.obliq/180*np.pi)*np.sin(ran_angle),np.cos(target.obliq/180*np.pi)]
    target.obliq = np.arccos(np.dot(normal,post_amomentum)/np.linalg.norm(post_amomentum))*180/np.pi

    delomegdrain =  omegdrainf(target,impactor)
    delomega = delomega + [0,0,delomegdrain]
    target.omega = target.omega + delomega
    print(f"New obliquity after impact is {target.obliq} and angular velocity is {target.omega}")
    myomega.append([impactor.impacttime,target.omega[2]])
    print("Omega after the collision", target.omega[2])
        
#%%
""" 
    Gives a diameter based on the number distribution of the asteroid. It selects first the bin in which the 
    the random number of asteroids resides. Then it finds the diameter by linear interpolation usig the log-log
    scale
    pos: This is the bin in which the required diameter exists. 
"""
def getdiaf(num,cumdistr):
    pos = next((i for i, val in enumerate(cumdistr[:, 1]) if val > num),None)
    low = cumdistr[pos]
    high = cumdistr[pos - 1]
    slope = math.log(high[1] / low[1]) / math.log(high[0] / low[0])
    di = high[0] * (num / high[1]) ** (1 / slope)
    return di

#%% 
""" Not currently relevant. Calculates angular momentum drained to ejecta."""

def omegdrainf(target,impactor):

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

#%% 
""" It calculates the values for catastrophic disruption. This uses the formula from the Holsapple paper"""

def qstarf(target, phi, vel):
    slope = 3 * target.mu / (1 - 2 *target.nsize)
    qstars = min((target.qconst1 * (target.d/2 ) ** slope), target.qconst1)
    qstarg = target.qconst2 * (target.d/2 / 5e5) ** (3 * target.mu)
    qstar = ((qstars  + qstarg ) * 
             (math.cos(phi) / math.cos(math.radians(45))) ** (-3 * target.mu) * 
             (vel / 5.5e+3) ** (2 - 3 * target.mu))
    
    massstar = 2 * qstar * target.M / vel ** 2
    dstar = ((6 / math.pi) * massstar / target.dens) ** (1 / 3)
    return [qstar, massstar, dstar]

#%% 
"""For calculating wobble. It is based on just adding wobble to the existing wobble. It also exist wobble decay."""

def wobblecalcf(target,impacttime,time):
    
    if target.omega[2]==0:
        print("No rotation")
        return
    
    wobble = np.arctan(math.sqrt(target.omega[0]**2 + target.omega[1]**2)/abs(target.omega[2]))
   
    if wobble < math.radians(1):
        target.omega[0] = math.sqrt(target.omega[0]**2 + target.omega[1]**2) * math.cos(np.random.uniform(0, 2 * math.pi))
        target.omega[1] = math.sqrt(target.omega[0]**2 + target.omega[1]**2) * math.sin(np.random.uniform(0, 2 * math.pi))
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


#%%  
""" Function for calcualting zeta. It is the efficiency factor for the spinup or it can be said to be the efficiency of angular momentum transferred"""

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
 

#%% verified
""" This creates the collision history. Using Poisson's distribution for determining the number of impactors and deducing th
    the collision time """

def Istuff(parameters,target,tmaxby,cumdistr):
    """
    Creates a list of impactors. 
    probi: is the intrinsic collisional probability
    tmaxby: the simulation time period
    explicit_cutoff: is the maximum limit on diameter that we considers for impact
    numgtd: The number of impactors greater than the explicit limit
    dialittle: Minimum size of impactors that we can consider
    velave: is the average velocity of the impactor
    energy_min: The mininum energy of the impactor that can cause global failure
    dexplicit: the minimum diameter for causing global landslide
    
    ----------

    Returns
    -------
    istuff : TYPE
        DESCRIPTION.

    """
    prob = probi * tmaxby * (target.d/2) ** 2
    explicit_cutoff = float(parameters['Explicit cutoff'])
    numgtd = round(astnum(explicit_cutoff*target.dstarave,cumdistr)[0])
    dialittle = cumdistr[-1][0]
    velave = float(parameters['Impactor velocity'])
    
    #Find minimum energy for the global failure
    energy_cons =  (math.pi* target.efficiency*velave**2*target.dens/12)*target.energy[-1]
    energy_min = target.cohesion_cons**2/(2*target.wave_speed**2*target.dens*np.tan(target.delta*math.pi/180)**2)
    
    #Find minimum explicit diameter
    dexplicit =  (energy_min/energy_cons)**(1/3)
    #Check whether the dexplicit lies inside the permissible range
    dexplicit = max(1.1*dialittle, dexplicit)
    
    #Find the bin in which the dexplicit lies
    binexplicit = astnum(dexplicit,cumdistr)[1]
    #Find number of impactors
    nexplicit = round(astnum(dexplicit,cumdistr)[0])
    
    #FInd expected number of impactors
    nexpimpactors = round(prob * (nexplicit-numgtd))
    #Number of impactors using Poisson's distribution
    nexpimpactors = poisson_from_exponential(nexpimpactors)
    print(f'Number of expected impactor is {nexpimpactors}')
    density = 1500
    
    #Not in use currently
    dimplicit = ((G**2*target.dens**3*target.d**5)/(9*target.efficiency*density*velave**2*target.f**2))**(1/3)*(
                np.exp(2*math.pi*target.f*target.d**2/(target.k_s*math.pi**2*target.Q)))
    binimplicit = astnum(dimplicit,cumdistr)[1]
    
    #File containing velocity distribution
    vel_dist = parameters['Velocity file']
    file =Output_File(parameters,'input',[vel_dist])

    #Create collisional history for the explicit impactors
    istuff=[]
    #  Added only to simulate manual collisional history
    #nexpimpactors = 20
    for j in range(nexpimpactors):
        istuff.append(Impactor(tmaxby=tmaxby, low=numgtd, high=nexplicit, cumdistr=cumdistr,file=file, explicit=True))

    #Create collisional history for the implicit impactors
    #implicit impactors are not required and are hence removed from further simulations
   # for j in range(binexplicit+1, min(binimplicit,len(cumdistr)-1)):
   #     istuff.append(Impactor(tmaxby=tmaxby, low=cumdistr[j,0], high=cumdistr[j+1,0], cumdistr=cumdistr,file=file, explicit=False))
    
    #Sort with the impact time
    istuff.sort(key=lambda x: x.impacttime)
    print(f"The minimum diameter for creating landslides is {dexplicit}")
    return istuff


def poisson_from_exponential(lambda_rate, max_time=1):
    """
    Simulate a Poisson-distributed random variable using exponential inter-arrival times.
    
    Parameters:
    - lambda_rate: The rate (λ) of the Poisson process.
    - max_time: The maximum time interval fGammaor the simulation (default is 1).
    
    Returns:
    - Number of events (Poisson-distributed) in the interval [0, max_time].
    """
    num_events = 0
    cumulative_time = 0
    if (lambda_rate == 0):
        print("No expected impact")
        return 0
    # Generate exponential random variables until the cumulative time exceeds max_time
    while cumulative_time <= max_time:
        # Generate the time between the next event (Exponential random variable)
        inter_arrival_time = np.random.exponential(1 / lambda_rate)
        
        # Update cumulative time
        cumulative_time += inter_arrival_time
        
        # If the event happens within the time interval, increase event count
        if cumulative_time <= max_time:
            num_events += 1

    return num_events
