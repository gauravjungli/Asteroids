#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 22:04:49 2024

@author: g
"""

import numpy as np
from scipy.special import jv, jvp,lpn  # Bessel function of first kind and its derivative and Legendre polynomial
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from gaurav import Output_File

# Define the equation
def bessel_eq(x, n):
    return 2 * x  * jvp(n + 0.5, x ) - jv(n + 0.5, x )

def dynamic_bounds_search(n, lower, step=1.0, max_iter=100):

    upper = lower + step
    iter_count = 0
    # Increment upper bound until sign changes
    while bessel_eq(lower, n) * bessel_eq(upper, n) > 0:
        lower = upper
        upper += step
        iter_count += 1
        if iter_count > max_iter:
            raise RuntimeError("Exceeded maximum iterations while searching for bounds.")
    return lower, upper


def max_time( t, theta, target ):
    n = target.roots.shape[0]
    m = target.roots.shape[1]
    legendre = lpn(n,theta)[0]
    k_d=np.pi*2*target.f/target.Q
    value = 3/2*k_d
    for i in range(0,n):
        for j in range(0,m):
            B = (2*i+1)/(1-i*(i+1)/(target.d/2*target.roots[i,j])**2)
            value += B*legendre[i]*np.exp(-target.roots[i,j]**2*target.k_s*t)*(k_d + target.roots[i,j]**2*target.k_s)
    print(value)
    return value
       



# Function to find the first n roots of the equation
def find_roots(n, R, num_roots=5, step=1.0):
    roots = []
    lower = max(0.1,np.sqrt(n*(n+1)))  # Start from a small value
    for _ in range(num_roots):
        # Dynamically find the next upper bound where the sign changes
        lower, upper = dynamic_bounds_search(n, lower, step)
        # Use the brentq method to find the root in the interval
        sol = root_scalar(bessel_eq, args=(n), bracket=[lower, upper], method='brentq')
        roots.append(sol.root/R)
        # Set new lower bound for the next root search
        lower = upper
    return np.array(roots)

def  energy(n,m,theta,t,roots,R=500/2,f=10.0,Q=1000,E=1e+6,k_s=0.3e+3):
    k_d=np.pi*2*f/Q
    G = 3/2 
    r = R
    legendre = lpn(n,theta)[0]
    for i in range(0,n+1):
        for j in range(m):
            B = (2*i+1)/(1 - i*(i+1)/(R * roots[i,j])**2 )
            G+=B*legendre[i]*np.exp(-roots[i,j]**2*k_s*t)*np.sqrt(R/r)*jv(i+1/2, r * roots[i,j])/jv(i+1/2, R * roots[i,j])
            
    G = G*E*np.exp(-k_d*t)/(2*np.pi*R**3)
    return G
       
        
#def Destabilize(parameters,target,impactor):
if __name__=="__main__":
    n=300
    m=50
    R=500/2

    t=np.exp(np.linspace(-2,7,1000))
    

    roots = np.zeros((n+1,m))

    # Find roots
    for i in range(0,n+1):
        
        roots[i,:] = find_roots(i, R, num_roots=m)

    # Display the roots
       # print(f"First {m} roots for {n}th order Bessel equation: {roots[i,:]}")

    
    theta = np.cos(np.linspace(np.pi/6, np.pi,20)) 
    G = energy(n=n,m=m,R=R,theta=np.cos(np.pi),t=t,roots=roots)
    #print(G)
    plt.plot(t,G,label='Spherical',linewidth = 2)
    plt.xscale('log')  # Set x-axis to logarithmic scale
    plt.xlabel('Time (t)',fontsize=16)
    plt.ylabel(r'$\epsilon_s(t)$',fontsize=16)
   # plt.title(r'$\epsilon_s(t)$ vs Time (log scale)',fontsize=16)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.legend()
    plt.show()
    