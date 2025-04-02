#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 06:53:52 2024

@author: g
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d
from scipy.linalg import eig
# Define the system of ODEs



def odes(X, Y,xi,B,R):
    if len(Y) != 2:
        raise ValueError("Input vector y must have length 2.")
    v, y = Y
    r=float(X)
    # Adding a small value to avoid division by zero
    epsilon = 1e-12
    if y+4/3*v<0:
        y=-4/3*v
        dy_dx=0
        dv_dx =  ((16*v**(3/2)*np.sqrt(6)*B)/(21*r)-alpha*r)/((2*np.sqrt(v)*np.sqrt(6)*B)/7)/(alpha/B/(6*np.sqrt(3)))
        print("Yield reached",r)
        return [dv_dx,dy_dx]
        
    denominator = B*((xi + 5/12)*v**2 + y*(xi + 7/6)*v + ((xi + 5/3)*y**2)/4)
    Fr=alpha*r/(alpha/B/(6*np.sqrt(3)))
    c = y +2*v
    if c <epsilon:
        print("Printing in ode solver",c,r)
        c = epsilon
    if abs(denominator.any())<epsilon:
        print("Denomintor 2 going to zero")
 #   if abs(denominator_1) < epsilon or abs(denominator_2) < epsilon:
 #       return [0, 0]  # Return zero derivatives to avoid divide by zero
    dv_dx = y/r - v/r
    #print(c)
    dy_dx = -(xi*(c/2)*Fr*r*np.sqrt(c)/3  -
              2*B*((xi + 3/4)*v**2 + y*(xi + 1)*v + y**2*(xi + 1)/4)*(-y + v))/denominator/r
    
        
    return [dv_dx, dy_dx]

# Initial conditions

# Define a function to compute the residual at y1(1) for a given guess of y1(0)
def shooting_residual(C,A,B,R,X_span,X_eval):
    # Initial conditions: y0(0) = 1 and y1(0) = y1_0_guess
    y_initial = [C*epsilon,C]
    
    # Solve the IVP from x=0 to x=1
    solution = solve_ivp(odes,X_span, y_initial, args=(xi,B,R), t_eval=X_eval,method='RK45')
    
    # Extract y1(1) from the solution
   # y1_at_end = (A+2*B)*solution.y[1,-1] + 2*A*solution.y[0,-1]
    u_at_end =solution.y[0,-1]
    y_at_end =solution.y[1,-1]
    # Compute the residual: we want y1(1) to be 0
    print(abs(u_at_end),C)
    return abs(u_at_end)+abs(y_at_end)-1e-4


xi=49/12
B=8.5e+9
R=250
rho=1200
Nx = 500
epsilon=1.9e-2
guess =  1e-6
alpha = R**2 * 4/3 * 6.67e-11 * (rho**2)*np.pi
# Time range for the solution
x_span = (guess,1-epsilon)  # Range of X values
x_eval = np.linspace(x_span[0], x_span[1], Nx)  # Points at which to store the solution
C =0.9337292
try:
    # Use a root-finding algorithm to adjust y1(0) until y1(1) = 0
    
   # result = root_scalar(shooting_residual, bracket=[1e-10,1e-6], args=(xi,B,R,X_span,X_eval), method='brentq')


  #  if not result.converged:
  #      print("The shooting method failed to converge.")
  #  else:
   #     print("The shooting method converged.")
        

     #   C = result.root
        #u_guess = -(np.sqrt(12)*alpha/B)**(2/3)*((-x_span[0]**(10/3) + 1))**(2/3)/(24*x_span[0]**(8/9))
        #u_guess = (28*np.sqrt(6)*alpha/B)**(2/3)*(x_span[0]**2 *(1-x_span[0]**2) )**(2/3)/16
        y_initial = [C,C]
        solution = solve_ivp(odes, x_span, y_initial, args=(xi,B,R), t_eval=x_eval, method='Radau',rtol=1e-12,atol=1e-12)
        
    # Extract the solution
  #  if solution.success:
        x = np.array(solution.t*R)
        v = solution.y[0]*(alpha/B/(6*np.sqrt(3)))**(2/3)
        y = solution.y[1]*(alpha/B/(6*np.sqrt(3)))**(2/3)
        alpha /= R**2
        x_end =  np.linspace(R*(1-epsilon),R*(1-1e-12), int(Nx/5))
        v_end = 28**(2/3)*(np.sqrt(6)*alpha*x_end**2*(R**2 - x_end**2)*B**2*R)**(2/3)/(16*B**2*R**2)
        y_end = -4/3*v_end
        x_tot =  np.concatenate((x,x_end))
        v_tot =  np.concatenate((v,v_end))
        y_tot =  np.concatenate((y,y_end))
        stress1 =  ((xi + 5/3)*y_tot**2 + 4*(xi + 1/6)*v_tot*y_tot + 4*(xi - 7/12)*v_tot**2)*B/(np.sqrt(y_tot + 2*v_tot)*xi)
        stress2 =  ((xi - 1/3)*y_tot**2 + 4*(xi - 1/3)*v_tot*y_tot + 4*(xi + 5/12)*v_tot**2)*B/(np.sqrt(y_tot + 2*v_tot)*xi)
        total_stress2 = (np.trapz(stress2*x_tot,x_tot) - 3*(R**2*(epsilon**2
                                             - 2*epsilon + 3/2)*(-2+epsilon)**2*epsilon**2*alpha)/12)/R**4
        total_stress = alpha/8
        print(total_stress2/total_stress)
        print(stress1[Nx-1]-alpha*(1/(2*(x[-1])**2)-1/(2*R**2))*(x[-1])**4)
        
        K =B*np.sqrt(y_tot+2*v_tot)*(1+1/3*(y_tot-v_tot)**2/(y_tot+2*v_tot)**2/xi)
        mu = B/xi*np.sqrt(y_tot+2*v_tot)
        speed1 = np.sqrt((K+4/3*mu)/rho)
        speed2 = np.sqrt(mu/rho) 
        plt.plot(x_tot,speed1,label =r'P-wave')
        plt.plot(x_tot,speed2,label = r'S-wave')
        plt.ylabel(r'Speed',fontsize=14)
        plt.xlabel("Radius",fontsize=14)
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        plt.grid()  
        plt.legend(fontsize=14)
        plt.show()
        plt.plot(x_tot,v_tot,label =r'$\epsilon_{\theta\theta}$')
        plt.plot(x_tot,y_tot,label = r'$\epsilon_{rr}$')
        plt.ylabel(r'Strain',fontsize=14)
        plt.xlabel("Radius",fontsize=14)
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        plt.grid()  
        plt.legend(fontsize=14)
        plt.show()
        time1 =  np.trapz(1/speed1,x_tot)
        avg_speed1 = R/time1
        time2 =  np.trapz(1/speed2,x_tot)
        avg_speed2 = R/time2  
        coeff_1 = np.polyfit(x_tot,stress1,2)
        poly_s1 = np.poly1d(coeff_1)


    # Evaluate the polynomial at specific points (e.g., for a smooth curve)
        x_fit = np.linspace(min(x), max(x), 100)
        s1_fit = poly_s1(x_fit)

    # Print the polynomial equation
        print(f"Fitted polynomial: {poly_s1}")
        # Plot the results
        #plt.figure(figsize=(12, 6))

        stress = 2/3*6.67e-11*(1200**2)*np.pi*(R**2-x_tot**2)

        plt.clf()
        plt.plot(x_tot, stress1, label=r'$\sigma_{rr}$')
        plt.plot(x_tot, stress2, label=r'$\sigma_{\theta\theta}$')
        plt.plot(x_tot, stress, label=r'$\sigma^{h}$')
        plt.ylabel(r'Stress',fontsize=14)
        plt.xlabel("Radius",fontsize=14)
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        plt.grid()  
        plt.legend(fontsize=14)

        plt.tight_layout()
        plt.show()
 #   else:
        print("The solver was not successful. Reason:", solution.message)
except ValueError as e:
    print("Error:", e)
 
    #%%


A = np.zeros((Nx-2, Nx-2))
M = np.zeros((Nx-2, Nx-2))

r =  x_eval
dr = abs(x_span[0]-x_span[1])/(Nx-1)
stress1/=B#sigma_mag
stress2/=B#sigma_mag
srp_y,srm_y,srp_v,srm_v,st_v = 0,0,0,0,0
for i in range(1, Nx-1):

    M[i-1,i-1] = ((r[i+1]+r[i-1])/2)**2
    srp_y = (stress1[i+1]-stress1[i])/(y[i+1]-y[i])
    srm_y = (stress1[i]-stress1[i-1])/(y[i]-y[i-1])
    srp_v = (stress1[i+1]-stress1[i])/(v[i+1]-v[i])
    srm_v = (stress1[i]-stress1[i-1])/(v[i]-v[i-1]) 
    st_v = (stress2[i+1]-stress2[i-1])/(v[i+1]-v[i-1])
    
    if i>1:
        A[i-1, i-2] = ((1/dr)**2 * srm_y *(r[i-1]**2 + r[i]**2)/2  
          +  ((1/dr)*(srp_v* (r[i+1]+r[i])/2 - srm_v* (r[i-1]+r[i])/2 ) - st_v)/2)
        
    A[i-1, i-1] = -((1/dr)**2 * srp_y *(r[i+1]**2 + r[i]**2)/2 + (1/dr)**2 * srm_y *(r[i-1]**2 + r[i]**2)/2  )

    if i<Nx-2:    
        A[i-1, i] =((1/dr)**2 * srp_y *(r[i+1]**2 + r[i]**2)/2 +
                    +  ((1/dr)*(srp_v* (r[i+1]+r[i])/2 - srm_v* (r[i-1]+r[i])/2 ) - st_v)/2)


A[-1, -1] = A[-1,-1] + ((1/dr)**2 * srp_y *(r[Nx-1]**2 + r[Nx-2]**2)/2
                        +  ((1/dr)*(srp_v* (r[Nx-1]+r[Nx-2])/2 - srm_v* (r[Nx-1]+r[Nx-2])/2 ) - st_v)/2)/(1+2/3*dr/r[Nx-1])*(1-2/3*dr/r[Nx-2])
# Solve the eigenvalue problem



eigenvalues, eigenvectors = eig(A,M)
neg_eigenvalues = eigenvalues[eigenvalues<0]
neg_eigenvectors = eigenvectors[:,eigenvalues<0]
frequencies = np.sqrt(np.abs(neg_eigenvalues) / rho / R**2 *B )
sorted_indices = np.argsort(frequencies)
freq_sorted = frequencies[sorted_indices]
vector_sorted = neg_eigenvectors[:,sorted_indices]
# Calculate natural frequencies (omega) from eigenvalues


#Plot the natural frequencies
plt.figure(figsize=(6, 6))
plt.plot(range(1,int( 11)), freq_sorted[0:10], 'bo-')
plt.xlabel('Mode Number',fontsize=16)
plt.ylabel('Frequency (rad/s)',fontsize=16)
#plt.title('Natural Frequencies of the Radial System (with Dirichlet BC at r=0 and Neumann BC at r=L)')
plt.grid(True)
plt.show()
print(min(frequencies),max(frequencies))
# Plot the first few mode shapes

L = len(freq_sorted)
plt.figure(figsize=(6, 6))
for i in range(0,5):  # Plot the first 3 mode shapes
    plt.clf()
    
    plt.plot(r[1:-1], vector_sorted[:, i], label=f'Mode {i+1}')

    plt.xlabel('Radius (r)',fontsize=16)
    plt.ylabel('Displacement (u_r)',fontsize=16)
    #plt.title('Mode Shapes (with Dirichlet BC at r=0 and Neumann BC at r=L)')
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.grid()  
    plt.legend(fontsize=16)
    plt.show()
    plt.pause(0.5)

#%%