#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 20:46:26 2024

@author: g
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0              # Length of the radial domain (maximum radius)
Nx = 100   
R= 250          # Number of spatial points
rho = 1200* R**2           # Density (assumed constant)
sigma_rr = 1.0       # Radial stress (assumed constant for simplicity)
dr = L / (Nx - 1)    # Radial step size

# Construct the matrix A
A = np.zeros((Nx-2, Nx-2))

# Fill in the finite difference coefficients for the radial part of the equation
for i in range(1, Nx-1):
    r_i = i * dr
    
    if i > 1:
        A[i-1, i-2] = sigma_rr / (dr**2) - sigma_rr / (r_i * 2 * dr)
        
    A[i-1, i-1] = -2 * sigma_rr / (dr**2)

    if i < Nx-2:
        A[i-1, i] = sigma_rr / (dr**2) + sigma_rr / (r_i * 2 * dr)

# Solve the eigenvalue problem
eigenvalues, eigenvectors = np.linalg.eigh(A)

# Calculate natural frequencies (omega) from eigenvalues
frequencies = np.sqrt(np.abs(eigenvalues) / rho)

# Plot the natural frequencies
plt.plot(range(1, len(frequencies) + 1), frequencies, 'bo-')
plt.xlabel('Mode Number')
plt.ylabel('Frequency (rad/s)')
plt.title('Natural Frequencies of the Radial System')
plt.grid(True)
plt.show()

# Plot the first few mode shapes
plt.figure()
for i in range(3):  # Plot the first 3 mode shapes
    plt.plot(np.linspace(0, L, Nx-2), eigenvectors[:, i], label=f'Mode {i+1}')

plt.xlabel('Radius (r)')
plt.ylabel('Displacement (u_r)')
plt.title('Mode Shapes')
plt.legend()
plt.grid(True)
plt.show()
