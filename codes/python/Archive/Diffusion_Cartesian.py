import numpy as np
import matplotlib.pyplot as plt

# Define constants (replace these with your actual values)
f = 1  # Frequency
Q = 1000  # Quality factor
Ks = 0.3e+3  # Constant Ks
R = 250
L = R*(4/3*np.pi)**(1/3)  # Constant L
E=1e+6
# Define time variable
t=np.exp(np.linspace(-2,9,1000))  # You can define an array for different time values if needed

# Define the number of terms for the summation
N = 1000

# Spatial coordinates
x0, y0, z0 = 0,0,0  # initial positions
x, y, z = 2*R/np.sqrt(3),2*R/np.sqrt(3),2*R/np.sqrt(3)     # current positions

# Function to compute the epsilon_s(t)
def epsilon_s():

    
    # Exponential factor
    decay_factor = np.exp(-2 * np.pi * f * t / Q)
    
    # Summation in x
    sum_x = sum(
        2 * np.cos(n * np.pi * x0 / L) * np.cos(n * np.pi * x / L) * np.exp(-Ks * n**2 * np.pi**2 * t / L**2)
        for n in range(1, N+1)
    )
    
    # Summation in y
    sum_y = sum(
        2 * np.cos(n * np.pi * y0 / L) * np.cos(n * np.pi * y / L) * np.exp(-Ks * n**2 *np.pi**2 * t / L**2)
        for n in range(1, N+1)
    )
    
    # Summation in z
    sum_z = sum(
        2 * np.cos(n * np.pi * z0 / L) * np.cos(n * np.pi * z / L) * np.exp(-Ks * n**2 * np.pi**2 * t / L**2)
        for n in range(1, N+1)
    )
    
    # Final result for epsilon_s
    epsilon_s = decay_factor * (1 + sum_x) * (1 + sum_y) * (1 + sum_z)*E/L**3

    
    return epsilon_s

# Evaluate epsilon_s for a given time
epsilon_value = epsilon_s()
plt.plot(t, epsilon_value, '-.', label=r'Cubical at $2R$', linewidth = 2 )
plt.legend()
