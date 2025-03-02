#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 08:52:43 2024

@author: g
"""
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def compute_power_series_coefficients(n, terms=10):
    """
    Compute the coefficients of the power series using the given recurrence relation.

    Parameters:
        n (float): Parameter `n` in the recurrence relation.
        terms (int): The number of terms in the series to compute.

    Returns:
        list: Coefficients of the power series.
    """
    # Initialize the coefficients
    a = [1]  # a_1 = 1 (initial coefficient)
    
    # Compute coefficients a_s for s = 3, 5, ..., up to the specified terms
    for s in range(1, 2 * terms, 2):
        if s + 2 >= 2 * terms:
            break
        next_a = -((n - s) * (n + s + 1)) / ((s + 2) *(s+1)) * a[-1]
        a.append(next_a)
    b = []
    for s in range(1, 2 * terms, 2):
            if s + 2 >= 2 * terms:
                break
            next_b = -((n - s) * (n + s + 1)) / ((s + 2)) * a[s//2]
            b.append(next_b)
    return b

def compute_u_coefficients(n, terms=10):
    """
    Compute the coefficients of the power series using the given recurrence relation.

    Parameters:
        n (float): Parameter `n` in the recurrence relation.
        terms (int): The number of terms in the series to compute.

    Returns:
        list: Coefficients of the power series.
    """
    # Initialize the coefficients
    a = [1]  # a_1 = 1 (initial coefficient)
    
    # Compute coefficients a_s for s = 3, 5, ..., up to the specified terms
    for s in range(1, 2 * terms, 2):
        if s + 2 >= 2 * terms:
            break
        next_a = -((n - s) * (n + s + 1)) / ((s + 2) *(s+1)) * a[-1]
        a.append(next_a)

    return a

def power_series(x, n, terms=10):
    """
    Evaluate the power series at a given x and n.

    Parameters:
        x (float): The point where the series is evaluated.
        n (float): Parameter `n` in the recurrence relation.
        terms (int): The number of terms in the series to compute.

    Returns:
        float: The value of the series at x.
    """
    coefficients = compute_power_series_coefficients(n, terms)
    return sum(coeff * (x ** (2 * i+1)) for i, coeff in enumerate(coefficients))


def u_series(x, n, terms=10):
    """
    Evaluate the power series at a given x and n.

    Parameters:
        x (float): The point where the series is evaluated.
        n (float): Parameter `n` in the recurrence relation.
        terms (int): The number of terms in the series to compute.

    Returns:
        float: The value of the series at x.
    """
    coefficients = compute_u_coefficients(n, terms)
    return sum(coeff * (x ** (2 * i)) for i, coeff in enumerate(coefficients))


def find_n_for_zero_series(x=1, terms=10, guess_range=(-10, 10), num_guesses=100):
    """
    Find the values of n for which the power series goes to zero at x = 1.

    Parameters:
        terms (int): The number of terms in the series to compute.
        guess_range (tuple): Range of initial guesses for n.
        num_guesses (int): Number of guesses to try within the range.

    Returns:
        list: Values of n for which the power series goes to zero at x = 1.
    """
    # Define the function for fsolve
    def series_func(n):
        #print(n, power_series_at_x1(n,terms))
        result = u_series(x,n, terms)
        if abs(result)<1e-9:
            print(n,result)
        return result
    
    # Generate initial guesses
    guesses = np.linspace(guess_range[0], guess_range[1], num_guesses)
    
    # Use fsolve to find roots
    roots = []
    for guess in guesses:
        root = fsolve(series_func, guess)[0]
        # Add the root if it's not already in the list (to handle duplicates)
        if not any(np.isclose(root, r, atol=1e-2) for r in roots):
            roots.append(root)
    
    return roots


def plot_power_series(n_values, terms=10, x_range=(-2, 2), x_points=100):
    """
    Plot the power series for the given values of n.

    Parameters:
        n_values (list): Values of n for which to plot the series.
        terms (int): Number of terms to use in the series.
        x_range (tuple): Range of x values for the plot.
        x_points (int): Number of points in the x range.
    """
    x = np.linspace(x_range[0], x_range[1], x_points)
    plt.figure(figsize=(10, 6))
    
    for n in n_values:
        y = [power_series(xi, n, terms) for xi in x]
        plt.plot(x, y, label=f'n = {n:.3f}')
    
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
    plt.title("Power Series for Roots of n")
    plt.xlabel("x")
    plt.ylabel("Power Series Value")
    plt.legend()
    plt.grid()
    plt.show()
    
    

def plot_u_series(n_values, terms=10, x_range=(-2, 2), x_points=100):
    """
    Plot the power series for the given values of n.

    Parameters:
        n_values (list): Values of n for which to plot the series.
        terms (int): Number of terms to use in the series.
        x_range (tuple): Range of x values for the plot.
        x_points (int): Number of points in the x range.
    """
    x = np.linspace(x_range[0], x_range[1], x_points)
    plt.figure(figsize=(10, 6))
    
    for n in n_values:
        y = [u_series(xi, n, terms) for xi in x]
        plt.plot(x, y, label=f'n = {n:.3f}')
    
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
    plt.title("Power Series for Roots of n")
    plt.xlabel("x")
    plt.ylabel("Power Series Value")
    plt.legend()
    plt.grid()
    plt.show()
    
    

# Example Usage
terms = 1000  # Number of terms in the series to compute
x = 0.99999999
roots = find_n_for_zero_series(x,terms, guess_range=(1,3), num_guesses=1)
print(f"Values of n for which the power series is zero at x = 1: {roots}")

#roots = np.linspace(3,3.01,10)
# Plot the power series for these values of n
plot_u_series(roots, terms, x_range=(0.99, 1.0001), x_points=100)