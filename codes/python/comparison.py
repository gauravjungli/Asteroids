import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import root_scalar, fsolve

# --- Core Mathematical Functions (from Scheeres 2015) ---

def r_sl(delta, omega):
    """Sea-level radius (Eq. 20)"""
    return np.sqrt((1 - omega**2) / (1 - omega**2 * np.cos(delta)**2))

def r_phi(delta, delta_1, omega, phi):
    """Radius profile for a constant slope angle phi (Eq. 45)"""
    return np.sqrt((1 - omega**2 * np.cos(delta_1)**2) / 
                   (1 - omega**2 * np.cos(delta)**2)) * np.exp((delta - delta_1) * np.tan(phi))

def r_H(delta, delta_0, omega):
    """Radius of zero-slope redistributed regolith (Eq. 47 & 50)"""
    H_plus_1 = np.sqrt((1 - omega**2 * np.cos(delta_0)**2) / (1 - omega**2))
    return H_plus_1 * np.sqrt((1 - omega**2) / (1 - omega**2 * np.cos(delta)**2))

def initial_failure_limits(omega, phi):
    """Finds where the apparent slope first exceeds the friction angle (Eq. 44)"""
    term1 = omega**2 / (2 * np.sqrt(1 - omega**2))
    if np.tan(phi) > term1:
        return None, None # No failure yet
    
    term2 = np.sqrt(omega**4 / (4 * (1 - omega**2)) - np.tan(phi)**2)
    tan_d1 = (term1 - term2) / np.tan(phi)
    tan_d2_approx = (term1 + term2) / np.tan(phi)
    
    return np.arctan(tan_d1), np.arctan(tan_d2_approx)

def roche_lobe(delta, omega):
    """Calculate Roche Lobe radius given latitude and spin rate (Eq. 30)"""
    # Solve (1/r)^3 - 1.5 * omega^(2/3) * (1/r)^2 + 0.5 * omega^2 * cos(delta)^2 = 0
    # Let u = 1/r
    r_rl = np.zeros_like(delta)
    for i, d in enumerate(delta):
        coeff = [1, -1.5 * omega**(2/3), 0, 0.5 * omega**2 * np.cos(d)**2]
        roots = np.roots(coeff)
        real_roots = np.real(roots[np.isreal(roots)])
        # u must be > 0 (so r > 0). Smallest positive real root gives largest r (exterior Roche Lobe)
        pos_roots = real_roots[real_roots > 0]
        if len(pos_roots) > 0:
            r_rl[i] = 1.0 / np.min(pos_roots)
        else:
            r_rl[i] = np.nan
    return r_rl

def limiting_radius(delta, omega):
    """Limiting radius where net outwards acceleration is zero, a_x = 0 (Eq. 13 limit)"""
    r_x = 1.0 / (omega**(2/3))
    return np.full_like(delta, r_x)

# --- Solvers for Failure Regimes ---

def solve_profile(omega, phi):
    """Determine the boundaries (delta_0, delta_1, delta_2) and return the radius profile function"""
    d1_init, _ = initial_failure_limits(omega, phi)
    
    if d1_init is None:
        return lambda d: np.ones_like(d), "Stable" # No failure
    
    # Function to find delta_2 given delta_1: r_phi(delta_2, delta_1) = 1
    def get_delta2(d1):
        res = root_scalar(lambda d2: r_phi(d2, d1, omega, phi) - 1.0, 
                          bracket=[d1, np.pi/2 - 1e-5], method='brentq')
        return res.root

    # Volume integral difference calculation
    def delta_V12(d1, d2):
        integrand = lambda d: (1.0 - r_phi(d, d1, omega, phi)**3) * np.cos(d)
        val, _ = quad(integrand, d1, d2)
        return val

    # Test Local Failure Scenario
    d2 = get_delta2(d1_init)
    dV12 = delta_V12(d1_init, d2)
    sin3_d0 = ((1 - omega**2) / omega**2) * dV12
    
    if sin3_d0 < 0:
        d0 = 0.0
    else:
        d0 = np.arcsin(sin3_d0**(1/3))

    if d0 <= d1_init:
        # LOCAL FAILURE Valid
        d1 = d1_init
        regime = "Local"
    else:
        # REGIONAL FAILURE (delta_0 = delta_1)
        regime = "Regional"
        def regional_residual(d1):
            d2_curr = get_delta2(d1)
            dV = delta_V12(d1, d2_curr)
            return ((1 - omega**2) / omega**2) * dV - np.sin(d1)**3
        
        # Upper bound: delta where tangent condition vanishes or reaches limits
        res = root_scalar(regional_residual, bracket=[1e-5, np.pi/2 - 0.01], method='brentq')
        d1 = res.root
        d0 = d1
        d2 = get_delta2(d1)
        
        # Global limit check (if delta_2 reaches pi/2)
        if d2 > np.pi/2 - 1e-4:
            d2 = np.pi/2
            regime = "Global"

    def profile(delta):
        delta = np.atleast_1d(delta)
        r = np.ones_like(delta)
        for i, d in enumerate(delta):
            if d <= d0:
                r[i] = r_H(d, d0, omega)
            elif d > d0 and d <= d1:
                r[i] = 1.0
            elif d > d1 and d <= d2:
                r[i] = r_phi(d, d1, omega, phi)
            else:
                r[i] = 1.0
        return r

    return profile, regime

# --- Plotting Scripts ---

def plot_fig9():
    """Reproduces Figure 9: Shape Profiles at specific spin rates for phi = 30, 35 deg."""
    plt.figure(figsize=(12, 5))
    
    cases = [
        (30, [0.816, 0.856, 0.923], ["Spin rate=0.82 (Local)", "Regional", "Global"]),
        (35, [0.854, 0.887, 0.954], ["Spin rate=0.85 (Local)", "Regional", "Global"])
    ]
    
    delta_range = np.linspace(0, np.pi/2, 500)
    
    for idx, (phi_deg, omegas, labels) in enumerate(cases):
        plt.subplot(1, 2, idx+1)
        phi = np.radians(phi_deg)
        
        for omega, label in zip(omegas, labels):
            profile_func, regime = solve_profile(omega, phi)
            r_vals = profile_func(delta_range)
            x = r_vals * np.cos(delta_range)
            z = r_vals * np.sin(delta_range)
            plt.plot(x, z, label=f"{label} ({regime})")
            
            # Plot the equilibrium point (r*)
            r_star = 1.0 / (omega**(2/3))
            plt.scatter([r_star], [0], zorder=5)

        # Baseline sphere
        plt.plot(np.cos(delta_range), np.sin(delta_range), 'k:', alpha=0.3, label="Radius=1")
        
        plt.title(f"Angle of Friction = {phi_deg}°")
        plt.xlabel("Equatorial Axis")
        plt.ylabel("Polar Axis")
        plt.legend(loc='upper right', fontsize='small')
        plt.xlim(0, 1.4)
        plt.ylim(0, 1.05)
        plt.gca().set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.show()

def plot_fig10():
    """Reproduces Figure 10: Shape Profiles for different phi at specific spin rates."""
    plt.figure(figsize=(12, 5))
    
    omegas = [0.90, 0.95]
    phis = [30, 35, 40, 45, 50]
    
    delta_range = np.linspace(0, np.pi/2, 500)
    
    for idx, omega in enumerate(omegas):
        plt.subplot(1, 2, idx+1)
        
        for phi_deg in phis:
            phi = np.radians(phi_deg)
            profile_func, regime = solve_profile(omega, phi)
            
            if regime == "Stable":
                continue # Skip if no failure occurs at this spin rate
                
            r_vals = profile_func(delta_range)
            x = r_vals * np.cos(delta_range)
            z = r_vals * np.sin(delta_range)
            plt.plot(x, z, label=f"Friction angle = {phi_deg}°")

        # Roche Lobe
        r_rl = roche_lobe(delta_range, omega)
        plt.plot(r_rl * np.cos(delta_range), r_rl * np.sin(delta_range), 
                 'c--', alpha=0.6, label="Roche Lobe")
        
        # Limiting radius a_x = 0
        r_lim = limiting_radius(delta_range, omega)
        # Only plot up to where it makes sense (usually a small arc near equator)
        mask = (np.sin(delta_range) < 0.3)
        plt.plot(r_lim[mask] * np.cos(delta_range[mask]), r_lim[mask] * np.sin(delta_range[mask]), 
                 'k-', linewidth=2, label="a_x=0 Limiting Radius")

        plt.title(f"Spin rate $\\bar{{\\omega}}$ = {omega}")
        plt.xlabel("Equatorial Axis")
        if idx == 0: plt.ylabel("Polar Axis")
        plt.xlim(0, 1.4)
        plt.ylim(0, 1.05)
        plt.legend(loc='upper right', fontsize='x-small')
        plt.gca().set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    plot_fig9()
    plot_fig10()