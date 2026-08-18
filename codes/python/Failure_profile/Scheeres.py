import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import root_scalar
import pdb
from collisions import G
import math
import os

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
    tan_phi = np.tan(phi)
    discriminant = omega**4 - 4 * (1 - omega**2) * tan_phi**2
    
    if discriminant < 0:
        return None, None 
        
    sqrt_disc = np.sqrt(discriminant)
    tan_d1 = (omega**2 - sqrt_disc) / (2 * tan_phi)
    tan_d2 = (omega**2 + sqrt_disc) / (2 * tan_phi)
    
    # d1 is the lower latitude (closer to equator), d2 is the upper latitude
    return np.arctan(tan_d1), np.arctan(tan_d2)

def get_delta2(d1,d2_init,omega,phi):
    """Find where the failure profile comes back up to meet the original sphere (r=1)"""
    upper_bound = np.pi/2 - 1e-5
    
    # Check if the curve reaches the pole before getting back up to r=1
    r_pole = r_phi(upper_bound, d1, omega, phi)
    if r_pole <= 1.0:
        return np.pi/2
    
    # The excavated radius dips below 1 and reaches its absolute minimum exactly at d2_init.
    # We start searching for the exit point from the minimum onwards to avoid numerical errors.
    if d2_init is None or d1 >= d2_init:
        return d1
        
    lower_bound = max(d1 + 1e-5, d2_init)
    
    try:
        res = root_scalar(lambda d2: r_phi(d2, d1, omega, phi) - 1.0, 
                          bracket=[lower_bound, upper_bound], method='brentq')
        return res.root
    except ValueError:
        return d1

def Scheeres(target,base):
    
    omega =target.omega[2]/(G * (4/3) * math.pi * target.dens)**0.5
    phi = math.radians(target.delta)
    res = target.res
    theta = base[0:int(res/2),0]
    """Determine the boundaries and return the radius profile function"""
    d1_init, d2_init = initial_failure_limits(omega, phi)
    
    if d1_init is None or d2_init is None:
        return np.zeros(res)
        
    d2_local = get_delta2(d1_init,d2_init,omega,phi)
    
    height = r_phi(theta, d1_init, omega, phi)
    
    out_of_range_mask = (theta < d1_init) | (theta > d2_local)

    height[out_of_range_mask] = 1
    
    height[:] = 1 - height[:] 
    
    height_rev = height[::-1]
    
    height = np.concatenate((height_rev,height))
    return height*target.d/2
    
    

def solve_profile(omega, phi):

    """Determine the boundaries and return the radius profile function"""
    d1_init, d2_init = initial_failure_limits(omega, phi)
    
    if d1_init is None:
        return lambda d: np.ones_like(d), "Stable"
    
    

    def delta_V12(d1, d2):
        """Calculate the excavated volume difference"""
        integrand = lambda d: (1.0 - r_phi(d, d1, omega, phi)**3) * np.cos(d)
        val, _ = quad(integrand, d1, d2)
        return val

    # Test Local Failure Scenario boundaries
    d1_local = d1_init
    d2_local = get_delta2(d1_local,d2_init,omega,phi)
    dV12 = delta_V12(d1_local, d2_local)
    
    sin3_d0 = ((1 - omega**2) / omega**2) * dV12
    d0_local = np.arcsin(sin3_d0**(1/3)) if sin3_d0 > 0 else 0.0

    if d0_local <= d1_local:
        # LOCAL FAILURE Valid: The equatorial bulge does not overlap the excavated region
        d0, d1, d2 = d0_local, d1_local, d2_local
        regime = "Local"
    else:
        # REGIONAL FAILURE: Bulge covers initial failure point. (delta_0 == delta_1)
        regime = "Regional"
        
        def regional_residual(d1_test):
            d2_test = get_delta2(d1_test,d2_init,omega,phi)
            dV = delta_V12(d1_test, d2_test)
            return ((1 - omega**2) / omega**2) * dV - np.sin(d1_test)**3
        
        try:
            # The intersection point must lie between the initial failure point and the pole
            res = root_scalar(regional_residual, bracket=[d1_init, np.pi/2 - 1e-5], method='brentq')
            d1 = res.root
        except ValueError:
            # Absolute fallback if limits are pushed completely
            d1 = d1_init 
            
        d0 = d1
        d2 = get_delta2(d1,d2_init,omega,phi)
        
        # If the excavation reaches the pole, transition to Global failure
        if d2 >= np.pi/2 - 1e-4:
            d2 = np.pi/2
            regime = "Global"
    print (d1,d2)
    def profile(delta):
        """Piecewise evaluation of the surface profile"""
        delta = np.atleast_1d(delta)
        r = np.ones_like(delta)
        for i, d in enumerate(delta):
            if d <= d0:
                # Equatorial Bulge (Zero-slope redistributed material)
                r[i] = r_H(d, d0, omega)
            elif d > d0 and d <= d1:
                # Unfailed original sphere separating bulge and excavation (Only in Local failure)
                r[i] = 1.0
            elif d > d1 and d <= d2:
                # Excavated constant-slope surface
                r[i] = r_phi(d, d1, omega, phi)
            else:
                # Unfailed original sphere near the pole
                r[i] = 1.0
        return r

    return profile, regime

def plot_fig9():
    """Reproduces Figure 9: Shape Profiles at specific spin rates for phi = 30, 35 deg."""
    plt.figure(figsize=(12, 5))
    
    # Cases formatted as: (Friction Angle, [Spin Rates], [Legend Labels])
    cases = [
        (30, [0.82, 0.856, 0.923], ["Spin rate=0.82", "Regional", "Global"]),
        (35, [0.86, 0.887, 0.954], ["Spin rate=0.86", "Regional", "Global"])
    ]
    
    delta_range = np.linspace(0, np.pi/2, 500)
    
    for idx, (phi_deg, omegas, labels) in enumerate(cases):
        plt.subplot(1, 2, idx+1)
        phi = np.radians(phi_deg)
        
        for omega, label in zip(omegas, labels):
            profile_func, regime = solve_profile(omega, phi)
            r_vals = profile_func(delta_range)
            
            # Map Polar coordinates to Cartesian (X = Equator, Z = Pole)
            x = r_vals * np.cos(delta_range)
            z = r_vals * np.sin(delta_range)
            
            plt.plot(x, z, label=label)
            
            # Plot the equilibrium point (Roche limit intersecting equator)
            r_star = 1.0 / (omega**(2/3))
            plt.scatter([r_star], [0], zorder=5)

        # Baseline undisturbed sphere
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

def compare(target):
    """Reproduces Figure 9: Shape Profiles at specific spin rates for phi = 30, 35 deg."""
    plt.figure(figsize=(12, 5))
    
    # Cases formatted as: (Friction Angle, [Spin Rates], [Legend Labels])
    cases = [(30, [ target.omega[2]]/(G * (4/3) * math.pi * target.dens)**0.5, [ "Scheeres 2015"])]
    
    delta_range = np.linspace(0, np.pi/2, 500)
    
    for idx, (phi_deg, omegas, labels) in enumerate(cases):
        phi = np.radians(phi_deg)
        
        for omega, label in zip(omegas, labels):
            profile_func, regime = solve_profile(omega, phi)
            r_vals = profile_func(delta_range)
            
            # Map Polar coordinates to Cartesian (X = Equator, Z = Pole)
            x = r_vals * np.cos(delta_range)
            z = r_vals * np.sin(delta_range)
            
            plt.plot(np.cos(delta_range), np.sin(delta_range), 'k:', alpha=0.3, label="Radius=1")
            plt.plot(x, z, label=label)
            

        # Baseline undisturbed sphere
        
        
        plt.title(f"Angle of Friction = {phi_deg}°")
        plt.xlabel("Equatorial Axis")
        plt.ylabel("Polar Axis")
        plt.legend(loc='upper right', fontsize='small')
        plt.xlim(0, 1.2)
        plt.ylim(0, 1.05)
        plt.gca().set_aspect('equal', adjustable='box')
        


        file1 = target.folder
        Res = target.res
        file2 = os.path.join(file1,'dia.txt')
        epsilon = np.loadtxt(file2,dtype=float,ndmin=2)[:,2]

        file = os.path.join(file1,'field_1.csv')

        if os.path.exists(file):
           
            w = np.loadtxt(file,delimiter=",",dtype='float')

            res = int(Res/2)
            theta = w[0:res,0]
            base = w[0:res,1]
            height = epsilon[0]*w[0:res,2]
            dbase = w[0:res,3]
        
            metric = base**2 + dbase**2
             
            rad = np.sqrt((2*base**2*np.sqrt(metric)*height + metric*(height**2 + base**2  ))/metric)

            z = np.cos(theta)*base + 1/np.sqrt(metric)*(height*dbase*np.sin(theta)+np.cos(theta)*base*height)
            theta_new = np.arccos(z/rad)
            
            x = rad*np.sin(theta_new)
            y = rad*np.cos(theta_new)

            plt.plot(x,y,color='red',label='Simulation')

            plt.axis('equal')



        else:
            break

    plt.tight_layout()
    plt.show()


def compare_new(target):
    """Reproduces Figure 9: Shape Profiles at specific spin rates for phi = 30, 35 deg."""
    plt.figure(figsize=(12, 5))
    
    # Cases formatted as: (Friction Angle, [Spin Rates], [Legend Labels])
    cases = [(30, [ target.omega[2]/(G * (4/3) * math.pi * target.dens)**0.5], [ "Scheeres 2015"])]
    
    delta_range = np.linspace(0, np.pi/2, 500)
    
    for idx, (phi_deg, omegas, labels) in enumerate(cases):
        phi = np.radians(phi_deg)
        
        for omega, label in zip(omegas, labels):
            profile_func, regime = solve_profile(omega, phi)
            r_vals = profile_func(delta_range)
            
            # Map Polar coordinates to Cartesian (X = Equator, Z = Pole)

            z =  np.where((1 - r_vals )*250>0,(1 - r_vals )*250,0 )
        
            plt.plot(delta_range, z, label=label)
            

        # Baseline undisturbed sphere
 
        


        file1 = "/home/g/Asteroids/output/" + target.folder +"/run1/data"
        Res = target.res
        file2 = os.path.join(file1,'dia.txt')

        file = os.path.join(file1,'field_1.csv')

        if os.path.exists(file):
           
            w = np.loadtxt(file,delimiter=",",dtype='float')

            res = int(Res/2)
            theta = w[0:res,0]
            base = 1
            dbase = 0
            
            rgrav = -1
            tgrav = 0
            rho =target.dens
            metric = np.sqrt(base**2 + dbase**2)
            J =  metric*base*np.sin(theta)
            kappa_phi_g = 1/J*(dbase*np.sin(theta)+base*np.cos(theta))
            kappa_phi_n = -1/J*(dbase*np.cos(theta)-base*np.sin(theta))
        
            normal_p = -rho*(rgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_n)
            tangential_p = rho*(tgrav + (target.omega[2]*base*np.sin(theta))**2*kappa_phi_g)
            
            c_0 = 0
            c_1 = 0
            
            mu = np.tan(phi)
            
            sign =np.sign(tangential_p)
           
            k = sign*tangential_p - mu*normal_p
            
            print(np.max(k))
            def mohr_coulomb(h,k):
                if k<0:
                    return 0
                return k*h - c_0*(np.exp(c_1*h)-1)
            
            
            height = np.ones(res)*1e-4
            
            for i in range(len(k)):
                    
                
                if k[i]<c_0*c_1:
                    
                    continue
                try:
                    sol = root_scalar(mohr_coulomb,bracket=[1e-4,50],args=(k[i]))
                    height[i] = max( sol.root, height[i])
                except: 
                    print(f"Could not find height for i = {i}")
                    
            plt.plot(theta,height,color='red',label='Simulation')



        else:
            break

    plt.tight_layout()
    plt.show()

#if __name__ == '__main__':
    #compare(parameters)