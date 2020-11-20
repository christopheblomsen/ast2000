"""Egen kode."""
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c
import os
import sys
import math

import load_orbit_sim as los
import tools

class orbit_sim:
    """
    A class to numericaly calculate the
    orbits from a system of stars
    if both semi major axes and radii
    is known
    """
    def __init__(self, seed=33382):
        """
        Need to figure this one out
        """
        self.system = SolarSystem(seed)          # Our famous system
        self.axes = self.system.semi_major_axes  # All the semi major axes for formulas
        self.e = self.system.eccentricities      # All the eccentricities for formulas
        self.G = c.G_sol                      # AU**3 yr**-2 SolarMass**-1
        self.M = self.system.star_mass           # Star mass in solar mass
        self.m = self.system.masses[5]           # Home planet mass
        self.omega = self.system.aphelion_angles # Aphelion angles

        self.r_numerical = []                    # List with the results of r from the numerical solution
        self.v_numerical = []                    # List with the results of v from the numerical solution
        self.r_analytical = []                   # List with the results of the analytical solution

        self.mu = self.G*(self.M+self.system.masses)     # Standard gravitational parameter
        self.timesteps_pr_orbit=10000            # Number of timesteps pr orbital period to calculate


    def euler_cromer(self, r0, v0, T, dt):
        """Euler Cromer."""
        G = self.G                               # For less writing
        N = int(T/dt)                            # Length of all our vectors
        t = np.zeros(N, float)                   # time array
        M = self.M                               # Star mass
        r = np.zeros((N, 2), float)              # Position vector
        v = np.zeros((N, 2), float)              # Velocity vector
        a = np.zeros((N, 2), float)              # Acceleration vector

        distance = np.zeros(N, float)            # Distance array

        distance[0] = np.linalg.norm(r0)         # Sets initial conditions
        r[0, :] = r0
        v[0, :] = v0
        t[0] = 0

        for i in range(N-1):
            a[i, :] = -G*M/(distance[0]**3) * r[i, :]
            v[i + 1, :] = v[i, :] + a[i, :]*dt
            r[i + 1, :] = r[i, :] + v[i + 1, :]*dt
            t[i + 1] = t[i] + dt

        return r, v, a, t

    def r(self,t):
        """Return the distance r for all planets at timestep t from the simulation."""
        r_vec = np.array(self.r_numerical)[:,t]
        return np.sqrt(r_vec[:,0]**2+r_vec[:,1]**2)

    def v(self,t):
        """Return the velocity for all planets at timestep t from the simulation."""
        v_vec = np.array(self.v_numerical)[:,t]
        return np.sqrt(v_vec[:,0]**2+v_vec[:,1]**2)

    def simulated_aphelion(self):
        """Find maximum distance from star for all planets."""
        r_vec = np.array(self.r_numerical)
        r = np.sqrt(r_vec[:,:,0]**2+r_vec[:,:,1]**2)
        return r[:,np.argmax(r,axis=1)]


    def mean_velocity(self):
        """Return the mean velocity of the planets."""
        v_vec = np.array(self.v_numerical)
        # Calculate velocity at all times
        velocities = np.sqrt(v_vec[:,:,0]**2+v_vec[:,:,1]**2)
        # Return the mean velocity for each planet
        return np.mean(velocities,axis=1)

    def dArea(self,t):
        """Calulate the area swept out at timestep t and over time self.dt."""
        return 0.5*self.r(t)[0]*(self.v(t)*self.dt)

    def area_swiped_during_period(self,t,days):
        """Calculate area swiped during a number of days."""
        dtPrDay = self.timesteps_pr_orbit/(self.hoth_period()*c.yr/c.day)
        n = int(dtPrDay*days)
        totalArea = 0
        for i in range(n):
            totalArea += self.dArea(t+i)

        return totalArea

    def distance_traveled_during_period(self,t,days):
        """Numerically calculate the distance traveled during # of days."""
        dtPrDay = self.timesteps_pr_orbit/(self.hoth_period()*c.yr/c.day)
        n = int(dtPrDay*days)
        totalDistance = 0
        for i in range(n):
            totalDistance += self.v(t)*self.dt

        return totalDistance

    def hoth_period(self):
        """Return the rotational periode of our planet in years."""
        return 2*np.pi*np.sqrt(self.axes[0]**3/self.mu[0])

    def hoth_perhelion_t(self):
        return int(self.hoth_period()/self.dt/2)

    def sim(self):
        """Simulate all the orbits."""
        print('Simulating orbits')

        mu = self.mu                                                # Standard gravitational parameter

        f = lambda r : -self.G*self.M/(np.linalg.norm(r)**3)

        N = len(self.system.masses)                                 # Length for for loop

        rotational_orbit_in_years = 2*np.pi*np.sqrt(self.axes**3/mu)      # Rotational orbit in years
        print(f'Orbit in years {rotational_orbit_in_years}')
        T = 45*rotational_orbit_in_years[0]                            # Total time of N rotations
        self.T = T
        self.dt = rotational_orbit_in_years[0]/self.timesteps_pr_orbit # Find dt based on timesteps pr year
        planet_pos = self.system.initial_positions                  # Initial planets positions
        planet_vel = self.system.initial_velocities                 # Initial planets velocities
        verification_r = np.zeros((2, N, int(T/self.dt)), float)
        #verification_t = []
        for i in range(N):
            print(f'Working on planet {i+1}')
            print(f'Orbit in years {2*np.pi*np.sqrt(self.axes[i]**3/mu[i])}')
            m = self.system.masses[i]                                   # Gets i'th planet mass

            r0 = planet_pos[:, i]                              # Gets i'th planet starting pos

            v0 = planet_vel[:, i]                              # --||--                    velocity
            r, v, a, t = tools.leapfrog(r0, v0, T, self.dt, f)    # runs leapfrog and returns
            self.r_numerical.append(r)
            self.v_numerical.append(v)
        self.analytical_solution()                      # analytical

    def center_mass(self, m, r):
        if(len(m) != len(r)):
            raise Exception("different number of bodies and positions")

        M = np.sum(m)
        bodies = len(m)
        dimensions = r.shape[1]
        print(f'Calculating CM for {bodies} bodies in {dimensions} D')
        R = np.zeros(dimensions)

        for i in range(bodies):
            R += m[i]*r[i]

        R = R/M

        return R
    def to_center_mass_positions(self,r,m1,m2):
        """
        Move bodies to center of mass coordinate system.

        r = distance vector between bodies
        m1 = mass of body 1
        m2 = mass of body 2

        Returns positions of the two bodies relative to the center of mass
        """
        mu = (m1*m2)/(m1+m2)
        r1cm = r*mu/m1
        r2cm = -r*mu/m2
        return r1cm, r2cm

    def solar_orbit(self, planet):
        """Simulate planet and star orbiting a common center of mass."""
        star_initial_pos = np.array([0, 0])
        star_initial_vel = np.array([0, 0])

        planet_pos = self.system.initial_positions[:, planet]
        planet_vel = self.system.initial_velocities[:, planet]

        M = self.M
        m = self.system.masses[planet]
        mu = (M * m)/(M + m)

        orbital_period = 2*np.pi*np.sqrt(self.axes[planet]**3/mu)      # One year
        dt = orbital_period/100000

        masses = np.array([m, M])

        r = np.array([planet_pos,star_initial_pos])

        center_of_mass_R = self.center_mass(masses, r)

        # Chaneg reference to center of mass
        r1cm, r2cm = self.to_center_mass_positions(r[0]-r[1],masses[0],masses[1])
        r1,r2,t = self.leapfrog2([r1cm,r2cm],[planet_vel,[0,0]],m,M,10,0.0001)

        plt.plot(r1[:,0],r1[:,1])
        plt.axis('equal')
        plt.plot(r2[:,0],r2[:,1])
        plt.axis('equal')
        plt.plot(0,0,"or")
        plt.show()

    def leapfrog2(self, r0, v0, m1, m2, T, dt):
        """Leapfrog integration for 2 body orbit. """
        G = self.G                               # For less writing
        N = int(T/dt)                            # Length of all our vectors
        t = np.zeros(N, float)                   # time array
        M = self.M                               # Star mass
        r1 = np.zeros((N, 2), float)              # Position vector
        r2 = np.zeros((N, 2), float)              # Position vector
        v1 = np.zeros((N, 2), float)              # Velocity vector
        v2 = np.zeros((N, 2), float)              # Velocity vector
        a1 = np.zeros((N, 2), float)              # Acceleration vector
        a2 = np.zeros((N, 2), float)              # Acceleration vector
        mu = (m1*m2)/(m1+m2)

        r = np.zeros((N,2),float)
        r[0] = r0[0]+r0[1]
        ar = math.sqrt(r0[0][0]**2+r0[0][1]**2)+math.sqrt(r0[1][0]**2+r0[1][1]**2)

        r1[0, :] = np.array(r0[0])
        r2[0, :] = np.array(r0[1])
        v1[0, :] = np.array(v0[0])
        v2[0, :] = np.array(v0[1])
        t[0] = 0

        a1[0, :] = -G*m2/(ar**2) * r[0]
        a2[0, :] = G*m1/(ar**2) * r[0]
        for i in range(N-1):
            """
            The actual leapfrog algorithm
            """
            r1[i + 1, :] = r1[i, :] + v1[i, :]*dt + 0.5*a1[i, :]*dt**2
            r2[i + 1, :] = r2[i, :] + v2[i, :]*dt + 0.5*a2[i, :]*dt**2
            r[i + 1] = r1[i+1,:]+r2[i+1,:]

            a1[i + 1, :] = -G*m2/(ar**2) * r[i + 1]
            a2[i + 1, :] = G*m1/(ar**2) * r[i + 1]
            v1[i + 1, :] = v1[i, :] + 0.5*(a1[i, :] + a1[i + 1, :])*dt
            v2[i + 1, :] = v2[i, :] + 0.5*(a2[i, :] + a2[i + 1, :])*dt
            t[i + 1] = t[i] + dt

        return r1, r2, t
    def cartesian_polar(self, r):
        """Convert to polar coordinates."""
        x = r[:, 0]                             # x values
        y = r[:, 1]                             # y values
        r_p = np.zeros(len(r))
        for i in range(len(r)):
            r_p[i] = np.linalg.norm(r[i])
        # r = np.sqrt(x**2 + y**2)              # distance from origo
        theta = np.arctan(y/x)                  # theta
        return r_p, theta

    def polar_cartesian(self, r, theta):
        """Convert to cartesian."""
        N = len(theta)
        x = np.zeros(N, float)
        y = np.zeros(N, float)
        for i in range(N):
            x[i] = r[i]*np.cos(theta[i])        # calculates x
            y[i] = r[i]*np.sin(theta[i])        # calculates y
        return x, y

    def analytical_solution(self):
        theta = np.linspace(0, 2*np.pi, 1000)   # array from 0 to 2pi

        def r(axes, e, theta, omega):
            """Analytical formula."""
            ans = (axes*(1 - e**2))/(1 + e*np.cos(theta - omega))
            return ans

        for i in range(len(self.axes)):
            x, y = self.polar_cartesian(r(self.axes[i], self.e[i], theta, self.omega[i]), theta)
            self.r_analytical.append([x, y])

    def plot(self):
        i = 0
        print(np.array(self.r_numerical).shape)
        for r in self.r_numerical:
            plt.plot(r[:, 0], r[:, 1],label=f'Planet {i+1} s')
            plt.axis('equal')
            i += 1

        i = 0
        print(np.array(self.r_analytical).shape)
        print(len(self.axes))
        for a in self.r_analytical:
            plt.plot(a[0],a[1],label=f'Planet {i+1} a')
            i += 1

        plt.xlabel('x in AU')
        plt.ylabel('y in AU')
        plt.legend(loc='lower right')
        plt.title('Hoth system')

        plt.show()

    def plot_planet(self,i):
        r = self.r_numerical[i]
        plt.plot(r[:, 0], r[:, 1],label=f'Planet {i} s')
        plt.axis('equal')
        a = self.r_analytical[i]
        plt.plot(a[0],a[1],label=f'Planet {i} a')
        i += 1

        plt.xlabel('x in AU')
        plt.ylabel('y in AU')
        plt.legend(loc='lower right')
        plt.title('Hoth system')

        plt.show()

    def least_squares(self, v0, v_end, v_r, P0, P_end, P, t0, t_end):
        """Least squares method for finding the most likely candidate."""
        v = np.linspace(v0, v_end, 1000)
        P = np.linspace(P0, P_end, 1000)
        t = np.linspace(t0, t_end, 1000)

        def f(v, P, t):
            return (v - v_r*np.cos(2*np.pi/P)*(t - t0))

        delta_min = f(v0, P0, t0)

        P_min = 0
        V_min = 0
        t_min = 0

        for i in range(len(v)):
            for j in range(len(P)):
                for k in range(len(t)):
                    delta = f(v[i], P[j], t[k])
                    if delta < delta_min:
                        delta_min = delta
                        if t[k] < t_min:
                            t_min = t[k]
                        if P[j] < P_min:
                            P_min = P[j]
                        if v[i] < V_min:
                            V_min = v[i]

        return delta_min, P_min, V_min, t_min



    def verify_planet_positions(self):
        planet_positions = np.moveaxis(np.array(self.r_numerical),[0,1,2],[1,2,0])
        self.system.verify_planet_positions(self.T,planet_positions,'planet_trajectories.npz')
        self.system.generate_system_snapshot('system_snapshot.xml')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-rs','--run-sim',action='store_true', help='Runs the simulation of orbits from scratch and saves the result')
    parser.add_argument('-d','--download',action='store_true', help='Downloads pickle file')
    args = parser.parse_args()

    filename = "simulated_orbits.pkl"
    orbit = los.orbit_sim_factory(filename,args)
    planet = 5
    m = np.array([orbit.system.masses[planet],orbit.system.star_mass])
    r = np.array([orbit.system.initial_positions[:,planet],np.array([0,0])])
    orbit.verify_planet_positions()
    print(f'Masses {m}')
    print(f'Positions {r}')
    R = orbit.center_mass(m,r)
    r1cm, r2cm = orbit.to_center_mass_positions(r[0]-r[1],m[0],m[1])
    print(f'CM positions {orbit.to_center_mass_positions(r[0]-r[1],m[0],m[1])}')
    plt.plot(r1cm[0],r1cm[1],"ob")
    plt.plot(r2cm[0],r2cm[1],"oy")
    plt.plot(0,0,".r")
    plt.show()

    orbit.solar_orbit(planet)
