# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c
import os
import sys

import argparse
import load_orbit_sim as los


parser = argparse.ArgumentParser()
parser.add_argument('-rs','--run-sim',action='store_true', help='Runs the simulation of orbits from scratch and saves the result')
parser.add_argument('-d','--download',action='store_true', help='Downloads pickle file')
args = parser.parse_args()

class orbit_sim:
    '''
    A class to numericaly calculate the
    orbits from a system of stars
    if both semi major axes and radii
    is known
    '''
    def __init__(self, seed=33382):
        '''
        Need to figure this one out
        '''
        self.system = SolarSystem(seed)          # Our famous system
        self.axes = self.system.semi_major_axes  # All the semi major axes for formulas
        self.e = self.system.eccentricities      # All the eccentricities for formulas
        self.G = c.G_sol                      # AU**3 yr**-2 SolarMass**-1
        self.M = self.system.star_mass           # Star mass in solar mass

        self.r_numerical = []                    # List with the results of r from the numerical solution
        self.v_numerical = []                    # List with the results of v from the numerical solution
        self.r_analytical = []                   # List with the results of the analytical solution

        self.mu = self.G*(self.M+self.system.masses)     # Standard gravitational parameter
        self.timesteps_pr_orbit=10000            # Number of timesteps pr orbital period to calculate


    def leapfrog(self, r0, v0, T, dt):
        '''
        Leapfrog integration
        '''
        G = self.G                               # For less writing
        N = int(T/dt)                            # Length of all our vectors
        print(f'{N} = {T}/{dt}')
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

        a[0, :] = -G*M/(distance[0]**3) * r[0, :]
        for i in range(N-1):
            '''
            The actual leapfrog algorithm
            '''
            r[i + 1, :] = r[i, :] + v[i, :]*dt + 0.5*a[i, :]*dt**2
            distance[i + 1] = np.linalg.norm(r[i + 1, :])

            a[i + 1, :] = -G*M/(distance[i + 1]**3) * r[i + 1, :]
            v[i + 1, :] = v[i, :] + 0.5*(a[i, :] + a[i + 1, :])*dt
            t[i + 1] = t[i] + dt

        return r, v, a, t

    def euler_cromer(self, r0, v0, T, dt):
        '''
        euler cromer
        '''
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
        '''
        Returns the distance r for all planets at timestep t from the simulation
        '''
        r_vec = np.array(self.r_numerical)[:,t]
        return np.sqrt(r_vec[:,0]**2+r_vec[:,1]**2)

    def v(self,t):
        '''
        Returns the velocity for all planets at timestep t from the simulation
        '''
        v_vec = np.array(self.v_numerical)[:,t]
        return np.sqrt(v_vec[:,0]**2+v_vec[:,1]**2)

    def simulated_aphelion(self):
        '''
        Finds maximum distance from star for all planets
        '''
        r_vec = np.array(self.r_numerical)
        r = np.sqrt(r_vec[:,:,0]**2+r_vec[:,:,1]**2)
        return r[:,np.argmax(r,axis=1)]


    def mean_velocity(self):
        '''
        Returns the mean velocity of the planets
        '''
        v_vec = np.array(self.v_numerical)
        # Calculate velocity at all times
        velocities = np.sqrt(v_vec[:,:,0]**2+v_vec[:,:,1]**2)
        # Return the mean velocity for each planet
        return np.mean(velocities,axis=1)

    def dArea(self,t):
        '''
        Calulate the area swept out at timestep t and over time self.dt
        '''
        return 0.5*self.r(t)[0]*(self.v(t)*self.dt)

    def area_swiped_during_period(self,t,days):
        '''
        Calculate area swiped during a number of days.
        '''
        dtPrDay = self.timesteps_pr_orbit/(self.hoth_period()*c.yr/c.day)
        n = int(dtPrDay*days)
        totalArea = 0
        for i in range(n):
            totalArea += self.dArea(t+i)

        return totalArea

    def distance_traveled_during_period(self,t,days):
        '''
        Numerically calculate the distance traveled during # of days
        '''
        dtPrDay = self.timesteps_pr_orbit/(self.hoth_period()*c.yr/c.day)
        n = int(dtPrDay*days)
        totalDistance = 0
        for i in range(n):
            totalDistance += self.v(t)*self.dt

        return totalDistance

    def hoth_period(self):
        '''
        The rotational periode of our planet in years
        '''
        return 2*np.pi*np.sqrt(self.axes[0]**3/self.mu[0])
    def hoth_perhelion_t(self):
        return int(self.hoth_period()/self.dt/2)

    def sim(self):
        '''
        Simulating all the orbits
        '''

        print('Simulating orbits')

        mu = self.mu                                                # Standard gravitational parameter


        N = len(self.system.masses)                                 # Length for for loop

        rotational_orbit_in_years = 2*np.pi*np.sqrt(self.axes**3/mu)      # Rotational orbit in years
        print(f'Orbit in years {rotational_orbit_in_years}')
        T = 45*rotational_orbit_in_years[0]                            # Total time of N rotations
        self.dt = rotational_orbit_in_years[0]/self.timesteps_pr_orbit # Find dt based on timesteps pr year
        planet_pos = self.system.initial_positions                  # Initial planets positions
        planet_vel = self.system.initial_velocities                 # Initial planets velocities
        for i in range(N):
            print(f'Working on planet {i}')
            print(f'Orbit in years {2*np.pi*np.sqrt(self.axes[i]**3/mu[i])}')
            m = self.system.masses[i]                                   # Gets i'th planet mass

            r0 = planet_pos[:, i]                              # Gets i'th planet starting pos

            v0 = planet_vel[:, i]                              # --||--                    velocity
            r, v, a, t = self.leapfrog(r0, v0, T, self.dt)    # runs leapfrog and returns
            self.r_numerical.append(r)
            self.v_numerical.append(v)
        self.analytical_solution()                      # analytical

    def cartesian_polar(self, r):
        '''
        Converts to polar coordinates
        '''
        x = r[:, 0]                             # x values
        y = r[:, 1]                             # y values
        r_p = np.zeros(len(r))
        for i in range(len(r)):
            r_p[i] = np.linalg.norm(r[i])
        # r = np.sqrt(x**2 + y**2)              # distance from origo
        theta = np.arctan(y/x)                  # theta
        return r_p, theta

    def polar_cartesian(self, r, theta):
        '''
        Converts to cartesian
        '''
        N = len(theta)
        x = np.zeros(N, float)
        y = np.zeros(N, float)
        for i in range(N):
            x[i] = r[i]*np.cos(theta[i])        # calculates x
            y[i] = r[i]*np.sin(theta[i])        # calculates y
        return x, y

    def analytical_solution(self):
        theta = np.linspace(0, 2*np.pi, 1000)   # array from 0 to 2pi

        def r(axes, e, theta):
            '''
            Analytical formula
            '''
            ans = (axes*(1 - e**2))/(1 + e*np.cos(theta))
            return ans

        for i in range(len(self.axes)):
            x, y = self.polar_cartesian(r(self.axes[i], self.e[i], theta), theta)
            self.r_analytical.append([y,x])

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

    def solar_orbit(self):
        '''
        Comment
        '''
        pass


if __name__ == '__main__':
    filename = "simulated_orbits.pkl"
    orbit = los.orbit_sim_factory(filename,args)

    orbit.plot()
