# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c
import os
import sys
try:
    import cPickle as pickle
except:
    import pickle
import argparse


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
        self.G = 4*np.pi**2                      # AU**3 yr**-2 SolarMass**-1
        self.M = self.system.star_mass           # Star mass in solar mass
        self.r_numerical = []                    # List with the results of the numerical solution
        self.r_analytical = []                   # List with the results of the analytical solution

    def leapfrog(self, r0, v0, T, dt):
        '''
        Leapfrog integration
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

        a[0, :] = -G*M/(distance[0]**3) * r[0, :]
        for i in range(N-1):
            '''
            The actual leapfrog algorithm
            '''
            r[i + 1, :] = r[i,:] + v[i, :]*dt + 0.5*a[i, :]*dt**2
            distance[i + 1] = np.linalg.norm(r[i + 1, :])

            a[i + 1, :] = -G*M/(distance[i + 1]**3) * r[i + 1, :]
            v[i + 1, :] = v[i, :] + 0.5*(a[i, :] + a[i + 1, :])*dt
            #a[i, :] = a[i + 1, :]
            t[i + 1] = t[i] + dt
            #print(f'F = {m*a[i, :]}')

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

    def sim(self):
        '''
        Simulating all the orbits
        '''
        print('Simulating orbits')
        mu = self.G*self.M                                  # Standard gravitational parameter

        N = len(self.system.masses)                         # Length for for loop

        planet_pos = self.system.initial_positions    # Initial planets positions
        planet_vel = self.system.initial_velocities   # Initial planets velocities
        for i in range(N):
            print(f'Working on planet {i}')
            orbital_period = 2*np.pi*np.sqrt(self.axes[i]**3/mu)      # One year
            T = 20*orbital_period                               # 20 years
            dt = orbital_period/10000                           # Timestep for 1 year

            m = self.system.masses[i]                       # Gets i'th planet mass

            r0 = planet_pos[:, i]                              # Gets i'th planet starting pos

            v0 = planet_vel[:, i]                              # --||--                    velocity
            r, v, a, t = self.leapfrog(r0, v0, T, dt)    # runs leapfrog and returns
            self.analytical_solution()                      # analytical
            self.r_numerical.append(r)

    def cartesian_polar(self, r):
        '''
        Converts to polar coordinates
        '''
        x = r[:, 0]                                         # x values
        y = r[:, 1]                                         # y values
        r_p = np.zeros(len(r))
        for i in range(len(r)):
            r_p[i] = np.linalg.norm(r[i])
        #r = np.sqrt(x**2 + y**2)                            # distance from origo
        theta = np.arctan(y/x)                              # theta
        return r_p, theta

    def polar_cartesian(self, r, theta):
        '''
        Converts to cartesian
        '''
        N = len(theta)
        x = np.zeros(N, float)
        y = np.zeros(N, float)
        for i in range(N):
            x[i] = r[i]*np.cos(theta[i])                                 # calculates x
            y[i] = r[i]*np.sin(theta[i])                                 # calculates y
        return x, y

    def analytical_solution(self):
        theta = np.linspace(0, 2*np.pi, 1000)               # array from 0 to 2pi

        def r(axes, e, theta):
            '''
            Analytical formula
            '''
            ans = (axes*(1 - e**2))/(1 + e*np.cos(theta))
            return ans

        for i in range(len(self.axes)):
            x, y = self.polar_cartesian(r(self.axes[i], self.e[i], theta), theta)
            self.r_analytical.append([x,y])

    def plot(self):
        for r in self.r_numerical:
            plt.plot(r[:, 0], r[:, 1])
            plt.axis('equal')

        for a in self.r_analytical:
            plt.plot(a[0],a[1])

        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Hoth system')

        plt.show()


if __name__ == '__main__':
    filename = "simulated_orbits.pkl"
    url = 'https://www.uio.no/studier/emner/matnat/astro/AST2000/h20/blogger/Flukten%20fra%20Hoth/data/simulated_orbits.pkl'
    if (os.path.exists(filename) == False or args.download==True):
        try:
            import requests
            r = requests.get(url, allow_redirects=True)

            open(filename, 'wb').write(r.content)
        except:
            print('You need to install requests to download file: pip install requests')
            
    if (os.path.exists(filename) == False or args.run_sim==True):
        orbit = orbit_sim()
        orbit.sim()

        with open(filename, "wb") as output:
            pickle.dump(orbit, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as input:
            orbit = pickle.load(input)

orbit.plot()
