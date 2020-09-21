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
        self.m = self.system.masses[5]           # Home planet mass
        self.omega = self.system.aphelion_angles # Aphelion angles

        self.r_numerical = []                    # List with the results of the numerical solution
        self.r_analytical = []                   # List with the results of the analytical solution

        self.mu = self.G*(self.M + self.m)       # Standard gravitational parameter


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

    def sim(self):
        '''
        Simulating all the orbits
        '''

        print('Simulating orbits')

        mu = self.mu                                                # Standard gravitational parameter


        N = len(self.system.masses)                                 # Length for for loop

        orbital_period = 2*np.pi*np.sqrt(self.axes[5]**3/mu)      # One year
        T = 20*orbital_period                                       # 20 years
        # dt = orbital_period/10000                           # Timestep for 1 year
        dt = 0.01

        planet_pos = self.system.initial_positions                  # Initial planets positions
        planet_vel = self.system.initial_velocities                 # Initial planets velocities
        verification_r = np.zeros((2, N, int(T/dt)), float)
        #verification_t = []
        for i in range(N):
            print(f'Working on planet {i+1}')

            r0 = planet_pos[:, i]                              # Gets i'th planet starting pos

            v0 = planet_vel[:, i]                              # --||--                    velocity
            r, v, a, t = self.leapfrog(r0, v0, T, dt)    # runs leapfrog and returns
            #print(np.shape(r))
            #print(np.shape(verification_r))
            verification_r[0, i, :] = r[:, 0]
            verification_r[1, i, :] = r[:, 1]
            self.analytical_solution()                      # analytical
            self.r_numerical.append(r)

        #verification_t = np.asarray(verification_t)
        #print(np.shape(verification_t[1, :]))
        verification_t = t.reshape(-1)
        print(np.shape(verification_t))
        return verification_r, verification_t, self.system

    def solar_orbit(self, planet):
        '''
        Comment
        '''
        star_initial_pos = np.array([0, 0])
        star_initial_vel = np.array([0, 0])

        planet_pos = self.system.initial_positions[:, planet]
        planet_vel = self.system.initial_velocities[:, planet]

        M = self.M
        m = self.system.masses[planet]
        mu = self.G*(M + self.system.masses[planet])

        orbital_period = 2*np.pi*np.sqrt(self.axes[planet]**3/mu)      # One year
        dt = orbital_period/100000

        masses = np.array([m, M])
        N = len(m)

        r = planet_pos - star_initial_pos

        center_of_mass_R = self.center_mass(masses, r)

    def sim_solar_orbit(self, r0, v0, dt):
        '''
        Simulating solar orbit
        '''
        G = self.G                               # For less writing
        N = int(T/dt)                            # Length of all our vectors
        t = np.zeros(N, float)                   # time array
        M = self.M                               # Star mass

        r_planet = np.zeros((N, 2), float)       # Position vector
        v_planet = np.zeros((N, 2), float)       # Velocity vector
        a_planet = np.zeros((N, 2), float)       # Acceleration vector

        r_sol = np.zeros((N, 2), float)          # Position vector
        v_sol = np.zeros((N, 2), float)          # Velocity vector
        a_sol = np.zeros((N, 2), float)          # Acceleration vector

        distance = np.zeros(N, float)            # Distance array

        distance[0] = np.linalg.norm(r0)         # Sets initial conditions
        r_planet[0, :] = r0
        v_planet[0, :] = v0
        t[0] = 0

        a_planet[0, :] = -G*M/(distance[0]**3) * r[0, :]

        r_sol[0, :] = r0
        v_sol[0, :] = v0
        a_sol[0, :] = -G*M/(distance[0]**3) * r[0, :]
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


    def center_mass(self, m, r):
        M = np.sum(m)
        R = np.array([0, 0])

        for i in range(len(r)):
            R = R + m[i] * r[i]

        self.R = R/M

        return self.R

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

        def r(axes, e, theta, omega):
            '''
            Analytical formula
            '''
            ans = (axes*(1 - e**2))/(1 + e*np.cos(theta - omega))
            return ans

        for i in range(len(self.axes)):
            x, y = self.polar_cartesian(r(self.axes[i], self.e[i], theta, self.omega[i]), theta)
            self.r_analytical.append([x, y])

    def plot(self):
        for r in self.r_numerical:
            plt.plot(r[:, 0], r[:, 1])
            plt.axis('equal')

        for a in self.r_analytical:
            plt.plot(a[0],a[1])

        plt.xlabel('x[AU]')
        plt.ylabel('y[AU]')
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
        r, T, system = orbit.sim()

        with open(filename, "wb") as output:
            pickle.dump(orbit, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as input:
            orbit = pickle.load(input)

orbit.plot()
system.verify_planet_positions(T, r, 'verification_of_planets')
system.generate_orbit_video(T, r)
