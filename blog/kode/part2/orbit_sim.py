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
        self.system = SolarSystem(seed)
        self.a = self.system.semi_major_axes
        self.e = self.system.eccentricities


    def leapfrog(self, r0, v0, T, dt, m):
        '''
        Leapfrog integration
        '''
        G = 4*np.pi**2                  # AU**3 yr**-2 SolarMass**-1
        N = int(T/dt)
        t = np.zeros(N, float)
        M = self.system.star_mass
        r = np.zeros((N, 2), float)
        v = np.zeros((N, 2), float)
        a = np.zeros((N, 2), float)
        R = np.zeros((N, 2), float)

        distance = np.zeros(N, float)

        distance[0] = np.linalg.norm(r0)
        r[0, :] = r0
        v[0, :] = v0
        t[0] = 0

        a[0, :] = -G*M/(distance[0]**3) * R[0, :]
        for i in range(N):
            '''
            The actual leapfrog algorithm
            '''
            r[i + 1, :] = v[i, :]*dt + 0.5*a[i, :]*dt**2
            distance[i + 1] = np.linalg.norm(r[i + 1, :])

            a[i + 1, :] = -G*M*m/(distance[i + 1]**2) * r[i + 1, :]
            v[i + 1, :] = 0.5*(a[i, :] + a[i + 1, :])*dt
            a[i, :] = a[i + 1, :]
            t[i + 1] = t[i] + dt

        return r, v, a, t

    def sim(self):
        '''
        Simulating all the orbits
        '''
        N = len(self.system.masses)
        dt = 1
        T = 365*81400                     # One year
        r0 = np.array([0, 0])
        planet_pos = self.system.initial_positions[0, :]
        planet_vel = self.system.initial_velocities[0, :]
        for i in range(N):
            m = self.system.masses[i]
            r0 = planet_pos[i]
            v0 = planet_vel[i]
            r, v, a , t = self.leapfrog(r0, v0, T, dt, m)
            plt.plot(r[:, 0], r[:, 1])
        plt.show()



if __name__ == "__main__":
    dt = 1
    filename = "simulated_orbits.pkl"

    if (os.path.exists(filename) == False):
        orbit = orbit_sim()
        orbit.sim()

        with open(filename, "wb") as output:
            pickle.dump(orbit, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as input:
            orbit = pickle.load(input)
