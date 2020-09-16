# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c

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

    def plot(self):
        '''
        For plotting
        '''
        pass

    def leapfrog(self, r0, v0, R0, T, dt):
        '''
        Leapfrog integration
        '''
        G = 4*np.pi**2                  # AU**3 yr**-2 SolarMass**-1
        N = int(T/dt)
        t = np.zeros(N, float)
        M = self.system.solarMass
        m = self.system.masses
        r = np.zeros((N, 2), float)
        v = np.zeros((N, 2), float)
        R = np.zeros((N, 2), float)

        distance = np.zeros(N, float)

        r[0, :] = r0
        v[0, :] = v0
        R[0] = R0
        t[0] = 0

        a[0, :] = -G*M*m/(distance[0]**2) * R[0, :]
        for i in range(N):
            '''
            The actual leapfrog algorithm
            '''
            r[i + 1, :] = v[i, :]*dt + 0.5*a[i, :]*dt**2
            distance[i + 1] = np.linalg.norm(r[i + 1, :] - [0, 0])

            a[i + 1, :] = -G*M*m/(distance[i + 1]**2) * R[i + 1, :]
            v[i + 1, :] = 0.5*(a[i, :] + a[i + 1, :])*dt
            a[i, :] = a[i + 1, :]
            t[i + 1] = t[i] + dt

        return r, v, a, t

    def orbit_sim(self, dt):
        '''
        Simulating all the orbits
        '''
        N = len(self.system.masses)
        T = 365                     # One year
        for i in range(N):
            r0 = planet_starting_pos[i]
            v0 = planet_velocity_at_t_0[i]
            R0 = r0 - [0, 0]        # Start distance
            r, v, a = self.leapfrog(r0, v0, R0, T, dt)
            plt.plot(r[:, 0], r[:, 1])
        plt.show()


