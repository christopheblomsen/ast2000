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

    def leapfrog(self, x0, v0, T, dt):
        '''
        Leapfrog integration
        '''
        G = c.G
        N = int(T/dt)
        t = np.linspace(0, T, N)

        x = np.zeros((N, 2), float)
        v = np.zeros((N, 2), float)
        R = np.zeros((N, 2), float)

        distance = np.zeros(N, float)

        x[0, :] = x0
        v[0, :] = v0

        for i in range(N):
            a[i, :] = -G*M*m/(distance[i]**2) * R[i, :]

            x[i + 1, :] = v[i, :]*dt + 0.5*a[i, :]*dt**2
            distance[i + 1] = np.linalg.norm(x[i + 1, :] - [0, 0])

            a[i + 1, :] = -G*M*m/(distance[i + 1]**2) * R[i + 1, :]
            v[i + 1, :] = 0.5*(a[i, :] + a[i + 1, :])*dt
            a[i, :] = a[i + 1, :]

        return x, v, a
