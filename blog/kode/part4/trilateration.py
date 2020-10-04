# Egen kode
import numpy as np
import matplotlib.pyplot
from ast2000tools.solar_system import SolarSystem
import ast2000tools.contstants as c
import orbit_sim as os


class trilateration:
    def __init__(self, seed=33382, dt=1000):
        '''
        Something
        '''
        self.system = SolarSystem(seed)
        self.orbit = os
        self.mu = c.G*(self.system.star_M + self.system.masses[0])
        self.dt = dt

        self.r_numerical = []
        self.v_numerical = []

    def just_calculation_orbit(self):
        mu = self.mu
        N = len(self.system.masses)

        rotational_orbit_in_years = 2*np.pi*np.sqrt(self.axes**3/mu)
        self.T = rotational_orbit_in_years[0]

        planet_pos = self.system.initial_positions
        planet_vel = self.system.initial_velocities

        for i in range(N):
            r0 = planet_pos[:, i]
            v0 = planet_vel[:, i]

            r, v, a, t, = self.orbit.leapfrog(r0, v0, self.T, self.dt)
            self.r_numerical.append(r)
            self.v_numerical.append(v)

    def comparrison(self, vec1, vec2, tol=1):
        '''
        Comparing 2 vectors
        '''
        for i in range(len(vec1)):
            tolerance_x = abs(vec1[i, 0] - vec2[i, 0])
            tolerance_y = abs(vec1[i, 1] - vec2[i, 1])
            if tolerance_x < tol and tolerance_y < tol:
                correct = vec1[i, :]

        return correct

    def circles(self, radii, a, b):
        '''
        calculating the x and y values around
        planets given from a radii from the radar
        a, and b are the planet positions in that
        timestep
        '''
        theta = np.linspace(0, 2*np.pi, 1000)
        N = len(radii)
        vec = np.zeros((N, 2), float)

        for i in range(N):
            '''
            x and y for itself is just to remove the
            E501 warning for PEP8
            '''
            x = (a[i]**2 + b[i]**2 + 2*b[i]*radii[i]*np.sin(theta))/(2*a[i])
            y = (2*a[i]*radii[i]*np.cos(theta) - a[i]**2 - b[i]**2)/(2*b[i])
            vec[i, 0] = x
            vec[i, 1] = y
        self.vec = vec

    def same_coordinates(self):
        '''
        Checks if the coordinates are correct
        '''
        self.circles()
        correct = []

        for i in range(len(self.vec)+1):
            '''
            Finds the correct x and y coordinates
            '''
            correct.append(float(self.comparrison(self.vec[i], self.vec[i+1])))

        return correct
