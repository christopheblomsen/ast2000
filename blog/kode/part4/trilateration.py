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

    def comparrison(self, vec1, vec2):
        '''
        Comparing 2 vectors
        '''
        for i in range(len(vec1)):
            if vec1[i, 0] == vec2[i, 0] and vec1[i, 1] == vec2[i, 1]:
                correct = i

        return correct

    def two(self, radii):
        '''
        r1, r2, r3 are the distances measured by the radar
        U is the distance between the star and the planet
        '''
        # x = (r1**2 - r2**2 + U**2)/(2*U)
        # y = (r1**2 - r3**2 + V**2 - 2*V[0]*x)/(2*V[1])
        theta = np.linspace(0, 2*np.pi, 1000)
        N = len(radii)
        vec = np.zeros((N, 2), float)

        for i in range(N):
            vec[i, 0] = radii[i]*np.sin(theta)
            vec[i, 1] = radii[i]*np.cos(theta)
        correct = []

        for i in range(len(vec)+1):
            correct.append(float(self.comparrison(vec[i], vec[i+1])
