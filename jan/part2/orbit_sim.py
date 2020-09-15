# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem


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
        a = self.system.semi_major_axes
        e = self.system.eccentricities

    def plot(self):
        '''
        For plotting
        '''
        pass

    def leapfrog(self, x0, v0, T, dt):
        '''
        Leapfrog integration
        '''
        def P(t):
            '''
            Insert stuff here
            '''
            pass

        def Q(t):
            '''
            Insert stuff here
            '''
            pass

        def R(t):
            '''
            Insert stuff here
            '''
            pass

        t = dt
        N = int(T/dt)

        x = x0
        v = v0

        a_i = R(t) - P(t)*v - Q(t)*x

        for i in range(N):
            t += dt
            x += v*dt + 0.5*a_i*dt**2
            a_i_pluss1 = R(t) - P(t)*v - Q(t)*x
            v += 0.5*(a_i + a_i_pluss1)*dt
            a_i = a_i_pluss1

        return x


