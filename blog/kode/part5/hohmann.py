import numpy as np
import ast2000tools.constants as c
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
import tools

class hohmann:
    def __init__(self, seed=33382):
        """
        Come back to this
        """
        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)

        self.G = c.G
        self.M = self.system.star_mass
        self.mu = self.G*self.M

        self.times, self.planet_pos = np.load('planet_trajectories.npz',
                                              allow_pickle=True)

    def theta_needed(self, P_2, T):
        """
        This will calculate the
        theta needed between 2 planets
        for the optimal launch
        """
        theta = np.pi*(P_2 - 2*T)/P_2

        return theta

    def time(self, a):
        """
        This will calculate the time it takes
        to complete the interplanetary travel
        """
        T = np.pi * np.sqrt(a/self.mu)
        return T

    def theta_now(self, planet_start, planet_goal):
        r1 = planet_start
        r2 = planet_goal
        N = len(r1)
        theta = np.zeros(N, float)
        for i in range(N):
            r1_norm = np.linalg.norm(r1[i])
            r2_norm = np.linalg.norm(r2[i])
            theta[i] = np.arccos(np.dot(r1[i], r2[i])/(r1_norm * r2_norm))
        return theta

    def new_axis(self, planet_start, planet_goal):
        """
        This will calculate the new semi major axis
        for the Hohmann transfer
        """
        r1 = planet_start
        r2 = planet_goal
        N = len(r1)
        a = np.zeros(N, float)

        for i in range(N):
            r1_norm = np.linalg.norm(r1[i])
            r2_norm = np.linalg.norm(r2[i])
            a[i] = (r1_norm + r2_norm)*.5

        return a

    def kepler_third(self, a):
        """
        Calculates te orbital period
        """
        P = np.sqrt(4*np.pi**2*a/self.mu)
        return P
    def v_per(self, a, P, 

