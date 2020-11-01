# Egen kode
import numpy as np
import ast2000tools.constants as c
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
import tools


class hohmann:
    def __init__(self, seed=33382):
        """
        All that will be used across the whole class
        """
        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)

        self.G = c.G
        self.M = self.system.star_mass
        self.mu = self.G*self.M

        self.times, self.planet_pos = np.load('planet_trajectories.npz',
                                              allow_pickle=True)

    def theta_needed(self, P_2, N):
        """
        This will calculate the
        theta needed between 2 planets
        for the optimal launch

        theta is in radians
        P_2 needs to be the orbital period for our destination
        N is time frame in days
        """
        theta = np.pi*(P_2 - 2*T)/P_2
        return theta

    def theta_now(self, planet_start, planet_goal, N):
        """
        Finds the theta between the planets

        planet_* are the planet id numbers from system
        N is time frame in days
        """
        r1 = planet_start
        r2 = planet_goal
        theta = np.zeros(N, float)
        for i in range(N):
            r1_norm = np.linalg.norm(r1[i])
            r2_norm = np.linalg.norm(r2[i])
            theta[i] = np.arccos(np.dot(r1[i], r2[i])/(r1_norm * r2_norm))
        return theta

    def find_the_launch(self, planet_start, planet_goal, N, tol=1):
        """
        This will compare the thetas and return the correct
        time of for transfer

        planet_* are the planet id numbers from system
        tol is the tolerance level of the theta, might need to be changed
        theta in radians
        N is how many years we will look at
        """
        a = self.new_axis(planet_start, planet_goal)
        a_goal = self.system.semi_major_axis[planet_goal]
        T = self.time(a)
        P_2 = self.kepler_third(a_goal)

        N = round(self.time_in_days(N))

        theta_needed = self.theta_needed(P_2, N)
        theta_now = self.theta(planet_start, planet_goal)

        launch_T = []  # Decided to have it as a list since we might get more than one
        for i in range(N):
            if abs(theta_now[i] - theta_needed[i]) < tol:
                launch_T.append(T[i])

        return launch_T

    def time(self, a):
        """
        This will calculate the time it takes
        to complete the interplanetary travel
        """
        T = np.pi * np.sqrt(a**3/self.mu)
        return T

    def time_in_seconds(self, T):
        """
        Converts years to seconds
        """
        T_new = self.time_in_days(T)*86400
        return T_new

    def time_in_days(self, T):
        """
        Converts years to days
        """
        T_new = T*365.2425
        return T_new


    def new_axis(self, planet_start, planet_goal, N):
        """
        This will calculate the new semi major axis
        for the Hohmann transfer

        a will then be an array of all the possible ellipses
        r1 is the posititon vector for start planet
        r2 is the posititon vector for goal planet
        N is the time frame in years
        """
        r1 = planet_start
        r2 = planet_goal
        N = round(time_in_days(N))
        a = np.zeros(N, float)

        for i in range(N):
            r1_norm = np.linalg.norm(r1[i])
            r2_norm = np.linalg.norm(r2[i])
            a[i] = (r1_norm + r2_norm)*.5

        return a

    def kepler_third(self, a):
        """
        Calculates the orbital period
        """
        P = np.sqrt(4*np.pi**2*a/self.mu)
        return P

    def planet_identifier(self, planet_id):
        """
        Gets all relavant info on the planet
        """
        a = self.system.semi_major_axis[planet_id]
        P = self.kepler_third(a)
        R = self.planet_pos[planet_id]
        R_norm = np.linalg.norm(R)

        return a, P, R, R_norm

    def v_at_planet(self, planet_id):
        """
        This is for calculating the v
        needed from periapsis to achieve
        the Hohmann transfer
        """
        a, P, R, R_norm = self.planet_identifier(planet_id)
        v = 2*np.pi*a/P * np.sqrt(2*a/R_norm - 1)

        return v

    def v_from_planet(self, planet_id):
        """
        Calculates the v from the planet
        seen from the sun
        """
        a, P, R, R_norm = self.planet_identifier(planet_id)

        v = 2*np.pi*R_norm/P
        return v

    def delta_v(self, planet_id):
        """
        Calculates the delta v from
        the v from planet and the initial
        velocity
        """
        v_planet = self.v_at_planet(planet_id)
        v = self.v_from_planet(planet_id)
        delta_v = v_planet - v
        return delta_v

    def check_if_close(self, planet_id):
        """
        Checks if spacecraft is close enough for
        injection maneuver
        """
        M_planet = self.system.masses[planet_id]
        M_sol = self.M
        current_positon = 1  # This needs to be an array from sun to spacecraft at that time
        l_test = current_positon*np.sqrt(M_planet/(10*M_sol))
        return l_test
