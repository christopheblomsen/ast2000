"""Egen kode."""
import numpy as np
import ast2000tools.constants as c
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from trajection_computer import TrajectionComputer
import matplotlib.pyplot as plt


class Hohmann:
    """Class for calculating a Hohmann transfer between to planets."""
    def __init__(self, seed=33382):
        """All that will be used across the whole class."""
        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)
        self.tr_comp = TrajectionComputer(1579)

        self.G = c.G_sol
        self.M = self.system.star_mass
        self.mu = self.G*self.M

        self.times = self.tr_comp.times
        self.planet_pos = self.tr_comp.planet_positions

    def theta_needed(self, P_2, T):
        """This will calculate the theta needed between 2 planets for the optimal launch.

        theta is in radians
        P_2 needs to be the orbital period for our destination
        T is time frame in days
        """
        theta = np.pi*(P_2 - 2*T)/P_2
        return theta

    def theta_now(self, planet_start, planet_goal, N):
        """Finds the theta between the planets.

        planet_* are the planet id numbers from system
        N is time frame in days
        """
        theta = np.zeros(N, float)
        for i in range(N):
            r1 = self.tr_comp.planet_position(planet_start, i*c.day/c.yr)
            r2 = self.tr_comp.planet_position(planet_goal, i*c.day/c.yr)
            r1_norm = np.linalg.norm(r1)
            r2_norm = np.linalg.norm(r2)
            theta[i] = np.arccos(np.dot(r1, r2)/(r1_norm * r2_norm))
        return theta

    def find_the_launch(self, planet_start, planet_goal, N, tol=0.01):
        """This will compare the thetas and return the correct time of for transfer.

        planet_* are the planet id numbers from system
        tol is the tolerance level of the theta, might need to be changed
        theta in radians
        N is how many years we will look at
        """
        print(f"Find launch time from planet {planet_start} to {planet_goal}")
        a = self.new_axis(planet_start, planet_goal, N)
        a_goal = self.system.semi_major_axes[planet_goal]
        a_start = self.system.semi_major_axes[planet_start]
        T = self.time(a)
        P_1 = self.kepler_third(a_start)
        P_2 = self.kepler_third(a_goal)

        print(f"Orbital period of start planet {planet_start} is {P_1} years")
        print(f"Orbital period of destination planet {planet_goal} is {P_2} years")

        theta_needed = self.theta_needed(P_2, T)
        N = round(self.time_in_days(N))
        theta_now = self.theta_now(planet_start, planet_goal, N)
        # print(f'Theta now {theta_now}')
        # print(f'Theta needed {theta_needed}')

        t = np.linspace(0, N, theta_now.shape[0])
        # plt.plot(t, T, label='Time needed')
        plt.plot(t, theta_now, label='Theta now')
        plt.plot(t, theta_needed, label='Theta needed')

        launch_T = []  # Decided to have it as a list since we might get more than one
        launch_tol = []
        launch_i = []
        found_theta = []
        for i in range(N):
            # print(f'theta_now[{i}]: {theta_now[i]}')
            if abs(theta_now[i] - theta_needed[i]) < tol:
                launch_T.append(T[i])
                launch_tol.append(abs(theta_now[i] - theta_needed[i]))
                launch_i.append(i)
                found_theta.append(theta_needed[i])

        return np.array(launch_T), np.array(launch_tol), np.array(launch_i), np.array(found_theta)

    def time(self, a):
        """This will calculate the time it takes to complete the interplanetary travel."""
        T = np.pi * np.sqrt(a**3/self.mu)
        return T

    def time_in_seconds(self, T):
        """Converts years to seconds."""
        T_new = self.time_in_days(T)*c.day
        return T_new

    def time_in_days(self, T):
        """Converts years to days."""
        T_new = T*c.yr/c.day
        return T_new

    def new_axis(self, planet_start, planet_goal, N):
        """This will calculate the new semi major axis for the Hohmann transfer.

        a will then be an array of all the possible ellipses
        N is the time frame in years
        """
        N = round(self.time_in_days(N))
        self.a = np.zeros(N, float)

        for i in range(N):
            r1_norm = np.linalg.norm(self.tr_comp.planet_position(planet_start, i*c.day/c.yr))
            r2_norm = np.linalg.norm(self.tr_comp.planet_position(planet_goal, i*c.day/c.yr))
            self.a[i] = (r1_norm + r2_norm)*.5

        return self.a

    def kepler_third(self, a):
        """Calculates the orbital period."""
        P = np.sqrt(4*np.pi**2*a**3/self.mu)
        return P

    def planet_identifier(self, planet_id, time):
        """Gets all relavant info on the planet."""
        a = self.system.semi_major_axes[planet_id]
        P = self.kepler_third(a)
        R = self.tr_comp.planet_position(planet_id, time)
        print(f'Planet {planet_id} position {R}')
        R_norm = np.linalg.norm(R)

        return a, P, R, R_norm

    def v_at_planet(self, planet_id, time):
        """This is for calculating the v needed from periapsis to achieve the Hohmann transfer."""
        a, P, R, R_norm = self.planet_identifier(planet_id, time)
        a = self.a[int(self.time_in_days(time))]
        P = self.kepler_third(a)
        print(f'v_at_planet: a {a}, P {P}, R {R}, R_norm {R_norm}')
        v = 2*np.pi*a/P * np.sqrt(2*a/R_norm - 1)
        print(f'v is {v}')

        return v

    def v_from_planet(self, planet_id, time):
        """Calculates the v from the planet seen from the sun."""
        a, P, R, R_norm = self.planet_identifier(planet_id, time)
        print(f'v_from_planet: a {a}, P {P}, R, {R} R_norm {R_norm}')
        v = 2*np.pi*R_norm/P
        print(f'v is {v}')
        return v

    def delta_v(self, planet_id, time):
        """Calculates the delta v from the v from planet and the initial velocity."""
        v_planet = self.v_at_planet(planet_id, time)
        v = self.v_from_planet(planet_id, time)
        delta_v = v_planet - v
        return delta_v

    def check_if_close(self, planet_id):
        """Checks if spacecraft is close enough for injection maneuver."""
        M_planet = self.system.masses[planet_id]
        M_sol = self.M
        current_positon = 1  # This needs to be an array from sun to spacecraft at that time
        l_test = current_positon*np.sqrt(M_planet/(10*M_sol))
        return l_test

    def delta_v_for_boost(self, position):
        """Return the delta velocity needed to boost the rocket."""
        theta = np.arctan(position[0]/position[1])
        dx = -np.sin(theta)
        dy = np.cos(theta)
        dv = [dx, dy]
        return np.array(dv)

    def plot_positions(self, time):
        hoth_pos = self.tr_comp.planet_position(0, time)
        pjuf_pos = self.tr_comp.planet_position(6, time)
        print(f'Hoth position at {time} is {hoth_pos[0]}, {hoth_pos[1]}')
        print(f'Hoth position at {time} is {pjuf_pos[0]}, {pjuf_pos[1]}')
        plt.scatter(hoth_pos[0], hoth_pos[1], label='Hoth')
        plt.scatter(pjuf_pos[0], pjuf_pos[1], label='Pjuf')
        plt.scatter(0, 0)
        plt.arrow(0, 0, hoth_pos[0], hoth_pos[1])
        plt.arrow(0, 0, pjuf_pos[0], pjuf_pos[1])
        plt.xlabel('AU')
        plt.ylabel('AU')
        plt.legend()
        plt.show()


if __name__ == "__main__":
    hohmann = Hohmann()

    launch_times, launch_tolerances, launch_i, found_theta = hohmann.find_the_launch(0, 6, 1)
    idx = np.argmin(launch_tolerances) + 1
    print(f'Found launch time {launch_i[idx]}')
    print(f'With tolerance {launch_tolerances}')
    print(f'At time {launch_i*c.day/c.yr} years')
    print(f'Angle is {np.rad2deg(found_theta[idx])}')
    time = launch_i[idx]*c.day/c.yr
    v_per = hohmann.v_at_planet(0, time)
    print(f'We need a velocity of {v_per} AU/yr to achieve the Hohmann transfer')
    print(f'Delta v to enter orbit {hohmann.delta_v(0, time)} AU/yr to achieve the Hohmann transfer')
    planet_pos = hohmann.tr_comp.planet_position(0, time)

    dvec = hohmann.delta_v_for_boost(planet_pos)
    print(f'Delta v vector to enter orbit [{dvec[0]},{dvec[1]}]')
    print('======================================')
    print(f'v at perhelion {v_per} AU/yr')
    print(f'Planet position is [{planet_pos[0]},{planet_pos[1]}]')
    cos_theta = planet_pos[0]/np.linalg.norm(planet_pos)
    sin_theta = planet_pos[1]/np.linalg.norm(planet_pos)
    print(f'cos theta = {cos_theta}')
    print(f'sin theta = {sin_theta}')
    v_theta = [-v_per*sin_theta, v_per*cos_theta]
    print(f'v theta = {v_theta}')

    plt.xlabel('Days')
    plt.ylabel('Radians')
    plt.legend()
    plt.show()

    hohmann.plot_positions(time)
"""
Kjøreeksempel:
$> python hohmann.py
Planet trajectories 8 planets loaded
Find launch time from planet 0 to 6
Orbital period of start planet 0 is 0.24747070147434883 years
Orbital period of destination planet 6 is 0.8067538643056904 years
Found launch time 91
With tolerance [0.00511095 0.00897038 0.00699234 0.00824632]
At time [0.02464116 0.24914954 0.38604489 0.73102117] years
Angle is 70.99142476279717
Planet 0 position [0.36550444 0.01480869]
v_at_planet: a 0.5612102290410104, P 0.48857271841706096, R [0.36550444 0.01480869], R_norm 0.3658043067868395
v is 10.379815420603816
We need a velocity of 10.379815420603816 AU/yr to achieve the Hohmann transfer
Planet 0 position [0.36550444 0.01480869]
v_at_planet: a 0.5612102290410104, P 0.48857271841706096, R [0.36550444 0.01480869], R_norm 0.3658043067868395
v is 10.379815420603816
Planet 0 position [0.36550444 0.01480869]
v_from_planet: a 0.35660598633414764, P 0.24747070147434883, R, [0.36550444 0.01480869] R_norm 0.3658043067868395
v is 9.287629735612652
Delta v to enter orbit 1.092185684991163 AU/yr to achieve the Hohmann transfer
Delta v vector to enter orbit [-0.9991802454195715,0.040482553813770444]
======================================
v at perhelion 10.379815420603816 AU/yr
Planet position is [0.36550443703081054,0.014808692534807199]
cos theta = 0.9991802454195715
sin theta = 0.040482553813770375
v theta = [-0.42020143634159757, 10.371306519368773]
Hoth position at 0.24914953763595418 is 0.36550443703081054, 0.014808692534807199
Hoth position at 0.24914953763595418 is 0.21076574649684443, 0.7266675997348173
"""
