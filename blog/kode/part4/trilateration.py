# Egen kode
import numpy as np
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c
from ast2000tools.space_mission import SpaceMission
import orbit_sim as os


class trilateration:
    def __init__(self, seed=33382, dt=1000):
        '''
        Something
        '''
        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)
        self.orbit = os
        self.G = c.G
        self.M = self.system.star_mass
        self.mu = self.G*(self.M + self.system.masses[0])
        self.dt = dt
        self.axes = self.system.semi_major_axes
        self.times, self.planet_pos = np.load('planet_trajectories.npz', allow_pickle=True)
        self.mission = self.mission.load('part1.bin', self.mission)

        self.r_numerical = []
        self.v_numerical = []

        self.just_calculation_orbit()

    def leapfrog(self, r0, v0, T, dt):
        '''
        Leapfrog integration
        couldn't use the one from orbit sim for
        some reason
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

            r, v, a, t, = self.leapfrog(r0, v0, self.T, self.dt)
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

    def circles_old(self, radii, a, b):
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
            Runs through all the planets
            '''
            if a == b and b == 0:
                x = radii[i]*np.cos(theta)
                y = radii[i]*np.sin(theta)
            elif b == 0:
                x = radii[i]*np.cos(theta) + a[i]
                y = radii[i]*np.sin(theta)
            elif a == 0:
                x = radii[i]*np.cos(theta)
                y = radii[i]*np.sin(theta) + b[i]
            else:
                x = (a[i]**2 + b[i]**2 - 2*b[i]*radii[i]*np.sin(theta))/(2*a[i])
                y = (a[i]**2 + b[i]**2 - 2*a[i]*radii[i]*np.cos(theta))/(2*b[i])
            vec[i, 0] = x
            vec[i, 1] = y
        self.vec = vec

    def circle(self, radii, a, b):
        '''
        New and improved circles
        radii is a list of distances to objects
        a[i] is that objects x coordinate
        b[i] is that objects y coordinate
        '''
        theta = np.linspace(0, 2*np.pi, 1000)
        vec = np.zeros(2, float)

        #for i in range(N):
        '''
        Runs through all the planets
        '''
        x = radii*np.cos(theta) + a
        y = radii*np.sin(theta) + b
        vec[0] = x
        vec[1] = y
        return vec

    def same_coordinates(self, radii, a, b):
        '''
        Checks if the coordinates are correct
        '''
        correct = []
        N = len(a)
        vec = np.zeros((N, 2), float)
        for i in range(N+1):
            '''
            Finds the correct x and y coordinates
            '''
            vec[i, :] = self.circle(radii, a, b)
            if i > 0:
                correct.append(float(self.comparrison(self.vec[i], self.vec[i+1])))

        return correct

    def tri_test(self):
        '''
        Test for t=0 at home planet
        '''
        radii = np.array(self.mission.measure_distances())
        times, planet_pos = self.times, self.planet_pos
        print(np.shape(planet_pos))
        for i in range(10):
            correct = self.same_coordinates(radii[i], planet_pos[0, i, 0], planet_pos[1, i, 0])
        return correct

    def radial_velocity(self):
        '''
        Calculates radial_velocity from doppler
        '''
        reference_wavelength = self.mission.reference_wavelength
        #reference_wavelength = 656.3E-9
        lambda_1, lambda_2 = self.mission.star_doppler_shifts_at_sun
        phi_1, phi_2 = self.mission.star_direction_angles

        radial_velocity = c.c*(lambda_1 + lambda_2)/reference_wavelength
        conts = 1/np.sin(phi_2 - phi_1)
        transformation = conts*np.array([[np.sin(phi_2), np.sin(phi_1)],
                                         [np.cos(phi_2), np.cos(phi_1)]])

        rocket_vel = np.dot(transformation, radial_velocity)
        return rocket_vel


tri = trilateration()
vel = tri.radial_velocity()
print(vel)
