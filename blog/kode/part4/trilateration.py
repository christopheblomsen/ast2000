# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c
from ast2000tools.space_mission import SpaceMission
import os as operatingsystem
import orbit_sim as os



class trilateration:
    def __init__(self, mission, seed=33382, dt=1000):
        '''
        Something
        '''
        self.system = SolarSystem(seed)
        self.mission = mission
        self.orbit = os
        self.G = c.G
        self.M = self.system.star_mass
        self.mu = self.G*(self.M + self.system.masses[0])
        self.dt = dt
        self.axes = self.system.semi_major_axes
        self.times, self.planet_pos = np.load('planet_trajectories.npz', allow_pickle=True)

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
        N = vec1.shape[1]
        mindiff = (tol,tol)
        minindex = 0
        for i in range(N):
            for j in range(N):
                diff_x = abs(vec1[0,i] - vec2[0,j])
                diff_y = abs(vec1[1,i] - vec2[1,j])
                if diff_x < tol and diff_y < tol and (mindiff[0] > diff_x and mindiff[1] > diff_y):
                    mindiff =(diff_x,diff_y)
                    minindex = i

        return vec1[:,minindex]
        #diff = np.abs(vec1[0,:] - vec2[:,1])
        #closest = diff.argmin()
        #print(f'Closest index {closest} with value {vec1[:,closest]}')
        #return vec1[closest]
        '''
        correct = []
        print(f'vec1 size {vec1.shape[1]}')
        for i in range(vec1.shape[1]):
            tolerance_x = abs(vec1[0,i] - vec2[0,i])
            tolerance_y = abs(vec1[1,i] - vec2[1,i])
            if tolerance_x < tol and tolerance_y < tol:
                correct.append(vec1[:,i])

        print(f'Candidate coordinates {correct}')
        return np.array(correct)
        '''
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
        N = len(a)
        theta = np.linspace(0, 2*np.pi, 5000)
        vec = np.zeros((2,N,5000), float)

        for i in range(N):
            '''
            Runs through all the planets
            '''
            x = radii[i]*np.cos(theta) + a[i]
            y = radii[i]*np.sin(theta) + b[i]
            vec[0,i] = x
            vec[1,i] = y
            plt.plot(vec[0,i],vec[1,i],label=f'Planet {i}')

        return vec

    def same_coordinates(self, radii, a, b):
        '''
        Checks if the coordinates are correct
        '''
        correct = []
        N = len(a)
        '''
        Finds the correct x and y coordinates
        '''
        vec = self.circle(radii, a, b)
        if(operatingsystem.path.exists('candidates.npy')):
            return np.load('candidates.npy',allow_pickle=True)

        print(f'vec.shape {vec.shape}')
        for i in range(N):
            for j in range(N):
                if j != i:
                    print(f'Comparing {i} and {j}')
                    correct.append(self.comparrison(vec[:,i], vec[:,j],tol=0.1))


        return np.array(correct,dtype=object)

    def find_closest_neighbours(self,x,y,tol):
        neighbours = {}
        for i in range(len(x)):
            neighbours[i] = []
            neighbours[i].append(i)
            for j in range(len(x)):
                if i != j:
                    distance = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
                    if distance < 0.01:
                        neighbours[i].append(j)
                        print(f'Length between {i} and {j} {distance}')

        return neighbours

    def tri_test(self, distances):
        '''
        Test for t=0 at home planet
        '''
        times, planet_pos = self.times, self.planet_pos

        # Find the three clostest planets, excluding our mother planet
        closest = np.argsort(distances)[1:4]
        print(f'Shortest distance from us {closest}')
        positions = np.concatenate((planet_pos[:,:,0],np.zeros((2,1))),axis=1)
        print(np.shape(positions))
        print(f'Planet positions: {positions[:,:]}')
        plt.scatter(positions[0,closest],positions[1,closest],marker='.')

        print(f'Positions {positions[:,closest]}')

        candidates = self.same_coordinates(distances[closest], positions[0,closest],positions[1,closest])

        if(operatingsystem.path.exists('candidates.npy') == False):
            np.save('candidates.npy',candidates)

        print(f'coordinate candidates {candidates.shape} ')
        x = candidates[:,0]
        y = candidates[:,1]
        neighbours = self.find_closest_neighbours(x,y,0.01)
        print(neighbours)
        for neighbour in neighbours:
            if(len(neighbours[neighbour]) == 3):
                print(f'x: {x[neighbours[neighbour]]}')
                print(f'y: {y[neighbours[neighbour]]}')
                centroid = [x[neighbours[neighbour]].sum()/3,y[neighbours[neighbour]].sum()/3]
        plt.scatter(candidates[:,0],candidates[:,1],c='r')
        plt.scatter(centroid[0],centroid[1],c='g')
        plt.xlabel('AU')
        plt.ylabel('AU')
        plt.legend()
        plt.show()

        return centroid

    def radial_velocity(self):
        '''
        Calculates radial_velocity from doppler
        '''
        #reference_wavelength = self.mission.reference_wavelength
        reference_wavelength = 656.3
        lambda_1_sun, lambda_2_sun = self.mission.star_doppler_shifts_at_sun
        lambda_1_rock, lambda_2_rock = self.mission.measure_star_doppler_shifts()
        phi = self.mission.star_direction_angles
        phi_1 = np.deg2rad(phi[0])
        phi_2 = np.deg2rad(phi[1])

        radial_velocity1 = c.c_AU_pr_yr*(lambda_1_sun - lambda_1_rock)/reference_wavelength
        radial_velocity2 = c.c_AU_pr_yr*(lambda_2_sun - lambda_2_rock)/reference_wavelength

        radial_velocity = np.array([radial_velocity1, radial_velocity2])
        conts = 1/np.sin(phi_2 - phi_1)
        transformation = conts*np.array([[np.sin(phi_2), -np.sin(phi_1)],
                                         [-np.cos(phi_2), np.cos(phi_1)]])

        rocket_vel = np.dot(transformation, radial_velocity)
        #rocket_vel_x = conts*(np.sin(phi_2) * radial_velocity1 - np.sin(phi_1) * radial_velocity2)
        #rocket_vel_y = conts*(-np.cos(phi_1) * radial_velocity1 + np.cos(phi_2) * radial_velocity2)
        #rocket_vel = np.array([rocket_vel_x, rocket_vel_y])
        rocket_vel_AU_y = rocket_vel
        print(f'The rocket velocity is {rocket_vel_AU_y}')
        return rocket_vel_AU_y

if __name__ == '__main__':
    mission = SpaceMission.load('part1.bin')
    tri = trilateration(mission)
    vel = tri.radial_velocity()
    print(vel)
    pos = tri.tri_test()
