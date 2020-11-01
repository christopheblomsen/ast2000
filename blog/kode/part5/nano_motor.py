"""Egen kode."""
import numpy as np
import math
import random
from ast2000tools import constants
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import sys
import os
try:
   import cPickle as pickle
except:
   import pickle


class nano_motor:
    """Class to simulate a gas driven rocket motor."""

    def __init__(self, L, N, T, dt):
        """
        Initialze a motor.

        L side length in meters
        N number of particles
        T temperature in Kelvin
        dt delta time in seconds for every step when calculating particle movement
        """
        self.L = L  # meters
        self.N = int(N)
        self.mu = 0
        self.T = T
        self.dt = dt
        self.sigma = math.sqrt((constants.k_B*self.T)/constants.m_H2)

        np.random.seed(90207617)
        # The x,y,z position of the particles
        self.pos = np.zeros((3, self.N))

        # The velocity in x,y,z direction of the particles
        self.v = np.zeros((3, self.N))

        # The exit box length, start and end points for x and y axis
        self.l = self.L/4
        self.start = (self.L - self.l)/2
        self.end = self.start + self.l

        # Pick a random distribution for the position
        self.pos[0] = np.random.uniform(0, L, self.N)
        self.pos[1] = np.random.uniform(0, L, self.N)
        self.pos[2] = np.random.uniform(0, L, self.N)

        # Pick a normal distribution for the velocity
        self.v[0] = np.random.normal(self.mu, self.sigma, self.N)
        self.v[1] = np.random.normal(self.mu, self.sigma, self.N)
        self.v[2] = np.random.normal(self.mu, self.sigma, self.N)

        # Accumulated momentum
        self.p = 0

        # Count number of steps
        self.steps = 0

        # Number of particles that have left the builing
        self.np = 0

    def step(self):
        """Move all particles one time step."""
        try:
            self.pos[0] += [self.dt*self.detect_collision(0, n) for n in range(self.N)]
            self.pos[1] += [self.dt*self.detect_collision(1, n) for n in range(self.N)]
            self.pos[2] += [self.dt*self.detect_collision(2, n) for n in range(self.N)]
            self.steps += 1
        finally:
            return

    def fill_new_particle(self, n):
        """
        Create a new particle at position n in the arrays.

        The particle is created randomly at the "top" of the box
        With a volcity pointing in some normally distributed direction downwards
        """
        self.pos[0, n] = np.random.uniform(0, self.L, 1)
        self.pos[0, n] = np.random.uniform(0, self.L, 1)
        self.pos[2, n] = self.L  # All particles enters from top
        self.v[0, n] = np.random.normal(self.mu, self.sigma)
        self.v[1, n] = np.random.normal(self.mu, self.sigma)
        # We assume the velocity in the z axis is negative because it is entering from "above"
        self.v[2, n] = abs(np.random.normal(self.mu, self.sigma))*-1

    def detect_collision(self, i, n):
        """Detect if particle goes beyond the box and returns updated velocity."""
        if self.pos[i, n] >= self.L or self.pos[i, n] <= 0:
            if(i == 2 and self.pos[i, n] <= 0 and self.detect_exit(n)):
                self.fill_new_particle(n)
                self.np = self.np + 1
            self.v[i, n] = self.v[i, n] * -1
        return self.v[i, n]

    def detect_exit(self, n):
        """Return True if particle has passed through exit area."""
        if self.pos[0, n] >= self.start and self.pos[0, n] <= self.end:
            if self.pos[1, n] >= self.start and self.pos[1, n] <= self.end:
                # Since exit is in the negative direction on the z axis we take the absolute value of velocity
                self.p += (abs(self.v[2, n])*constants.m_H2)  # Calculate and update momentum p = m*v
                return True
        return False

    #
    def fuel_consumption(self):
        """
        Return the fuel consumption in kg/s.

        n particles with mass / time we have run the simulation
        """
        return (self.np * constants.m_H2)/(self.dt*self.steps)

    #
    def calculated_thrust(self):
        """
        Thrust = P/dt in N where P is momentum and dt is time.

        Momentum / time
        """
        return self.p/(self.dt*self.steps)


    def plot_velocity(self, i, label, color='blue'):
        """
        Plot velocity histogram for given direction.

        0 = x
        1 = y
        2 = z
        """
        plt.hist(self.v[i, :], bins=50, density=True, label=label, color=color)
        plt.xlabel('Velocity m/s')
        plt.ylabel('Probability')

    def plot_absolute_velocity(self, label):
        """Plot the absollute velocity."""
        abs_v = np.sqrt(self.v[0, :]**2+self.v[1, :]**2+self.v[2, :]**2)
        plt.hist(abs_v, bins=50, density=True, label=label)
        plt.xlabel('Velocity m/s')
        plt.ylabel('Probability')

    def plot_position(self):
        """Plot the position of the particles."""
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim3d(0, self.L*1.1)
        ax.set_ylim3d(0, self.L*1.1)
        ax.set_zlim3d(0, self.L*1.1)
        ax.set_xlabel('X meters')
        ax.set_ylabel('Y meters')
        ax.set_zlabel('Z meters')
        ax.scatter(self.pos[0], self.pos[1], self.pos[2], marker='o')

        self.plot_box(ax)

        self.plot_exit(ax)

    def plot_exit(self,ax):
        """Plot the exit area as 25% of the side area."""
        ax.plot([self.start, self.end], [self.start, self.start], [0, 0], color='red')
        ax.plot([self.start, self.start], [self.start, self.end], [0, 0], color='red')
        ax.plot([self.start, self.end], [self.end, self.end], [0, 0], color='red')
        ax.plot([self.end, self.end], [self.start, self.end], [0, 0], color='red')


    def plot_box(self, ax):
        """Plot the position of the box."""
        # Plot the box
        ax.plot([0,self.L] ,[0,0],[0,0],color='black')
        ax.plot([0,0] ,[0,self.L],[0,0],color='black')
        ax.plot([self.L,self.L] ,[0,self.L],[0,0],color='black')
        ax.plot([0,self.L] ,[self.L,self.L],[0,0],color='black')

        ax.plot([0,0] ,[0,0],[0,self.L],color='black')
        ax.plot([self.L,self.L] ,[0,0],[0,self.L],color='black')
        ax.plot([self.L,self.L] ,[self.L,0],[self.L,self.L],color='black')

        ax.plot([self.L,self.L] ,[self.L,self.L],[0,self.L],color='black')
        ax.plot([0,self.L] ,[0,0],[self.L,self.L],color='black')
        ax.plot([0,0] ,[0,self.L],[self.L,self.L],color='black')
        ax.plot([0,self.L] ,[self.L,self.L],[self.L,self.L],color='black')
        ax.plot([0,0] ,[self.L,self.L],[self.L,0],color='black')

    def __str__(self):
        """Print the object."""
        if(self.steps == 0):
            return 'No engine simulation has been run yet. Run the step method n times to produce results'

        line1 = f'Accumulated momentum {self.p} kg m/s from {self.np} particles in {self.dt*self.steps} sec\n'
        line2 = f'Mass loss rate is {self.np*constants.m_H2/(self.dt*self.steps)} kg/s\n'
        line3 = f'Thrust is {self.calculated_thrust()} N\n'
        line4 = f'There are {self.np/(self.dt*self.steps):g} exiting the box pr second'
        return line1 + line2 + line3 + line4

    def __repr__(self):
        """Print the object."""
        return self.__str__()


if __name__ == "__main__":
    dt = 1e-12  # 15**-12
    steps = 1000
    N = [8**5, 9**5, 9.8**5, 11**5, 12**5]  # Number of particles to simulate
    T = [2800, 3000, 3200, 3500, 3600]  # Temperature
    L = 10**-6  # Box side length in m
    filename = 'nano_motor.pkl'

    # if(os.path.exists(filename) is False):
    '''    motor = nano_motor(L, N, T, dt)
        print('Arguments')
        print('-------------------------------------')
        print(f'Temperature: {T} K')
        print(f'# particles: {N}')
        print(f' Box length: {L} cm')
        print(f'         dt: {dt} sec')
        print('-------------------------------------')
        print(f'Running engine particle simulation for {steps} steps')
    '''
    mass_loss = []
    thrust = []
#    for i in range(5):
#        for j in range(5):
    i = 2
    j = 3
    print(f'{i},{j}# particles {N[i]}, temp {T[j]} K')
    motor = nano_motor(L, N[i], T[j], dt)
    for k in range(steps):
        print(f'{k:4d}\b\b\b\b', end='', flush=True)
        motor.step()
    print(f'Trhust: {motor.calculated_thrust()}')
    print(f'Mass loss {motor.fuel_consumption()}')
    thrust.append(motor.calculated_thrust())
    mass_loss.append(motor.fuel_consumption())

    print(f'N: {N}')
    print(f'T: {T}')
    print(f'Thrust: {thrust}')
    print(f'Mass loss: {mass_loss}')

    with open('nano_motor.pkl', 'wb') as output:
        pickle.dump(motor, output, pickle.HIGHEST_PROTOCOL)
#    else:
#        with open(filename, 'rb') as input:
#            motor = pickle.load(input)

    print(motor)

"""
Thrust:
[8.95302246034942e-11, 9.267258677858525e-11, 9.571183606491261e-11, 9.719582899619638e-11, 1.3555623342856187e-10,
5.36991093045451e-11, 7.845256015628057e-11, 8.102545572069602e-11, 8.228173925350814e-11, 8.594047378325178e-11,
3.814410691054312e-11, 3.948290170625775e-11, 4.077776553820819e-11, 4.141001666095551e-11, 4.3251351800553175e-11,
1.631220062463307e-10, 1.688473177221318e-10, 1.7438475989056238e-10, 1.770885608167686e-10, 1.8496296938130674e-10,
3.5188446897624916e-10, 3.6423502936197395e-10, 4.0793887469565186e-10, 4.1426388559643146e-10, 4.707979262565577e-10]
Mass loss:
[2.1715952686146868e-14, 2.1715952686146868e-14, 2.1715952686146868e-14, 2.1715952686146868e-14, 2.6059143223376244e-14,
 1.7372762148917496e-14, 2.1715952686146868e-14, 2.1715952686146868e-14, 2.1715952686146868e-14, 2.1715952686146868e-14,
 1.3029571611688122e-14, 1.3029571611688122e-14, 1.3029571611688122e-14, 1.3029571611688122e-14, 1.3029571611688122e-14,
 3.908871483506437e-14, 3.908871483506437e-14, 3.908871483506437e-14, 3.908871483506437e-14, 3.908871483506437e-14,
9.120700128181685e-14, 9.120700128181685e-14, 9.989338235627559e-14, 9.989338235627559e-14, 1.0423657289350497e-13]
A
"""
