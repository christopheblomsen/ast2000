# Egen kode
import numpy as np
import math
import random
from ast2000tools import constants
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import sys


class nano_motor:

    #
    # Initialze a motor
    # L side length in meters
    # N number of particles
    # T temperature in Kelvin
    # dt delta time in seconds for every step when calculating particle movement
    #
    def __init__(self,L,N,T,dt):
        self.L = L # meters
        self.N = N
        self.mu = 0
        self.T = T
        self.dt = dt
        self.sigma = math.sqrt((constants.k_B*T)/constants.m_H2)

        np.random.seed(90207617)
        # The x,y,z position of the particles
        self.pos = np.zeros((3,N))

        # The velocity in x,y,z direction of the particles
        self.v = np.zeros((3,N))

        # The exit box length, start and end points for x and y axis
        self.l = self.L/4
        self.start = (self.L - self.l)/2
        self.end = self.start + self.l

        # Pick a random distribution for the position
        self.pos[0] = np.random.uniform(0,L,N)
        self.pos[1] = np.random.uniform(0,L,N)
        self.pos[2] = np.random.uniform(0,L,N)

        # Pick a normal distribution for the velocity
        self.v[0] = np.random.normal(self.mu,self.sigma,N)
        self.v[1] = np.random.normal(self.mu,self.sigma,N)
        self.v[2] = np.random.normal(self.mu,self.sigma,N)

        # Accumulated momentum
        self.p = 0

        # Count number of steps
        self.steps = 0

        # Number of particles that have left the builing
        self.np = 0

    # Move all particles one time step
    def step(self):
        try:
            self.pos[0] += [self.dt*self.detect_collision(0,n) for n in range(self.N)]
            self.pos[1] += [self.dt*self.detect_collision(1,n) for n in range(self.N)]
            self.pos[2] += [self.dt*self.detect_collision(2,n) for n in range(self.N)]
            self.steps = self.steps + 1
        #print(f'particle 0: ({self.pos[0,0]},{self.pos[1,0]},{self.pos[2,0]})')
        finally:
            return

    # Create a new particle at position n in the arrays
    # The particle is created randomly at the "top" of the box
    # With a volcity pointing in some normally distributed direction downwards
    def fill_new_particle(self,n):
        self.pos[0,n] = np.random.uniform(0,self.L,1)
        self.pos[0,n] = np.random.uniform(0,self.L,1)
        self.pos[2,n] = self.L # All particles enters from top
        self.v[0,n] = np.random.normal(self.mu,self.sigma)
        self.v[1,n] = np.random.normal(self.mu,self.sigma)
        # We assume the velocity in the z axis is negative because it is entering from "above"
        self.v[2,n] = abs(np.random.normal(self.mu,self.sigma))*-1

    # Detect if particle goes beyond the box
    def detect_collision(self,i, n):
        if self.pos[i,n] >= self.L or self.pos[i,n] <= 0:
            #print(f'COLLISION for particle {n} on axis {i}')
            #print(f'particle {n}: ({self.pos[0,n]},{self.pos[1,n]},{self.pos[2,n]})')
            if(i == 2 and self.pos[i,n] <= 0 and self.detect_exit(n)):
                self.fill_new_particle(n)
                self.np = self.np + 1
                #self.Trhust = self.Thrust +
                #sys.exit()
            self.v[i,n] = self.v[i,n] * -1
        return self.v[i,n]

    # Returns True if particle has passed through exit area
    def detect_exit(self,n):
        #print(f'Particle position for {n} ({self.pos[0,n]},{self.pos[1,n]},{self.pos[2,n]}) ({self.start},{self.end})')
        if self.pos[0,n] >= self.start and self.pos[0,n] <= self.end:
            if self.pos[1,n] >= self.start and self.pos[1,n] <= self.end:
                # Since exit is in the negative direction on the z axis we take the absolute value of velocity
                self.p = self.p + abs(self.v[2,n])*constants.m_H2 #Calculate and update momentum.
                return True
        return False

    def calculated_thrust(self):
        # Momentum / time
        return self.p/(self.dt*self.steps)
    #
    # Plots velocity histogram for given direction:
    # 0 = x
    # 1 = y
    # 2 = z
    #
    def plot_velocity(self,i,label):
        plt.hist(self.v[i,:],bins=50,density=True,label=label)

    #
    # Plots the position of the particles
    #
    def plot_position(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim3d(0,self.L*1.1)
        ax.set_ylim3d(0,self.L*1.1)
        ax.set_zlim3d(0,self.L*1.1)
        ax.set_xlabel('X meters')
        ax.set_ylabel('Y meters')
        ax.set_zlabel('Z meters')
        ax.scatter(self.pos[0], self.pos[1], self.pos[2], marker='o')

        self.plot_box(ax)

        self.plot_exit(ax)

    # Plots the exit area as 25% of the side area
    def plot_exit(self,ax):
        ax.plot([self.start,self.end] ,[self.start,self.start],[0,0],color='red')
        ax.plot([self.start,self.start] ,[self.start,self.end],[0,0],color='red')
        ax.plot([self.start,self.end] ,[self.end,self.end],[0,0],color='red')
        ax.plot([self.end,self.end] ,[self.start,self.end],[0,0],color='red')

    # Plots the position of the box
    def plot_box(self,ax):
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

if __name__ == "__main__":
    motor = nano_motor(10**-6,10**5,3000,1e-12)

    motor.plot_velocity(0,'V_x')
    motor.plot_velocity(1,'V_y')
    motor.plot_velocity(2,'V_z')
    #motor.plot_position()
    for i in range(1000):
        motor.step()
    #motor.plot_position()
    plt.legend()
    plt.show()
    print(f'Thrust {motor.calculated_thrust()} N from {motor.np} particles')
