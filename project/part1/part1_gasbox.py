import numpy as np
import math
from ast2000tools import constants
'''
N = 10**5               Number of particles
L = 10**-6              Box Length in meters
T = 3000                Temperature in Kelvin
A = (L/4)**2            Box hole area
'''
def simulate(N,L,T,A):
    #CONSTANTS
    k = 1.38064852e-23          #Boltzmann's Constant

    #PHYSICAL VARIABLES


    #
    m = 3.3474472e-27           #Mass of individual particle in kg

    #INITIAL CONDITIONS
    sigma = math.sqrt((constants.k_B*T)/constants.m_H2)      #The standard deviation of our particle velocities

    x =  np.random.uniform(0,L, size = (int(N), 3))                     # Position vector, fill in uniform distribution in all 3 dimensions
    v =  np.random.normal(0,sigma, size = (int(N), 3))                       # Velocity vector, fill in correct distribution
    #print(x)
    #print(v)
    sidemax = np.full_like(x,L) # Array with maximum side length to detect collision
    sidemin = np.zeros_like(x) # Array with minimum side length to detect collision
    exit = np.full_like(x,math.sqrt(A)) # Our exit box
    exit[:,2] = 0 # Set z to 0
    #print(exit)
    '''

    An array with 10 particles (such that N = 10) would look like this:

                          x =  [[x0, y0, z0],
                                [x1, y1, z1],
                                [x2, y2, z2],
                                [x3, y3, z3],
                                [x4, y4, z4],
                                [x5, y5, z5],
                                [x6, y6, z6],
                                [x7, y7, z7],
                                [x8, y8, z8],
                                [x9, y9, z9]]
    '''

    #SIMULATION VARIABLES
    time = 1e-9                       #Simulation Runtime in Seconds
    steps = 1000                 #Number of Steps Taken in Simulation
    dt = time/steps                  #Simulation Step Length in Seconds

    #PREPARING THE SIMULATION
    exiting = 0         #The total number of particles that have exited the gas box
    f = 0                   #Used to calculate Force/Second later in the simulation

    def fill_new_particles(v,x,n):
        x[n] =  np.random.uniform(0,L, size = (len(n), 3))                     # Position vector, fill in uniform distribution in all 3 dimensions
        v[n] =  np.random.normal(0,sigma, size = (len(n), 3))

    #RUNNING THE CONDITIONAL INTEGRATION LOOP
    for i in range(int(steps)):
        '''
        Each loop must begin by allowing each particle to move depending on its velocity
        '''
        x += v*dt #fill in
        #print(f'x {x}')

        '''
        Now you need to check which particles are colliding with a wall and make them bounce:
        you may do this the slow but easy way: (1) looping over particles and using if tests:

        for j in range(int(N)):

        OR YOU CAN check for collisions (2) the fast and elegant way using vectorization with masking:
        Look at some examples of masking in the beginning of the Numpy section in the Numerical Compendium
        in addition to more detailed explanations inside the Numpy section.
        '''
        collision_points = np.greater_equal(x,sidemax)#np.logical_or(np.greater_equal(x,sidemax),np.less_equal(x,sidemin)) #
        collisions_indices = np.where(collision_points == True)
        v[collisions_indices] *= -1
        collision_points = np.less_equal(x,sidemin)#np.logical_or(np.greater_equal(x,sidemax),np.less_equal(x,sidemin)) #
        collisions_indices = np.where(collision_points == True)
        v[collisions_indices] *= -1

        '''
        To check that these conditions are fulfilled for your particles, you can use
        NumPy's logical operators, which return an array of booleans giving an
        elementwise evaluation of these conditions for an entire matrix.  You may
        need:

            (a)     np.logical_or(array1, array2)            or
            (b)     np.logical_and(array1, array2)          and
            (c)     np.less(array1, array2)                   <
            (d)     np.greater(array1, array2)                >
            (e)     np.less_equal(array1, array2)            <=
            (f)     np.greater_equal(array1, array2)         >=

        '''

        '''Now, in the same way as for the collisions, you need to count the particles
        escaping, again you can use the slow way or the fast way.
        For each particle escaping, make sure to replace it!
        Then update the total force from each of the exiting particles:
        '''
        exit_points = np.less_equal(x,exit)
        exit_indices = np.where(exit_points == True)
        exit_unique,counts = np.unique(exit_indices[0],return_counts=True)
        exited = np.where(np.equal(counts,3) == True)
        #
        #if(exit_indices != False):
        '''
        print(f'exit_points\n{exit_points}')
        print(f'exit_unique\n{exit_unique}')
        print(f'counts\n{counts}')
        print(f'x[exit_unique]\n{x[exit_unique]}')
        print(f'exited\n{exited}')

        print(f'Particle indices that have left the building\n{exit_unique[exited]}')
        print(f'x[exit_unique][exited]\n{x[exit_unique[exited]]}')
        '''
        exiting += len(exit_unique[exited])
        '''
        print(f'exiting {exiting}')
        print(f'V_z {v[exit_unique[exited],2]}')
        '''
        f += np.sum((v[exit_unique[exited],2]*constants.m_H2)/dt)
        fill_new_particles(v,x,exit_unique[exited])
        #print(f'f {f}')



    particles_per_second = exiting/time  #The number of particles exiting per second
    mean_force = f/steps                             #The box force averaged over all time steps
    box_mass = particles_per_second*m                     #The total fuel loss per second

    print('There are {:g} particles exiting the gas box per second.'.format(particles_per_second))
    print('The gas box exerts a thrust of {:g} N.'.format(mean_force))
    print('The box has lost a mass of {:g} kg/s.'.format(box_mass))

    return mean_force, box_mass

if __name__ == "__main__":
    N = 10**5              # Number of particles
    L = 10**-6             # Box Length in meters
    T = 3000               # Temperature in Kelvin
    A = (L/4)**2           # Box hole area
    simulate(N,L,T,A)
