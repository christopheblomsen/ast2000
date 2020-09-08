# egen kode
import nano_motor as nm
from ast2000tools.space_mission import SpaceMission
import matplotlib.pyplot as plt
try:
   import cPickle as pickle
except:
   import pickle

import os
import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--run-steps', help='foo help')
args = parser.parse_args()

print(args.run_steps)

class rocket_engine:

    # Initiates the rocket engine with N nano_motors of type motor
    def __init__(self,N, motor):
        # Number of nano_motors that makes up the engine
        self.N = N

        # Save a reference to the engine prototype. We don't need N of these
        # since we can just multiply the results by N
        self.motor = motor

    # Returns the thrust of the engine in Newtons
    def thrust(self):
        return self.motor.calculated_thrust() * self.N

    # Returns fuel consumption in kg/s
    def fuel_consumption(self):
        return (self.N * self.motor.fuel_consumption())

    # Boost the rocket of spacecraft_mass in kg with delta v m/s
    # Returns fuel consumed in kg
    def boost(self, spacecraft_mass, dv):
        force_needed = spacecraft_mass * dv
        # calculate the time we need to run the engine to produce the
        # thrust we need to change the velocity
        self.time_needed = force_needed/self.thrust()
        return self.time_needed * self.fuel_consumption()
if __name__ == "__main__":
    dt = 1e-12
    fuel_mass = 2500
    dv = 5
    steps = 1000
    motors = 10**13
    filename = 'nano_motor.pkl'
    mission = SpaceMission(33382)

    # Check if previously saved file of nano_motor exists
    # If not simulate motor performance and save it to file
    # Run with paramterer --run_steps true if you want to simulate again
    if(os.path.exists(filename) == False or args.run_steps == 'true'):
        motor = nm.nano_motor(10**-6,10**5,3000,dt)

        for i in range(steps):
            print(f'{i:4d}\b\b\b\b', end = '',flush = True)
            motor.step()

        with open('nano_motor.pkl', 'wb') as output:
            pickle.dump(motor, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as input:
            motor = pickle.load(input)


    print('Single motor performance:')
    print(motor)

    engine = rocket_engine(motors,motor)
    fuel_consumed = engine.boost(mission.spacecraft_mass+fuel_mass,dv)

    print(f'ENIGNE with {motors:g} motors:')
    print(f'          Thrust {engine.thrust():g} N')
    print(f'Fuel consumption {engine.fuel_consumption():g} kg/s')
    print(f'Fuel consumed after boost to {dv} m/s in {engine.time_needed} sec is {fuel_consumed:.2f} kg')

    #motor.plot_absolute_velocity('Velocity')
    #motor.plot_velocity(1,'Vy')
    #motor.plot_velocity(2,'Vz')
    #plt.legend()
    #plt.show()
