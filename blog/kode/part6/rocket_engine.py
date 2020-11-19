"""Egen kode."""
import nano_motor
from ast2000tools.space_mission import SpaceMission

try:
    import cPickle as pickle
except Exception:
    import pickle

import os


class RocketEngine:
    """The rockcet engine to take us of this planet."""

    def __init__(self, N, motor):
        """Initiate the rocket engine with N nano_motors of type motor."""
        self.N = N

        # Save a reference to the engine prototype. We don't need N of these
        # since we can just multiply the results by N
        self.motor = motor

    def thrust(self):
        """Return the thrust of the engine in Newtons."""
        return self.motor.calculated_thrust() * self.N

    def fuel_consumption(self):
        """Return fuel consumption in kg/s."""
        return (self.N * self.motor.fuel_consumption())

    def boost(self, spacecraft_mass, dv):
        """Boost the rocket of spacecraft_mass in kg with delta v m/s.

        Return fuel consumed in kg
        """
        force_needed = spacecraft_mass * dv
        # calculate the time we need to run the engine to produce the
        # thrust we need to change the velocity
        self.time_needed = force_needed/self.thrust()
        return self.time_needed * self.fuel_consumption()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-steps', help='foo help')
    args = parser.parse_args()

    print(args.run_steps)

    dt = 1e-12
    fuel_mass = 2500
    dv = 5
    steps = 1000

    motors = 11.0699895**15
    filename = 'nano_motor.pkl'
    mission = SpaceMission(33382)
    N = 10**5  # Number of particles to simulate
    T = 3500  # Temperature
    L = 10**-6  # Box side length in m
    # Check if previously saved file of nano_motor exists
    # If not simulate motor performance and save it to file
    # Run with paramterer --run_steps true if you want to simulate again
    if(os.path.exists(filename) is False or args.run_steps == 'true'):
        motor = nano_motor.nano_motor(L, N, T, dt)

        for i in range(steps):
            print(f'{i:4d}\b\b\b\b', end='', flush=True)
            motor.step()

        with open(filename, 'wb') as output:
            pickle.dump(motor, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as input:
            motor = pickle.load(input)

    print('Single motor performance:')
    print(motor)

    engine = RocketEngine(motors, motor)
    fuel_consumed = engine.boost(mission.spacecraft_mass+fuel_mass, dv)

    print(f'ENGINE with {motors:g} motors:')
    print(f'          Thrust {engine.thrust():g} N')
    print(f'Fuel consumption {engine.fuel_consumption():g} kg/s')
    print(f'Fuel consumed after boost to {dv} m/s in {engine.time_needed} sec \
    is {fuel_consumed:.2f} kg')