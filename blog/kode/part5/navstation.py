"""Egen kode."""
# import numpy as np

# import matplotlib.pyplot as plt
# from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission

import spacecraft
# from astrogation_computer import AstrogationComputer
# from trajection_computer import TrajectionComputer
from commands import Launch


class NavStation:
    """
    This is the Navigation controller that performs the commands to navigate.

    Add commands in the order they should be executed.
    """

    def __init__(self):
        """Initialize the NavStation with an empyt commands list."""
        self.commands = []
        pass

    def AddCommand(self, cmd):
        """Add command to commands list."""
        self.commands.append(cmd)

    def Execute(self):
        """Execute all the commands in command list."""
        for command in self.commands:
            command.execute(None)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=float, default=0.0, help='Time in years')
    parser.add_argument('-la', '--launch-angle', type=float, default=0.0,
                        help='Launch angle along equator')
    parser.add_argument('--run-steps', help='Run engine simulation')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print debug statements')
    args = parser.parse_args()

    mission = SpaceMission(33382)

    steps = 1000
    motors = 11.0699895**15
    max_launch_time = 1000
    rocket_filename = "nano_motor.pkl"
    temperature = 3500
    fuel_mass = 43000  # Fuel mass in kg
    particles = 9.8**5
    t = args.t       # Time in years after t=0 to launch at
    angle = args.launch_angle  # Angle along equator to launch from in radians
    dt = 0.01

    nav_station = NavStation()

    rocket_engine = spacecraft.rocket_engine_factory(rocket_filename, motors,
                                                     temperature, steps,
                                                     particles, args)

    space_craft = spacecraft.Spacecraft(mission, fuel_mass, rocket_engine,
                                        dt, args.verbose)

    nav_station.AddCommand(Launch(0, 0, space_craft, mission))

    nav_station.Execute()
