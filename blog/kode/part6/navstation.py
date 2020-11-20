"""Egen kode."""

import spacecraft
from nano_motor import nano_motor
from commands import Launch, Boost, Coast


class NavStation:
    """This is the Navigation controller that performs the commands to navigate.

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
                                                     particles, dt, args)

    space_craft = spacecraft.Spacecraft(fuel_mass, rocket_engine,
                                        dt, 33382, args.verbose)

    nav_station.AddCommand(Launch(t, angle, space_craft))

    nav_station.AddCommand(Boost([-1.6712756907642885, -0.25203293684429795], space_craft))

    nav_station.AddCommand(Coast(0.24469793490545463, spacecraft))

    nav_station.Execute()

    cur_time, cur_pos, cur_vel = space_craft.interplanetary_travel.orient()
    print(f'Position {cur_pos} and velocity {cur_vel} at {cur_time}')

    dest_pos = space_craft.planet_position(6, cur_time)
    print(f'Destination pos {dest_pos}')
