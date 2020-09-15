# -*- coding: utf-8 -*-
"""
Template for the main solution code in Part 2 of the project.

The objective of the template files is to give you an idea of
the most important functions that you have to implement, what
input they will need and what output they should produce.

To make things work in practice you will have to add more
functionalities than the ones outlined, and you might have
to adapt the interfaces to fit your specific approach.

You are of course free (and encouraged) to structure the code
in any way you want!
"""
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission


def simulate_orbits(system, duration, n_time_steps_per_year):
    """
    Here you can implement the orbital simulation for challenge A
    of Part 2 of the project.
    """

    # Insert awesome code here...

    # You will probably also need these quantities:
    # const.G_sol
    # system.number_of_planets
    # system.initial_positions
    # system.initial_velocities
    # system.star_mass

    return times, planet_positions


# Prevent the following code from executing when calling `import part_2`
if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available
    utils.check_for_newer_version()

    # Construct SpaceMission instance for my mission
    seed = utils.get_seed('my_username')
    mission = SpaceMission(seed)

    # Extract associated SolarSystem object
    system = mission.system

    # Simulate orbits
    times, planet_positions = simulate_orbits(system,
                                              duration,
                                              n_time_steps_per_year)

    # Verify simulated orbits and write exact orbits to planet_trajectories.npy
    mission.verify_planet_positions(duration,
                                    planet_positions,
                                    filename='planet_trajectories.npy')

    # Generate orbit video to view in SSView
    mission.generate_orbit_video(times,
                                 planet_positions,
                                 filename='orbit_video.xml')
