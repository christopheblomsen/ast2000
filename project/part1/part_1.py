# -*- coding: utf-8 -*-
"""
Template for the main solution code in Part 1 of the project.

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


L = 10**(-6)

def simulate_engine_performance(number_of_particles_in_box,
                                box_side_length,
                                temperature,
                                box_hole_area):
    """
    Here you can implement the particle simulation for challenges B
    and C of Part 1 of the project.
    """

    # Insert awesome code here...

    # You will probably also need these constants:
    # const.k_B
    # const.m_H2

    return thrust_per_box, \
           mass_loss_rate_per_box


def compute_fuel_mass_needed_for_boost(mission,
                                       rocket_thrust,
                                       rocket_mass_loss_rate,
                                       initial_fuel_mass,
                                       target_delta_v):
    """
    Here you can implement the function for challenge D of Part 1
    of the project.
    """

    # Insert awesome code here...

    # You will probably also need this quantity:
    # mission.spacecraft_mass

    return fuel_mass_used


def simulate_rocket_launch(mission,
                           rocket_thrust,
                           rocket_mass_loss_rate,
                           initial_fuel_mass):
    """
    Here you can implement the rocket launch simulation for challenge E
    of Part 1 of the project.
    """

    # Insert awesome code here...

    # You will probably also need these quantities:
    # const.G
    # const.m_sun
    # mission.spacecraft_mass
    # system.masses
    # system.radii

    return final_height_above_surface, \
           final_upward_speed,         \
           fuel_mass_after_launch,     \
           launch_duration


def convert_to_solar_system_frame_simple(system,
                                         final_height_above_surface,
                                         final_upward_speed,
                                         launch_duration):
    """
    Here you can implement the coordinate system conversion for
    challenge F of Part 1 of the project.
    """

    # Insert awesome code here...

    # You will probably also need these quantities:
    # const.AU
    # const.day
    # const.yr
    # system.radii
    # system.initial_positions
    # system.initial_velocities
    # system.rotational_periods

    # Return quantities in the solar system units and frame of reference
    return rocket_position_before_launch, \
           rocket_position_after_launch,  \
           rocket_velocity_after_launch,  \
           time_after_launch


def run_engine_and_launch_simulation(mission,
                                     number_of_boxes,
                                     number_of_particles_in_box,
                                     box_side_length,
                                     temperature,
                                     box_hole_area,
                                     initial_fuel_mass):
    """
    This function executes the simulation code for the rocket engine and launch,
    and returns the results.
    """

    # Simulate engine
    thrust_per_box,        \
    mass_loss_rate_per_box \
      = simulate_engine_performance(number_of_particles_in_box,
                                    box_side_length,
                                    temperature,
                                    box_hole_area)

    # Compute rocket_thrust and rocket_mass_loss_rate here...

    # Simulate launch
    final_height_above_surface, \
    final_upward_speed,         \
    fuel_mass_after_launch,     \
    launch_duration             \
      = simulate_rocket_launch(mission,
                               rocket_thrust,
                               rocket_mass_loss_rate,
                               initial_fuel_mass)

    return thrust_per_box,             \
           mass_loss_rate_per_box,     \
           final_height_above_surface, \
           final_upward_speed,         \
           fuel_mass_after_launch,     \
           launch_duration


# Prevent the following code from executing when calling `import part_1`
if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available
    utils.check_for_newer_version()

    # Construct SpaceMission instance for my mission
    seed = utils.get_seed('my_username')
    mission = SpaceMission(seed)

    # Extract associated SolarSystem object
    system = mission.system

    # Run engine and launch simulations and get results relative to planet
    thrust_per_box,             \
    mass_loss_rate_per_box,     \
    final_height_above_surface, \
    final_upward_speed,         \
    fuel_mass_after_launch,     \
    launch_duration             \
      = run_engine_and_launch_simulation(mission,
                                         number_of_boxes,
                                         number_of_particles_in_box,
                                         box_side_length,
                                         temperature,
                                         box_hole_area,
                                         initial_fuel_mass)

    # Convert simulated launch results to the solar system frame,
    # assuming launch along the x-direction at time = 0 years
    rocket_position_before_launch, \
    rocket_position_after_launch,  \
    rocket_velocity_after_launch,  \
    time_after_launch              \
      = convert_to_solar_system_frame_simple(system,
                                             final_height_above_surface,
                                             final_upward_speed,
                                             launch_duration)
    time_of_launch = 0.0

    # Perform real launch
    mission.set_launch_parameters(thrust_per_box*number_of_boxes,
                                  mass_loss_rate_per_box*number_of_boxes,
                                  initial_fuel_mass,
                                  launch_duration,
                                  rocket_position_before_launch,
                                  time_of_launch)
    mission.launch_rocket()

    # Verify simulated launch results
    mission.verify_launch_result(rocket_position_after_launch)