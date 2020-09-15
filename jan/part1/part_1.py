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
import nano_motor as nm
import part1_gasbox as gasbox
import numpy as np
import math

def simulate_engine_performance(number_of_particles_in_box,
                                box_side_length,
                                temperature,
                                box_hole_area):
    """
    Here you can implement the particle simulation for challenges B
    and C of Part 1 of the project.
    """

    #return thrust_per_box, \
    #       mass_loss_rate_per_box
    return gasbox.simulate(number_of_particles_in_box,box_side_length,temperature,box_hole_area)

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
    force_needed = (mission.spacecraft_mass + initial_fuel_mass) * target_delta_v
    # calculate the time we need to run the engine to produce the
    # thrust we need to change the velocity
    time_needed = force_needed/rocket_thrust
    return time_needed * rocket_mass_loss_rate
    #return fuel_mass_used


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
    # Calculate escape velocity
    v_esc = np.sqrt((2*const.G*mission.system.masses[0]*const.m_sun)/(mission.system.radii[0]*1000**3))

    print(f'Escape velocity {v_esc} km/s')
    print(f'Thrust {rocket_thrust} N')
    v = 0
    x = 0
    dt = 1

    r = mission.system.radii[0]*1000
    spc_mass = mission.spacecraft_mass + initial_fuel_mass
    loss_rate = rocket_mass_loss_rate*dt

    rp_sec = mission.system.rotational_periods[0]*86400 #rotational period in sec
    rv = 2*np.pi*r/rp_sec # Rotational velocity

    for i in range(86400):
        F_G = (const.G*mission.system.masses[0]*spc_mass)/((r+x)**2)
        a = (rocket_thrust - F_G)/spc_mass

        v_e = np.sqrt(v**2+rv**2)/1000
        if v_e >= v_esc:
            launch_duration = dt*i+12
            break

        v = v + a*dt
        x = x + v*dt
        spc_mass -= loss_rate*dt

    final_height_above_surface = x
    final_upward_speed = v
    fuel_mass_after_launch = mission.spacecraft_mass + initial_fuel_mass - spc_mass

    print(f'Launch duration is {launch_duration} s')
    print(f'Velocity {v/1000} km/s')
    print(f'Height {x/1000} km')
    print(f'Fuel consumed {initial_fuel_mass - fuel_mass_after_launch}')
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
    r_m = mission.system.radii[0]*1000 # Planet radius in m
    rp_sec = mission.system.rotational_periods[0]*86400 #rotational period in sec

    # Orbital Velocity
    ov = mission.system.initial_velocities[:,0]
    ov_ms = (ov*const.AU)/const.yr
    p_vr = 2*np.pi*r_m/rp_sec # Rotational velocity
    p_vo = (np.sqrt(ov_ms[0]**2+ov_ms[1]**2)*launch_duration)/const.AU # Orbital velocity


    start_x = mission.system.initial_positions[0,0]+(r_m/const.AU)
    end_x = start_x + math.sqrt(final_height_above_surface**2+p_vr**2+p_vo**2)/const.AU
    start_y = 0
    end_y = start_y + (p_vo*launch_duration)/const.AU

    rocket_position_before_launch = [start_x,start_y]
    rocket_position_after_launch = [end_x,end_y]
    rocket_velocity_after_launch = [0,0]
    time_after_launch = launch_duration
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
    rocket_thrust = thrust_per_box * number_of_boxes
    rocket_mass_loss_rate = mass_loss_rate_per_box * number_of_boxes

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
    seed = utils.get_seed('janman')
    mission = SpaceMission(seed)

    # Extract associated SolarSystem object
    system = mission.system

    number_of_boxes = 11.005**15
    number_of_particles_in_box = 10**5
    box_side_length = 10**-6
    temperature = 3000
    box_hole_area = (box_side_length/4)**2
    initial_fuel_mass = 50000
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
                                  launch_duration+500,
                                  rocket_position_before_launch,
                                  time_of_launch)
    mission.launch_rocket()

    # Verify simulated launch results
    mission.verify_launch_result(rocket_position_after_launch)
