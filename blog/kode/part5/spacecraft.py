"""Egen kode."""
import numpy as np
import rocket_engine as re
from nano_motor import nano_motor
from astrogation_computer import AstrogationComputer
from trajection_computer import TrajectionComputer
from trilateration import Trilateration
import ast2000tools.constants as c
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
import matplotlib.pyplot as plt
import os
import tools
try:
    import cPickle as pickle
except Exception:
    import pickle


import orbit_sim as orbits


class Spacecraft:
    """Class for launching the rocket_thrust.

    Takes mission, initial_mass of rocket (including fuel) and one nano_enigne
    as arguments to the constructor

    Call launch_process to actually launch the rocket.

    - mission       The mission object from AST2000tools
    - fuel_mass     Fuel mass in kg
    - Engine        The engine object to launch with
    - launch_angle  The angle on planet surface from 0 to 2 pi to launch from
    - times         The time in years at which to launch

    """

    def __init__(self, fuel_mass, engine, dt, seed=33382, verbose=False):
        """Initialize Spacecraft."""
        self.mission = SpaceMission(seed)
        self.fuel_mass = fuel_mass
        self.initial_mass = self.mission.spacecraft_mass + fuel_mass
        self.system = SolarSystem(seed)
        self.engine = engine

        # Number of time steps to try
        self.N = 75000
        self.verbose = verbose
        self.dt = dt

        # Positions at time step i in m. Not necessary to store all unless we
        # want to plot afterwards
        self.pos = np.zeros((self.N, 2), float)
        # Velocity at time step i in m/s
        self.vel = np.zeros((self.N, 2), float)

        self.net_force = np.zeros((self.N, 2), float)
        self.F_g = np.zeros((self.N, 2), float)
        self.a = np.zeros((self.N, 2), float)

        # Planet mass in kg
        self.M = self.system.masses[0]*c.m_sun

        # Planet radius in m
        self.r = self.system.radii[0]*1000

        # The rotational period of our planet in seconds
        self.rp_sec = self.system.rotational_periods[0]*86400

        # Calculate planet rotational vector
        self.vpr = np.array([0, 2*np.pi*self.r/self.rp_sec])

        # Set initial rocket velocity
        self.vel[0] = self.vpr

        # Planet Orbital Velocity
        self.vpo = self.system.initial_velocities[:, 0]
        self.vpo_ms = (self.vpo*c.AU)/c.yr

        self.fuel_consumption = self.engine.fuel_consumption()

        # Engine thrust in N
        self.thrust = np.array([self.engine.thrust(), 0])

        # Escape velocity
        self.v_esc = np.sqrt((2*c.G*self.system.masses[0]*c.m_sun/(self.system.radii[0]*1000**3)))

        self.orbit_sim = orbits.orbit_sim()
        self.orbit_sim.analytical_solution()
        self.analytical = False  # Flag for using analytically calculated orbits
        self.load_planet_trajectories()

        # This objet is initialized after successfull launch
        self.interplanetary_travel = None

        self.astrogation_computer = AstrogationComputer()

        self.trilateration = Trilateration(self.mission)

    def __str__(self):
        """Return string representation."""
        str = "Launch parameters:\n"
        str = str + "                      Name | Value \n"
        str = str + "=========================================\n"
        str = str + f"               Planet mass | {self.M:g} kg\n"
        str = str + f"             Planet radius | {self.r:.2f} m\n"
        str = str + f"  Planet rotational period | {self.system.rotational_periods[0]:.2f} days\n"
        str = str + f"  Planet rotational period | {self.rp_sec:.2f} s\n"
        str = str + f"      Planet circumference | {2*np.pi*self.r:.2f} m\n"
        str = str + f"Planet rotational velocity | {np.linalg.norm(self.vpr):.2f} m/s\n"
        str = str + f"    Planet escape velocity | {self.v_esc:.2f} km/s\n"
        str = str + f"   Planet orbital velocity | {self.vpo_ms} m/s\n"
        orb_vel_vec = self.planet_orbital_velocity(0, self.time)*c.AU/c.yr
        str = str + f"   Planet orbital velocity | {orb_vel_vec} {np.linalg.norm(orb_vel_vec)} m/s\n"
        str = str + f"     Planet aphelion angle | {self.system.aphelion_angles[0]} r\n"
        str = str + f"      Planet orbital angle | {self.system._initial_orbital_angles[0]} r\n"
        str = str + f"          Fuel consumption | {self.fuel_consumption:.2f} kg/s\n"
        str = str + f"                    Thrust | {self.thrust[0]/1000:.2f} kN\n"
        str = str + f"                 Fuel mass | {self.fuel_mass} kg\n"
        str = str + f"   Cartesian launch coord. | {self.get_launch_pos(self.time)} AU\n"
        str = str + f"       Polar launch coord. | {tools.cartesian_polar(self.get_launch_pos(self.time))} AU, rad.\n"
        return str

    def launch_process(self, launch_angle, time):
        """Actually launches the rocket and calculates exit velocity and position of rocket.

        dt is delta t in sec. for the Euler Cromer calculation og the ascent

        Verifies the launch process and begins the interplanetary travel

        """
        self.launch_angle = launch_angle
        self.time = time
        print(f"               Launch time | {self.time:g} year")
        print(f"              Launch angle | {self.launch_angle:g} radians")
        print(f"   Escape velocity of Hoth | {self.v_esc} km/s")

        current_mass = self.initial_mass

        for i in range(self.N):
            self.i = i
            v = np.linalg.norm(self.vel[i])/1000
            if v >= self.v_esc:
                print(f"Achieved escape velocity {v:.2f} km/s")
                print(f"Our position is {self.pos[i]} relative to launch site and it took {i*self.dt:.2f} seconds")
                self.rocket_mass = current_mass
                print(f"It took {self.initial_mass-current_mass:.2f} kg of fuel")
                print(f'Spacecraft mass is now {current_mass} kg')
                self.__verify_launch_result()
                break
            else:
                current_mass = current_mass - self.fuel_consumption * self.dt
                if(current_mass <= self.mission.spacecraft_mass):
                    print(f"We ran out of fuel after {i*dt:.2f} seconds")
                    break

                self.F_g[i] = current_mass*self.g(self.M, self.r + np.linalg.norm(self.pos[i]))
                self.net_force[i] = self.thrust - self.F_g[i]
                # Accelleration m/s^2
                self.a[i] = self.net_force[i]/current_mass

                # Euler-cromer
                if(np.linalg.norm(self.a[i]) > 0):
                    self.vel[i+1] = self.vel[i] + self.a[i]*self.dt
                    self.pos[i+1] = self.pos[i] + self.vel[i+1]*self.dt

        print(f"Final velocity is {v:.6f} km/s")

    def __verify_launch_result(self):
        self.mission.set_launch_parameters(self.engine.thrust(),
                                           self.engine.fuel_consumption(),
                                           self.fuel_mass,
                                           1000,
                                           self.get_launch_pos(self.time),
                                           self.time)

        self.mission.launch_rocket()

        self.plot_orbit(0)

        self.mission.verify_launch_result(self.final_position())

        self.manual_orientation(self.time)

        self.interplanetary_travel = self.mission.begin_interplanetary_travel()

    def g(self, m, r):
        """Calculate g at distance r with mass m."""
        return np.array([(c.G*m)/r**2, 0])

    def final_position(self):
        """Return the rocket position relative to the star after launch."""
        self.escape_time = self.time + (self.i*self.dt/c.yr)

        planet_position = self.planet_position(0, self.escape_time)

        print(f'orbital_movement {np.sqrt(self.orbital_movement[0]**2+self.orbital_movement[1]**2)*c.AU} m')

        self.launch_site_position = planet_position + self.get_planet_launch_pos()

        if(self.verbose):
            print(f'Launch site position {self.launch_site_position} AU')

        escape_pos = self.get_planet_escape_pos()  # +self.orbital_movement

        # Add escape position vector to planet position vector
        rotated_escape_pos_x, rotated_escape_pos_y = tools.polar_cartesian(escape_pos[0], escape_pos[1])
        if(self.verbose):
            print(f'Rotated escape pos [{rotated_escape_pos_x} {rotated_escape_pos_y}] m')
        final_pos = self.launch_site_position + np.array([rotated_escape_pos_x, rotated_escape_pos_y])/c.AU

        print(f'Final escape pos [{final_pos[0]}, {final_pos[1]}] AU')

        return final_pos

    def final_velocity(self):
        """Get final velocity after launch."""
        return (self.vel[self.i]/c.AU*c.yr) + self.planet_orbital_velocity(0, self.escape_time)

    def get_planet_orbital_angel(self, t):
        """Get planet rotational angel at time t."""
        pos = np.array(tools.cartesian_polar(self.planet_position(0, t)))
        if(self.verbose):
            print(f'Planet orbital angle at {t} is {pos[1]}')
        return pos[1]

    def get_planet_rotational_angel(self, t):
        """Gets the angle that the planet has rotated to at time t."""
        rot_per_yr = c.yr/self.rp_sec  # Rotational periods pr year
        rotations_at_t = rot_per_yr*t
        last_rotation_fraction = rotations_at_t % 1
        angle = last_rotation_fraction * 2*np.pi
        if(self.verbose):
            print(f'Rotations at time {t} is {rotations_at_t}')
            print(f'Rotation fraction is {last_rotation_fraction}')
            print(f'Planet rotations pr year {rot_per_yr}')
            print(f'Rotational angle is {angle} at time {t}')
        return angle

    def get_planet_escape_pos(self):
        """Rotate escape position."""
        escape_pos = np.array(tools.cartesian_polar(self.pos[self.i]))
        escape_pos[1] += self.launch_angle
        return escape_pos

    def get_planet_launch_pos(self):
        """Get launch pos relative to planet center."""
        return tools.polar_cartesian(self.r/c.AU, self.launch_angle)

    def get_launch_pos(self, t):
        """Calculate initial launch position relative to star."""
        pos = self.planet_position(0, t)  # Initial planet center position

        ramp_pos = self.get_planet_launch_pos()  # Initial rocket position relative to planet center
        return pos + ramp_pos

    def load_planet_trajectories(self):
        """Load the planet trajectories from file."""
        trajectories_filename = "planet_trajectories.npz"

        url = 'https://www.uio.no/studier/emner/matnat/astro/AST2000/h20/blogger/Flukten%20fra%20Hoth/\
                data/planet_trajectories.npz'
        if (os.path.exists(trajectories_filename) is False):
            try:
                import requests
                r = requests.get(url, allow_redirects=True)

                open(trajectories_filename, 'wb').write(r.content)
            except Exception:
                print('You need to install requests to download file: pip install requests')

        # Load the saved planet trajectories
        trajectories = np.load(trajectories_filename, allow_pickle=True)
        self.times = trajectories[0]
        self.planet_positions = trajectories[1]

    def planet_position_analytical(self, n, t):
        """Calculate the analytical position og planet n at time t."""
        if(self.verbose):
            print(f'Orbital periode {self.orbit_sim.hoth_period()}')
        times = np.arange(0, self.orbit_sim.hoth_period(), self.orbit_sim.hoth_period()/1000)
        if(self.verbose):
            print(f'times {times}')
        idx = (np.abs(times - t)).argmin()
        if(self.verbose):
            print(f'idx {idx}')
        return np.array(self.orbit_sim.r_analytical)[n, :, idx]

    def planet_position(self, n, t):
        """Find position of planet n at time t."""
        if(self.analytical is True):
            return self.planet_position_analytical(n, t)

        # Find index of time closest to our t
        idx = (np.abs(self.times - t)).argmin()
        print(f'Time is {t}')
        print(f'Planet position is {self.planet_positions[:,n,idx]} at index {idx} at time {self.times[idx]}')
        print(f'Planet position is {self.planet_positions[:,n,idx+1]} at index {idx+1} at time {self.times[idx+1]}')

        # Get planet position at time t
        planet_position = self.planet_positions[:, n, idx]

        # Interpolate position to exact time
        dt = t - self.times[idx]
        # Calculate orbital movement
        orbital_velocity = self.planet_orbital_velocity(n, self.time)

        print(f'Difference in time between position got and actual time {dt}')
        print(f'Orbital velocity is  {orbital_velocity}')
        self.orbital_movement = orbital_velocity*dt

        return planet_position + self.orbital_movement

    def planet_orbital_velocity(self, n, t):
        """Calculate orbital velocity of the planet n at time t in years."""
        # t = 0 is special case
        if t == 0:
            return self.system.initial_velocities[:, n]

        # Find index of time closest to our t
        idx = (np.abs(self.times - self.time)).argmin()

        # Calculate delta t in years
        deltaT = (self.times[idx]-self.times[idx-1])
        if(self.verbose):
            print(f'deltaT {t} {idx} {deltaT} = {self.times[idx]} - {self.times[idx-1]}')
        # Calculate delta s in AU
        deltaS = (self.planet_positions[:, n, idx] - self.planet_positions[:, n, idx-1])
        # Planet orbital velocity vector in AU/year
        planet_velocity = deltaS/deltaT
        if(self.verbose):
            print(f'Calculated orbital velocity {planet_velocity} = {deltaS}/{deltaT} :\
             {np.linalg.norm(planet_velocity*c.AU/c.yr)} m/s')

        return planet_velocity

    def plot_launch_position(self):
        """Plot launch position."""
        fig, ax = plt.subplots()
        r = self.r/c.AU
        planet = plt.Circle((0, 0), r, color='b', fill=True)
        ax.set(xlim=(-r*1.1, r*1.1), ylim=(-r*1.1, r*1.1))
        ax.add_artist(planet)
        rocket_pos = tools.polar_cartesian(r, self.launch_angle)
        rocket = plt.Circle((rocket_pos[0], rocket_pos[1]), r/40, color='r', fill=True)
        ax.add_artist(rocket)
        plt.show()

    def get_trajectories(self, n):
        """Return the trajectories."""
        if(self.analytical is True):
            return np.array(self.orbit_sim.r_analytical)[n, :, :]

        return self.planet_positions[:, n, :100000:10]

    def plot_orbit(self, n, zoom_to_planet=True):
        """Plot orbit of planet n and mark position at time self.time."""
        rocket_position = self.final_position()
        planet_position = self.planet_position(0, self.escape_time)
        trajectories = self.get_trajectories(n)
        # Plot the planet orbit around the star
        fig, ax = plt.subplots()
        ax.set_xlabel('AU')
        ax.set_ylabel('AU')
        if(zoom_to_planet is True):
            ax.set_xlim([planet_position[0] - self.mission.system.radii[n]*3*1000/c.AU,
                        planet_position[0] + self.mission.system.radii[n]*3*1000/c.AU])
            ax.set_ylim([planet_position[1]-self.mission.system.radii[n]*3*1000/c.AU,
                        planet_position[1]+self.mission.system.radii[n]*3*1000/c.AU])
        else:
            ax.set_xlim([trajectories[0, :].min()*1.05, trajectories[0, :].max()*1.05])
            ax.set_ylim([trajectories[1, :].min()*1.05, trajectories[0, :].max()*1.05])

        # Plot the star maginified 10 times
        star = plt.Circle((0, 0), 10*self.mission.system.star_radius*1000/c.AU, color='y', fill=True)

        # Plot the planet
        planet = plt.Circle((planet_position), self.mission.system.radii[n]*1000/c.AU, color='b', fill=True)

        ax.plot(self.launch_site_position[0], self.launch_site_position[1], 'r.')
        ax.plot(rocket_position[0], rocket_position[1], 'g.')
        ax.add_artist(star)
        ax.add_artist(planet)
        # Plot the center of the planet
        ax.plot(planet_position[0], planet_position[1], 'k.')

        ax.plot(trajectories[0], trajectories[1], lw=0.1, label='Planet trajectory')

        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.show()

    def boost(self, dv):
        """Boost the engine with vector dv."""
        self.interplanetary_travel.boost(dv)

    def coast(self, duration):
        """Coast through space for duration years."""
        self.interplanetary_travel.coast(duration)

    def manual_orientation(self, time):
        """Does a manual orientation of the spacecraft."""
        self.mission.take_picture()

        dopler_shifts = self.mission.measure_star_doppler_shifts()

        print(f'Doppler shifts {dopler_shifts}')

        distances = np.array(self.mission.measure_distances())

        print(f'Planet and star distances\n{distances}')

        angle = self.astrogation_computer.find_orientation_angle('sky_picture.png')
        print(f'Spacecraft is pointing in direction {angle} degrees')

        velocity = self.trilateration.radial_velocity()

        print(f'The spacecraft velocity is [{velocity[0]},{velocity[1]}]')

        trajection_computer = TrajectionComputer(self.rocket_mass)

        position = self.trilateration.tri_test(distances, trajection_computer.planet_positions_at_t(time))

        print(f'Spacecraft position is {position}')

        self.mission.verify_manual_orientation(position, velocity, angle)

        SpaceMission.save('after_launch.bin', self.mission)


def rocket_engine_factory(filename, number_of_motors, temperature, steps, particles, dt, args):
    """Cosntruct, or read from file, a rocket enigine with given number of motors."""
    if(os.path.exists(filename) is False or args.run_steps == 'true'):
        motor = nano_motor(10**-6, particles, temperature, dt)
        print(f'Running engine particle simulation for {steps} steps')
        for i in range(steps):
            print(f"{i:4d}\b\b\b\b", end="", flush=True)
            motor.step()

        with open(filename, "wb") as output:
            pickle.dump(motor, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, "rb") as input:
            motor = pickle.load(input)

    # Construct the engine
    engine = re.RocketEngine(number_of_motors, motor)

    return engine


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=float, default=0.0, help='Time in years')
    parser.add_argument('-la', '--launch-angle', type=float, default=0.0, help='Launch angle along equator')
    parser.add_argument('--run-steps', help='Run engine simulation')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print debug statements')
    args = parser.parse_args()

    dt = 10**-12
    fuel_mass = 43000
    dv = 1000
    steps = 1000
    motors = 11.0699895**15
    max_launch_time = 1000
    rocket_filename = "nano_motor.pkl"
    temperature = 3500
    particles = 9.8**5
    t = args.t       # Time in years after t=0 to launch at
    angle = args.launch_angle  # Angle along equator to launch from in radians

    # Construct the engine
    engine = rocket_engine_factory(rocket_filename, motors, temperature, steps, particles, dt, args)

    # Construct the launch object
    dt = 0.01
    space_craft = Spacecraft(fuel_mass, engine, dt, 33382, args.verbose)

    space_craft.time = t
    space_craft.launch_angle = angle
    # Prints details of the rocket
    print(space_craft)

    # Launch the rocket
    space_craft.launch_process(angle, t,)

    space_craft.plot_orbit(0)

    SpaceMission.save('part5.bin', space_craft.mission)


"""
Example running the code:

janmagneandersen$  python spacecraft.py -t 0.24914953763595418 -la 0.8                (ast2000) 52.5s
Time is 0.24914953763595418
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time -3.899717475663866e-06
Orbital velocity is  [-0.36480495  8.81593169]
Time is 0.24914953763595418
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time -3.899717475663866e-06
Orbital velocity is  [-0.36480495  8.81593169]
Launch parameters:
                      Name | Value
=========================================
               Planet mass | 7.3112e+24 kg
             Planet radius | 6622416.98 m
  Planet rotational period | 1.03 days
  Planet rotational period | 89145.88 s
      Planet circumference | 41609873.06 m
Planet rotational velocity | 466.76 m/s
    Planet escape velocity | 12.14 km/s
   Planet orbital velocity | [    0.         41827.38053534] m/s
   Planet orbital velocity | [-1729.3826138 41792.5219251] 41828.28771399426 m/s
     Planet aphelion angle | 0.0 r
      Planet orbital angle | 0.0 r
          Fuel consumption | 123.03 kg/s
                    Thrust | 582.93 kN
                 Fuel mass | 43000 kg
   Cartesian launch coord. | [0.36553528 0.01484045] AU
       Polar launch coord. | (0.3658364102366809, 0.040576939163003735) AU, rad.

               Launch time | 0.24915 year
              Launch angle | 0.8 radians
   Escape velocity of Hoth | 12.139587606214251 km/s
Achieved escape velocity 12.14 km/s
Our position is [787378.05040441 161312.81260592] relative to launch site and it took 345.60 seconds
It took 42519.25 kg of fuel
Spacecraft mass is now 1580.7467320600986 kg
Time is 0.24914953763595418
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time -3.899717475663866e-06
Orbital velocity is  [-0.36480495  8.81593169]
Rocket was moved up by 963.224 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 343.79 s.
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
orbital_movement 9308333.4147471 m
Final escape pos [0.3655341771636786, 0.014941524270976511] AU
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
orbital_movement 9308333.4147471 m
Final escape pos [0.3655341771636786, 0.014941524270976511] AU
Your spacecraft position was satisfyingly calculated. Well done!
*** Achievement unlocked: No free launch! ***
Picture written to sky_picture.png.
Doppler shifts (-0.10960464349253882, 0.06693221889417213)
Planet and star distances
[4.93755691e-05 9.50555542e-01 2.89686137e+00 1.89849876e+00
 1.06073918e+00 4.72010385e+00 7.28395216e-01 4.32597359e-01
 3.65839283e-01]
Comparing to angle  359
Spacecraft is pointing in direction 37 degrees
The spacecraft velocity is [1.251074254422691,10.62333945621307]
Planet trajectories 8 planets loaded
Shortest distance from us [8 7 6]
planet_pos.shape (8, 2)
planet_pos [[ 0.36550444  0.01480869]
 [-0.58050164 -0.07766207]
 [ 2.7420217   1.67147304]
 [ 1.95497118  1.05319478]
 [ 1.00992576 -0.82764346]
 [ 0.07087376 -4.69595837]
 [ 0.21076575  0.7266676 ]
 [-0.0418933   0.15998709]]
(2, 9)
Planet positions: [[ 0.36550444 -0.58050164  2.7420217   1.95497118  1.00992576  0.07087376
   0.21076575 -0.0418933   0.        ]
 [ 0.01480869 -0.07766207  1.67147304  1.05319478 -0.82764346 -4.69595837
   0.7266676   0.15998709  0.        ]]
Positions [[ 0.         -0.0418933   0.21076575]
 [ 0.          0.15998709  0.7266676 ]]
coordinate candidates (6, 2)
Length between 1 and 2 0.004583222114380304
Length between 1 and 4 0.012812755010802033
Length between 1 and 5 0.0170099430032587
Length between 2 and 1 0.004583222114380304
Length between 2 and 4 0.012012368114612678
Length between 2 and 5 0.015554114949931563
Length between 4 and 1 0.012812755010802033
Length between 4 and 2 0.012012368114612678
Length between 4 and 5 0.027461658924070962
Length between 5 and 1 0.0170099430032587
Length between 5 and 2 0.015554114949931563
Length between 5 and 4 0.027461658924070962
{0: [0], 1: [1, 2, 4, 5], 2: [2, 1, 4, 5], 3: [3], 4: [4, 1, 2, 5], 5: [5, 1, 2, 4]}
x: [0.36572371854579994 0.36523789913991483 0.3533069816178027
 0.38013158952085263]
y: [0.009194725207393683 0.0137521262710055 0.012355632518237014
 0.018236374113447473]
x: [0.36523789913991483 0.36572371854579994 0.3533069816178027
 0.38013158952085263]
y: [0.0137521262710055 0.009194725207393683 0.012355632518237014
 0.018236374113447473]
x: [0.3533069816178027 0.36572371854579994 0.36523789913991483
 0.38013158952085263]
y: [0.012355632518237014 0.009194725207393683 0.0137521262710055
 0.018236374113447473]
x: [0.38013158952085263 0.36572371854579994 0.36523789913991483
 0.3533069816178027]
y: [0.018236374113447473 0.009194725207393683 0.0137521262710055
 0.012355632518237014]
Spacecraft position is [0.3661000472060925, 0.013384714527520918]
Pointing angle after launch correctly calculated. Well done!
Velocity after launch correctly calculated. Well done!
Position after launch correctly calculated. Well done!
Your manually inferred orientation was satisfyingly calculated. Well done!
*** Achievement unlocked: Well-oriented! ***
SpaceMission instance saved as after_launch.bin.
Final velocity is 12.140128 km/s
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
orbital_movement 9308333.4147471 m
Final escape pos [0.3655341771636786, 0.014941524270976511] AU
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
SpaceMission instance saved as part5.bin.
"""
