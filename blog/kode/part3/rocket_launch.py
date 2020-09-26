# Egen kode
import numpy as np
import nano_motor as nm
import rocket_engine as re
import ast2000tools.utils as utils
import ast2000tools.constants as c
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
import matplotlib.pyplot as plt
import sys
import os
import tools
try:
    import cPickle as pickle
except:
    import pickle
import argparse
from ast2000tools.shortcuts import SpaceMissionShortcuts

import orbit_sim as orbits
import load_orbit_sim as los

parser = argparse.ArgumentParser()
parser.add_argument('--run-steps', help='foo help')
args = parser.parse_args()

'''
Class for launching the rocket_thrust

Takes mission, initial_mass of rocket (including fuel) and one nano_enigne
as arguments to the constructor

Call launch_process to actually launch the rocket.

- launch_pos is angle on planet surface from 0 to 2 pi
'''
class rocket_launch:
    # A class for calculating the launch of our rocket
    def __init__(self, mission, initial_mass, engine, launch_pos):
        self.mission = mission
        self.initial_mass = initial_mass
        self.system = SolarSystem(33382)
        self.engine = engine
        self.launch_pos = launch_pos
        # Number of time steps to try
        self.N = 50000

        # Positions at time step i in m. Not necessary to store all unless we want to plot afterwards
        self.pos = np.zeros((self.N, 2), float)
        # Velocity at time step i in m/s
        self.vel = np.zeros((self.N, 2), float)

        self.net_force = np.zeros((self.N,2),float)
        self.F_g = np.zeros((self.N,2),float)
        self.a = np.zeros((self.N,2),float)

        # Planet mass in kg
        self.M = self.system.masses[0]*c.m_sun

        # Planet radius in m
        self.r = self.system.radii[0]*1000

        # The rotational period of our planet in seconds
        self.rp_sec = self.system.rotational_periods[0]*86400

        # Calculate planet rotational vector at launch position
        self.vpr = np.array([0,2*np.pi*self.r/self.rp_sec]) # Planet rotational velocity

        # Set initial rocket velocity
        self.vel[0] = self.vpr

        # Planet Orbital Velocity
        self.vpo = self.system.initial_velocities[:,0]
        self.vpo_ms = (self.vpo*c.AU)/c.yr

        self.fuel_consumption = self.engine.fuel_consumption()

        # Engine trhust i N
        self.thrust = np.array([self.engine.thrust(),0])

        # Escape velocity
        self.v_esc = np.sqrt((2*c.G*self.system.masses[0]*c.m_sun)/(self.system.radii[0]*1000**3))
    def __str__(self):
        str = f"Launch parameters:\n"
        str = str + f"                      Name | Value \n"
        str = str + f"=========================================\n"
        str = str + f"               Planet mass | {self.M:g} kg\n"
        str = str + f"             Planet radius | {self.r:.2f} m\n"
        str = str + f"  Planet rotational period | {self.system.rotational_periods[0]:.2f} days\n"
        str = str + f"  Planet rotational period | {self.rp_sec:.2f} s\n"
        str = str + f"      Planet circumference | {2*np.pi*self.r:.2f} m\n"
        str = str + f"Planet rotational velocity | {np.linalg.norm(self.vpr):.2f} m/s\n"
        str = str + f"    Planet escape velocity | {self.v_esc:.2f} km/s\n"
        str = str + f"   Planet orbital velocity | {self.vpo_ms} m/s\n"
        str = str + f"     Planet aphelion angle | {self.system.aphelion_angles[0]} r\n"
        str = str + f"      Planet orbital angle | {self.system._initial_orbital_angles[0]} r\n"
        str = str + f"          Fuel consumption | {self.fuel_consumption:.2f} kg/s\n"
        str = str + f"                    Thrust | {self.thrust[1]/1000:.2f} kN\n"
        str = str + f"   Cartesian launch coord. | {self.get_launch_pos()} AU\n"
        str = str + f"       Polar launch coord. | {tools.cartesian_polar(self.get_launch_pos())} AU, radians\n"
        return str

    '''
    Actually launches the rocket and calculates exit velocity and position of rocket.

    dt is delta t in sec. for the Euler Cromer calculation og the ascent

    '''
    def launch_process(self,dt):

        current_mass = self.initial_mass
        self.dt = dt

        print(f"Escape velocity of Hoth is {self.v_esc} km/s")

        for i in range(self.N):
            self.i = i
            v = np.linalg.norm(self.vel[i])/1000
            if v >= self.v_esc:
                print(f"Achieved escape velocity {v:.2f} km/s")
                print(f"Our position is {self.pos[i]} relative to launch site and it took {i*dt:.2f} seconds")
                print(f"It took {self.initial_mass-current_mass:.2f} kg of fuel")
                self.time = i*dt
                break
            else:
                current_mass = current_mass - self.fuel_consumption * dt
                if(current_mass <= self.mission.spacecraft_mass):
                    print(f"We ran out of fuel after {i*dt:.2f} seconds")
                    break

                self.F_g[i] = current_mass*self.g(self.M,self.r + np.linalg.norm(self.pos[i]))
                self.net_force[i] = self.thrust - self.F_g[i]
                # Accelleration m/s^2
                self.a[i] = self.net_force[i]/current_mass
                #print(f"Pos {self.pos[i]} net force {net_force} F_G {F_G} kg/ms^2 a {a} m/s^2")
                # Euler-cromer
                if(np.linalg.norm(self.a[i]) > 0):
                    self.vel[i+1] = self.vel[i] + self.a[i]*dt
                    self.pos[i+1] = self.pos[i] + self.vel[i+1]*dt


        print(f"Final velocity is {v:.6f} km/s")

    def g(self,m,r):
        '''
        Calculate g at distance r with mass m
        '''
        return np.array([(c.G*m)/r**2,0])

    def final_position(self,launch_pos):
        '''
        Returns the rocket position relative to the star after launch
        '''
        # Add planet orbital movement to launch position
        planet_pos = launch_pos + (self.vpo_ms/c.AU)*self.i*self.dt

        # Rotate escape position
        escape_pos = np.array(tools.cartesian_polar(self.pos[self.i]))
        escape_pos[1] += self.launch_pos

        # Add escape position vector to planet position vector
        rotated_escape_pos_x, rotated_escape_pos_y  = tools.polar_cartesian(escape_pos[0],escape_pos[1])
        final_pos = planet_pos + np.array([rotated_escape_pos_x,rotated_escape_pos_y])/c.AU

        return final_pos

    def get_launch_pos(self):
        '''
        Calculate initial launch position relative to star
        '''
        pos = np.array(self.mission.system.initial_positions[:,0])                  # Initial planet center position
        rocket_pos = tools.polar_cartesian(self.r/c.AU,self.launch_pos)   # Initial rocket position relative to planet center
        return pos + rocket_pos

def rocket_engine_factory(filename, number_of_motors, args):
    '''
        Cosntruct, or read from file, a rocket enigine with given number of motors
    '''
    if(os.path.exists(filename) == False or args.run_steps == 'true'):
        motor = nm.nano_motor(10**-6, 10**5, temperature, dt)

        for i in range(steps):
            print(f"{i:4d}\b\b\b\b", end = "", flush = True)
            motor.step()

        with open(filename, "wb") as output:
            pickle.dump(motor, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, "rb") as input:
            motor = pickle.load(input)

    # Construct the engine
    engine = re.rocket_engine(number_of_motors, motor)

    return engine


if __name__ == '__main__':
    mission = SpaceMission(33382)
    dt = 10**-12
    fuel_mass = 41000
    dv = 1000
    steps = 1000
    motors = 11.0699895**15
    max_launch_time = 1000
    filename = "test_launch.pkl"
    temperature = 3000

    # Construct the engine
    engine = rocket_engine_factory(filename,motors,args)

    # Calculate fuel consumed based on boost
    fuel_consumed = engine.boost(mission.spacecraft_mass+fuel_mass, dv)

    print(f"Fuel mass {fuel_mass} kg")

    # Construct the launch object
    launch = rocket_launch(mission,mission.spacecraft_mass+fuel_mass,engine,2)

    # Prints details of the rocket
    print(launch)

    # Launch the rocket
    launch.launch_process(0.01)

    launch_pos = launch.get_launch_pos()
    escape_pos = launch.final_position(launch_pos)
    print(f"Start coordinates {launch_pos} AU")
    print(f"  End coordinates {escape_pos} AU")
    '''
    fig, ax1 = plt.subplots()

    # Plot Trajectory
    plt.ylabel('km')
    plt.xlabel('km')
    plt.plot(launch.pos[:launch.i,0]/1000,launch.pos[:launch.i,1]/1000,label='Trajectory')
    plt.legend()
    plt.show()

    # Plot velocity as a function of time
    plt.ylabel('km/s')
    plt.xlabel('s')
    plt.plot(np.linspace(0,launch.time,launch.i),np.linalg.norm(launch.vel[:launch.i],axis=1)/100,'r',label='Velocity')
    plt.legend()
    #plt.show()
    '''
    # Launch the rocket using AST2000tools
    mission.set_launch_parameters(engine.thrust(),engine.fuel_consumption(),fuel_mass,max_launch_time,launch_pos,0)

    mission.launch_rocket()

    mission.verify_launch_result(escape_pos)

    SpaceMission.save('part1.bin',mission)

'''
Example running the code:

python rocket_launch.py --run-steps true
true
Fuel mass 41000 kg
Launch parameters:
                      Name | Value
=========================================
               Planet mass | 7311200101015753339174912.00 kg
             Planet radius | 6622.42 km
  Planet rotational period | 89145.88 s
    Planet escape velocity | 12.14 km/s
Planet rotational velocity | 0.47 km/s
   Planet orbital velocity | [ 0.         41.82738054] km/s
          Fuel consumption | 126.49 kg/s
                    Thrust | 542.27 kN

Escape velocity of Hoth is 12.139587606214251 km/s
Achieved escape velocity 12.14 km/s
Our position is 677361.3055238673 meters relative to launch site and it took 324.08 seconds
It took 40993.03 kg of fuel
Final velocity is 12.143855 km/s
Start coordinates [0.36585653049782685, 0.0] AU
  End coordinates [0.36586105837848687, 2.8271887793580364e-16] AU
Rocket was moved down by 3.29129e-06 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 322.83 s.
Your spacecraft position deviates too much from the correct position.
The deviation is approximately 9.12702e-05 AU.
Make sure you have included the rotation and orbital velocity of your home planet.
Note that units are AU and relative the the reference system of the star.
Traceback (most recent call last):
  File "rocket_launch.py", line 190, in <module>
    mission.verify_launch_result(escape_pos)
  File "build/bdist.linux-x86_64/egg/ast2000tools/space_mission.py", line 348, in verify_launch_result
RuntimeError: Incorrect spacecraft position after launch.
'''
