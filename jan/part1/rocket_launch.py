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
try:
    import cPickle as pickle
except:
    import pickle
import argparse
from ast2000tools.shortcuts import SpaceMissionShortcuts
parser = argparse.ArgumentParser()
parser.add_argument('--run-steps', help='foo help')
args = parser.parse_args()

'''
Class for launching the rocket_thrust

Takes mission, initial_mass of rocket (including fuel) and one nano_enigne
as arguments to the constructor

Call launch_process to actually launch the rocket.

'''
class rocket_launch:
    # A class for calculating the launch of our rocket
    def __init__(self, mission, initial_mass, engine):
        self.mission = mission
        self.initial_mass = initial_mass
        self.system = SolarSystem(33382)
        self.engine = engine
        # Number of time steps to try
        self.N = 50000

        # Positions at time step i in m. Not necessary to store all unless we want to plot afterwards
        self.pos = np.zeros(self.N, float)
        # Velocity at time step i in m/s
        self.vel = np.zeros(self.N, float)

        # Planet mass in kg
        self.M = self.system.masses[0]*c.m_sun

        # Planet radius in m
        self.r = self.system.radii[0]*1000

        # The rotational period of our planet in seconds
        self.rp_sec = self.system.rotational_periods[0]*86400

        # Rotational velocity in m/s
        self.rv = 2*np.pi*self.r/self.rp_sec

        # Orbital Velocity
        self.ov = self.system.initial_velocities[:,0]
        self.ov_ms = (self.ov*c.AU)/c.yr

        self.fuel_consumption = self.engine.fuel_consumption()

        # Engine trhust i N
        self.thrust = self.engine.thrust()

        # Escape velocity
        self.v_esc = np.sqrt((2*c.G*self.system.masses[0]*c.m_sun)/(self.system.radii[0]*1000**3))
    def __str__(self):
        str = f"Launch parameters:\n"
        str = str + f"                      Name | Value \n"
        str = str + f"=========================================\n"
        str = str + f"               Planet mass | {self.M:.2f} kg\n"
        str = str + f"             Planet radius | {self.r/1000:.2f} km\n"
        str = str + f"  Planet rotational period | {self.rp_sec:.2f} s\n"
        str = str + f"    Planet escape velocity | {self.v_esc:.2f} km/s\n"
        str = str + f"Planet rotational velocity | {self.rv/1000:.2f} km/s\n"
        str = str + f"   Planet orbital velocity | {self.ov_ms/1000} km/s\n"
        str = str + f"          Fuel consumption | {self.fuel_consumption:.2f} kg/s\n"
        str = str + f"                    Thrust | {self.thrust/1000:.2f} kN\n"
        return str

    '''
    Actually launches the rocket and calculates exit velocity and position of rocket.

    dt is delta t in sec. for the Euler Cromer calculation og the ascent

    '''
    def launch_process(self,dt):

        current_mass = self.initial_mass

        print(f"Escape velocity of Hoth is {self.v_esc} km/s")

        for i in range(self.N):
            self.i = i
            v = np.sqrt(self.vel[i]**2+self.rv**2+np.sqrt(self.ov_ms[0]**2+self.ov_ms[1]**2))/1000
            if v >= self.v_esc:
                print(f"Achieved escape velocity {v:.2f} km/s")
                print(f"Our position is {self.pos[i]} meters relative to launch site and it took {i*dt:.2f} seconds")
                print(f"It took {self.initial_mass-current_mass:.2f} kg of fuel")
                self.time = i*dt
                break
            else:
                current_mass = current_mass - self.fuel_consumption * dt
                if(current_mass <= self.mission.spacecraft_mass):
                    print(f"We ran out of fuel after {i*dt:.2f} seconds")
                    break
                # The Force of gravity GM/r^2 in kg/ms^2
                F_G = (c.G*(current_mass*self.M))/((self.r+self.pos[i])**2)
                net_force = self.thrust - F_G
                # Accelleration m/s^2
                a = net_force/current_mass
                #print(f"Pos {self.pos[i]} net force {net_force} F_G {F_G} kg/ms^2 a {a} m/s^2")
                # Euler-cromer
                if(a > 0):
                    self.vel[i+1] = self.vel[i] + a*dt
                    self.pos[i+1] = self.pos[i] + self.vel[i+1]*dt


        print(f"Final velocity is {v:.6f} km/s")
mission = SpaceMission(33382)
dt = 10**-12
fuel_mass = 41000
dv = 1000
steps = 1000
motors = 11.0699895**15
max_launch_time = 1000
filename = "test_launch.pkl"
temperature = 3000

'''
    Check if we should run simulation from scracth or read from file
'''
if(os.path.exists(filename) == False or args.run_steps == 'true'):
    motor = nm.nano_motor(10**-6, 10**5, temperature, dt)

    for i in range(steps):
        print(f"{i:4d}\b\b\b\b", end = "", flush = True)
        motor.step()

    with open("test_launch.pkl", "wb") as output:
        pickle.dump(motor, output, pickle.HIGHEST_PROTOCOL)
else:
    with open(filename, "rb") as input:
        motor = pickle.load(input)

# Construct the engine
engine = re.rocket_engine(motors, motor)

# Calculate fuel consumed based on boost
fuel_consumed = engine.boost(mission.spacecraft_mass+fuel_mass, dv)

print(f"Fuel mass {fuel_mass} kg")

# Launch the rocket
launch = rocket_launch(mission,mission.spacecraft_mass+fuel_mass,engine)

# Prints details of the rocket
print(launch)

# Launch the rocket
launch.launch_process(0.01)

# Calculate coordinates in AU after launch
xp = mission.system.initial_positions[0,0]
yp = mission.system.initial_positions[1,0]
pr = (mission.system.radii[0]*1000)/c.AU
dradiell=(launch.rv*(dt*steps))/c.AU
dorbitell=(np.sqrt(launch.ov_ms[0]**2+launch.ov_ms[1]**2)*(dt*steps))/c.AU

launch_pos = [xp+pr,yp]
escape_pos = [xp+pr+(launch.pos[launch.i]/c.AU),yp+dradiell+dorbitell]
print(f"Start coordinates {launch_pos} AU")
print(f"  End coordinates {escape_pos} AU")

# Plot heigth as a function of time
plt.plot(np.linspace(0,launch.time,launch.i),launch.pos[:launch.i])
plt.xlabel('Time in sec')
plt.ylabel('Height in meters')
plt.show()
# Plot velocity as a function of time
plt.plot(np.linspace(0,launch.time,launch.i),launch.vel[:launch.i])
plt.xlabel('Time in sec')
plt.ylabel('Velocity m/s')
plt.show()

# Launch the rocket using AST2000tools
mission.set_launch_parameters(engine.thrust(),engine.fuel_consumption(),fuel_mass,max_launch_time,launch_pos,0)

mission.launch_rocket()

# Appluy shortcut
shortcuts = SpaceMissionShortcuts(mission, [38433])
consumed_fuel_mass, final_time, final_position, final_velocity = shortcuts.get_launch_results()
print(f'Consumed fuel {consumed_fuel_mass}')
print(f'Final time {final_time}')
print(f'Final position {final_position}')
print(f'Final velocity {final_velocity}')
mission.verify_launch_result(final_position)

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
