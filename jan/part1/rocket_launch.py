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
parser = argparse.ArgumentParser()
parser.add_argument('--run-steps', help='foo help')
args = parser.parse_args()

class rocket_launch:
    # A class for calculating the launch of our rocket
    def __init__(self, mission, initial_mass, engine):
        self.mission = mission
        self.initial_mass = initial_mass
        self.system = SolarSystem(33382)
        self.engine = engine

    def launch_process(self,dt):
        # Escape velocity
        v_esc = np.sqrt((2*c.G*self.system.masses[0]*c.m_sun)/(self.system.radii[0]*1000**3))

        print(f"Escape velocity of Hoth is {v_esc} km/s")
        # Number of time steps to try
        N = 50000
        # Positions at time step i in m. Not necessary to store all unless we want to plot afterwards
        self.pos = np.zeros(N, float)
        # Velocity at time step i in m/s
        self.vel = np.zeros(N, float)

        # Planet mass in kg
        M = self.system.masses[0]*c.m_sun
        # Planet radius in m
        r = self.system.radii[0]*1000

        fuel_consumption = self.engine.fuel_consumption()
        current_mass = self.initial_mass
        # Engine trhust i N
        thrust = self.engine.thrust()

        for i in range(N):
            self.i = i
            if self.vel[i]/1000 >= v_esc:
                print(f"Achieved escape velocity {self.vel[i]:.2f} km/s")
                print(f"Our height is {self.pos[i]*1000:.2f} meters and it took {i*dt:.2f} seconds")
                print(f"It took {self.initial_mass-current_mass:.2f} kg of fuel")
                break
            else:
                current_mass = current_mass - fuel_consumption * dt
                if(current_mass <= self.mission.spacecraft_mass):
                    print(f"We ran out of fuel after {i*dt:.2f} seconds")
                    break
                # The Force of gravity GM/r^2 in kg/ms^2
                F_G = (c.G*(current_mass+M))/((r+self.pos[i])**2)
                net_force = thrust - F_G * current_mass
                # Accelleration m/s^2
                a = net_force/current_mass
                #print(f"Pos {self.pos[i]} net force {net_force} F_G {F_G} kg/ms^2 a {a} m/s^2")
                # Euler-cromer
                if(a > 0):
                    self.vel[i+1] = self.vel[i] + a*dt
                    self.pos[i+1] = self.pos[i] + self.vel[i+1]*dt

        print(f"Final velocity is {self.vel[i]/1000} km/s")
mission = SpaceMission(33382)
dt = 1e-12
fuel_mass = 35000
dv = 1000 #np.sqrt((2*c.G*mission.massses[0])/mission.radii[0]*1000**3)
steps = 1000
motors = 11**15
max_launch_time = 1000
filename = "test_launch.pkl"
temperature = 3000

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

engine = re.rocket_engine(motors, motor)
fuel_consumed = engine.boost(mission.spacecraft_mass+fuel_mass, dv)
print(fuel_consumed)

launch = rocket_launch(mission,mission.spacecraft_mass+fuel_mass,engine)

launch.launch_process(0.01)

start_pos = [mission.system.initial_positions[0,0]+((mission.system.radii[0]*1000)/c.AU),mission.system.initial_positions[1,0]]
print(f"Start coordinates {start_pos} AU")
print(f"Engine thrust {engine.thrust()} N")
print(f"Fuel consumption {engine.fuel_consumption()} kg/s")
print(f"Fuel mass {fuel_mass} kg")
print(f"Planet radius {mission.system.radii[0]} km")
print(f"Planet mass {mission.system.masses[0]*c.m_sun} kg")
mission.set_launch_parameters(engine.thrust(),engine.fuel_consumption(),fuel_mass,max_launch_time,start_pos,0)
mission.launch_rocket()

# Plot heigth as a function of time
plt.plot(np.linspace(0,launch.i,launch.i),launch.pos[:launch.i])

# Plot velocity as a function of time
plt.plot(np.linspace(0,launch.i,launch.i),launch.vel[:launch.i])
plt.show()
