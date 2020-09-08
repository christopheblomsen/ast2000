# Egen kode
import numpy as np
import nano_motor as nm
import rocket_engine as re
import ast2000tools.utils as utils
import ast2000tools.constants as c
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
import sys
import os
try:
    import cPickle as pickle
except:
    import pickle


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

        print(f"Escape velocity is {v_esc} m/s")
        N = 15000
        # Positions at time step i. Not necessary to store all unless we want to plot afterwards
        pos = np.zeros(N, float)
        # Velocity at time step i
        vel = np.zeros(N, float)
        # Momentum at time step i
        p = np.zeros(N, float)

        pos[0] = 1
        pos0 = np.array([0, 0, 0])
        v = np.array([0, 0, 0])

        h = 1

        time = np.linspace(0, N, N)

        M = self.system.masses[0]

        fuel_consumption = self.engine.fuel_consumption()
        current_mass = self.initial_mass
        thrust = self.engine.thrust()
        for i in range(len(time)-1):
            if vel[i] >= v_esc:
                print(f"Achieved escape velocity {vel[i]:.2f} km/s")
                print(f"Our height is {pos[i]:.2f} meters and it took {time[i]*dt:.2f} seconds")
                print(f"It took {self.initial_mass-current_mass:.2f} kg of fuel")
                break
            else:
                current_mass = current_mass - fuel_consumption * dt
                if(current_mass <= self.mission.spacecraft_mass):
                    print(f"We ran out of fuel after {time[i]*dt:.2f} seconds")
                    sys.exit()
                # The Force of gravity
                F_G = (c.G*current_mass*M)/(self.system.radii[0]+pos[i])**2
                # Accelleration
                a = (thrust-F_G)/current_mass

                # Euler-cromer
                vel[i+1] = vel[i] + a*dt
                pos[i+1] = pos[i] + vel[i+1]*dt

        print(f"Final velocity is {vel[i]} km/s")
mission = SpaceMission(33382)
dt = 1e-12
fuel_mass = 200
dv = 1000 #np.sqrt((2*c.G*mission.massses[0])/mission.radii[0]*1000**3)
steps = 1000
motors = 500
max_launch_time = 500
filename = "test_launch.pkl"

if(os.path.exists(filename) == False):
    motor = nm.nano_motor(10**-6, 10**5, 3000, dt)

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
print(f"Engine thrust {engine.thrust()}")
print(f"Fuel consumption {engine.fuel_consumption()}")
print(f"Fuel mass {fuel_mass}")
print(f"Planet radius {mission.system.radii[0]}")
print(f"Planet mass {mission.system.masses[0]}")
mission.set_launch_parameters(engine.thrust(),engine.fuel_consumption(),fuel_mass,max_launch_time,start_pos,0)
mission.launch_rocket()
