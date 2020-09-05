# Egen kode
import numpy as np
import nano_motor as nm
import rocket_engine as re
import ast2000tools.utils as utils
import ast2000tools.constants as c
from ast2000tools.solar_system import SolarSystem
import sys
import os
try:
    import cPickle as pickle
except:
    import pickle

"""
class rocket_launch:
    # A class for calculating the launch of our rocket
    def __init__(self, mission, initial_mass):
        self.mission = mission
        self.initial_mass = initial_mass
        self.system = SolarSystem(33382)


    def launch_process(self, L, N, T, dt):
        self.motor = nm(L, N, T, dt)

        v_esc = np.sqrt((2*c.G*self.system.masses[0]*c.m_sun)/(self.system.radii[0]*1000**3))

        pos = np.zeros(1000, float)
        vel = np.zeros(1000, float)
        p = np.zeros(1000, float)

        pos0 = np.array([0, 0, 0])
        vel = np.array([0, 0, 0])

        h = 1

        time = np.linspace(0, 1000, 1000)

        M = self.system.masses[0]

        for i in range(len(time)):
            if v[0] >= v_esc:
                print("Achieved escape velocity")
                print(f"Our height is {pos[1,i]} meters and it took {time[i]} seconds")
                print(f"It took {fuel_consumption}")
                sys.exit()
            else:
                fuel_consumption = self.motor.fuel_consumption()
                p[i + 1] = self.motor.calculated_thrust()
                current_mass = self.initial_mass - fuel_consumption
                try:
                    F_G = -(c.G*self.initial_mass*M)/(pos[i]**3)
                    a = (p[i + 1] - p[i])/(m*(time[i + 1] - time[i])) - (c.G*M)/(pos[i]**3)
                except ZeroDivisionError:
                    F_G = -(c.G*self.initial_mass*M/h**3)
                    a = (p[i + 1] - p[i])/(m*(time[i + 1] - time[i])) - (c.G*M)/(h**3)

                # Euler-cromer
                vel[i+1] = vel[i] + a*dt
                pos[i+1] = pos[i] + vel[i+1]*dt

"""
mission = SolarSystem(33382)
dt = 1e-12
fuel_mass = 2500
dv = 1000 #np.sqrt((2*c.G*mission.massses[0])/mission.radii[0]*1000**3)
steps = 1000
motors = 10**17
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

engine = re(motors, motor)
fuel_consumed = engine.boost(mission.spacecraft_mass+fuel_mass, dv)
