# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from PIL import Image
import logging
import sys
import os
import astrogation_computer
from ast2000tools.shortcuts import SpaceMissionShortcuts
import trilateration

ast_comp = astrogation_computer.astrogation_computer()

#mission = SpaceMission(33382)

mission = SpaceMission.load('part1.bin')

mission.take_picture()

dopler_shifts = mission.measure_star_doppler_shifts()

print(f'dopler_shifts {dopler_shifts}')

distances = np.array(mission.measure_distances())

print(f'Planet and star distances\n{distances}')

angle = ast_comp.find_orientation_angle('sky_picture.png')
print(f'Spacecraft is pointing in direction {angle} degrees')

trilateration = trilateration.trilateration(mission)

velocity = trilateration.radial_velocity()

print(f'The spacecraft velocity is {velocity}')

position = trilateration.tri_test(distances)

print(f'Spacecraft position is {position}')

mission.verify_manual_orientation(position,velocity,angle)
