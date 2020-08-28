import numpy as np
import ast2000tools.utils as utils
import ast2000tools.constants as conts

seed = utils.get_seed('chriskbl')

from ast2000tools.solar_system import SolarSystem
system = SolarSystem(seed)

print('My system has a {:g} solar mass star with a radius of {:g} kilometers.'
      .format(system.star_mass, system.star_radius))

for planet_idx in range(system.number_of_planets):
    print('Planet {:d} is a {}  planet with a semi-major axis of {:g} AU.'
          .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx]))

"""
times, planet_positions = seed
system.generate_orbit_video(times, planet_positions, filename='orbit_video.xml')
"""
from ast2000tools.space_mission import SpaceMission
mission = SpaceMission(seed)

home_planet_idx = 0 # The home planet always has index 0
print('My mission starts on planet {:d}, which has a radius of {:g} kilometers.'
      .format(home_planet_idx, mission.system.radii[home_planet_idx]))

print('My spacecraft has a mass of {:g} kg and a cross-sectional area of {:g} m^2.'
      .format(mission.spacecraft_mass, mission.spacecraft_area))

print(f'My home planet has a mass of {system.masses[0]:g}')

system.print_info()

v_esc = np.sqrt(2*system.masses[0]*9.223e18*conts.G/system.radii[0])
print(f"Escape velocity is {v_esc}")
