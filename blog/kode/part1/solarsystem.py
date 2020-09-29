import ast2000tools.utils as utils
import ast2000tools.constants as c
import math
import numpy as np

from ast2000tools.solar_system import SolarSystem
system = SolarSystem(33382)

print('My system has a {:g} solar mass star with a radius of {:g} kilometers.'
      .format(system.star_mass, system.star_radius))

for planet_idx in range(system.number_of_planets):
    print('Planet {:g} is a {} planet with a semi-major axis of {:g} AU. and a distance {:g}'
          .format(planet_idx + 1, system.types[planet_idx], system.semi_major_axes[planet_idx], np.linalg.norm(system.initial_positions[:, planet_idx])))
    if(planet_idx == 0):
        print('My planet has Mass {:.2g} kg,  Radius {:.2f} km, Escape velocity {:.4f} km/s'
          .format( system.masses[planet_idx]*c.m_sun,system.radii[planet_idx],math.sqrt((2*c.G*system.masses[planet_idx]*c.m_sun)/(system.radii[planet_idx]*1000**3))))
        print('')

print(f'Earth escape velocity {math.sqrt((2*c.G*5.972*10**24)/(6371*1000**3)):g} km/s')
#times, planet_positions = ... # Your own orbit simulation code
#system.generate_orbit_video(times, planet_positions, filename='orbit_video.xml')
