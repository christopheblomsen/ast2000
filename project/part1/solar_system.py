import ast2000tools.utils as utils
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
