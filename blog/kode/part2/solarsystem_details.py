# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem

system = SolarSystem(33382)

print('<table>')
print(f'<tr><th>Planet #</th><th>Radius m</th><th>Mass kg</th><th>Semi major axes AU</th><th>Eccentrisity AU</th></tr>')
for i in range(len(system.radii)):
    '''
    List planet details
    '''
    print(f'<tr><td>{i:02}</td><td>{system.radii[i]:g}</td><td>{system.masses[i]:g}</td><td>{system.semi_major_axes[i]:g}</td><td>{system.eccentricities[i]:g}</td></tr>')

print('</table>')

print(f'Stjerne masse {system.star_mass:g} Sol masser, Overflate temperatur {system.star_temperature:g} K, radius {system.star_radius:g} km')
