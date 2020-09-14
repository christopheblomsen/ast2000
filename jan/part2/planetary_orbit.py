# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem

system = SolarSystem(33382)


e = 0
a = []

theta = np.linspace(0, 2*np.pi, 1000)
r = lambda a, e, theta : (a*(1 - e**2))/(1 + e*np.cos(theta))
radii = system.radii

for planet_idx in range(system.number_of_planets):
    a.append(system.semi_major_axes[planet_idx])

for i in range(len(a)):
    plt.polar(theta, r(a[i], e, theta)) #, label=f'Planet number {i}')

# plt.legend()
plt.show()
