import numpy as np
from ast2000tools.star_population import StarPopulation
import matplotlib.pyplot as plt


stars = StarPopulation()

T = stars.temperatures # [K]

L = stars.luminosities # [L_sun]

r = stars.radii        # [R_sun]


c = stars.colors

s = np.maximum(1e3*(r - r.min())/(r.max() - r.min()), 1.0) # Make point areas proportional to star radii


fig, ax = plt.subplots()

ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)


ax.set_xlabel('Temperature [K]')

ax.invert_xaxis()

ax.set_xscale('log')

ax.set_xticks([35000, 18000, 10000, 6000, 4000, 3000])

ax.set_xticklabels(list(map(str, ax.get_xticks())))

ax.set_xlim(40000, 2000)

ax.minorticks_off()


ax.set_ylabel(r'Luminosity [$L/L_\odot$]')

ax.set_yscale('log')

ax.set_ylim(1e-4, 1e6)


plt.savefig('HR_diagram.png')
