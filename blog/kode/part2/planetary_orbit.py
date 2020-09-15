# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem

system = SolarSystem(33382)

theta = np.linspace(0, 2*np.pi, 1000)                           # Angle running full circle with 1000 steps
r = lambda a, e, theta : (a*(1 - e**2))/(1 + e*np.cos(theta))   # Function for position

a = system.semi_major_axes                                      # assigns the semi major axes
e = system.eccentricities                                       # assigns the eccentricities

for i in range(len(a)):
    '''
    Plots all the planets
    '''
    plt.polar(theta, r(a[i], e[i], theta))

plt.title('Hoth system')
plt.show()


