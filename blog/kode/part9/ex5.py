# Egen kode
import numpy as np
from ast2000tools.relativity import RelativityExperiments
import ast2000tools.constants as c

seed = 33382

experiments = RelativityExperiments(seed)
planet_idx = 4  # I want to perform the experiment near planet 4

#experiments.gps(planet_idx)

MASS = 1.840138835198e23        # kg
RADIUS = 1929.1055309           # km
R = RADIUS

t_earth = 11.9326414            # s

t_sat1 = 11.9145299
x_sat1 = -3869.124
y_sat1 = -2902.911

t_sat2 = 11.9119655
x_sat2 = -1899.305
y_sat2 = -4448.557

radius_1 = np.sqrt(x_sat1**2 + y_sat1**2)
radius_2 = np.sqrt(x_sat2**2 + y_sat2**2)


R_1 = radius_1-RADIUS
R_2 = radius_2-RADIUS
print(f'{R_1:g}, {R_2:g}')
print(f'{R_1}, {R_2}')


"""
(2pi R/v_tan)**2 = (4pi**2 R**3)/(GM)

v_tan = np.sqrt(G*M/R)
"""

v_tan_1 = np.sqrt(c.G*MASS/(R_1))
v_tan_2 = np.sqrt(c.G*MASS/(R_2))

print(v_tan_1, v_tan_2)
conts = (c.c**2*(t_earth-t_sat1)**2 - R_1**2 - R**2)/(2*R*R_1)
print(conts)
alpha = np.arccos(conts)
print(alpha)
