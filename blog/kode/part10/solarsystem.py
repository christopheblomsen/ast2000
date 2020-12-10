from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c
system = SolarSystem(33382)

print('My system has a {:g} solar mass star with a radius of {:g} kilometers, and a surface temperatur of {:g} K.'
      .format(system.star_mass, system.star_radius, system.star_temperature))

R = system.star_radius*1000
T = system.star_temperature
M = system.star_mass

L = 4*c.pi*R**2*c.sigma*T**4
print('Luminosity is {:g} Watts' .format(L))

print(f'Spectral height {L/c.L_sun}')
