from ast2000tools.solar_system import SolarSystem
system = SolarSystem(33382)

print('My system has a {:g} solar mass star with a radius of {:g} kilometers, and a surface temperatur of {:g} K.'
      .format(system.star_mass, system.star_radius, system.star_temperature))

R = system.star_radius
