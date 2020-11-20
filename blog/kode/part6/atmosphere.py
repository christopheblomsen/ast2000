"""Egen kode."""

from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c


class Atmosphere:
    """Models the atmosphere of a planet."""
    def __init__(self, planet, gamma=1.4, seed=33382):
        """Initialize the class."""
        self.system = SolarSystem(seed)
        self.surface_density = self.system.atmospheric_densities[planet]
        self.surface_temperature = 334.85
        self.R = self.system.radii[planet]
        self.g = -c.G * (self.system.masses[planet] * c.m_sun) / self.R**2
        self.k = self.surface_density**(1-gamma) * self.surface_temperature**gamma
        self.min_temp = self.surface_temperature / 2

    def T(self, r):
        """Calculates the temperature at height r meters."""
        t = self.k*self.density(r)**(self.gamma-1)
        T = t**(1/self.k)
        if T < self.min_temp:
            return T

        return self.min_temp

    def density(self, r):
        """Calculates the density at heigth r meters."""
        pass


if __name__ == '__main__':
    atmosphere = Atmosphere(4)
