# Egen kode
import numpy as np
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission


class new_orbit:
    """
    This class is for the finding the new orbit
    it will need the bin file for times and
    positions and velocity in polar coordinates
    """
    def __init__(self, planet_id, bin_file, seed=33382):
        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)
        self.times, self.r, self.v = np.load(bin_file, allow_picle=True)
        self.m_p = self.system.masses[planet_id]
        self.m_r = self.mission.rocket_mass
        self.v_tangetial = self.v[1]

        self.find_energy()
        self.find_a()
        self.h = self.r*self.v_tangetial  # insert spin her
        self.find_p()
        self.find_b()
        self.r_max()
        self.r_min()

    def find_a(self):
        """
        Finds the new semi major axis
        """
        a = - (self.m_r*self.m_p)/(2*self.find_energy())
        self.a = a

    def find_energy(self):
        """
        Finds the new energy in the system
        """
        kin = .5*self.mu*self.v**2
        pot = -self.G(self.m_r*self.m_p)/(np.linalg.norm(self.r))
        self.E = kin + pot

    def find_e(self):
        """
        Finds the new ecentrisity
        """
        b = self.b
        a = self.a
        e = np.sqrt(1 - (b/a)**2)
        self.e = e

    def find_p(self):
        h = self.h
        p = h**2/(self.G(self.m_p + self.m_r))
        self.p = p

    def find_b(self):
        """
        Finds the new minor axis
        """
        b = np.sqrt(self.p*self.a)
        self.b = b

    def r_analytical(self, f):
        """
        analytical form of the ellipses
        """
        a = self.a
        e = self.r
        r = (a*(1 - e**2))/(1 + e*np.cos(f))
        return r

    def r_min(self):
        """
        Finds the periapsis
        """
        r_min = self.r_analytical(np.pi*.5)
        self.r_min = r_min

    def r_max(self):
        """
        Finds the apoapsis
        """
        r_max = self.r_analytical(0)
        self.r_max = r_max
