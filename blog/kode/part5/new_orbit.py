# Egen kode
import numpy as np


class new_orbit:
    def __init__(self):
        self.find_energy()
        self.find_a()
        self.h = h # insert spin her
        self.find_p()
        self.find_b()
        self.r_max()
        self.r_min()

    def find_a(self):
        a = - (m_r*m_p)/(2*self.find_energy())
        self.a = a

    def find_energy(self):
        kin = .5*self.mu*v**2
        pot = -self.G(m_r*m_p)/(np.linalg.norm(r))
        self.E = kin + pot

    def find_e(self):
        b = self.b
        a = self.a
        e = np.sqrt(1 - (b/a)**2)
        self.e = e

    def find_p(self):
        h = self.h
        p = h**2/(self.G(m_p + m_r))
        self.p = p

    def find_b(self):
        b = np.sqrt(self.p*self.a)
        self.b = b

    def r_analytical(self, f):
        a = self.a
        e = self.r
        r = (a*(1 - e**2))/(1 + e*np.cos(f))
        return r

    def r_min(self):
        r_min = self.r_analytical(np.pi*.5)
        self.r_min = r_min

    def r_max(self):
        r_max = self.r_analytical(0)
        self.r_max = r_max
