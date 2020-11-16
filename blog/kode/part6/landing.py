"""Egen kode."""
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from ast2000tools.shortcuts import SpaceMissionShortcuts

class landing:
    def __init__(self, seed=3382,
                A_para=1, A_rocket=1,
                C_para=1, C_rocket=1,
                codes=[92725]):
        """
        Initilize
        """
        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)
        self.C_rocket = C_rocket
        self.C_para = C_para
        self.A_rocket = A_rocket
        self.A_para = A_para
        self.drag_coeff_full()
        self.shortcut = SpaceMissionShortcuts(self.mission, codes)

    def cheat_stable_orbit(self):
        self.shortcut.place_spacecraft_in_stable_orbit(0, 30000, 0, 4)
        landing = self.mission.begin_landing_sequence()
        self.landing = landing
        landing.start_video()
        landing.finish_video(filename='landing_video.xml',
                             number_of_frames=100)

    def drag_same(self, h):
        drag = 0.5*atmosphear(h)*v**2
        return drag

    def drag_coeff_rocket(self):
        self.coeff_rocket = self.C_rocket*self.A_rocket
        return self.coeff_rocket

    def drag_coeff_full(self):
        self.coeff_full = self.drag_coeff_rocket() + self.C_para*self.A_para

    def drag_force_without(self, h):
        F_D = self.drag_same(h)*self.coeff_rocket
        return F_D

    def drag_force_with(self, h):
        F_D = self.drag_same(h)*self.coeff_full
        return F_D

    def theta_pos(self, x, y, omega, t):
        theta = np.arctan(y/x) + omega*t
    
    def t_used(self):
        start_pos = self.start_pos
        start_vel = self.start_vel

if __name__ == '__main__':
    test = landing()
    test.cheat_stable_orbit()
