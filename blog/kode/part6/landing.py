"""Egen kode."""
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from ast2000tools.shortcuts import SpaceMissionShortcuts

class landing:
    def __init__(self, seed=33382,
                A_para=1, A_rocket=1,
                C_para=1, C_rocket=1,
                codes=[92725]):
        """
        Initilize
        """
        position_after_launch = np.load('part1.bin', allow_pickle=True)
        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)
        self.C_rocket = C_rocket
        self.C_para = C_para
        self.A_rocket = A_rocket
        self.A_para = A_para
        self.drag_coeff_full()
        self.shortcut = SpaceMissionShortcuts(self.mission, codes)
        self.mission.verify_launch_result(position_after_launch)

    def cheat_stable_orbit(self):
        self.shortcut.place_spacecraft_in_stable_orbit(0, 30000, 0, 4)
        landing = self.mission.begin_landing_sequence()
        self.landing = landing
        landing.start_video()
        landing.finish_video(filename='landing_video.xml',
                             number_of_frames=100)

    def drag_same(self, h):
        '''
        The drag force that is the same for both instances
        '''
        drag = 0.5*atmosphear(h)*v**2
        return drag

    def drag_coeff_rocket(self):
        '''
        The drag coefficient of just the rocket
        '''
        self.coeff_rocket = self.C_rocket*self.A_rocket
        return self.coeff_rocket

    def drag_coeff_full(self):
        '''
        Drag coeff with the parachute
        '''
        self.coeff_full = self.drag_coeff_rocket() + self.C_para*self.A_para

    def drag_force_without(self, h):
        '''
        Drag force without the parachute
        '''
        F_D = self.drag_same(h)*self.coeff_rocket
        return F_D

    def drag_force_with(self, h):
        '''
        Drag force with parachute
        '''
        F_D = self.drag_same(h)*self.coeff_full
        return F_D

    def theta_pos(self, x, y, omega, t):
        '''
        Calculates the theta pos for landing
        '''
        theta = np.arctan(y/x) + omega*t
        return theta

    def t_used(self):
        '''
        Should find the t used to get to the ground
        '''
        start_pos = self.start_pos
        start_vel = self.start_vel

    def find_omega(self, T):
        '''
        Finds the rotational speed
        '''
        omega = 2*np.pi/T
        return omega

    def coordintates(self, planet_idx, original_coord, phi):
        '''
        Finds the coordinates we need to aim for
        '''
        x, y, z = original_coord[0], original_coord[1], original_coord[2]
        T = self.system.rotational_periods[planet_idx]

        t = self.t_used()

        omega = self.find_omega(T)
        theta = self.theta_pos(x, y, omega, t)

        r = self.system.radii[planet_idx]
        x = r*np.sin(phi)*np.cos(theta)
        y = r*np.sin(phi)*np.sin(theta)
        z = r*np.sin(phi)
        coord = [x, y, z]
        coord = np.asarray(coord)

        return coord


if __name__ == '__main__':
    test = landing()
    test.cheat_stable_orbit()
