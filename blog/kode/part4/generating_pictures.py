# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from PIL import Image

class generate_pic(self, seed=33382):
    def __init__(self):
        '''
        Need to figure this out
        '''
        self.system = SolarSystem(seed)

    def open_picture(self, picture):
        '''
        Will open pictures
        '''

    def field_of_view(self):
        '''
        Caclulates the x_min/max
        and the y_min/max
        '''
        x_numerator = 2*np.sin(alpha_phi/2)
        x_denomerator = 1 + np.cos(alpha_phi/2)

        y_numerator = 2*np.sin(alpha_theta/2)
        y_denomerator = 1 + np.cos(alpha_theta/2)

        x_max = x_numerator/x_denomerator
        x_min = - x_max

        y_max = y_numerator/y_denomerator
        y_min = -y_max




