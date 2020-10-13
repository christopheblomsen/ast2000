# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from PIL import Image


class generate_pic:
    def __init__(self, seed=33382, image='himmelkule.npy'):
        '''
        Need to figure this out
        '''
        self.system = SolarSystem(seed)
        self.sky_image = np.load(image)

    def field_of_view(self, alpha_phi=70, alpha_theta=70):
        '''
        Caclulates the x_min/max
        and the y_min/max
        '''
        x_numerator = 2*np.sin(alpha_phi/2)
        x_denomerator = 1 + np.cos(alpha_phi/2)

        y_numerator = 2*np.sin(alpha_theta/2)
        y_denomerator = 1 + np.cos(alpha_theta/2)

        self.x_max = x_numerator/x_denomerator
        self.x_min = - self.x_max

        self.y_max = y_numerator/y_denomerator
        self.y_min = -self.y_max

    def sample(self, image='sample0000.png'):
        img = Image.open(image)
        pixels = np.array(img)
        width = len(pixels[0, :])
        height = len(pixels[:, 0])
        return width, height

    def full_range(self):
        X, Y = self.sample()
        x, y = np.meshgrid(X, Y)
        plt.plot(x, y)

pic = generate_pic()
print(pic.sample())
pic.full_range()
plt.show()
