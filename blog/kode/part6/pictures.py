# This is only for generating pictures for the blog since I can't draw

import numpy as np
import matplotlib.pyplot as plt

def r_ellipticali(p, a, c, theta):
    return p/(a - c*np.cos(theta))

def r_circle(r, theta):
    return r*np.cos(theta)

theta = np.linspace(0, 2*np.pi, 1000)
r = 1
p = 1
a = 1
c = 1

plt.polar(r_ellipticali(p, a, c, theta))
plt.show()

