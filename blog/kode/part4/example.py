import numpy as np
import matplotlib.pyplot as plt

a = 0
b = 1
r = 2

theta = np.linspace(0, 2*np.pi, 100)

def g(a, b, r):
    if a == 0 and b == 0:
        x = r*np.cos(theta)
        y = r*np.sin(theta)
    elif b == 0:
        x = r*np.cos(theta) + a
        y = r*np.sin(theta)
    elif a == 0:
        x = r*np.cos(theta)
        y = r*np.sin(theta) + b
    else:
        x = (a**2 + b**2 - 2*b*r*np.sin(theta))/(2*a)
        y = (- 2*a*r*np.cos(theta) + a**2 + b**2)/(2*b)
    plt.plot(x, y, label=f'Circle center = ({a}, {b})')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.plot(a, b, 'o')
    plt.grid()

g(0, 0, 1)
g(1, 0, 1)
g(0, -1, 1)
plt.legend()
plt.show()
