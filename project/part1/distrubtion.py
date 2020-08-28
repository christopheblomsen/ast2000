import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as conts

v_x = np.linspace(-2.5e4, 2.5e4, 1000)
m = conts.m_H2
k = conts.k_B
T = 3000  # K
N = 10E5


def probability_density(v):
    ans = ((m)/(2*np.pi*k*T))**(0.5) * np.exp(-0.5*((m*v**2)/(k*T)))
    return ans


def integration(f, a, b, N):
    dx = (b - a)/N
    s = 0
    x = a
    for i in range(N):
        s += (f(x) + f(x+dx))/2*dx
        x += dx
    return s


Pv_x = probability_density(v_x)

plt.plot(v_x, Pv_x)
plt.show()
plt.close()

integrated_Pv_x = integration(probability_density, -5e3, 30e3, 1000)
print(f'{integrated_Pv_x:g}')


def velocity_distrubtion(v):
    ans = (m/(2*np.pi*k*T))**(3/2) * np.exp(-0.5*(m*v**2)/(k*T))*4*np.pi*v**2
    return ans


v = np.linspace(0, 3e4, 1000)

Pv = velocity_distrubtion(v)

plt.plot(v, Pv)
plt.show()
