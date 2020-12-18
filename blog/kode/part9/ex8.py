# Egen kode
import numpy as np
import matplotlib.pyplot as plt

def V(r):
    return np.sqrt((1 - 2*M/r)/(r**2))

M = 1
r = np.linspace(0, 40*M, 1000)

plt.plot(r, V(r))
plt.title('M = 1')
plt.xlabel('r')
plt.ylabel('V(r)')
plt.show()
