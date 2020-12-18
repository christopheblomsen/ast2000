import numpy as np
import matplotlib.pyplot as plt

def V(r):
    return np.sqrt((1 - (2*M)/r)/np.sqrt(1-r**2)-r**2)

E = 1
M = 1
m = 1

r = np.linspace(0, 202, 1000)

plt.plot(r, V(r))
plt.show()
