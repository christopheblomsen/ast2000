import numpy as np

a = -1
b = 1
sigma = 1
mu = 1


def f(x):
    ans = 1/(np.sqrt(2*np.pi) * sigma) * np.exp(-0.5*((x-mu)/sigma)**2)
    return ans

integral = np.trapz(f, a, b)

