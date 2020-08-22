import numpy as np


def trapezoidal(f, a, b, n):
    h = (b-a)/(n)
    s = 0.5*(f(mu, sigma, a) + f(mu, sigma, b))
    for i in range(1, n, 1):
        s = s + f(a + i*h)
    return h*s


def f(mu, sigma, x):
    return 1/(np.sqrt(2*np.pi)*sigma) * np.exp(-0.5*((x-mu)/sigma)**2)


mu, sigma, x = 1, 1, 1

a = - sigma
b = sigma

P = trapezoidal(f(mu, sigma, x), a, b, 10)
print(P)
