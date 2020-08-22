import numpy as np
import matplotlib.pyplot as plt


def integral(g, a, x, N=20):
    index_set = range(N+1)
    x = np.linspace(a, x, N+1)
    g_ = np.zeros_like(x)
    f = np.zeros_like(x)
    g_[0] = g(x[0])
    f[0] = 0

    for n in index_set[1:]:
        g_[n] = g(x[n])
        f[n] = f[n-1] + 0.5*(x[n] - x[n-1])*(g_[n-1] + g_[n])
    return x, f


def gauss():
    def g(t):
        return 1.0/np.sqrt(2*np.pi)*np.exp(-0.5*t**2)
    x, f = integral(g, a=-3, x=3, N=200)
    integrand = g(x)
    plt.plot(x, f)
    plt.show()

gauss()

mu, sigma, x0 = 1, 1, 0

u, t = runge_kutta2(f, x0, 20, 5)

plt.plot(t, u)
plt.show()
