import numpy as np


class integral:
    def __init__(self, Sigma, mu, a, b):
        self.Sigma = Sigma
        self.mu = mu
        self.a = a
        self.b = b

    def f(self, x):
        self.x = x
        ans = 1/(np.sqrt(2*np.pi)*self.Sigma)*np.exp(-1/2*((x - self.mu)/self.Sigma)**2)

        return ans

    def P(self):
        """
        integrasjon midtpunktsmetoden
        """
        N = 1000                # Antall tidssteg
        dx = (self.b-self.a)/N  # tidsintevall
        s = 0
        x = self.a              # Start
        for i in range(N):
            s += (self.f(x)+self.f((x+dx)))/2*dx
            x += dx
        return s
