import numpy as np
import unittest


class integral(unittest.TestCase):
    def __init__(self, sigma, mu, a, b):
        self.sigma = sigma
        self.mu = mu
        self.a = a
        self.b = b

    def f(self, x):
        self.x = x
        ans = 1/(np.sqrt(2*np.pi)*self.sigma)*np.exp(-1/2*((x - self.mu)/self.sigma)**2)

        return ans

    def P(self, N):
        """
        integrasjon midtpunktsmetoden
        """
        dx = (self.b-self.a)/N  # tidsintevall
        s = 0
        x = self.a              # Start
        for i in range(N):
            s += (self.f(x)+self.f((x+dx)))/2*dx
            x += dx
        return s


class test(unittest.TestCase):
    def test_P(self):
        expected1 = 0.68
        expected2 = 0.95
        expected3 = 0.997

        N = 1000

        p1 = integral(1, 0, -1, 1)
        computed1 = p1.P(N)
        msg1 = f'expected {expected1} but got {computed1}'

        p2 = integral(1, 0, -2, 2)
        computed2 = p2.P(N)
        msg2 = f'expected {expected2} but got {computed2}'

        p3 = integral(1, 0, -3, 3)
        computed3 = p3.P(N)
        msg3 = f'expected {expected3} but got {computed3}'

        self.assertTrue(np.allclose(expected1, computed1, rtol=1e-02), msg1)
        self.assertTrue(np.allclose(expected2, computed2, rtol=1e-02), msg2)
        self.assertTrue(np.allclose(expected3, computed3, rtol=1e-02), msg3)


if __name__ == '__main__':
    unittest.main()
