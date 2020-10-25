"""Egen kode."""
import numpy as np
import math
import unittest
import matplotlib.pyplot as plt


def cartesian_polar(r):
    """Convert to polar coordinates."""
    x = r[0]                             # x values
    y = r[1]                             # y values
    r_p = np.linalg.norm(r)

    theta = np.arctan2(y, x)              # theta
    return r_p, theta


def polar_cartesian(r, theta):
    """Convert to cartesian."""
    x = r*math.cos(theta)        # calculates x
    y = r*math.sin(theta)        # calculates y
    return x, y


def leapfrog(r0, v0, T, dt, f):
    """Leapfrog integration."""
    N = int(T//dt)                            # Length of all our vectors
    print(f'{N} = {T}/{dt}')
    t = np.zeros(N, float)                   # time array
    r = np.zeros((N, 2), float)              # Position vector
    v = np.zeros((N, 2), float)              # Velocity vector
    a = np.zeros((N, 2), float)              # Acceleration vector

    r[0, :] = r0
    v[0, :] = v0
    t[0] = 0

    a[0, :] = f(r[0, :], 0)  # * r[0, :]
    for i in range(N-1):
        '''
        The actual leapfrog algorithm
        '''
        r[i + 1, :] = r[i, :] + v[i, :]*dt + 0.5*a[i, :]*dt**2

        a[i + 1, :] = f(r[i + 1, :], t[i])  # * r[i + 1, :]

        v[i + 1, :] = v[i, :] + 0.5*(a[i, :] + a[i + 1, :])*dt

        t[i + 1] = t[i] + dt

    return r, v, a, t


class TestCoordinateConversion(unittest.TestCase):
    """Class to test module."""

    def test_polar_cartesian(self):
        """Test method for polar to cartesioan conversion."""
        x, y = polar_cartesian(1, 0)
        self.assertEqual(x, 1)
        self.assertEqual(y, 0)
        x, y = polar_cartesian(1, math.pi/2)
        self.assertEqual(round(x, 4), 0)
        self.assertEqual(round(y, 4), 1)
        x, y = polar_cartesian(1, math.pi)
        self.assertEqual(round(x, 4), -1)
        self.assertEqual(round(y, 4), 0)
        x, y = polar_cartesian(1, 3*math.pi/2)
        self.assertEqual(round(x, 4), 0)
        self.assertEqual(round(y, 4), -1)

    def test_cartesian_polar(self):
        """Test method for cartesian to polar conversion."""
        r, t = cartesian_polar(np.array([1, 0]))
        self.assertEqual(r, 1)
        self.assertEqual(t, 0)
        r, t = cartesian_polar(np.array([0, 1]))
        self.assertEqual(r, 1)
        self.assertEqual(t, math.pi/2)
        r, t = cartesian_polar(np.array([-1, 0]))
        self.assertEqual(r, 1)
        self.assertEqual(t, math.pi)
        r, t = cartesian_polar(np.array([0, -1]))
        self.assertEqual(r, 1)
        self.assertEqual(t, 3*math.pi/2)


class TestLeapFrog(unittest.TestCase):
    """Class for testing Leafrog integration."""

    def test_leapfrog(self):
        """Test ethod to test Leapfrog."""
        def f(r, t):
            return np.array([-9.8, 0.3/(t+1)])

        r, v, a, t = leapfrog(np.array([1000, 0]), 0, 10, 1, f)
        plt.scatter(r[:, 1], r[:, 0])
        plt.show()
        print(r)
        print('')
        print(t)


if __name__ == '__main__':
    unittest.main()
