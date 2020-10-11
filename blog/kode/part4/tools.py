# Egen kode
import numpy as np
import math
import unittest

def cartesian_polar(r):
    '''
    Converts to polar coordinates
    '''
    x = r[0]                             # x values
    y = r[1]                             # y values
    r_p = np.linalg.norm(r)

    theta = np.arctan2(y,x)              # theta
    return r_p, theta

def polar_cartesian( r, theta):
    '''
    Converts to cartesian
    '''
    x = r*math.cos(theta)        # calculates x
    y = r*math.sin(theta)        # calculates y
    return x, y

class TestCoordinateConversion(unittest.TestCase):

    def test_polar_cartesian(self):
        x,y = polar_cartesian(1,0)
        self.assertEqual(x,1)
        self.assertEqual(y,0)
        x,y = polar_cartesian(1,math.pi/2)
        self.assertEqual(round(x,4),0)
        self.assertEqual(round(y,4),1)
        x,y = polar_cartesian(1,math.pi)
        self.assertEqual(round(x,4),-1)
        self.assertEqual(round(y,4),0)
        x,y = polar_cartesian(1,3*math.pi/2)
        self.assertEqual(round(x,4),0)
        self.assertEqual(round(y,4),-1)

    def test_cartesian_polar(self):
        r,t = cartesian_polar(np.array([1,0]))
        self.assertEqual(r,1)
        self.assertEqual(t,0)
        r,t = cartesian_polar(np.array([0,1]))
        self.assertEqual(r,1)
        self.assertEqual(t,math.pi/2)
        r,t = cartesian_polar(np.array([-1,0]))
        self.assertEqual(r,1)
        self.assertEqual(t,math.pi)
        r,t = cartesian_polar(np.array([0,-1]))
        self.assertEqual(r,1)
        self.assertEqual(t,3*math.pi/2)


if __name__ == '__main__':
    unittest.main()
