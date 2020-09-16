# egen kode
from scipy.integrate import quad
import numpy as np
import math
import ast2000tools.constants as c

def gaussian_distibution(v, mean, stdd):
    return (1/(math.sqrt(2*math.pi)*stdd))*np.exp(-0.5*((v-mean)**2/(2*stdd**2)))

def maxwell_boltzmann_distribution(v,m,T):
    return (m/(2*math.pi*c.k_B*T))**(3/2)*np.exp(-0.5*(m*v**2)/(c.k_B*T))

def integrate(func,a,b,c,d):
    return quad(func,a,b, args=(c,d))
