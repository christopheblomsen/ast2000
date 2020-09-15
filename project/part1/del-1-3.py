# egen kode
import numpy as np
import math
import ast2000tools.constants as c
import matplotlib.pyplot as plt
import distribution
import sys

a = int(-2.5*10**4)
b = int(2.5*10**4)
d = int((b-a)/500)
T = 3000

mean = 0
stdd = math.sqrt((c.k_B*T)/c.m_H2) # standardavviket gitt temperatur og masse
print(f'Standdardavviket for fordelingen av hastigheten i gassen er {stdd}')
vx = []
bins = []

# Oppgave 3.1
res, err = distribution.integrate(lambda x,m,T: x*distribution.maxwell_boltzmann_distribution(x,m,T),0,float('inf'),c.m_H2,T)
print(f'Den absolutte middelhastigheten i en en H2 gass med temperatur {T}K er {res:g} m/s')

sys.exit()
