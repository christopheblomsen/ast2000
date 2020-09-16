# egen kode
import numpy as np
import math
import ast2000tools.constants as c
import matplotlib.pyplot as plt
import distribution

a = int(-2.5*10**4)
b = int(2.5*10**4)
d = int((b-a)/500)

mean = 0
stdd = math.sqrt((c.k_B*3000)/c.m_H2) # standardavviket gitt temperatur og masse
print(f'Standdardavviket for fordelingen av hastigheten i gassen er {stdd}')
vx = []
bins = []

# Oppgave 2.1
# Finner sansynlighetsfordeling for hastigheten vx over intervallet a - b, med delta lik d
for v in range(a,b,d):
    res, err = distribution.integrate(distribution.gaussian_distibution,v,v+d,mean,stdd)
    vx.append(res)
    bins.append((v + v +d)/2)

plt.bar(bins,vx,width=d)
plt.show()

# Oppgave 2.2 sannsynligheten for å finne vx i området [5,30]*10^3
res, err = distribution.integrate(distribution.gaussian_distibution,5*10**3,30*10**3,mean,stdd)

print(f'Sannsynligheten for vx i området [5,30]*10^3 er {res:.4f}')

a = 0
b = int(3*10**4)
d = int(b/500)
vx = []
bins = []

# Oppgave 2.3
for v in range(a,b,d):
    res, err = distribution.integrate(distribution.gaussian_distibution,v,v+d,mean,stdd)
    vx.append(res)
    bins.append((v + v +d)/2)

plt.bar(bins,vx,width=d)
plt.show()
