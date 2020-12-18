# Egen kode
import numpy as np
import matplotlib.pyplot as plt


def ty_mark(ty, xy, v):
    return ty - xy*v


L0 = 200  # lyr

v = np.linspace(0.99, 0.0, 1000)
xy = np.linspace(L0, 2*L0, 1000)
ty = np.linspace(202, 296, 1000)

ty_m = ty_mark(ty, xy, v)
plt.plot(ty_m, ty, label="$t_{Y'} = t_y - x_Y v$")
plt.xlabel("$t_{Y'}[yr]$")
plt.ylabel("$t_Y$[yr]")
plt.legend()
plt.show()
plt.show()
