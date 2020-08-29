import ast2000tools.constants as const
import numpy as np
import random

k = const.k_B
m = const.m_H2

box_side_length = 1e-05             # m
temperature = 3e3                   # K
number_of_particles_in_box = 100
timestep = 1000
time = np.linspace(0, 10e-09, timestep)

random.seed(95)

mean = 0
standard_deviation = np.sqrt(k*temperature/m)


def particles_in_box(N, L, meaa, standard_deviation):
    particle_positions = np.zeros((N, 3), float)   # Placeholder array for pos
    particle_velocities = np.zeros((N, 3), float)  # Placeholder array of vel

    for i in range(N):
        particle_positions[i, :] = random.uniform(0, L)
        particle_velocities[i, :] = random.gauss(mean, standard_deviation)

    return particle_positions, particle_velocities


def gasses_in_box(number_of_particles_in_box,
                  box_side_length,
                  temperature):
