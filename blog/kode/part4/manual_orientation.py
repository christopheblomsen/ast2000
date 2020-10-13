# Egen kode
import numpy as np
from ast2000tools.space_mission import SpaceMission
import astrogation_computer
import trilateration

ast_comp = astrogation_computer.astrogation_computer()

mission = SpaceMission.load('part1.bin')

mission.take_picture()

dopler_shifts = mission.measure_star_doppler_shifts()

print(f'Doppler shifts {dopler_shifts}')

distances = np.array(mission.measure_distances())

print(f'Planet and star distances\n{distances}')

angle = ast_comp.find_orientation_angle('sky_picture.png')
print(f'Spacecraft is pointing in direction {angle} degrees')

trilateration = trilateration.trilateration(mission)

velocity = trilateration.radial_velocity()

print(f'The spacecraft velocity is {velocity}')

position = trilateration.tri_test(distances)

print(f'Spacecraft position is {position}')

mission.verify_manual_orientation(position, velocity, angle)

'''
Test run
$ python manual_orientation.py
SpaceMission instance loaded from part1.bin.
Picture written to sky_picture.png.
Doppler shifts (-0.0994931407937419, 0.046595090974288794)
Planet and star distances
[4.87204793e-05 2.49507657e-01 2.85664654e+00 1.84465903e+00
 1.37842272e+00 4.74243099e+00 5.88339176e-01 4.16691335e-01
 3.65860973e-01]
Comparing to angle  359
Spacecraft is pointing in direction 37 degrees
The spacecraft velocity is [2.43963171 8.92174606]
Shortest distance from us [1 8 7]
(2, 9)
Planet positions: [[ 0.36581226  0.57485571  3.04552524  2.20274777 -0.0909102  -0.5333506
   0.51367347 -0.01680852  0.        ]
 [ 0.          0.13633947  0.98994901  0.16923892 -1.30046735 -4.65631179
  -0.56941104  0.1647003   0.        ]]
Positions [[ 0.57485571  0.         -0.01680852]
 [ 0.13633947  0.          0.1647003 ]]
coordinate candidates (6, 2)
Length between 0 and 2 0.00834484457537991
Length between 0 and 5 0.007044370148718203
Length between 2 and 0 0.00834484457537991
Length between 2 and 5 0.003475441017506861
Length between 5 and 0 0.007044370148718203
Length between 5 and 2 0.003475441017506861
{0: [0, 2, 5], 1: [1], 2: [2, 0, 5], 3: [3], 4: [4], 5: [5, 0, 2]}
x: [0.3705584800760562 0.3658609725371775 0.3645506144236]
y: [-0.006897090249476601 0.0 -0.003218951984784918]
x: [0.3658609725371775 0.3705584800760562 0.3645506144236]
y: [0.0 -0.006897090249476601 -0.003218951984784918]
x: [0.3645506144236 0.3705584800760562 0.3658609725371775]
y: [-0.003218951984784918 -0.006897090249476601 0.0]
Spacecraft position is [0.36699002234561123, -0.003372014078087173]
Pointing angle after launch correctly calculated. Well done!
Velocity after launch correctly calculated. Well done!
Position after launch correctly calculated. Well done!
Your manually inferred orientation was satisfyingly calculated. Well done!
*** Achievement unlocked: Well-oriented! ***

'''

