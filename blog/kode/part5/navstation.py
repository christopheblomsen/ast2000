"""Egen kode."""

import spacecraft
from nano_motor import nano_motor
from commands import Launch, Boost, Coast


class NavStation:
    """This is the Navigation controller that performs the commands to navigate.

    Add commands in the order they should be executed.
    """

    def __init__(self):
        """Initialize the NavStation with an empyt commands list."""
        self.commands = []
        pass

    def AddCommand(self, cmd):
        """Add command to commands list."""
        self.commands.append(cmd)

    def Execute(self):
        """Execute all the commands in command list."""
        for command in self.commands:
            command.execute(None)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=float, default=0.0, help='Time in years')
    parser.add_argument('-la', '--launch-angle', type=float, default=0.0,
                        help='Launch angle along equator')
    parser.add_argument('--run-steps', help='Run engine simulation')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print debug statements')
    args = parser.parse_args()

    steps = 1000
    motors = 11.0699895**15
    max_launch_time = 1000
    rocket_filename = "nano_motor.pkl"
    temperature = 3500
    fuel_mass = 43000  # Fuel mass in kg
    particles = 9.8**5
    t = args.t       # Time in years after t=0 to launch at
    angle = args.launch_angle  # Angle along equator to launch from in radians
    dt = 0.01

    nav_station = NavStation()

    rocket_engine = spacecraft.rocket_engine_factory(rocket_filename, motors,
                                                     temperature, steps,
                                                     particles, dt, args)

    space_craft = spacecraft.Spacecraft(fuel_mass, rocket_engine,
                                        dt, 33382, args.verbose)

    nav_station.AddCommand(Launch(t, angle, space_craft))

    nav_station.AddCommand(Boost([-1.6712756907642885, -0.25203293684429795], space_craft))

    nav_station.AddCommand(Coast(0.24469793490545463, spacecraft))

    nav_station.Execute()

    cur_time, cur_pos, cur_vel = space_craft.interplanetary_travel.orient()
    print(f'Position {cur_pos} and velocity {cur_vel} at {cur_time}')

    dest_pos = space_craft.planet_position(6, cur_time)
    print(f'Destination pos {dest_pos}')

"""
Kjøreeksempel:
$>  python navstation.py -t 0.24914953763595418 -la 0.8
               Launch time | 0.24915 year
              Launch angle | 0.8 radians
   Escape velocity of Hoth | 12.139587606214251 km/s
Achieved escape velocity 12.14 km/s
Our position is [787378.05040441 161312.81260592] relative to launch site and it took 345.60 seconds
It took 42519.25 kg of fuel
Spacecraft mass is now 1580.7467320600986 kg
Time is 0.24914953763595418
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time -3.899717475663866e-06
Orbital velocity is  [-0.36480495  8.81593169]
Rocket was moved up by 963.224 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 343.79 s.
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
orbital_movement 9308333.4147471 m
Final escape pos [0.3655341771636786, 0.014941524270976511] AU
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
Time is 0.24916048926398213
Planet position is [0.36550301 0.01484307] at index 10068 at time 0.24915343735342985
Planet position is [0.36549385 0.01506124] at index 10069 at time 0.24917818441713202
Difference in time between position got and actual time 7.051910552280383e-06
Orbital velocity is  [-0.36480495  8.81593169]
orbital_movement 9308333.4147471 m
Final escape pos [0.3655341771636786, 0.014941524270976511] AU
Your spacecraft position was satisfyingly calculated. Well done!
*** Achievement unlocked: No free launch! ***
Picture written to sky_picture.png.
Doppler shifts (-0.10960464349253882, 0.06693221889417213)
Planet and star distances
[4.93755691e-05 9.50555542e-01 2.89686137e+00 1.89849876e+00
 1.06073918e+00 4.72010385e+00 7.28395216e-01 4.32597359e-01
 3.65839283e-01]
Comparing to angle  359
Spacecraft is pointing in direction 37 degrees
The spacecraft velocity is [1.251074254422691,10.62333945621307]
Planet trajectories 8 planets loaded
Shortest distance from us [8 7 6]
planet_pos.shape (8, 2)
planet_pos [[ 0.36550444  0.01480869]
 [-0.58050164 -0.07766207]
 [ 2.7420217   1.67147304]
 [ 1.95497118  1.05319478]
 [ 1.00992576 -0.82764346]
 [ 0.07087376 -4.69595837]
 [ 0.21076575  0.7266676 ]
 [-0.0418933   0.15998709]]
(2, 9)
Planet positions: [[ 0.36550444 -0.58050164  2.7420217   1.95497118  1.00992576  0.07087376
   0.21076575 -0.0418933   0.        ]
 [ 0.01480869 -0.07766207  1.67147304  1.05319478 -0.82764346 -4.69595837
   0.7266676   0.15998709  0.        ]]
Positions [[ 0.         -0.0418933   0.21076575]
 [ 0.          0.15998709  0.7266676 ]]
coordinate candidates (6, 2)
Length between 1 and 2 0.004583222114380304
Length between 1 and 4 0.012812755010802033
Length between 1 and 5 0.0170099430032587
Length between 2 and 1 0.004583222114380304
Length between 2 and 4 0.012012368114612678
Length between 2 and 5 0.015554114949931563
Length between 4 and 1 0.012812755010802033
Length between 4 and 2 0.012012368114612678
Length between 4 and 5 0.027461658924070962
Length between 5 and 1 0.0170099430032587
Length between 5 and 2 0.015554114949931563
Length between 5 and 4 0.027461658924070962
{0: [0], 1: [1, 2, 4, 5], 2: [2, 1, 4, 5], 3: [3], 4: [4, 1, 2, 5], 5: [5, 1, 2, 4]}
x: [0.36572371854579994 0.36523789913991483 0.3533069816178027
 0.38013158952085263]
y: [0.009194725207393683 0.0137521262710055 0.012355632518237014
 0.018236374113447473]
x: [0.36523789913991483 0.36572371854579994 0.3533069816178027
 0.38013158952085263]
y: [0.0137521262710055 0.009194725207393683 0.012355632518237014
 0.018236374113447473]
x: [0.3533069816178027 0.36572371854579994 0.36523789913991483
 0.38013158952085263]
y: [0.012355632518237014 0.009194725207393683 0.0137521262710055
 0.018236374113447473]
x: [0.38013158952085263 0.36572371854579994 0.36523789913991483
 0.3533069816178027]
y: [0.018236374113447473 0.009194725207393683 0.0137521262710055
 0.012355632518237014]
Spacecraft position is [0.3661000472060925, 0.013384714527520918]
Pointing angle after launch correctly calculated. Well done!
Velocity after launch correctly calculated. Well done!
Position after launch correctly calculated. Well done!
Your manually inferred orientation was satisfyingly calculated. Well done!
*** Achievement unlocked: Well-oriented! ***
SpaceMission instance saved as after_launch.bin.
Final velocity is 12.140128 km/s
Spacecraft boosted with delta-v (-1.67128, -0.252033) AU/yr (1471.02 kg of fuel was used).
You ran out of fuel. Too bad!
Traceback (most recent call last):
  File "navstation.py", line 67, in <module>
    nav_station.Execute()
  File "navstation.py", line 26, in Execute
    command.execute(None)
  File "/Users/janmagneandersen/Workspace/UiO/AST2000/prosjekt/ast2000/blog/kode/part5/commands.py", line 40, in execute
    self.spacecraft.boost(self.dv)
  File "/Users/janmagneandersen/Workspace/UiO/AST2000/prosjekt/ast2000/blog/kode/part5/spacecraft.py", line 396, in boost
    self.interplanetary_travel.boost(dv)
  File "build/bdist.linux-x86_64/egg/ast2000tools/space_mission.py", line 873, in boost
RuntimeError: Out of fuel.
"""
