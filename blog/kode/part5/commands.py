# Egen kode
from ast2000tools.space_mission import InterplanetaryTravel
from spacecraft import Spacecraft
'''
These classes are commands that can be executed to control the journey of our spacecraft.
There is a superclass that defines the interface.
This is based on the Command pattern from the gang of four:
https://archive.org/details/designpatternsel00gamm/page/232/mode/2up
'''
class Command:
    '''
    Super class defining command interface
    '''
    def __init__(self):
        pass
    def execute(self, args):
        return False

class Boost(Command):
    '''
    Command to bost engine with velocity vector dv
    '''
    def __init__(self, dv, interplanetary_travel):
        self.dv = dv
        self.int_travel = interplanetary_travel

    def execute(self, args):
        self.int_travel.boost(self.dv)
        return True

class Coast(Command):
    '''
    Coast the spacecraft for a given timeperiode
    '''
    def __init__(self, duration, interplanetary_travel):
        self.duration = duration
        self.int_travel = interplanetary_travel

    def execute(self, args):
        self.int_travel.caost(self.duration)
        return True

class Launch(Command):
    '''
    Execute the launch sequence at time t and from the equatorial position angle.
    '''
    def __init__(self,t,angle, spacecraft):
        self.t = t
        self.angle = angle
        self.spacecraft = spacecraft

    def execute(self, args):
        self.spacecraft.launch_process(self.t, self.angle)
        pass

class CorrectionalBoost(Command):
    '''
    Perform a correctional boost based on the error vector supplied.
    The error vector is a vector pointing from the current position to
    the planned position.
    '''
    def __init__(self, error, interplanetary_travel):
        self.error = error
        self.int_travel = interplanetary_travel

    def execute(self, args):
        pass
