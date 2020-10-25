"""Egen kode."""
# from ast2000tools.space_mission import InterplanetaryTravel
# from spacecraft import Spacecraft
'''
These classes are commands that can be executed to control
the journey of our spacecraft.
There is a superclass that defines the interface.
This is based on the Command pattern from the gang of four:
https://archive.org/details/designpatternsel00gamm/page/232/mode/2up
'''


class Command:
    """Abstract Super class defining command interface."""

    def __init__(self):
        """Initialize the Command class."""
        pass

    def execute(self, args):
        """Execute the command."""
        return False


class Boost(Command):
    """Command to bost engine with velocity vector dv."""

    def __init__(self, dv, spacecraft):
        """Construct the Boost command.

        dv is the vector for change in velocity
        interplanetary_travel is the interplanetary_travel object from the
        AST2000tools package
        """
        self.dv = dv
        self.spacecraft = spacecraft

    def execute(self, args):
        """Execute the boost command on the interplanetary_travel object."""
        self.spacecraft.boost(self.dv)
        return True


class Coast(Command):
    """Coast the spacecraft for a given timeperiode."""

    def __init__(self, duration, spacecraft):
        """Initialize the Coast class.

        duration is time in years to coast.
        interplanetary_travel is the receiver of the command.
        """
        self.duration = duration
        self.spacecraft = spacecraft

    def execute(self, args):
        """Execute the boost command on the interplanetary_travel object."""
        self.spacecraft.caost(self.duration)
        return True


class Launch(Command):
    """Execute the launch sequence."""

    def __init__(self, t, angle, spacecraft):
        """Initialize the Launch class.

        t is the time in years
        angle is the equatorial angle on the planet to launch from.
        spacecraft is an instance of the Spacecraft class.
        mission is the mission object from ast2000tools
        """
        self.t = t
        self.angle = angle
        self.spacecraft = spacecraft
        self.mission = spacecraft.mission
        self.max_launch_time = 1000

    def execute(self, args):
        """Execute launch command on spacecraft."""
        self.spacecraft.launch_process(self.t, self.angle)

        return True


class CorrectionalBoost(Command):
    """Perform a correctional boost based on the error vector supplied.

    The error vector is a vector pointing from the current position to
    the planned position.
    """

    def __init__(self, error, spacecraft):
        """Initialize CorrectionalBoost."""
        self.error = error
        self.spacecraft = spacecraft

    def execute(self, args):
        """Execute command, Not implemented yet."""
        pass
