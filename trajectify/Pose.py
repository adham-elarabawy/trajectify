import math


class Pose:
    """
    A 2D Pose that contains x,y displacement with heading

    ...

    Attributes
    ----------
    x : double
        x position in inches
    y : double
        y position in inches
    theta :
        angle/heading in radians
    dydx :
        angle/heading in slope form

    Methods
    -------
    """

    def __init__(self, x, y, theta):
        self.x = x
        self.y = y
        self.theta = math.radians(theta)
        self.dydx = math.tan(math.radians(theta))

    @staticmethod
    def lerp(pose0, pose1, t):
        end_minus_start = [pose1.x - pose0.x, pose1.y - pose0.y]
        times = [element * t for element in end_minus_start]

        return [pose0.x + times[0], pose0.y + times[1]]
