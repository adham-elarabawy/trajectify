class State:

    def __init__(self, t, time, distance, pose, velocity, acceleration, curvature):
        self.t = t
        self.time = time
        self.distance = distance
        self.pose = pose
        self.velocity = velocity
        self.acceleration = acceleration
        self.curvature = curvature

    def __str__(self):
        return str([self.t, self.time, self.distance, self.pose, self.velocity, self.acceleration, self.curvature])
