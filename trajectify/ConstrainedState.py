class ConstrainedState:

    def __init__(self, t, distance, max_velocity, min_acceleration, max_acceleration):
        self.t = t
        self.distance = distance
        self.max_velocity = max_velocity
        self.min_acceleration = min_acceleration
        self.max_acceleration = max_acceleration

    def __str__(self):
        return str([self.t, self.distance, self.max_velocity, self.min_acceleration, self.max_acceleration])
