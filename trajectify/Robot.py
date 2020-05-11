class Robot:

    def __init__(self, track_width, max_velocity, max_acceleration):
        self.track_width = track_width
        self.max_velocity = max_velocity
        self.max_acceleration = max_acceleration
        self.min_acceleration = -max_acceleration

    def get_wheel_speeds_from_state(self, state):
        # get curvature: w = v(in/s) * curvature(rad/in) = rad/s
        angular_velocity = state.velocity * state.curvature

        # get velocities:
        right_velocity = state.velocity + self.track_width * angular_velocity
        left_velocity = state.velocity - self.track_width * angular_velocity

        return left_velocity, right_velocity
