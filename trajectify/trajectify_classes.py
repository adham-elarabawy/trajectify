import math

import numpy as np


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


class QuinticSpline:
    """
    An individual quintic hermite spline

    ...

    Attributes
    ----------
    pose0 : Pose
        2D Pose for the 0th point in the spline
    pose1 : Pose
        2D Pose for the 1th point in the spline
    x_control_vector : numpy array
        vector (length 6) describing: initial x pos, initial x vel, initial x accel, final x pos, final x vel, final x accel
    y_control_vector : numpy array
        vector (length 6) describing: initial y pos, initial y vel, initial y accel, final y pos, final y vel, final y accel

    Methods
    -------
    """

    hermite_basis = np.array([[-6, 15, -10, 0, 0, 1],
                              [-3, 8, -6, 0, 1, 0],
                              [-0.5, 1.5, -1.5, 0.5, 0, 0],
                              [6, -15, 10, 0, 0, 0],
                              [-3, 7, -4, 0, 0, 0],
                              [0.5, -1, 0.5, 0, 0, 0]])

    hermite_basis_d = np.array([[0, -30, 60, -30, 0, 0],
                                [0, -15, 32, -18, 0, 1],
                                [0, -2.5, 6, -4.5, 1, 0],
                                [0, 30, -60, 30, 0, 0],
                                [0, -15, 28, -12, 0, 0],
                                [0, 2.5, -4, 1.5, 0, 0]])

    hermite_basis_dd = np.array([[0, 0, -120, 180, -60, 0],
                                 [0, 0, -60, 96, -36, 0],
                                 [0, 0, -10, 18, -9, 1],
                                 [0, 0, 120, -180, 60, 0],
                                 [0, 0, -60, 84, -24, 0],
                                 [0, 0, 10, -12, 3, 0]])

    def __init__(self, pose0, pose1, safety_scaling=1.3):
        self.pose0 = pose0
        self.pose1 = pose1
        self.safety_scaling = 1

        euclidian_distance = safety_scaling * \
            math.sqrt((pose1.x - pose0.x)**2 + (pose1.y - pose0.y)**2)

        vx0 = math.cos(pose0.theta) * euclidian_distance
        vx1 = math.cos(pose1.theta) * euclidian_distance
        ax0 = 0
        ax1 = 0

        self.x_control_vector = np.array(
            [pose0.x, vx0, ax0, pose1.x, vx1, ax1])

        vy0 = math.sin(pose0.theta) * euclidian_distance
        vy1 = math.sin(pose1.theta) * euclidian_distance
        ay0 = 0
        ay1 = 0

        self.y_control_vector = np.array(
            [pose0.y, vy0, ay0, pose1.y, vy1, ay1])

    @staticmethod
    def get_hermite_vector(t, d=0):
        """returns the hermite vector of length 6: [h0(t), h1(t), h2(t), h3(t), h4(t), h5(t)] with each element evaluated at t"""
        assert ((d >= 0) and (
            d <= 2)), "Attempted to evaluate a derivative greater than available hermite basis (or a negative derivative)"
        assert ((t >= 0) and (t <= 1)
                ), "Attempted to extrapolate out of the region of spline"
        t_vector = np.array([t**5, t**4, t**3, t**2, t, 1])
        if d == 0:
            return QuinticSpline.hermite_basis.dot(t_vector)
        if d == 1:
            return QuinticSpline.hermite_basis_d.dot(t_vector)
        if d == 2:
            return QuinticSpline.hermite_basis_dd.dot(t_vector)

    def evaluate(self, t, d=0):
        """returns the point on the trajectory by evaluating x(t) and y(t) at provided t parameter value (0<=t<=1)"""
        assert ((d >= 0) and (
            d <= 2)), "Attempted to evaluate a derivative greater than available hermite basis (or a negative derivative)"
        assert ((t >= 0) and (t <= 1)
                ), "Attempted to extrapolate out of the region of spline"
        hermite_vector = QuinticSpline.get_hermite_vector(t, d)
        return np.array([hermite_vector.dot(self.x_control_vector), hermite_vector.dot(self.y_control_vector)])

    def compute_curvature(self, t):
        return ((self.evaluate(t, 1)[0] * self.evaluate(t, 2)[1]) - (self.evaluate(t, 2)[0] * self.evaluate(t, 1)[1])) / (math.sqrt((self.evaluate(t, 1)[0]**2 + self.evaluate(t, 1)[1]**2)**3))


class Path:

    def __init__(self, waypoints):
        assert len(
            waypoints) > 1, "Path cannot be generated with only one waypoint."
        self.waypoints = waypoints
        self.num_waypoints = len(waypoints)

        self.splines = []

        for i, waypoint in enumerate(waypoints):
            if (i < self.num_waypoints - 1):
                self.splines.append(QuinticSpline(
                    waypoints[i], waypoints[i + 1]))

    def map_parameter(self, t):
        return t * (len(self.splines))

    def get_spline(self, t):
        assert ((t >= 0) and (t <= 1)), "Attempted to extrapolate out of the Path"
        normalized_t = self.map_parameter(t)
        spline_index = int(normalized_t)
        spline_local_t = normalized_t - spline_index

        if spline_index == len(self.splines):
            spline_index = len(self.splines) - 1
            spline_local_t = 1

        return self.splines[spline_index], spline_local_t

    def evaluate(self, t, d=0):
        assert ((t >= 0) and (t <= 1)), "Attempted to extrapolate out of the Path"

        spline, local_t = self.get_spline(t)
        return spline.evaluate(local_t, d)

    def compute_curvature(self, t):
        assert ((t >= 0) and (t <= 1)), "Attempted to extrapolate out of the Path"

        spline, local_t = self.get_spline(t)
        return spline.compute_curvature(local_t)

    def theta(self, t):
        """returns radians"""
        path_deriv = self.evaluate(t, 1)
        dydt = path_deriv[1]
        dxdt = path_deriv[0]
        slope = dydt / dxdt

        return math.atan(slope)

    def get_plot_values(self, d=0, resolution=100):
        t = np.linspace(0, 1, num=resolution)
        x, y = [], []
        for step in t:
            point = self.evaluate(step, d)
            x.append(point[0])
            y.append(point[1])
        return x, y

    @staticmethod
    def get_distance_between(point0, point1):
        return math.sqrt((point0[0] - point1[0])**2 + (point0[1] - point1[1])**2)

    @staticmethod
    def transform(pose0, pose1):
        initial_translation = [pose0.x, pose0.y]
        last_translation = [pose1.x, pose1.y]
        initial_rotation = [math.cos(pose0.theta), math.sin(pose0.theta)]
        last_rotation = [math.cos(pose1.theta), math.sin(pose1.theta)]

        initial_unary = [math.cos(math.radians(-math.degrees(pose0.theta))),
                         math.sin(math.radians(-math.degrees(pose0.theta)))]

        matrix0 = [last_translation[0] - initial_translation[0],
                   last_translation[1] - initial_translation[1]]

        m_translation = [matrix0[0] * initial_unary[0] - matrix0[1] * initial_unary[1],
                         matrix0[0] * initial_unary[1] + matrix0[1] * initial_unary[0]]
        m_rotation = [last_rotation[0] * initial_unary[0] - last_rotation[1] * initial_unary[1],
                      last_rotation[1] * initial_unary[0] + last_rotation[0] * initial_unary[1]]

        # normalize rotation matrix
        magnitude = math.sqrt(m_rotation[0]**2 + m_rotation[1]**2)
        if magnitude > 10**-9:
            m_rotation[0] /= magnitude
            m_rotation[1] /= magnitude
        else:
            m_rotation[0] = 1
            m_rotation[1] = 0

        return m_translation, m_rotation

    @staticmethod
    def twistify(pose0, pose1):
        transform_translation, transform_rotation = Path.transform(
            pose0, pose1)
        dtheta = math.atan2(transform_rotation[1], transform_rotation[0])

        half_dtheta = dtheta / 2
        cos_minus_one = transform_rotation[0] - 1

        if (abs(cos_minus_one) < 10**-9):
            half_theta_by_tan_of_half_dtheta = 1 - 1 / 12 * dtheta * dtheta
        else:
            half_theta_by_tan_of_half_dtheta = - \
                (half_dtheta * transform_rotation[1]) / cos_minus_one

        # rotation

        rotate_by = [half_theta_by_tan_of_half_dtheta, -half_dtheta]
        times_by = math.sqrt(
            half_theta_by_tan_of_half_dtheta**2 + half_dtheta**2)

        rotated = [transform_translation[0] * rotate_by[0] - transform_translation[1] * rotate_by[1],
                   transform_translation[0] * rotate_by[1] + transform_translation[1] * rotate_by[0]]
        final = [rotated[0] * times_by, rotated[1] * times_by]

        return final[0], final[1], dtheta


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


class ConstrainedState:

    def __init__(self, t, distance, max_velocity, min_acceleration, max_acceleration):
        self.t = t
        self.distance = distance
        self.max_velocity = max_velocity
        self.min_acceleration = min_acceleration
        self.max_acceleration = max_acceleration

    def __str__(self):
        return str([self.t, self.distance, self.max_velocity, self.min_acceleration, self.max_acceleration])


class Trajectory:

    max_dx = 0.127 * 39.3701
    max_dy = 0.00127 * 39.3701
    max_dtheta = 0.0872

    kEpsilon = 10**-6

    def integrand(self, t):
        deriv_point = self.path.evaluate(t, 1)
        dx = deriv_point[0]
        dy = deriv_point[1]
        return math.sqrt((dx)**2 + (dy)**2)

    def __init__(self, robot, path, v_initial=0, v_final=0, a_initial=0, spline_resolution=1000, max_trajectory_time=10, min_trajectory_time=1, optimization_dt=0.1):
        self.robot = robot
        self.path = path
        self.v_initial = v_initial
        self.v_final = v_final
        self.a_initial = a_initial
        self.step_size = 1 / (spline_resolution * len(self.path.splines))
        self.max_trajectory_time = max_trajectory_time
        self.min_trajectory_time = min_trajectory_time
        self.optimization_dt = optimization_dt
        # self.total_arc_length = quad(self.integrand, 0, 1)[0]
        self.control_points = []
        self.trajectory = []

        done = False
        t0 = 0
        t1 = 1
        count = 0

        # recursive subdivision
        self.control_points.append(t0)
        while not done:
            point0 = self.path.evaluate(t0)
            point1 = self.path.evaluate(t1)
            theta0 = self.path.theta(t0)
            theta1 = self.path.theta(t1)

            dx, dy, dtheta = Path.twistify(
                Pose(point0[0], point0[1], theta0), Pose(point1[0], point1[1], theta1))

            if (abs(dx) <= Trajectory.max_dx) and (abs(dy) <= Trajectory.max_dy) and (abs(dtheta) <= Trajectory.max_dtheta):
                self.control_points.append(t1)
                if t1 >= 1:
                    done = True
                t0 = t1
                t1 = 1
                count += 1
                print(t0, end='\r')
            else:
                t1 = (t1 + t0) / 2

        print('\n Found all %f control points, moving on to time-parameterization.\n' %
              len(self.control_points))

        states = []  # ConstrainedState: t, distance, max_velocity, min_acceleration, max_acceleration

        # STEP 1: forward pass
        for i, t in enumerate(self.control_points):
            states.append(ConstrainedState(t, 0, 0, 0, 0))

            if t == 0:
                ds = 0
                states[i].max_velocity = v_initial
            else:
                ds = Path.get_distance_between(self.path.evaluate(
                    states[i].t), self.path.evaluate(states[i - 1].t))
                states[i].distance = states[i - 1].distance + ds

            while True:
                states[i].min_acceleration = self.robot.min_acceleration
                states[i].max_acceleration = self.robot.max_acceleration

                if (ds < Trajectory.kEpsilon):
                    break

                # for state in states:
                #     print(state)

                # vf = sqrt(vi^2 + 2ad)
                temp_velocity = math.sqrt(
                    states[i - 1].max_velocity**2 + (states[i - 1].max_acceleration * 2 * ds))
                states[i].max_velocity = min(
                    self.robot.max_velocity, temp_velocity)

                # if the actual acceleration for this state is higher than the max acceleration that we applied, then we need to reduce the maximum acceleration of the predecessor and try again.
                actual_accel = (states[i].max_velocity**2 -
                                states[i - 1].max_velocity**2) / (2 * ds)

                # if we violate the max acceleration constraint, let's modify the predecessor
                if(states[i].max_acceleration < (actual_accel - Trajectory.kEpsilon)):
                    states[i - 1].max_acceleration = states[i].max_acceleration
                else:
                    if (actual_accel > states[i - 1].min_acceleration):
                        states[i - 1].max_acceleration = actual_accel
                    break

        # ConstrainedState: t, distance, max_velocity, min_acceleration, max_acceleration
        # STEP 2: backward pass
        for i, (reversed_i, control_point_t) in enumerate(reversed(list(enumerate(self.control_points)))):
            if i == 0:
                states[reversed_i].max_velocity = v_final
                ds = 0
            else:
                ds = states[reversed_i].distance - \
                    states[reversed_i + 1].distance  # negative

            while True:
                if ds > -Trajectory.kEpsilon:
                    break

                # enforce max velocity limit (reverse): vf = sqrt(vi^2 + 2*a*d), where vi is the control point after this one (chronologically on the trajectory)
                new_max_velocity = math.sqrt(
                    states[reversed_i + 1].max_velocity**2 + states[reversed_i + 1].min_acceleration * 2 * ds)

                # this state can be finalized if this is true
                if (new_max_velocity >= states[reversed_i].max_velocity):
                    break

                states[reversed_i].max_velocity = new_max_velocity

                # if the actual acceleration for this state is lower than the minimum acceleration, we need to lower the minimum acceleration of the control point after this one and try again
                actual_accel = (states[reversed_i].max_velocity**2 -
                                states[reversed_i + 1].max_velocity**2) / (2 * ds)

                if (states[reversed_i].min_acceleration > actual_accel + Trajectory.kEpsilon):
                    states[reversed_i +
                           1].min_acceleration = states[reversed_i].min_acceleration
                else:
                    states[reversed_i + 1].min_acceleration = actual_accel
                    break

        # STEP 3: integrate the constrained states forward in time to get the trajectory states
        # t, distance, velocity, acceleration, curvature
        time = 0  # seconds
        distance = 0  # inches
        velocity = 0  # inches/second

        for i, state in enumerate(states):

            # calculate the change in position between the current state and the previous state
            ds = state.distance - distance

            accel = 0

            # calcualte dt
            dt = 0
            if (i > 0):
                # calculate the acceleration between the current state and the previous state
                accel = (state.max_velocity**2 - velocity**2) / (2 * ds)
                self.trajectory[i - 1].acceleration = accel
                if abs(accel) > Trajectory.kEpsilon:
                    # vf = v0 + a * t
                    dt = (state.max_velocity - velocity) / accel
                elif abs(velocity) > Trajectory.kEpsilon:
                    # delta_x = v * t
                    dt = ds / velocity
                else:
                    raise RuntimeError

            velocity = state.max_velocity
            distance = state.distance
            time += dt

            self.trajectory.append(State(state.t, time, distance, Pose(self.path.evaluate(state.t)[0], self.path.evaluate(
                state.t)[1], self.path.theta(state.t)), velocity, accel, self.path.compute_curvature(state.t)))
        self.total_time = self.trajectory[-1].time

    @staticmethod
    def lerp(start_value, end_value, t):
        """
        linearly interpolates between two values

        start_value : double
            the start value
        end_value : double
            the end value
        t : double
            the fraction for interpolation
        """
        return start_value + (end_value - start_value) * t

    def interpolate(self, prev_state, end_state, i):
        """
        prev_state : State
            ith state
        end_state : State
            i+1 th state
        i :
            interpolant (fraction)
        """

        # find new t value
        # this might be using the wrong interpolant
        new_t = Trajectory.lerp(prev_state.t, end_state.t, i)

        # find new time value
        new_time = Trajectory.lerp(prev_state.time, end_state.time, i)

        # find the delta time between the current state and the interpolated state
        delta_t = new_time - prev_state.time

        # if delta_time is negative, flip the order of interpolation
        if delta_t < 0:
            return self.interpolate(end_state, prev_state, 1 - i)

        # calculate the new velocity: vf = v0 + a*t
        new_v = prev_state.velocity + (prev_state.acceleration * delta_t)

        # calculate the chane in position
        # delta_s = v0 * t + 0.5 * a * t^2
        new_s = prev_state.velocity * delta_t + \
            0.5 * prev_state.acceleration * delta_t**2

        # to find the new position for the new state, we need to interpolate between the two endpoint poses. The freaction for interpolation is the chagne in position (delta_s) divided by the total distance between the two endpoints.
        interpolation_frac = new_s / \
            (math.sqrt((prev_state.pose.x - end_state.pose.x)
                       ** 2 + (prev_state.pose.y - end_state.pose.y)**2))

        new_pos = Pose.lerp(
            prev_state.pose, end_state.pose, interpolation_frac)
        new_pose = Pose(new_pos[0], new_pos[1], self.path.theta(new_t))
        new_curvature = Trajectory.lerp(
            prev_state.curvature, end_state.curvature, interpolation_frac)

        return State(new_t, new_time, prev_state.distance + new_s, new_pose, new_v, prev_state.acceleration, new_curvature)

    def sample(self, time):
        # if the sample time is smaller than 0, then return the first state
        if (time <= self.trajectory[0].time):
            return self.trajectory[0]

        # if the sample time is larger than the total trajectory time, then return the last state
        if (time >= self.total_time):
            return self.trajectory[-1]

        # to get the element we want, we will use a binary search algorithm instead of iterating using a for loop.

        # we start at 1 ebcause we use the previous state later on for interpolation

        low = 1
        high = len(self.trajectory) - 1

        while not low == high:
            mid = int((low + high) / 2)
            if self.trajectory[mid].time < time:
                low = mid + 1
            else:
                high = mid

        # if we reach this point, high and low must be the same

        # the sample's timestamp is now >= to the requested timestep. If it is greater, we need to interpolate between the previous state and the current state to get the exact state that we want
        sample = self.trajectory[low]
        prev_sample = self.trajectory[low - 1]

        # if the difference between the states is below kEpsilon, then we are pretty much spot on and can return this state
        if abs(sample.time - prev_sample.time) < 10**-9:
            return sample

        interpolated_sample = self.interpolate(
            prev_sample, sample, (time - prev_sample.time) / (sample.time - prev_sample.time))
        return interpolated_sample
