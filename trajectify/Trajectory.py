import math

from trajectify.ConstrainedState import ConstrainedState
from trajectify.Path import Path
from trajectify.Pose import Pose
from trajectify.State import State


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
