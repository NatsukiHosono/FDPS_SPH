from math import pi, sqrt, cos, sin, atan, asin
import pandas as pd
import numpy as np
from copy import copy
import matplotlib.pyplot as plt

def center_of_mass(x_coords, y_coords, z_coords, masses, particle_tags=None, target_iron=False, target_iron_id=1):
    if target_iron:
        masses = np.array([i for index, i in enumerate(masses) if particle_tags[index] == target_iron_id],
                          dtype=np.float32)
        x_coords = np.array([i for index, i in enumerate(x_coords) if particle_tags[index] == target_iron_id],
                            dtype=np.float32)
        y_coords = np.array([i for index, i in enumerate(y_coords) if particle_tags[index] == target_iron_id],
                            dtype=np.float32)
        z_coords = np.array([i for index, i in enumerate(z_coords) if particle_tags[index] == target_iron_id],
                            dtype=np.float32)
    else:
        masses = np.array(masses, dtype=np.float32)
    total_mass = float(np.sum(masses))
    x_center = sum([a * b for a, b in zip(x_coords, masses)]) / total_mass
    y_center = sum([a * b for a, b in zip(y_coords, masses)]) / total_mass
    z_center = sum([a * b for a, b in zip(z_coords, masses)]) / total_mass
    return x_center, y_center, z_center

class ReverseTime:
    """
    Reversing the initial position of the impactor backwards in time.
    """

    def __init__(self, target_file_path, impactor_file_path, impact_parameter, dt=-0.1, center_target=False,
                 v_esc_multiple=1.0):
        self.G = 6.674 * 10 ** -11
        self.target_df = pd.read_csv(target_file_path, sep='\t', skiprows=2, header=None).dropna()
        self.impactor_df = pd.read_csv(impactor_file_path, sep='\t', skiprows=2, header=None).dropna()
        self.center_target = center_target
        self.target_mass = sum([float(i) for i in self.target_df[2]])
        self.impactor_mass = sum([float(i) for i in self.impactor_df[2]])
        self.total_mass = self.target_mass + self.impactor_mass
        self.imp_total_mass_ratio = self.impactor_mass / self.total_mass
        self.com_target = self.center_of_mass(self.target_df)
        self.com_impactor = self.center_of_mass(self.impactor_df)
        self.target_centered = self.center_body(
            x=self.target_df[3],
            y=self.target_df[4],
            z=self.target_df[5],
            com=self.com_target
        )
        self.impactor_centered = self.center_body(
            x=self.impactor_df[3],
            y=self.impactor_df[4],
            z=self.impactor_df[5],
            com=self.com_impactor
        )
        self.radius_target = self.radius_body(centered_coords=self.target_centered)
        self.radius_impactor = self.radius_body(centered_coords=self.impactor_centered)

        self.v_esc = v_esc_multiple * sqrt((2 * self.G * self.total_mass) / (self.radius_target + self.radius_impactor))
        self.v_target_x = (self.impactor_mass / self.total_mass) * self.v_esc
        self.v_impactor_x = - (self.target_mass / self.total_mass) * self.v_esc
        self.v_target_y = 0.0
        self.v_impactor_y = -0.0
        self.v_target_z = 0.0
        self.v_impactor_z = -0.0

        self.angle = self.impact_angle(impact_parameter=impact_parameter)

        self.time = 0
        self.dt = dt

        self.current_target_position = copy(self.target_centered)
        self.current_impactor_position = self.initial_impactor_position(centered_coords=self.impactor_centered)
        self.com_target = center_of_mass(
            x_coords=self.current_target_position[0],
            y_coords=self.current_target_position[1],
            z_coords=self.current_target_position[2],
            masses=self.target_df[2]
        )
        self.com_impactor = center_of_mass(
            x_coords=self.current_impactor_position[0],
            y_coords=self.current_impactor_position[1],
            z_coords=self.current_impactor_position[2],
            masses=self.impactor_df[2]
        )
        self.calculate_distance()

        self.v_target_history = [[self.v_target_x, self.v_target_y, self.v_target_z]]
        self.v_impactor_history = [[self.v_impactor_x, self.v_impactor_y, self.v_impactor_z]]
        self.position_target_history = [self.com_target]
        self.position_impactor_history = [self.com_impactor]
        self.time_history = [0]


    def center_of_mass(self, df):
	df[2]=df[2].astype(float)
        total_mass = sum(df[2])
        x = sum(df[2] * df[3]) / total_mass
        y = sum(df[2] * df[4]) / total_mass
        z = sum(df[2] * df[5]) / total_mass
        return x, y, z

    def center_body(self, x, y, z, com):
        x = [i - com[0] for i in x]
        y = [i - com[1] for i in y]
        z = [i - com[2] for i in z]
        return [x, y, z]

    def radius_body(self, centered_coords):
        reformatted_coords = [[centered_coords[0][index], centered_coords[1][index], centered_coords[2][index]]
                              for index, i in enumerate(centered_coords[0])]
        return max([np.linalg.norm(i) for i in reformatted_coords])

    def impact_angle(self, impact_parameter):
        return asin(impact_parameter) * (180 / pi)

    def initial_impactor_position(self, centered_coords):
        x = [i + ((self.radius_target + self.radius_impactor) * cos(self.angle * (pi / 180))) for i in
             centered_coords[0]]
        y = [i + ((self.radius_target + self.radius_impactor) * sin(self.angle * (pi / 180))) for i in
             centered_coords[1]]
        z = centered_coords[0]
        return [x, y, z]

    def plot_current_position(self):
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot(111)
        ax.scatter(
            self.current_target_position[0],
            self.current_target_position[1],
            marker='+',
            color='blue',
            label="TARGET"
        )
        ax.scatter(
            self.current_impactor_position[0],
            self.current_impactor_position[1],
            marker='+',
            color='red',
            label="IMPACTOR"
        )
        ax.plot(
            [i[0] for i in self.position_target_history],
            [i[1] for i in self.position_target_history],
            linewidth=2.0,
            color='yellow',
            label="PATH TARGET"
        )
        ax.plot(
            [i[0] for i in self.position_impactor_history],
            [i[1] for i in self.position_impactor_history],
            linewidth=2.0,
            color='green',
            label="PATH IMPACTOR"
        )
        ax.plot(
            [self.com_target[0], self.com_impactor[0]],
            [self.com_target[1], self.com_impactor[1]],
            linewidth=2.0,
            color='black',
            label="Angle: {}".format(
                round(
                    atan((self.com_impactor[1] - self.com_target[1]) / (self.com_impactor[0] - self.com_target[0])) * (
                            180 / pi), 4)
            )
        )
        ax.plot(
            [self.com_target[0], self.com_impactor[0]],
            [self.com_target[1], self.com_target[1]],
            linestyle="--",
            color='black'
        )
        ax.plot(
            [self.com_impactor[0], self.com_impactor[0]],
            [self.com_target[1], self.com_impactor[1]],
            linestyle="--",
            color='black'
        )
        # ax.plot(
        #     [self.com_target[0], self.com_target[0] + self.radius_target],
        #     [self.com_target[1], self.com_target[1]],
        #     color='green',
        #     linewidth=2.0
        # )
        # ax.plot(
        #     [self.com_impactor[0], self.com_impactor[0] + self.radius_impactor],
        #     [self.com_impactor[1], self.com_impactor[1]],
        #     color='green',
        #     linewidth=2.0,
        # )
        ax.set_xlim(-2.1e7, 2.1e7)
        ax.set_ylim(-2.1e7, 2.1e7)
        ax.set_ylim()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Time Before Impact: {} s".format(self.time))
        ax.grid()
        ax.legend(loc='upper left')
        return fig

    def plot_velocity_history(self):
        fig = plt.figure(figsize=(16, 9))
        ax_x = fig.add_subplot(311)
        ax_y = fig.add_subplot(312)
        ax_z = fig.add_subplot(313)
        ax_x.plot(
            self.time_history,
            [i[0] for i in self.v_target_history],
            linewidth=2.0,
            color="blue",
            label="TARGET"
        )
        ax_x.plot(
            self.time_history,
            [i[0] for i in self.v_impactor_history],
            linewidth=2.0,
            color="red",
            label="IMPACTOR"
        )
        ax_y.plot(
            self.time_history,
            [i[1] for i in self.v_target_history],
            linewidth=2.0,
            color="blue",
            label="TARGET"
        )
        ax_y.plot(
            self.time_history,
            [i[1] for i in self.v_impactor_history],
            linewidth=2.0,
            color="red",
            label="IMPACTOR"
        )
        ax_z.plot(
            self.time_history,
            [i[2] for i in self.v_target_history],
            linewidth=2.0,
            color="blue",
            label="TARGET"
        )
        ax_z.plot(
            self.time_history,
            [i[2] for i in self.v_impactor_history],
            linewidth=2.0,
            color="red",
            label="IMPACTOR"
        )
        ax_x.set_xlabel("Time (s)")
        ax_y.set_xlabel("Time (s)")
        ax_z.set_xlabel("Time (s)")
        ax_x.set_ylabel("v_x")
        ax_y.set_ylabel("v_y")
        ax_z.set_ylabel("v_z")
        ax_x.grid()
        ax_y.grid()
        ax_z.grid()
        ax_x.legend()

        return fig

    def calculate_distance(self):
        """
        Distance between target and impactor COM.
        :return:
        """
        self.x_distance = self.com_impactor[0] - self.com_target[0]
        self.y_distance = self.com_impactor[1] - self.com_target[1]
        self.z_distance = self.com_impactor[2] - self.com_target[2]
        self.distance = sqrt(
            self.x_distance ** 2 + self.y_distance ** 2 + self.z_distance ** 2
        )

    def delta_v(self):
        # grav acceleration exerted by impactor on target, calculate velocity change, v = a * dt
        self.target_delta_v_x = ((self.G * self.impactor_mass) / (self.distance ** 3)) * self.x_distance * self.dt
        self.target_delta_v_y = ((self.G * self.impactor_mass) / (self.distance ** 3)) * self.y_distance * self.dt
        self.target_delta_v_z = ((self.G * self.impactor_mass) / (self.distance ** 3)) * self.z_distance * self.dt
        self.impactor_delta_v_x = (self.target_mass / self.impactor_mass) * (-self.target_delta_v_x)
        self.impactor_delta_v_y = (self.target_mass / self.impactor_mass) * (-self.target_delta_v_y)
        self.impactor_delta_v_z = (self.target_mass / self.impactor_mass) * (-self.target_delta_v_z)

        self.v_target_x += self.target_delta_v_x
        self.v_target_y += self.target_delta_v_y
        self.v_target_z += self.target_delta_v_z
        self.v_impactor_x += self.impactor_delta_v_x
        self.v_impactor_y += self.impactor_delta_v_y
        self.v_impactor_z += self.impactor_delta_v_z

        self.v_target_history.append([self.v_target_x, self.v_target_y, self.v_target_z])
        self.v_impactor_history.append([self.v_impactor_x, self.v_impactor_y, self.v_impactor_z])

    def move(self):
        # x = v_x * dt
        self.current_target_position[0] = [i + (self.v_target_x * self.dt) for i in self.current_target_position[0]]
        self.current_target_position[1] = [i + (self.v_target_y * self.dt) for i in self.current_target_position[1]]
        self.current_target_position[2] = [i + (self.v_target_z * self.dt) for i in self.current_target_position[2]]

        self.current_impactor_position[0] = [i + (self.v_impactor_x * self.dt) for i in
                                             self.current_impactor_position[0]]
        self.current_impactor_position[1] = [i + (self.v_impactor_y * self.dt) for i in
                                             self.current_impactor_position[1]]
        self.current_impactor_position[2] = [i + (self.v_impactor_z * self.dt) for i in
                                             self.current_impactor_position[2]]

        self.com_target = center_of_mass(
            x_coords=self.current_target_position[0],
            y_coords=self.current_target_position[1],
            z_coords=self.current_target_position[2],
            masses=self.target_df[2]
        )
        self.com_impactor = center_of_mass(
            x_coords=self.current_impactor_position[0],
            y_coords=self.current_impactor_position[1],
            z_coords=self.current_impactor_position[2],
            masses=self.impactor_df[2]
        )

        if self.center_target:
            self.current_target_position[0] = [i - self.com_target[0] for i in self.current_target_position[0]]
            self.current_target_position[1] = [i - self.com_target[1] for i in self.current_target_position[1]]
            self.current_target_position[2] = [i - self.com_target[2] for i in self.current_target_position[2]]

            self.current_impactor_position[0] = [i - self.com_target[0] for i in self.current_impactor_position[0]]
            self.current_impactor_position[1] = [i - self.com_target[1] for i in self.current_impactor_position[1]]
            self.current_impactor_position[2] = [i - self.com_target[2] for i in self.current_impactor_position[2]]

            self.com_impactor = [i - j for i, j in zip(self.com_impactor, self.com_target)]
            self.com_target = [i - j for i, j in zip(self.com_target, self.com_target)]

        self.position_target_history.append(self.com_target)
        self.position_impactor_history.append(self.com_impactor)

    def reverse(self):
        # recalulate target-impactor distance
        self.calculate_distance()
        # recalculate velocities of target and impactor
        self.delta_v()
        # reposition target and impactor
        self.move()
        # plot
        # self.plot_current_position()
        # change time
        self.time += self.dt
        self.time_history.append(self.time)
