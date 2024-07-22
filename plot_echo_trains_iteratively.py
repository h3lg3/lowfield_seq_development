from enum import Enum
from math import pi

import numpy as np
import pypulseq as pp

from console.interfaces.interface_acquisition_parameter import Dimensions
from console.utilities.sequences.system_settings import raster, system

from matplotlib import pyplot as plt
class Trajectory(Enum):
    """Trajectory type enum."""

    INOUT = 1
    LINEAR = 2

echo_time = 15e-3
repetition_time = 600e-3
etl = 7
dummies = 0
rf_duration = 400e-6
ramp_duration= 200e-6
gradient_correction = 0.
ro_bandwidth = 20e3
echo_shift= 0.0
trajectory: Trajectory = Trajectory.INOUT
excitation_angle = pi / 2
excitation_phase = 0.
refocussing_angle = pi
refocussing_phase = pi / 2
channel_ro, channel_pe1, channel_pe2 = 'x', 'y', 'z'

# %% 
fov_ro, fov_pe1, fov_pe2 = 220e-3, 220e-3, 220e-3
n_enc_ro, n_enc_pe1, n_enc_pe2 = 16, 16, 32
system.rf_ringdown_time = 0
seq = pp.Sequence(system)
seq.set_definition("Name", "tse_3d")

# Calculate center out trajectory
pe1 = np.arange(n_enc_pe1) - (n_enc_pe1 - 1) / 2
pe2 = np.arange(n_enc_pe2) - (n_enc_pe2 - 1) / 2

pe0_pos = np.arange(n_enc_pe1)
pe1_pos = np.arange(n_enc_pe2)

pe_points = np.stack([grid.flatten() for grid in np.meshgrid(pe1, pe2)], axis=-1)
pe_positions = np.stack([grid.flatten() for grid in np.meshgrid(pe0_pos, pe1_pos)], axis=-1)

pe_mag = np.sum(np.square(pe_points), axis=-1)  # calculate magnitude of all gradient combinations
pe_mag_sorted = np.argsort(pe_mag)

if trajectory is Trajectory.INOUT:
    pe_mag_sorted = np.flip(pe_mag_sorted)

pe_traj = pe_points[pe_mag_sorted, :]  # sort the points based on magnitude
pe_order = pe_positions[pe_mag_sorted, :]  # kspace position for each of the gradients

if trajectory is Trajectory.LINEAR:
    center_pos = 1 / 2  # where the center of kspace should be in the echo train
    num_points = np.size(pe_mag_sorted)
    linear_pos = np.zeros(num_points, dtype=int) - 10
    center_point = int(np.round(np.size(pe_mag) * center_pos))
    odd_indices = 1
    even_indices = 1
    linear_pos[center_point] = pe_mag_sorted[0]

    for idx in range(1, num_points):
        # check if its in bounds first
        if center_point - (idx + 1) / 2 >= 0 and idx % 2:
            k_idx = center_point - odd_indices
            odd_indices += 1
        elif center_point + idx / 2 < num_points and idx % 2 == 0:
            k_idx = center_point + even_indices
            even_indices += 1
        elif center_point - (idx + 1) / 2 < 0 and idx % 2:
            k_idx = center_point + even_indices
            even_indices += 1
        elif center_point + idx / 2 >= num_points and idx % 2 == 0:
            k_idx = center_point - odd_indices
            odd_indices += 1
        else:
            print("Sorting error")
        linear_pos[k_idx] = pe_mag_sorted[idx]

    pe_traj = pe_points[linear_pos, :]  # sort the points based on magnitude
    pe_order = pe_positions[linear_pos, :]  # kspace position for each of the gradients

# calculate the required gradient area for each k-point
pe_traj[:, 0] /= fov_pe1
pe_traj[:, 1] /= fov_pe2

# Divide all PE steps into echo trains
num_trains = int(np.ceil(pe_traj.shape[0] / etl))
trains = [pe_traj[k::num_trains, :] for k in range(num_trains)]

# Create a list with the kspace location of every line of kspace acquired, in the order it is acquired
trains_pos = [pe_order[k::num_trains, :] for k in range(num_trains)]
acq_pos = []
for train_pos in trains_pos:
    acq_pos.extend(train_pos)
    
# Plot k-space acquisition order
plt.ion()  # Turn on interactive mode
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

for array in trains[0:10]:
    #ax1.clear()
    ax2.clear()
    
    # Plot the whole array on the first subplot
    ax1.scatter(array[:, 0], array[:, 1])
    ax1.set_title('Array-wise Plot')
    
    ax1.set_xlim(min(pe_traj[:, 0]), max(pe_traj[:, 0]))
    ax1.set_ylim(min(pe_traj[:, 1]), max(pe_traj[:, 1]))
    ax2.set_xlim(min(pe_traj[:, 0]), max(pe_traj[:, 0]))
    ax2.set_ylim(min(pe_traj[:, 1]), max(pe_traj[:, 1]))
    # Plot each element of the array iteratively on the second subplot
    for point in array:
        ax2.scatter(point[0], point[1], color='r')
        ax2.set_title('Element-wise Plot')
        plt.draw()
        plt.pause(0.2)  # Pause for half a second

plt.ioff()  # Turn off interactive mode
plt.show()    