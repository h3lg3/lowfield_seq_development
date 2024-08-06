import numpy as np
from enum import Enum
import math

class Trajectory(Enum):
    """Trajectory type enum."""

    OUTIN = 1   # sampling pattern is in fact OUTIT and should be renamed accordingly
    LINEAR = 2
    INOUT = 3
    SYMMETRIC = 4
    
    # Define trajectories corresponding to https://colab.research.google.com/github/pulseq/MR-Physics-with-Pulseq/blob/main/tutorials/03_k_space_sampling/notebooks/01_cartesian_ordering_and_undersampling_solution.ipynb#scrollTo=zYC6t2eOCt_L


def constructor(
    n_enc_pe1: int = 128,
    n_enc_pe2: int = 128,
    etl: int = 8,
    trajectory: Trajectory = Trajectory.INOUT,
) -> tuple[list, list]:


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
        pe_traj = pe_points[pe_mag_sorted, :]  # sort the points based on magnitude
        pe_order = pe_positions[pe_mag_sorted, :]  # kspace position for each of the gradients
    elif trajectory is Trajectory.OUTIN:
        pe_mag_sorted = np.flip(pe_mag_sorted)
        pe_traj = pe_points[pe_mag_sorted, :]  # sort the points based on magnitude
        pe_order = pe_positions[pe_mag_sorted, :]  # kspace position for each of the gradients
    elif trajectory is Trajectory.LINEAR:
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
    elif trajectory is Trajectory.SYMMETRIC:    # just example for symmetric encoding given n_pe2 = 1, why is scheme so different from linear?
        # PE dir 1
        n_ex = math.floor(n_enc_pe1 / etl)
        pe_steps = np.arange(1, etl * n_ex + 1) - 0.5 * etl * n_ex - 1
        if divmod(etl, 2)[1] == 0:
            pe_steps = np.roll(pe_steps, [0, int(-np.round(n_ex / 2))])
        pe_traj = np.array([[pe_steps[i],0.0] for i in np.arange(etl * n_ex)])      
        
    return (pe_traj, pe_order)