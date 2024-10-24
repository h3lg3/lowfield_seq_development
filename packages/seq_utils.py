import numpy as np
from enum import Enum
import math
# from console.interfaces.dimensions import Dimensions
from dataclasses import dataclass

"""Nexus Console Functions"""
"""Interface class for dimensions."""

@dataclass(frozen=True)
class Dimensions:
    """Dataclass for definition of dimensional parameters."""

    x: float | int  # pylint: disable=invalid-name
    """X dimension."""

    y: float | int  # pylint: disable=invalid-name
    """Y dimension."""

    z: float | int  # pylint: disable=invalid-name
    """Z dimension."""

# Helper function
def raster(val: float, precision: float) -> float:
    """Fit value to gradient raster.

    Parameters
    ----------
    val
        Time value to be aligned on the raster.
    precision
        Raster precision, e.g. system.grad_raster_time or system.adc_raster_time

    Returns
    -------
        Value wih given time/raster precision
    """
    # return np.round(val / precision) * precision
    gridded_val = round(val / precision) * precision
    return gridded_val
    # decimals = abs(Decimal(str(precision)).as_tuple().exponent)
    # return round(gridded_val, ndigits=decimals)
    

class Trajectory(Enum):
    """Trajectory type enum."""

    OUTIN = 1   # sampling pattern is in fact OUTIT and should be renamed accordingly
    LINEAR = 2
    INOUT = 3
    SYMMETRIC = 4
    # Define trajectories corresponding to https://colab.research.google.com/github/pulseq/MR-Physics-with-Pulseq/blob/main/tutorials/03_k_space_sampling/notebooks/01_cartesian_ordering_and_undersampling_solution.ipynb#scrollTo=zYC6t2eOCt_L

def get_traj(
    n_enc_pe1: int = 128,
    n_enc_pe2: int = 128,
    etl: int = 8,
    trajectory: Trajectory = Trajectory.INOUT,
    ) -> list:

    # Calculate center out trajectory
    pe1 = np.arange(n_enc_pe1) - (n_enc_pe1 - 1) / 2
    pe2 = np.arange(n_enc_pe2) - (n_enc_pe2 - 1) / 2
    pe_points = np.stack([grid.flatten() for grid in np.meshgrid(pe1, pe2)], axis=-1)

    pe_mag = np.sum(np.square(pe_points), axis=-1)  # calculate magnitude of all gradient combinations
    pe_mag_sorted = np.argsort(pe_mag)

    if trajectory is Trajectory.INOUT:
        pe_traj = pe_points[pe_mag_sorted, :]  # sort the points based on magnitude
    elif trajectory is Trajectory.OUTIN:
        pe_mag_sorted = np.flip(pe_mag_sorted)
        pe_traj = pe_points[pe_mag_sorted, :]  # sort the points based on magnitude
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
    elif trajectory is Trajectory.SYMMETRIC:    # symmetric encoding, why is scheme so different from linear?
        assert n_enc_pe1 % etl == 0
        # PE dir 1
        n_ex = math.floor(n_enc_pe1 / etl)
        pe_traj = np.arange(1, etl * n_ex + 1) - 0.5 * etl * n_ex - 1
      
    return pe_traj

def get_trains(
    pe_traj:list,
    n_enc_pe1:int,
    n_enc_pe2:int,
    fov_pe1:float,
    fov_pe2:float,
    etl:int,
    trajectory:Trajectory,
    ) -> list:
    
    # Divide all PE steps into echo trains
    if trajectory.name == 'SYMMETRIC':
        num_trains = int(np.ceil(n_enc_pe1/etl))    # due to slice wise acquisition, only first pe direction is divided into trains
        temp = [pe_traj[k::num_trains] for k in range(num_trains)]
        trains = []
        for k in np.arange(n_enc_pe2):
            for i in np.arange(num_trains):
                for j in np.arange(etl):
                    trains.append([temp[i][j]/fov_pe1, k/fov_pe2])
                    
        trains = [trains[k*etl:(k+1)*etl] for k in range(num_trains*n_enc_pe2)]  
    else:
        num_trains = int(np.ceil(n_enc_pe1*n_enc_pe2/etl)) # both pe directions are divided into trains
        pe_traj[:, 0] /= fov_pe1         # calculate the required gradient area for each k-point
        pe_traj[:, 1] /= fov_pe2
        trains = [pe_traj[k::num_trains, :] for k in range(num_trains)]        
                
    return trains

# Print min_esp (echo spacing) and max_etl (echo train length) and recommended etl for given PE steps
def get_esp_etl(
    tau_1: float = 3e-3,
    tau_2: float = 4e-3,
    tau_3: float = 5e-3,
    echo_time: float = 25e-3,
    T2: int = 100,              # T2 value of main tissue of interest
    n_enc_pe1: int = 64,
    ) -> dict: 
    
    # minimum echo spacing to accomodate gradients
    min_esp = -(min(tau_1, tau_2, tau_3) - echo_time/2)*2   
    min_esp = np.ceil(min_esp*1e3)*1e-3    # round to 1 ms -> new echo time
    
    # sampling duration [ms] till signal drops to 20%
    max_sampling = -math.log(0.2)*T2*1e-3  
               
    # maximum numbers of 180° echoes fitting in sampling duration     
    max_etl = int(np.floor(max_sampling/min_esp))
    
    # ETL that is multiple of n_enc_pe1 and closest to max_etl
    cc = [(i, n_enc_pe1%i, np.abs(max_etl-i)) for i in np.arange(1, 2*max_etl)]
    # Filter the list to find tuples where the second element is zero
    filtered_cc = [item for item in cc if item[1] == 0]
    # Closest ETL value that is multiple of n_enc_pe1
    rec_etl = int(min(filtered_cc, key=lambda x: x[2])[0])     

    return dict({'min_esp':min_esp, 'max_etl':max_etl, 'rec_etl':rec_etl})

def validate_inputs(
    x: str, 
    y: str, 
    z: str
) -> tuple[str, str, str]:   
    """
    Validates the input strings x, y, and z.
    This function performs the following validations:
    1. Converts each input to lowercase.
    2. Ensures each variable is of length one.
    3. Ensures each variable is one of 'x', 'y', or 'z'.
    4. Ensures all variables are unique.
    Parameters:
    x (str): The first input string.
    y (str): The second input string.
    z (str): The third input string.
    Returns:
    tuple[str, str, str]: A tuple containing the validated and converted input strings.
    Raises:
    AssertionError: If any of the validation checks fail.
    Example:
    >>> validate_inputs("x", "y", "z")
    ('x', 'y', 'z')
    >>> validate_inputs("x", "y", "y")
    AssertionError: y must be unique
    """
    # Convert each input to lowercase
    x = x.lower()
    y = y.lower()
    z = z.lower()
        
    # Ensure each variable is of length one
    assert len(x) == 1, "x must be of length 1"
    assert len(y) == 1, "y must be of length 1"
    assert len(z) == 1, "z must be of length 1"
    
    # Ensure each variable is one of 'a', 'b', or 'c'
    allowed_values = {"x", "y", "z"}
    assert x in allowed_values, "x must be one of 'x', 'y', or 'z'"
    assert y in allowed_values, "y must be one of 'x', 'y', or 'z'"
    assert z in allowed_values, "z must be one of 'x', 'y', or 'z'"
    
    # Ensure all variables are unique
    assert len({x, y, z}) == 3, "x, y, and z must be unique"

    print("All inputs are valid.")
    
    return x, y, z

    # # Example usage:
    # validate_inputs("x", "y", "z")  # This will pass validation.
    # validate_inputs("x", "y", "y")  # This will raise an AssertionError.
    
def calculate_enc_fov_order(
    channel_ro:str, 
    channel_pe1:str, 
    n_enc:Dimensions, 
    fov:Dimensions
    )-> tuple[int,float,int,float,int,float]:
    # Dynamically access the n_enc and fov attributes using getattr
    n_enc_ro = getattr(n_enc, channel_ro)
    fov_ro = getattr(fov, channel_ro)

    # Determine the remaining channel for pe2
    remaining_channels = {"x", "y", "z"} - {channel_ro, channel_pe1}
    channel_pe2 = remaining_channels.pop()

    # Assign the values for pe1 and pe2
    n_enc_pe1 = getattr(n_enc, channel_pe1)
    fov_pe1 = getattr(fov, channel_pe1)
    n_enc_pe2 = getattr(n_enc, channel_pe2)
    fov_pe2 = getattr(fov, channel_pe2)

    return (n_enc_ro, fov_ro, n_enc_pe1, fov_pe1, n_enc_pe2, fov_pe2)

    # # Example usage:
    # n_enc = Dimensions(x=64, y=80, z=100)
    # fov = Dimensions(x=200, y=220, z=240)

    # channel_ro = "x"
    # channel_pe1 = "y"
    # n_enc_ro, fov_ro, n_enc_pe1, fov_pe1, n_enc_pe2, fov_pe2 = calculate_enc_fov(channel_ro, channel_pe1, n_enc, fov)
    # print(n_enc_ro, fov_ro, n_enc_pe1, fov_pe1, n_enc_pe2, fov_pe2)    
    
def calculate_acq_pos(
    k_traj_adc:np.array,
    n_enc_ro:int, 
    channel_pe1:str, 
    fov_pe1:float, 
    channel_pe2:str, 
    fov_pe2:float
    )->list:
    xyz = ["x", "y", "z"]
    n_adc = n_enc_ro*2
    
    # RO samples
    # idx_ro = xyz.index(channel_ro)
    # k_traj_ro = -np.round(k_traj_adc[idx_ro][0:n_adc]*fov_ro*10)/10 
    # k_traj_ro = k_traj_ro - min(k_traj_ro)   
    
    # PE1 samples
    idx_pe1 = xyz.index(channel_pe1)
    k_traj_pe1 = -np.round(k_traj_adc[idx_pe1][0::n_adc]*fov_pe1*10)/10 
    k_traj_pe1 = k_traj_pe1 - min(k_traj_pe1)
    
    # PE2 samples
    idx_pe2 = xyz.index(channel_pe2)
    k_traj_pe2 = -np.round(k_traj_adc[idx_pe2][0::n_adc]*fov_pe2*10)/10  
    k_traj_pe2 = k_traj_pe2 - min(k_traj_pe2)
    acq_pos = [np.array([int(x), int(y)]) for x,y in zip(k_traj_pe1, k_traj_pe2)]
    # if cast to integer not necessary: acq_pos = [list([x, int(y)]) for x,y in zip(k_traj_pe1, k_traj_pe2)]
    
    return acq_pos

