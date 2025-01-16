# Reconstruct Siemens data
import pypulseq as pp
import numpy as np
import glob
import nibabel as nib

import packages.utils as utils  # from MR-Physics-with-Pulseq\Tutorials\utils
from packages.mr_systems import lumina as system
from packages.seq_utils import plot_3d

data_path = r".\siemens_data\241218"
seq_name = "tse_3d_mte_lumina_64_64_32_TR300"

seq_file = f"{data_path}\\{seq_name}.seq"
data_pattern = f"{data_path}\\meas_*_{seq_name}.dat"

# Use glob to find the actual data file matching the pattern
data_files = glob.glob(data_pattern)
if not data_files:
    raise FileNotFoundError(f"No files matching pattern {data_pattern}")
data_file = data_files[0]  # Use the first matching file

# import sequence and data
seq = pp.Sequence(system=system)
seq.read(seq_file, detect_rf_use = True)
 
kdata_unsorted = utils.read_raw_data(data_file)
kdata_sorted = utils.sort_data_labels(kdata_unsorted, seq)  # order: 'COIL', 'LIN', 'PAR', 'ADC' or order: 'COIL', 'REP', 'LIN', 'PAR', 'ADC'

# Reconstruct image
image = utils.ifft_3d(kdata_sorted)
image_sos = np.sqrt((abs(image)**2).sum(axis=0))
if image_sos.ndim == 4:
    image_sos =  image_sos.transpose(3, 1, 2, 0)                # order: 'ADC', 'LIN', 'PAR', 'REP'
    plot_3d(image_sos)
else:
    image_sos = image_sos.transpose(2, 0, 1)                # order: 'ADC', 'LIN', 'PAR'
    # Save image_sos to data_path with the name seq_name + "image.dat"
    output_file = f"{data_path}\\{seq_name}_image.dat"
    np.save(output_file, image_sos)

    # Plot image
    # plot_3d(image_sos)

    # Save image_sos as a nifti file
    dx = seq.definitions['FOV'][0]/seq.definitions['number_of_readouts']*1000
    dy = seq.definitions['FOV'][1]/seq.definitions['k_space_encoding1']*1000
    dz = seq.definitions['FOV'][2]/seq.definitions['k_space_encoding2']*1000

    voxel_size = [dx, dy, dz]
    nifti_img = nib.Nifti1Image(image_sos, np.diag(voxel_size + [1]))
    nifti_output_file = f"{data_path}\\{seq_name}_image.nii"

    # # Add voxel size information
    nifti_img.header.set_zooms(tuple(voxel_size))

    # Add field-of-view information
    nifti_img.header['dim'][1:4] = tuple(seq.definitions['FOV']*1000)

    nib.save(nifti_img, nifti_output_file)