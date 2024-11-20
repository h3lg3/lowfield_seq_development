# %% Import packages
import numpy as np
import matplotlib.pyplot as plt
import pypulseq as pp
import MRzeroCore as mr0
import pickle
import os 
import torch
import nibabel as nib

from packages.mr_systems import lumina as system



# choose flags
FLAG_PLOTS = True
FLAG_SAVE = False

# define sequence filename
seq_filename = "tse_3d_lumina"
seq_path = "./sequences/"

# define qMR data 
qMR_path = "E:/MATLAB/Projekte/Phantome/Results/Phantom-Aga_20241113-174629/"
qMR_mat = "qMR_data.mat"
qMR_t2star = "T2star_map.nii"

# define output filename
sim_filename = "agar_tse_3d_lumina"
sim_path: str = "./simulation/"

# define simulation parameters
fov = (250e-3, 250e-3, 250e-3)
n_enc = (64, 64, 1)


# %% Create phantom, simulate sequence, reconstruct image

seq_file = seq_path + seq_filename + ".seq"
sim_name = "sim_" + seq_filename

seq = pp.Sequence(system=system)
seq.read(seq_file, detect_rf_use = True)
seq0 = mr0.Sequence.import_file(seq_file)

# load qMR data from mat file
obj_p = mr0.VoxelGridPhantom.load_mat(qMR_path + qMR_mat) 

# load the T2star_map from nifti file
nifti_file = qMR_path + qMR_t2star
nifti_img = nib.load(nifti_file)
nifti_data = nifti_img.get_fdata()
nifti_tensor = torch.tensor(nifti_data, dtype=torch.float32)
obj_p.T2dash = nifti_tensor

# set diffusion to 0
obj_p.D[:] = 0

if n_enc[2] == 1:
    center_slice = obj_p.PD.shape[2]//2
    obj_p = obj_p.slices([center_slice])    # select center slice
    obj_p = obj_p.interpolate(n_enc[0], n_enc[1], 1)  # interpolate     
else:
    range_slices = tuple(np.linspace(40, 90, n_enc[2], dtype=int))
    obj_p = obj_p.slices(range_slices)      # select slices within the range
    obj_p = obj_p.interpolate(n_enc[0], n_enc[1], n_enc[2])      # interpolate

obj_sim = obj_p.build(PD_threshold=0.005)

# SIMULATE the external.seq file and add acquired signal to ADC plot
graph=mr0.compute_graph(seq0, obj_sim, 200, 1e-3)
signal=mr0.execute_graph(graph, seq0, obj_sim)
reco = mr0.reco_adjoint(signal, seq0.get_kspace(), resolution=n_enc, FOV=fov) # Recommended: RECO has same Reso and FOV as sequence
if FLAG_PLOTS:
    # plt.imshow(reco.abs(), cmap='viridis', origin="lower")
    # plt.show()  
    if reco.shape[2] > 1:
        # Plot reco.T.abs() along the third dimension in one figure
        fig, axes = plt.subplots(1, reco.shape[2], figsize=(15, 5))
        
        for i in range(reco.shape[2]):
            axes[i].imshow(reco[:, :, i].abs(), cmap='viridis', origin="lower")
            axes[i].set_title(f'Slice {i+1}')
            axes[i].axis('off')  # Hide axes
        
        # Adjust layout to prevent overlap
        plt.tight_layout()
        plt.show()
    else:
        plt.figure()
        plt.subplot(121)
        plt.title("Magnitude")
        plt.imshow(reco[:, :, 0].abs(), origin="lower")
        plt.colorbar()
        plt.subplot(122)
        plt.title("Phase")
        plt.imshow(reco[:, :, 0].angle(), origin="lower", vmin=-np.pi, vmax=np.pi)
        plt.colorbar()
        plt.show()

# %% save results
if FLAG_SAVE:
    # Check if directory exists
    if not os.path.exists(sim_path):
        # Create the directory
        os.makedirs(sim_path)
        print(f"Directory '{sim_path}' created.")

    with open(sim_path+ sim_name + '_obj_p.pkl', 'wb') as file:
        pickle.dump(obj_p, file)

    np.save(sim_path + sim_name + '_signal.npy', signal)

    with open(sim_path + sim_name + '_reco.pkl', 'wb') as file:
        pickle.dump(reco, file)
        
    seq_file = sim_name + '_seq.seq'
    seq.write(sim_path + seq_file)

    with open(sim_path + sim_name + '_seq0.pkl', 'wb') as file:
        pickle.dump(seq0, file)
