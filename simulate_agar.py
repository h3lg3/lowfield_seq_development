# %% Import packages
import numpy as np
import matplotlib.pyplot as plt
import pypulseq as pp
import MRzeroCore as mr0
import pickle
import os 

import torch
from torch.nn.functional import normalize

import nibabel as nib

from packages.seq_utils import plot_3d
from packages.mr_systems import lumina as system
from packages.write_seq_definitions import custom_seq_definitons

# choose flags
FLAG_PLOTS = True
FLAG_SAVE = True

# define sequence filename
seq_filename = "tse_3d_lumina_32_32_16_160mm"
seq_path = "./sequences/"

# %%# define qMR data 
qMR_path = "E:/MATLAB/Projekte/Phantome/Results/Phantom_Agar_coronal_160mm/"
qMR_mat = "qMR_data.mat"
qMR_T2dash = "T2strich_map.nii"

# define output filename
sim_filename = "agar_tse_3d_lumina_160mm"
sim_path: str = "./simulation/"

# for multi-slice simulation
n_slices = tuple(range(12, 29, 1))
# %% Create phantom, simulate sequence, reconstruct image
seq_file = seq_path + seq_filename + ".seq"
sim_name = "sim_" + seq_filename

seq = pp.Sequence(system=system)
seq.read(seq_file, detect_rf_use = True)

n_enc = (seq.definitions['number_of_readouts'], seq.definitions['k_space_encoding1'], seq.definitions['k_space_encoding2'])
n_enc = tuple(map(int, n_enc))

fov = tuple(seq.definitions['FOV'])

# Remove definitions from the sequence because they cause import error in mr0.Sequence.import_file
temp_seq_file = seq_path + seq_filename + '_temp.seq'
for definition in custom_seq_definitons:
    seq.definitions.pop(definition, None)   
seq.write(temp_seq_file)
# import the sequence with MRzero
seq0 = mr0.Sequence.import_file(temp_seq_file)
# Delete the temporary sequence file
os.remove(temp_seq_file)

# load qMR data from mat file
obj_p = mr0.VoxelGridPhantom.load_mat(qMR_path + qMR_mat, size = [0.16, 0.16, 0.16]) 

# # # load the T2dash_map from nifti file
nifti_file = qMR_path + qMR_T2dash
nifti_img = nib.load(nifti_file)
nifti_data = nifti_img.get_fdata()
nifti_tensor = torch.tensor(nifti_data, dtype=torch.float32)
obj_p.T2dash = nifti_tensor

# set diffusion to 0
obj_p.D[:] = 0

if n_enc[2] == 1:
    center_slice = obj_p.PD.shape[2]//2
    obj_p = obj_p.slices([4])    # select center slice
    obj_p = obj_p.interpolate(n_enc[0], n_enc[1], 1)  # interpolate     
else:
    obj_p = obj_p.slices(n_slices)                              # select slices within the range
    obj_p = obj_p.interpolate(n_enc[0], n_enc[1], n_enc[2])     # interpolate

# plot_3d(obj_p.PD)
# plot_3d(obj_p.T1)
# plot_3d(obj_p.T2)
# plot_3d(obj_p.T2dash)

obj_sim = obj_p.build(PD_threshold=0.005)

# SIMULATE the external.seq file and add acquired signal to ADC plot
graph=mr0.compute_graph(seq0, obj_sim, 200, 1e-3)
signal=mr0.execute_graph(graph, seq0, obj_sim)
reco = mr0.reco_adjoint(signal, seq0.get_kspace(), resolution=n_enc, FOV=fov) # Recommended: RECO has same Reso and FOV as sequence

# Normalize reco by its max value
reco = normalize(reco)

if FLAG_PLOTS:
    plot_3d(reco)

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
