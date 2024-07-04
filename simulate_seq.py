# %% Import packages
#import math
#import warnings

import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt

#import utils # several helper functions for simulation and recon
import MRzeroCore as mr0

# %% VISUALIZATION
seq_file = "./Examples/epi_se_pypulseq.seq" 
seq = pp.Sequence()
seq.read(seq_file)
seq0 = mr0.Sequence.import_file(seq_file)

# % SIMULATION, SETUP SPIN SYSTEM/object on which we can run the MR sequence external.seq from above
# title 1D SE in a pixel phantom - simulation
plot_phantom = False 
pixel_phantom = True 
sz = [64, 64]

print('load phantom')
if pixel_phantom:
    obj_p = mr0.CustomVoxelPhantom(
        pos=[[0., 0., 0]], # [0., 0., 0], [-0.25, -0.25, 0]
        PD=[1.0],
        T1=[3.0],
        T2=[0.5],
        T2dash=[30e-3],
        D=[0.0],
        B0=0,
        voxel_size=0.1,
        voxel_shape="box"
    )
    
else:
    obj_p = mr0.VoxelGridPhantom.load_mat('numerical_brain_cropped.mat')
    brain_phantom_res = 64 #@param {type:"slider", min:16, max:128, step:16}
    obj_p = obj_p.interpolate(brain_phantom_res, brain_phantom_res, 1)
    obj_p.B0[:] = 0
    # obj_p.D[:] = 0

if plot_phantom:
    obj_p.plot()
    
obj_p = obj_p.build()

# SIMULATE  the external.seq file and add acquired signal to ADC plot
# Simulate the sequence
graph=mr0.compute_graph(seq0, obj_p, 200, 1e-3)
signal=mr0.execute_graph(graph, seq0, obj_p)
# %% Plot sequence and signal
sp_adc, t_adc = mr0.util.pulseq_plot(seq=seq,signal=signal.numpy())

# Unfortunately, we need to limit the resolution as reco_adjoint is very RAM-hungy
print('reconstruct and plot')
seq0.plot_kspace_trajectory()

reco = mr0.reco_adjoint(signal, seq0.get_kspace(), resolution=(64, 64, 1), FOV=(0.22, 0.22, 1))
plt.figure()
plt.subplot(121)
plt.title("Magnitude")
plt.imshow(reco[:, :, 0].T.abs(), origin="lower")
plt.colorbar()
plt.subplot(122)
plt.title("Phase")
plt.imshow(reco[:, :, 0].T.angle(), origin="lower", vmin=-np.pi, vmax=np.pi)
plt.colorbar()
plt.show()
# # %% Numerical brain phantom
# sz = [64, 64]
# obj_p = mr0.VoxelGridPhantom.load_mat('numerical_brain_cropped.mat')
# brain_phantom_res = 64 #@param {type:"slider", min:16, max:128, step:16}
# obj_p = obj_p.interpolate(brain_phantom_res, brain_phantom_res, 1)
# obj_p.B0[:] = 0
# plot_phantom = True #@param {type:"boolean"}
# if plot_phantom: obj_p.plot()

# obj_p = obj_p.build()

# %%
