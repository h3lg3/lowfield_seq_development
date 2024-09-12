# %% load results
import pickle
import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt
import MRzeroCore as mr0
import torch

def plot_sim(plot: dict, seq_filename: str, sim_path: str = "./simulation/"):
    sim_name = 'sim_' + seq_filename
    signal = np.load(sim_path + sim_name + '_signal.npy')

    with open(sim_path + sim_name + '_obj_p.pkl', 'rb') as file:
        obj_p = pickle.load(file)
        
    with open(sim_path + sim_name + '_reco.pkl', 'rb') as file:
        reco = pickle.load(file)
        
    seq = pp.Sequence()
    seq.read(sim_path + sim_name + '_seq.seq')

    with open(sim_path + sim_name + '_seq0.pkl', 'rb') as file:
        seq0 = pickle.load(file)
    # %%
    if plot["phantom"]:
        obj_p.plot()
        
    if plot["seq"]:
        try:
            sp_adc, t_adc = mr0.util.pulseq_plot(seq=seq,signal=signal)
        except:
            seq.plot()
            
        # plt.plot(np.real(signal))
        # plt.title('real part of RF signal')
        # plt.xlabel('sampling points')
        # plt.ylabel('RF signal')
        # plt.show()
        
    if plot["kspace"]:
        seq0.plot_kspace_trajectory()
            
    if plot["reco"]:
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
    
if __name__ == "__main__":
    plot_sim(plot={
        "phantom": True,
        "seq": True,
        "kspace": True,
        "reco": True
        },
        seq_filename='tse_pypulseq')        
# %% 3D plots for reco

# # Determine the shortest dimension
#     shortest_dim = torch.argmin(torch.tensor(reco.shape)).item()

#     # Slice the tensor along the shortest dimension and plot each 2D slice
#     num_slices = reco.shape[shortest_dim]

#     # Prepare subplots
#     fig, axes = plt.subplots(2, num_slices, figsize=(15, 5))
#     # Handle case where there's only one subplot (axes is not an array)
#     if num_slices == 1:
#         axes = [axes]
        
#     for i in range(num_slices):
#         if shortest_dim == 0:
#             slice_2d = reco[i, :, :].T.abs()
#         elif shortest_dim == 1:
#             slice_2d = reco[:, i, :].T.abs()
#         else:
#             slice_2d = reco[:, :, i].T.abs()
            
#         # Plot the 2D slice
#         im1 = axes[i][0].imshow(slice_2d, cmap='viridis', origin="lower")
#         axes[i][0].set_title(f'Magnitude Slice {i+1}')
#         axes[i][0].axis('off')  # Hide axes
#         im2 = axes[i][1].imshow(slice_2d, cmap='viridis', origin="lower")
#         axes[i][1].set_title(f'Magnitude Slice {i+1}')
#         axes[i][1].axis('off')  # Hide axes
        
#     fig.colorbar(im1, ax=axes[i][0])
#     fig.colorbar(im2, ax=axes[i][1])

#     # Adjust layout to prevent overlap
#     plt.tight_layout()
#     plt.show()  

