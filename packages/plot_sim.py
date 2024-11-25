# %% load results
import pickle
import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

import MRzeroCore as mr0

def plot_sim(plot: dict, seq_filename: str, system:pp.Opts, sim_path: str = "./simulation/"):
    sim_name = 'sim_' + seq_filename
    signal = np.load(sim_path + sim_name + '_signal.npy')

    with open(sim_path + sim_name + '_obj_p.pkl', 'rb') as file:
        obj_p = pickle.load(file)
        
    with open(sim_path + sim_name + '_reco.pkl', 'rb') as file:
        reco = pickle.load(file)
        
    seq = pp.Sequence(system=system)
    seq.read(sim_path + sim_name + '_seq.seq', detect_rf_use = True)

    with open(sim_path + sim_name + '_seq0.pkl', 'rb') as file:
        seq0 = pickle.load(file)
    # %%
    if plot["phantom"]:
        obj_p.plot()
        plt.pause(1)
    if plot["seq"]:
        try:
            sp_adc, t_adc = mr0.util.pulseq_plot(seq=seq,signal=signal)
        except:
            seq.plot()
            plt.pause(1)
        # plt.plot(np.real(signal))
        # plt.title('real part of RF signal')
        # plt.xlabel('sampling points')
        # plt.ylabel('RF signal')
        # plt.show()
        
    if plot["kspace"]:
        seq0.plot_kspace_trajectory()
            
    if plot["reco"]:  
        if reco.shape[2] > 1:
            # Create a figure and axis
            fig, ax = plt.subplots()
            plt.subplots_adjust(left=0.25, bottom=0.25)

            # Initial plot
            slice_index = 0
            img = ax.imshow(reco[:, :, slice_index].T.abs(), cmap='viridis', origin="lower")
            ax.set_title(f'Slice {slice_index + 1}')
            plt.colorbar(img, ax=ax)

            # Create a slider axis and slider
            ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])
            slider = Slider(ax_slider, 'Slice', 0, reco.shape[2] - 1, valinit=slice_index, valstep=1)

            # Update function for the slider
            def update(val):
                slice_index = int(slider.val)
                img.set_data(reco[:, :, slice_index].T.abs())
                ax.set_title(f'Slice {slice_index + 1}')
                fig.canvas.draw_idle()

            # Attach the update function to the slider
            slider.on_changed(update)
            # Create a slider axis and slider for contrast adjustment
            ax_contrast_slider = plt.axes([0.25, 0.15, 0.65, 0.03])
            contrast_slider = Slider(ax_contrast_slider, 'Contrast', 0.1, 100.0, valinit=1.0, valstep=0.1)

            # Update function for the contrast slider
            def update_contrast(val):
                contrast = contrast_slider.val
                img.set_clim(vmin=0, vmax=reco[:, :, slice_index].T.abs().max() * contrast)
                fig.canvas.draw_idle()

            # Attach the update function to the contrast slider
            contrast_slider.on_changed(update_contrast)
            plt.show()

        else:
            # Add adjustable contrast for magnitude and phase images
            fig, (ax1, ax2) = plt.subplots(1, 2)
            plt.subplots_adjust(left=0.25, bottom=0.25)

            # Initial plots
            img1 = ax1.imshow(reco[:, :, 0].abs(), origin="lower")
            ax1.set_title("Magnitude")
            plt.colorbar(img1, ax=ax1)

            img2 = ax2.imshow(reco[:, :, 0].angle(), origin="lower", vmin=-np.pi, vmax=np.pi)
            ax2.set_title("Phase")
            plt.colorbar(img2, ax=ax2)

            # Create a slider axis and slider for contrast adjustment
            ax_contrast_slider = plt.axes([0.25, 0.15, 0.65, 0.03])
            contrast_slider = Slider(ax_contrast_slider, 'Contrast', 0.1, 100.0, valinit=1.0, valstep=0.1)

            # Update function for the contrast slider
            def update_contrast(val):
                contrast = contrast_slider.val
                img1.set_clim(vmin=0, vmax=reco[:, :, 0].abs().max() * contrast)
                img2.set_clim(vmin=-np.pi * contrast, vmax=np.pi * contrast)
                fig.canvas.draw_idle()

            # Attach the update function to the contrast slider
            contrast_slider.on_changed(update_contrast)
            plt.show()         
        
if __name__ == "__main__":
    plot_sim(plot={
        "phantom": True,
        "seq": True,
        "kspace": True,
        "reco": True
        },
        seq_filename='tse_pypulseq')        