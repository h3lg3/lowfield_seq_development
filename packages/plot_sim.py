# %% load results
import pickle
import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt
from packages.seq_utils import plot_3d
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
            plt.pause(0.1)
        # plt.plot(np.real(signal))
        # plt.title('real part of RF signal')
        # plt.xlabel('sampling points')
        # plt.ylabel('RF signal')
        # plt.show()
    if plot["kspace"]:
        seq0.plot_kspace_trajectory()
            
    if plot["reco"]:   
        plot_3d(reco)
        #plot_3d(np.transpose(np.squeeze(reco[0, :, :, :]), (0, 2, 1)))

    
        
if __name__ == "__main__":
    plot_sim(plot={
        "phantom": True,
        "seq": True,
        "kspace": True,
        "reco": True
        },
        seq_filename='tse_pypulseq')        