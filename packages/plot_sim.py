# %% load results
import pickle
import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt
from packages.seq_utils import t2_fit, plot_3d
import MRzeroCore as mr0

def plot_sim(plot: dict, seq_filename: str, system:pp.Opts, sim_path: str = "./simulation/"):
    sim_name = 'sim_' + seq_filename

    if plot["seq"]:
        seq = pp.Sequence(system=system)
        seq.read(sim_path + sim_name + '_seq.seq', detect_rf_use = True)
        try:
            signal = np.load(sim_path + sim_name + '_signal.npy')
            sp_adc, t_adc = mr0.util.pulseq_plot(seq=seq,signal=signal)
        except:
            seq.plot()
            plt.pause(0.1)
            # seq.plot(time_range=(0, 2*tr), plot_now=False)
            # import mplcursors
            # mplcursors.cursor()
            # plt.show()
        # plt.plot(np.real(signal))
        # plt.title('real part of RF signal')
        # plt.xlabel('sampling points')
        # plt.ylabel('RF signal')
        # plt.show()

    if plot["kspace"]:        
        with open(sim_path + sim_name + '_seq0.pkl', 'rb') as file:
            seq0 = pickle.load(file)
        seq0.plot_kspace_trajectory()

    if plot["phantom"]:
        with open(sim_path + sim_name + '_obj_p.pkl', 'rb') as file:
            obj_p = pickle.load(file)
        obj_p.plot()
        plt.pause(1)

    if plot["reco"]: 
        with open(sim_path + sim_name + '_reco.pkl', 'rb') as file:
            reco = pickle.load(file)
        plot_3d(reco)
        #plot_3d(np.transpose(np.squeeze(reco[0, :, :, :]), (0, 2, 1)))

    if seq.definitions['name'] == 'tse_3d_mte':
        with open(sim_path + sim_name + '_reco.pkl', 'rb') as file:
            reco = pickle.load(file)
        a_map, T2_map = t2_fit(reco, seq) 
        # plot_3d(a_map)
        plot_3d(T2_map)
    
        
if __name__ == "__main__":
    plot_sim(plot={
        "phantom": True,
        "seq": True,
        "kspace": True,
        "reco": True
        },
        seq_filename='tse_pypulseq')        