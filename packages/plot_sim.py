# %% load results
import pickle
import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt
import MRzeroCore as mr0

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
            
        plt.plot(np.real(signal))
        plt.title('real part of RF signal')
        plt.xlabel('sampling points')
        plt.ylabel('RF signal')
        plt.show()
        
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
# %%
