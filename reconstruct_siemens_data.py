# Reconstruct Siemens data
# %%
import packages.utils as utils  # from MR-Physics-with-Pulseq\Tutorials\utils
import pypulseq as pp
import matplotlib.pyplot as plt
# %%
data_file = r".\siemens_data\meas_MID00108_FID109501_tse_pulseq.dat"
seq_file =  r".\siemens_data\tse_pypulseq.seq"

seq = pp.Sequence()
seq.read(seq_file)
 
kdata_unsorted = utils.read_raw_data(data_file)
# kdata_sorted = utils.sort_data_implicit(kdata = kdata_unsorted, seq=seq)
# %% Reconstruct
sos1 = utils.recon_cartesian_2d(kdata = kdata_unsorted, 
                         seq = seq,
                         use_labels = False)
utils.plot_nd(rec = sos1)
plt.show()