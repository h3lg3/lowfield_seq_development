# Reconstruct Siemens data
# %%
import packages.utils as utils  # from MR-Physics-with-Pulseq\Tutorials\utils
import pypulseq as pp
import matplotlib.pyplot as plt
import glob
# %%
data_path = r".\siemens_data"
seq_name = "tse_3d_lumina"

seq_file = f"{data_path}\\{seq_name}.seq"
data_pattern = f"{data_path}\\meas_*_{seq_name}.dat"

# Use glob to find the actual data file matching the pattern
data_files = glob.glob(data_pattern)
if not data_files:
    raise FileNotFoundError(f"No files matching pattern {data_pattern}")
data_file = data_files[0]  # Use the first matching file

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