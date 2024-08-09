"""3D turbo spin echo sequence."""
# %%
# imports
import sys
sys.path.append("../")
import logging

import matplotlib.pyplot as plt
import numpy as np

import console.spcm_control.globals as glob
# import console.utilities.sequences as sequences
from console.interfaces.interface_acquisition_data import AcquisitionData
from console.interfaces.interface_acquisition_parameter import Dimensions
from console.spcm_control.acquisition_control import AcquisitionControl
from packages.tse_trajectory import Trajectory
from packages import tse_3d

# %%
# Create acquisition control instance
configuration = "device_config_ptb.yaml"
acq = AcquisitionControl(configuration_file=configuration, console_log_level=logging.INFO, file_log_level=logging.DEBUG)


# %%
# Create sequence
params = dict(
    # echo_time=14e-3,
    # repetition_time=600e-3,
    # noise_scan=False,
    # etl=18,
    # # dummies=3,
    # trajectory=sequences.tse_3d.Trajectory.INOUT,
    # gradient_correction=80e-6,
    # rf_duration=100e-6,
    # fov=Dimensions(x=140e-3, y=140e-3, z=140e-3),
    # channel_ro="z",
    # channel_pe1="y",
    # channel_pe2="x",
    # ro_bandwidth=20e3,
    # n_enc=Dimensions(x=1, y=120, z=120), # 2D, y = RO
    echo_time=28e-3,
    repetition_time=2000e-3, 
    etl=8, # define max sampling period (tmax = 200ms?), etl_max = round(tmax/esp), nr. of pe1 steps should be multiple of etl
    # dummies=2, 
    ro_bandwidth=20e3, 
    fov=Dimensions(x=140e-3, y=140e-3, z=140e-3), 
    n_enc=Dimensions(x=1, y=120, z=120),
    trajectory=Trajectory.LINEAR,
    excitation_phase=np.pi/2,
    refocussing_phase=0,
    channel_ro="z",
    channel_pe1="y",
    channel_pe2="x"
)
# seq, _ = sequences.tse.tse_3d.constructor(**params)
seq, traj, dim = tse_3d.constructor(**params)

seq.plot()
                

# Calculate decimation:
decimation = int(acq.rx_card.sample_rate * 1e6 / params["ro_bandwidth"])
glob.update_parameters(decimation=decimation)

# %%
# Unroll equence
acq.set_sequence(sequence=seq)

# %%
# Run sequence
acq_data: AcquisitionData = acq.run()

# Sort k-space
# ksp = sequences.tse.tse_3d.sort_kspace(acq_data.raw, seq=seq)
ksp = tse_3d.sort_kspace(acq_data.raw, trajectory=traj, kdims=dim)
ksp = np.mean(ksp, axis=0).squeeze()
img = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(ksp)))


 # %%
# Plot 2D k-space and image
fig, ax = plt.subplots(1, 2, figsize=(8, 4))
ax[0].imshow(np.abs(ksp), cmap="gray")
ax[1].imshow(np.abs(img), cmap="gray")
plt.show()

# %%
# Plot 3D images
num_slices = img.shape[0]
num_cols = int(np.ceil(np.sqrt(num_slices)))
num_rows = int(np.ceil(num_slices/num_cols))
fig, ax = plt.subplots(num_rows, num_cols, figsize=(10, 10))
ax = ax.ravel()
total_max = np.amax(np.abs(img))
total_min = 0   # np.amin(np.abs(img))
for k, x in enumerate(img[:, ...]):
    ax[k].imshow(np.abs(x), vmin=total_min, vmax=total_max, cmap="gray")
    ax[k].axis("off")
_ = [a.remove() for a in ax[k+1:]]
fig.set_tight_layout(True)
fig.set_facecolor("black")

# %%
# Save data
acq_data.add_info(dict(
    subject="",
    note=""
))

acq_data.add_info(params)
acq_data.add_data({
    "kspace": ksp,
    "image": img
})

acq_data.save(
    # overwrite=True
    # save_unprocessed=True, 
    # user_path="/home/guest/charite/data/"
)

# %%
del acq

# %%
