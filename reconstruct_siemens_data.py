# Reconstruct Siemens data
# %%
import packages.utils as utils  # from MR-Physics-with-Pulseq\Tutorials\utils
import pypulseq as pp
import matplotlib.pyplot as plt
import glob
from matplotlib.widgets import Slider
import numpy as np
# %%
data_path = r".\siemens_data\241021"
seq_name = "tse_3d_lumina_no_spoiler"

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
kdata_cropped = kdata_unsorted[:, 10:54, 10:54] 
# %%
axes = (-2, -1)
#fft_data_unsorted = np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(kdata_unsorted, axes=axes), axes=axes), axes=axes)

fft_data_unsorted = np.fft.fftshift(np.fft.ifftn(kdata_cropped, axes=axes), axes=axes)
image = np.sqrt((abs(fft_data_unsorted)**2).sum(axis=0))

# # Plot the first slice of fft_data_unsorted
plt.figure(figsize=(10, 8))
plt.imshow(image, cmap='gray')
plt.title('FFT Data Unsorted - Slice 0')
plt.axis('off')
plt.show()
# %%
data2plot = fft_data_unsorted

# Plot abs(kdata_unsorted) as 2D image for kdata_unsorted[0:17, :, :]

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 8))
plt.subplots_adjust(left=0.25, bottom=0.25)

# Initial contrast limits
vmin = 0
vmax = abs(data2plot).max()

# Plot the initial image
im = ax.imshow(abs(data2plot[0, :, :]), cmap='gray', vmin=vmin, vmax=vmax)
ax.set_title('Slice 0')
ax.axis('off')

# Create sliders for contrast adjustment
axcolor = 'lightgoldenrodyellow'
ax_vmin = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_vmax = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

s_vmin = Slider(ax_vmin, 'Min', 0, abs(data2plot).max(), valinit=vmin)
s_vmax = Slider(ax_vmax, 'Max', 0, abs(data2plot).max(), valinit=vmax)

# Update function for sliders
def update(val):
    im.set_clim([s_vmin.val, s_vmax.val])
    fig.canvas.draw_idle()

s_vmin.on_changed(update)
s_vmax.on_changed(update)

# Create a slider for slice selection
ax_slice = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
s_slice = Slider(ax_slice, 'Slice', 0, data2plot.shape[0] - 1, valinit=0, valstep=1)

# Update function for slice selection
def update_slice(val):
    slice_idx = int(s_slice.val)
    im.set_data(abs(data2plot[slice_idx, :, :]))
    ax.set_title(f'Slice {slice_idx}')
    fig.canvas.draw_idle()

s_slice.on_changed(update_slice)

plt.show()


# plt.figure(figsize=(10, 8))
# for i in range(18):
#     plt.subplot(3, 6, i + 1)
#     plt.imshow(abs(kdata_unsorted[i, :, :]), cmap='gray')
#     plt.title(f'Slice {i}')
#     plt.axis('off')
# plt.tight_layout()
# plt.show()


#kdata_sorted = utils.sort_data_implicit(kdata = kdata_unsorted, seq=seq)
# %% Reconstruct
# sos1 = utils.recon_cartesian_2d(kdata = kdata_unsorted, 
#                          seq = seq,
#                          use_labels = False)
# utils.plot_nd(rec = sos1) # , vmin=0, vmax=1e-5
# plt.show()
# %%
