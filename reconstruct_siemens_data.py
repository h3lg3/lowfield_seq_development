# Reconstruct Siemens data
# %%
import packages.utils as utils  # from MR-Physics-with-Pulseq\Tutorials\utils

import pypulseq as pp

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import glob
import numpy as np
# %%
data_path = r".\siemens_data\241029"
seq_name = "tse_3d_lumina_alternate180"

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
kdata_sorted = utils.sort_data_labels(kdata_unsorted, seq)
kdata_sorted = kdata_sorted.transpose(0,1,3,2)

# %%
axes = (-3, -2, -1)
#fft_data_unsorted = np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(kdata_unsorted, axes=axes), axes=axes), axes=axes)

image = np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(kdata_sorted, axes=axes), axes=axes), axes=axes)
image_sos = np.sqrt((abs(image)**2).sum(axis=0))
# Save image_sos to data_path with the name seq_name + "image.dat"
output_file = f"{data_path}\\{seq_name}_image.dat"
np.save(output_file, image_sos)


# # # # Plot the first slice of fft_data_unsorted
# plt.figure(figsize=(10, 8))
# plt.imshow(image_sos[:, 4, :], cmap='gray')
# plt.title('FFT Data Unsorted - Slice 0')
# plt.axis('off')
# plt.show()
# plt.pause(0.1)
# %%
data2plot = image_sos

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 8))
plt.subplots_adjust(left=0.25, bottom=0.25)

# Initial contrast limits
vmin = 0
vmax = abs(data2plot).max()

# Plot the initial image
im = ax.imshow(abs(data2plot[:, :, 0]), cmap='gray', vmin=vmin, vmax=vmax)
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
s_slice = Slider(ax_slice, 'Slice', 0, data2plot.shape[2] - 1, valinit=0, valstep=1)

# Update function for slice selection
def update_slice(val):
    slice_idx = int(s_slice.val)
    im.set_data(abs(data2plot[:, :, slice_idx]))
    ax.set_title(f'Slice {slice_idx}')
    fig.canvas.draw_idle()

s_slice.on_changed(update_slice)

plt.show()