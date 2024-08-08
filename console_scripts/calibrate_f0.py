"""Spin-echo spectrum."""
# %%
import logging

import matplotlib.pyplot as plt
import numpy as np

import console.spcm_control.globals as glob
from console.interfaces.interface_acquisition_data import AcquisitionData
from console.spcm_control.acquisition_control import AcquisitionControl
from console.utilities.sequences.spectrometry import se_spectrum
from console.utilities.snr import signal_to_noise_ratio


# import pypulseq as pp
# %%
# Create acquisition control instance
acq = AcquisitionControl(configuration_file="device_config_ptb.yaml", console_log_level=logging.INFO, file_log_level=logging.DEBUG)

# %%
# Set SE parameter

params = dict(
    echo_time=12e-3,
    # rf_duration=400e-6,
    rf_duration=200e-6,
    time_bw_product=0,
    use_sinc=False,
    
    adc_duration=50e-3,
    use_fid=True
    
    # adc_duration=6e-3,
    # use_fid=False
)
seq = se_spectrum.constructor(**params)

# Adjust b1-scaling
# glob.update_parameters(b1_scaling=4.7) # ACR phantom
# glob.update_parameters(b1_scaling=2.1) # small sphere


# %%
#acquire data
current_f0 = glob.parameter.larmor_frequency

acq.set_sequence(sequence=seq)
acq_data: AcquisitionData = acq.run()

# FFT
data = np.mean(acq_data.raw, axis=0)[0].squeeze()
data_fft = np.fft.fftshift(np.fft.fft(np.fft.fftshift(data)))
fft_freq = np.fft.fftshift(np.fft.fftfreq(data.size, acq_data.dwell_time))

max_spec = np.max(np.abs(data_fft))
f_0_offset = fft_freq[np.argmax(np.abs(data_fft))]

#update global larmor frequency to measured f0
if abs(f_0_offset) <= 6e3:
    glob.update_parameters(larmor_frequency=current_f0-f_0_offset)

snr = signal_to_noise_ratio(data_fft, dwell_time=acq_data.dwell_time)


# Add information to acquisition data
acq_data.add_info({
    "true f0": current_f0 - f_0_offset,
    "magnitude spectrum max": max_spec,
    "snr dB": snr,
})
# acq_data.add_info(params)

print(f"Frequency offset [Hz]: {f_0_offset}\nNew frequency f0 [Hz]: {current_f0 - f_0_offset}")
print(f"Frequency spectrum max.: {max_spec}")
print("Acquisition data shape: ", acq_data.raw.shape)
print("SNR [dB]: ", snr)

# Plot spectrum
time_axis = np.arange(data.size)*acq_data.dwell_time*1e3
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(time_axis, np.abs(data), label="Abs")
ax[0].plot(time_axis, np.real(data), label="Re")
ax[0].plot(time_axis, np.imag(data), label="Im")
ax[0].set_xlabel("Time [ms]")
ax[0].set_ylabel("RX signal [mV]")
# ax[0].set_xlim((0, np.min([np.max(time_axis), 6])))
ax[0].legend(loc="upper right")
ax[1].plot(fft_freq, np.abs(data_fft))
# ax[1].set_xlim([-10e3, 10e3])
ax[1].set_ylim([0, max_spec * 1.05])
ax[1].set_ylabel("Abs. FFT Spectrum [a.u.]")
ax[1].set_xlabel("Frequency [Hz]")
plt.show()

# %%
# Save acquisition data

acq_data.add_info({
    "note": ""
})

acq_data.save(save_unprocessed=True)
# %%
