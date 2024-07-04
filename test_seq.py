# Test seq timing

# %% Import packages
import math
import warnings

import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt

# import utils # several helper functions for simulation and recon
import MRzeroCore as mr0
# %% Define sequence events
system = pp.Opts(
    max_grad=32,
    grad_unit="mT/m",
    max_slew=130,
    slew_unit="T/m/s",
    rf_ringdown_time=100e-6,
    rf_dead_time=100e-6,
    adc_dead_time=10e-6,
)

seq = pp.Sequence(system)  # Create a new sequence object
fov = 256e-3  # Define FOV and resolution
Nx, Ny = 64, 64
n_slices = 1
slice_thickness = 5e-3
dwell_time = 50e-6
TE = 20e-3  # Echo time
TR = 100e-3  # Repetition time

rf_ex, gz, _ = pp.make_sinc_pulse(
    flip_angle=90*np.pi/180,
    system=system,
    duration=1e-3,
    slice_thickness=slice_thickness,
    apodization=0.5,
    time_bw_product=4,
    phase_offset=0,
    return_gz=True,
)

rf_ref, gz, _ = pp.make_sinc_pulse(
    flip_angle=180*np.pi/180,
    system=system,
    duration=1e-3,
    slice_thickness=slice_thickness,
    apodization=0.5,
    time_bw_product=4,
    phase_offset=0,
    use="refocusing",
    return_gz=True,
)

adc = pp.make_adc(
    num_samples=Nx, duration=Nx*dwell_time, phase_offset=180*np.pi/180, delay=system.adc_dead_time, system=system
)

TE1_delay = TE/2-rf_ex.shape_dur/2-rf_ex.ringdown_time-rf_ref.shape_dur/2-rf_ref.delay
TE2_delay = TE/2-rf_ref.shape_dur/2-rf_ref.ringdown_time-adc.delay-(adc.num_samples*adc.dwell)/2

#@title Calculate delays
# ct_ex=pp.calc_rf_center(rf_ex)    # rf center time returns time and index of the center of the pulse
# ct_ref=pp.calc_rf_center(rf_ex)    # rf center time returns time and index of the center of the pulse
# minTE2 = ct_ex[0]+rf_ex.ringdown_time+rf_ref.delay+ct_ref[0]   #  center pulse to center pulse
# d1=TE/2-minTE2
#       # TE/2   - half ref pulse                 - half adc      -adc.delay
# d2 = minTE2+d1 - (ct_ref[0]+rf_ex.ringdown_time)  - Nx*dwell_time/2 -adc.delay

# % CONSTRUCT SEQUENCE
seq.add_block(rf_ex, gz) # , gz
seq.add_block(pp.make_delay(TE1_delay))
seq.add_block(rf_ref, gz) #, gz
seq.add_block(pp.make_delay(TE2_delay))
seq.add_block(adc)

## Check whether the timing of the sequence is compatible with the scanner
(
    ok,
    error_report,
) = seq.check_timing()  # Check whether the timing of the sequence is correct
if ok:
    print("Timing check passed successfully")
else:
    print("Timing check failed. Error listing follows:")
    [print(e) for e in error_report]

# Prepare the sequence output for the scanner
seq.set_definition('Name', 'se')
seq.write('test_seq.seq')

# % SIMULATION
# title 1D SE in a pixel phantom - simulation
dB0 = 0
# S4: SETUP SPIN SYSTEM/object on which we can run the MR sequence external.seq from above
# set phantom  manually to a pixel phantom. Coordinate system is [-0.5, 0.5]^3
obj_p = mr0.CustomVoxelPhantom(
        pos=[[0., 0., 0]],
        PD=[1.0],
        T1=[3.0],
        T2=[0.5],
        T2dash=[30e-3],
        D=[0.0],
        B0=0,
        voxel_size=0.1,
        voxel_shape="box"
    )
# Manipulate loaded data
obj_p.B0+=dB0
obj_p.D*=0
obj_p.plot()
# Convert Phantom into simulation data
obj_p=obj_p.build()

# SIMULATE  the external.seq file and add acquired signal to ADC plot
# Read in the sequence
seq0 = mr0.Sequence.from_seq_file("test_seq.seq")
#seq0.plot_kspace_trajectory()
# Simulate the sequence
graph=mr0.compute_graph(seq0, obj_p, 200, 1e-3)
signal=mr0.execute_graph(graph, seq0, obj_p)
# PLOT sequence with signal in the ADC subplot
seq.plot(plot_now=False)
mr0.util.insert_signal_plot(seq=seq, signal =signal.numpy())
plt.show()

# %%
