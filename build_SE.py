# %% Import packages
#import math
#import warnings

import numpy as np
import pypulseq as pp
from matplotlib import pyplot as plt
import math

#import packages.utils as utils# several helper functions for simulation and recon
import MRzeroCore as mr0
# %% Define sequence events
######################## SEQUENCE SETTINGS ########################
system = pp.Opts(
    max_grad=32,
    grad_unit="mT/m",
    max_slew=130,
    slew_unit="T/m/s",
    rf_ringdown_time=100e-6,
    rf_dead_time=100e-6,
    adc_dead_time=10e-6,
    adc_raster_time=2e-6,
    grad_raster_time=20e-6,
    B0=3
)
seq = pp.Sequence(system)  # Create a new sequence object

fov = 256e-3  # Define FOV and resolution
delta_k = 1/fov

n_read, n_phase = 32, 32
n_slices = 1
slice_thickness = 10e-3

dwell_time = 100e-6
ro_bandwidth = 1/dwell_time
ro_duration = np.ceil(n_read*dwell_time/system.grad_raster_time)*system.grad_raster_time
rf_duration = 1e-3

TE = 30e-3  # Echo time
TR = 10  # Repetition time

######################## RF PULSES ########################
rf_ex, gz_ex, gzr_ex = pp.make_sinc_pulse(
    flip_angle=90*np.pi/180,
    system=system,
    duration=rf_duration,
    slice_thickness=slice_thickness,
    apodization=0.5,
    time_bw_product=4,
    phase_offset=0,
    return_gz=True,
)

## Selective or non-selective rf pulse?
rf_ref, gz_ref, _ = pp.make_sinc_pulse(
    flip_angle=180*np.pi/180,
    system=system,
    duration=rf_duration,
    slice_thickness=slice_thickness,
    apodization=0.5,
    time_bw_product=4,
    phase_offset=0,
    use="refocusing",
    return_gz=True,
)

# rf_ref = pp.make_block_pulse(
#     flip_angle=np.pi, system=system, duration=500e-6, use="refocusing"
# )
# gz_spoil = pp.make_trapezoid(
#     channel="z", system=system, area=gz_ex.area * 2
# )

######################## ADC ######################## 
adc = pp.make_adc(
    num_samples=n_read, duration=ro_duration, phase_offset=180*np.pi/180, delay=gx.rise_time, system=system
    )

######################## X GRADIENTS ########################    
## Calculate readout gradient, define area or flat_area?
gx = pp.make_trapezoid(
    channel='x', flat_area=delta_k*n_read, flat_time=ro_duration, system=system
    )
## Account for adc samples taken at dwell time
gx_pre = pp.make_trapezoid(
    channel='x', area=gx.area / 2 - delta_k / 2, duration=ro_duration/2, system=system
    )
#gx_pre = pp.make_trapezoid(channel='x', flat_area=delta_k*n_read/2, flat_time=ro_duration/2, system=system)

## Spoiler gradient on x (used three times: before excitation, before refocusing, after refocusing) 
gx_sp = pp.make_trapezoid(
    channel='x', area=2*gx.area, duration=ro_duration, system=system
    )

######################## Y GRADIENTS ######################## 
## Calculate phase encode gradient, define area or flat_area?
max_gy = pp.make_trapezoid(
    channel='y', area=delta_k*n_phase/2, duration=pp.calc_duration(gx_pre), system=system
    )
# max_gy = pp.make_trapezoid(channel='y', flat_area=delta_k*n_phase/2, flat_time=pe_duration, system=system)
scale_gy = np.arange(n_phase) - n_phase//2
scale_gy = scale_gy/n_phase
## Phase encode rewinder, is it necessary?
max_gy_rew = pp.make_trapezoid(
    channel='y', area=delta_k*n_phase/2, duration=pp.calc_duration(gx_pre), system=system
    )

## Allow for arbitrary phase encoding scheme
## https://github.com/pulseq/MR-Physics-with-Pulseq/blob/main/slides/D1-0930-PulSeqCourse_IntroductionMRI_Michael_Bock.pdf
## Linear
# phase_encode_table = [i for i in range(n_phase)]

## Center-out order
assert n_phase%2 == 0, 'Number of shots must be even'
phase_encode_table = []
for i in range(n_phase):
        phase_encode_table.append(n_phase//2 + (i//2) * (1 if (i)%2==0 else -1) - (0 if i%2 == 0 else 1))

######################## Z GRADIENTS ######################## 
## Calculate z-crusher gradient
gz_sp_pre = pp.make_trapezoid(
    channel='z', area=2*gx_sp.area,  system=system
    )
gz_sp_pre.delay = pp.calc_duration(gx_sp) - pp.calc_duration(gz_sp_pre)
gz_sp_post = pp.make_trapezoid(
    channel='z', area=2*gx_sp.area,  system=system
    )

## Set gradients=0 and nphase=1 to test for echo symmetry
# n_phase = 1
# max_gy.amplitude = 0
# max_gy_rew.aplitude = 0
# gx.amplitude = 0
# gx_pre.amplitude = 0
# gx_sp.amplitude = 0

######################## TIMING ########################
te_fill1 = (
    math.ceil(
    (
        TE/2 - 
        (
        pp.calc_duration(gz_ex)/2 
        + pp.calc_duration(max_gy, gzr_ex, gx_pre) 
        + pp.calc_duration(gx_sp, gz_sp_pre) 
        + pp.calc_duration(gz_ref)/2
        )
    ) / seq.grad_raster_time
    ) * seq.grad_raster_time
)
te_fill2 = (
    math.ceil(
    (
        TE/2 - 
        (
        pp.calc_duration(gz_ref)/2 
        + pp.calc_duration(gx_sp, gz_sp_post) 
        + pp.calc_duration(gx)/2
        )
    ) / seq.grad_raster_time
    ) * seq.grad_raster_time
)

tr_fill = (
    math.ceil(
    (
        TR - (
        pp.calc_duration(gx_sp, max_gy_rew) 
        + pp.calc_duration(gz_ex)/2 
        + TE 
        + pp.calc_duration(gx)/2
        )
    ) / seq.grad_raster_time
    ) * seq.grad_raster_time
)

assert(te_fill1 >= 0)
assert(te_fill2 >= 0)
assert(tr_fill >= 0)

######################## CONSTRUCT SEQUENCE ########################
pe_index = 0
for i in range(n_phase): 
    ## Scale PE gradients 
    gy = pp.scale_grad(max_gy, scale_gy[phase_encode_table[pe_index]])
    gy_rew = pp.scale_grad(max_gy_rew, scale_gy[phase_encode_table[pe_index]])    
    ## Vary RF phase quasi-randomly, is it necessary for SE?
    # rand_phase = np.mod(117*(pe_index^2 + pe_index + 2), 360)*np.pi/180
    # rf_ex.phase_offset = rand_phase
    # rf_ref.phase_offset = rand_phase
    # adc.phase_offset = rand_phase
    
    seq.add_block(rf_ex, gz_ex)
    seq.add_block(gx_pre, gy, gzr_ex)
    
    seq.add_block(pp.make_delay(te_fill1))
    
## Use 180 slice select rf
    seq.add_block(gx_sp,gz_sp_pre)
    seq.add_block(rf_ref, gz_ref)
    seq.add_block(gx_sp, gz_sp_post)
## Use non-select rf    
    # seq.add_block(gx_sp,gz_sp_pre)
    # seq.add_block(rf_ref)
    # seq.add_block(gx_sp, gz_sp_post)

    seq.add_block(pp.make_delay(te_fill2))
    
    seq.add_block(adc, gx)
    seq.add_block(gx_sp, gy_rew)
    
    seq.add_block(pp.make_delay(tr_fill))
    
    pe_index += 1
    
## Check whether the timing of the sequence is compatible with the scanner
(ok,error_report,) = seq.check_timing()  # Check whether the timing of the sequence is correct
if ok:
    print("Timing check passed successfully")
else:
    print("Timing check failed. Error listing follows:")
    [print(e) for e in error_report]

## Prepare the sequence output for the scanner
seq.set_definition('Name', 'se')
seq_file = 'SE_test.seq'
seq.write('./sequences/' + seq_file)

# %% Create phantom, simulate sequence, reconstruct image
seq = pp.Sequence()
seq.read('./sequences/' + seq_file)
seq0 = mr0.Sequence.import_file('./sequences/' + seq_file)

## Setup spin system/object on which we can run the MR sequence
sel_phantom = 'invivo' # select phantom type: invivo, simbrain, pixel
reso = (64, 64, 1)

print('load phantom')
if sel_phantom == 'pixel':
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
elif sel_phantom == 'simbrain':
    obj_p = mr0.VoxelGridPhantom.load_mat('./data/numerical_brain_cropped.mat')
    obj_p = obj_p.interpolate(reso[0], reso[1], reso[2])
    obj_p.B0[:] = 0
    obj_p.D[:] = 0
elif sel_phantom == 'invivo':
    obj_p = mr0.VoxelGridPhantom.brainweb("./data/subject05.npz")
    obj_p = obj_p.interpolate(reso[0], reso[1], 32).slices([16])
else:
    print('Select proper phantom')
obj_sim = obj_p.build()

## Simulate the external.seq file and add acquired signal to ADC plot
graph=mr0.compute_graph(seq0, obj_sim, 200, 1e-3)
signal=mr0.execute_graph(graph, seq0, obj_sim)
reco = mr0.reco_adjoint(signal, seq0.get_kspace(), resolution=reso, FOV=(0.22, 0.22, 1))
# %% Reconstruct and plot phantom, sequence and signal
plot_seq_zoom = True
plot_phantom = False 
plot_seq = False
plot_kspace = True 
plot_reco = True

if plot_seq_zoom: 
    seq.plot(time_range=(0, TE+ro_duration+pp.calc_duration(gz_ex)/2+pp.calc_duration(gx_sp)))
    plt.figure()
    plt.plot(signal.numpy())
    plt.title('ADC signal')
    
if plot_phantom:
    obj_p.plot()
    
if plot_seq:
    sp_adc, t_adc = mr0.util.pulseq_plot(seq=seq,signal=signal.numpy())

if plot_kspace:
    seq0.plot_kspace_trajectory()

if plot_reco:
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
# %%
