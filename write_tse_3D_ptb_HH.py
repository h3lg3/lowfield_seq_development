# %%
from enum import Enum
from math import pi

import numpy as np
import pypulseq as pp

from console.interfaces.interface_acquisition_parameter import Dimensions
from console.utilities.sequences.system_settings import raster, system

from matplotlib import pyplot as plt

plot_animation = False
plot_kspace = False
plot_seq = False

disable_pe = False

write_seq = True
seq_file = 'tse_3D_ptb_HH_dummies.seq'

class Trajectory(Enum):
    """Trajectory type enum."""

    INOUT = 1
    LINEAR = 2
    OUTIN = 3

echo_time = 28e-3   # 12, 28
repetition_time = 2000e-3
etl = 8             # 16, 8
dummies = 5
rf_duration = 400e-6
ramp_duration= 200e-6
gradient_correction = 0.
ro_bandwidth = 10e3
echo_shift= 0.0
trajectory: Trajectory = Trajectory.LINEAR
excitation_angle = pi / 2
excitation_phase = 0.
refocussing_angle = pi
refocussing_phase = pi / 2
channel_ro, channel_pe1, channel_pe2 = 'x', 'y', 'z'

# %% 
fov_ro, fov_pe1, fov_pe2 = 220e-3, 220e-3, 220e-3
n_enc_ro, n_enc_pe1, n_enc_pe2 = 64, 64, 1  # 100, 100, 1 / 64, 64, 1
system.rf_ringdown_time = 0
system.max_slew = 100 * system.gamma

# adjust system to siemens scanner
system.adc_raster_time=1e-7
system.block_duration_raster=1e-5
system.grad_raster_time=1e-5
system.rf_raster_time=1e-6
system.rf_ringdown_time=100e-6
system.rf_dead_time=100e-6
system.adc_dead_time=10e-6

seq = pp.Sequence(system)

# Calculate center out trajectory
pe1 = np.arange(n_enc_pe1) - (n_enc_pe1 - 1) / 2
pe2 = np.arange(n_enc_pe2) - (n_enc_pe2 - 1) / 2

pe0_pos = np.arange(n_enc_pe1)
pe1_pos = np.arange(n_enc_pe2)

pe_points = np.stack([grid.flatten() for grid in np.meshgrid(pe1, pe2)], axis=-1)
pe_positions = np.stack([grid.flatten() for grid in np.meshgrid(pe0_pos, pe1_pos)], axis=-1)

pe_mag = np.sum(np.square(pe_points), axis=-1)  # calculate magnitude of all gradient combinations
pe_mag_sorted = np.argsort(pe_mag)

if trajectory is Trajectory.OUTIN:
    pe_mag_sorted = np.flip(pe_mag_sorted)

pe_traj = pe_points[pe_mag_sorted, :]  # sort the points based on magnitude
pe_order = pe_positions[pe_mag_sorted, :]  # kspace position for each of the gradients

if trajectory is Trajectory.LINEAR:
    center_pos = 1 / 2  # where the center of kspace should be in the echo train
    num_points = np.size(pe_mag_sorted)
    linear_pos = np.zeros(num_points, dtype=int) - 10
    center_point = int(np.round(np.size(pe_mag) * center_pos))
    odd_indices = 1
    even_indices = 1
    linear_pos[center_point] = pe_mag_sorted[0]

    for idx in range(1, num_points):
        # check if its in bounds first
        if center_point - (idx + 1) / 2 >= 0 and idx % 2:
            k_idx = center_point - odd_indices
            odd_indices += 1
        elif center_point + idx / 2 < num_points and idx % 2 == 0:
            k_idx = center_point + even_indices
            even_indices += 1
        elif center_point - (idx + 1) / 2 < 0 and idx % 2:
            k_idx = center_point + even_indices
            even_indices += 1
        elif center_point + idx / 2 >= num_points and idx % 2 == 0:
            k_idx = center_point - odd_indices
            odd_indices += 1
        else:
            print("Sorting error")
        linear_pos[k_idx] = pe_mag_sorted[idx]

    pe_traj = pe_points[linear_pos, :]  # sort the points based on magnitude
    pe_order = pe_positions[linear_pos, :]  # kspace position for each of the gradients

# calculate the required gradient area for each k-point
pe_traj[:, 0] /= fov_pe1
pe_traj[:, 1] /= fov_pe2

# Divide all PE steps into echo trains
num_trains = int(np.ceil(pe_traj.shape[0] / etl))
trains = [pe_traj[k::num_trains, :] for k in range(num_trains)]

# Create a list with the kspace location of every line of kspace acquired, in the order it is acquired
trains_pos = [pe_order[k::num_trains, :] for k in range(num_trains)]
acq_pos = []
for train_pos in trains_pos:
    acq_pos.extend(train_pos)
    
# # # Plot k-space acquisition order
if plot_animation:
    plt.ion()  # Turn on interactive mode
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    for array in trains_pos[0:10]:
        #ax1.clear()
        ax2.clear()
        
        # Plot the whole array on the first subplot
        ax1.scatter(array[:, 0], array[:, 1])
        ax1.set_title('Array wise Phase Encode Steps ky per Echotrain')
        ax1.set_xlabel('ky')
        ax1.set_xlabel('kz')
        ax1.set_xlim(min(pe_traj[:, 0]), max(pe_traj[:, 0]))
        ax1.set_ylim(min(pe_traj[:, 1]), max(pe_traj[:, 1]))
        ax2.set_xlim(min(pe_traj[:, 0]), max(pe_traj[:, 0]))
        ax2.set_ylim(min(pe_traj[:, 1]), max(pe_traj[:, 1]))
        # Plot each element of the array iteratively on the second subplot
        for point in array:
            ax2.scatter(point[0], point[1], color='r')
            ax2.set_title('Element-wise Phase Encode Steps ky per Echotrain')
            ax1.set_xlabel('ky')
            ax1.set_xlabel('kz')
            plt.draw()
            plt.pause(0.2)  # Pause for half a second

    plt.ioff()  # Turn off interactive mode
    plt.show()    

    plt.ion()  # Turn on interactive mode
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    plt.ion()  # Turn on interactive mode
    fig, ax1= plt.subplots(figsize=(10, 10))
    for array in trains_pos[0:10]:
        ax1.set_title('Array wise Phase Encode Steps ky per Echotrain')
        ax1.set_xlabel('# PE')
        ax1.set_xlabel('ky')
        # Plot the whole array on the first subplot
        ax1.plot(array[:, 0])
        plt.draw()
        plt.pause(0.5)  # Pause for half a second
    plt.ioff()  # Turn off interactive mode
    plt.show()          
# %%

# Definition of RF pulses
rf_90 = pp.make_block_pulse(
    system=system,
    flip_angle=excitation_angle,
    phase_offset=excitation_phase,
    duration=rf_duration,
    use="excitation"
)
rf_180 = pp.make_block_pulse(
    system=system,
    flip_angle=refocussing_angle,
    phase_offset=refocussing_phase,
    duration=rf_duration,
    use="refocusing"
)

# ADC duration
adc_duration = n_enc_ro / ro_bandwidth

# Define readout gradient and prewinder
grad_ro = pp.make_trapezoid(
    channel=channel_ro,
    system=system,
    flat_area=n_enc_ro / fov_ro,
    rise_time=ramp_duration,
    fall_time=ramp_duration,
    # Add gradient correction time and ADC correction time
    flat_time=raster(adc_duration, precision=system.grad_raster_time),
)
grad_ro = pp.make_trapezoid(  # using the previous calculation for the amplitde, hacky, should find a better way
    channel=channel_ro,
    system=system,
    amplitude=grad_ro.amplitude,
    rise_time=ramp_duration,
    fall_time=ramp_duration,
    # Add gradient correction time
    flat_time=raster(adc_duration + 2 * gradient_correction, precision=system.grad_raster_time),
)

# Calculate readout prephaser without correction times
ro_pre_duration = pp.calc_duration(grad_ro) / 2

grad_ro_pre = pp.make_trapezoid(
    channel=channel_ro,
    system=system,
    area=grad_ro.area / 2,
    rise_time=ramp_duration,
    fall_time=ramp_duration,
    duration=raster(ro_pre_duration, precision=system.grad_raster_time),
)

## Spoiler gradient on x (used three times: before excitation (or after ADC), before refocusing, after refocusing) 
grad_ro_sp = pp.make_trapezoid(
    channel='x', area=1*grad_ro.area, duration=pp.calc_duration(grad_ro), system=system
    )


adc = pp.make_adc(
    system=system,
    num_samples = n_enc_ro,
    #num_samples=int((adc_duration) / system.adc_raster_time), 64000 ADC samples?
    duration=raster(val=adc_duration, precision=system.adc_raster_time),
    # Add gradient correction time and ADC correction time
    delay=raster(val=2 * gradient_correction + grad_ro.rise_time, precision=system.adc_raster_time)
)

# Calculate delays
# Note: RF dead-time is contained in RF delay
# Delay duration between RO prephaser after initial 90 degree RF and 180 degree RF pulse

# Muss hier die rf-duration nicht durch 2 geteilt werden?  Vielleicht weil 180 genauso lang wie 90 ist.
tau_1 = echo_time / 2 - rf_duration - rf_90.ringdown_time - rf_180.delay - ro_pre_duration - pp.calc_duration(grad_ro_sp)
# Delay duration between Gy, Gz prephaser and readout
tau_2 = (echo_time - rf_duration - adc_duration) / 2 - 2 * gradient_correction \
    - ramp_duration - rf_180.ringdown_time - ro_pre_duration + echo_shift - pp.calc_duration(grad_ro_sp)
# Delay duration between readout and Gy, Gz gradient rephaser
tau_3 = (echo_time - rf_duration - adc_duration) / 2 - ramp_duration - rf_180.delay - pp.calc_duration(grad_ro_pre) - echo_shift

for dummy in range(dummies):
    seq.add_block(rf_90)
    seq.add_block(pp.make_delay(raster(val=echo_time / 2 - rf_duration, precision=system.grad_raster_time)))
    for idx in range(etl):
        seq.add_block(rf_180)
        seq.add_block(pp.make_delay(raster(
            val=echo_time - rf_duration,
            precision=system.grad_raster_time
        )))
    seq.add_block(pp.make_delay(raster(
        val=repetition_time - (etl + 0.5) * echo_time - rf_duration,
        precision=system.grad_raster_time
    )))

for train in trains:
    seq.add_block(rf_90)
    seq.add_block(grad_ro_pre)
    seq.add_block(pp.make_delay(raster(val=tau_1, precision=system.grad_raster_time)))

    for echo in train:
        pe_1, pe_2 = echo
        if disable_pe:
            pe_1 = 0
            pe_2 = 0
        seq.add_block(grad_ro_sp)
        seq.add_block(rf_180)

        seq.add_block(
            pp.make_trapezoid(
                channel=channel_pe1,
                area=-pe_1,
                duration=ro_pre_duration, # warum werden diese nicht auf Rasterzeit angepasst? 
                system=system,
                rise_time=ramp_duration,
                fall_time=ramp_duration
            ),
            pp.make_trapezoid(
                channel=channel_pe2,
                area=-pe_2,
                duration=ro_pre_duration,
                system=system,
                rise_time=ramp_duration,
                fall_time=ramp_duration
            ),
            grad_ro_sp
        )

        seq.add_block(pp.make_delay(raster(val=tau_2, precision=system.grad_raster_time)))

        seq.add_block(grad_ro, adc)

        seq.add_block(
            pp.make_trapezoid(
                channel=channel_pe1,
                area=pe_1,
                duration=ro_pre_duration,
                system=system,
                rise_time=ramp_duration,
                fall_time=ramp_duration
            ),
            pp.make_trapezoid(
                channel=channel_pe2,
                area=pe_2,
                duration=ro_pre_duration,
                system=system,
                rise_time=ramp_duration,
                fall_time=ramp_duration
            )
        )

        seq.add_block(pp.make_delay(raster(val=tau_3, precision=system.grad_raster_time)))
        
    seq.add_block(grad_ro_sp) # add spoiler after last 180 pulse in echo train
    
    # recalculate TR each train because train length is not guaranteed to be constant
    tr_delay = repetition_time - echo_time * len(train) - adc_duration / 2 - ro_pre_duration \
        - tau_3 - rf_90.delay - rf_duration / 2 - ramp_duration - pp.calc_duration(grad_ro_sp)

    seq.add_block(pp.make_delay(raster(
        val=tr_delay,
        precision=system.block_duration_raster
    )))

## Check whether the timing of the sequence is compatible with the scanner
(ok,error_report,) = seq.check_timing()  # Check whether the timing of the sequence is correct
if ok:
    print("Timing check passed successfully")
else:
    print("Timing check failed. Error listing follows:")
    [print(e) for e in error_report]

## kspace trajectory: PLOT true and calculated kspace locations
# Plot k-space trajectory
if plot_kspace:
    k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

    plt.figure()
    plt.plot(k_traj[0],k_traj[1])
    plt.plot(k_traj_adc[0],k_traj_adc[1],'.')
    plt.show()
    
    # plot acquired and true k space locations
    # plt.figure()
    # plt.plot(k_traj[0],k_traj[1], 'x')
    # #plt.xlim(-238, -200)
    # plt.xlim(-234.1, -233.9)
    # plt.ylim(-226, -200)
    # plt.plot(k_traj_adc[0],k_traj_adc[1],'.')
    # plt.plot(np.ones(len(pe_traj[:, 0]))*-234,pe_traj[:, 0], 'o')
    # plt.title('Precalculated vs actual locations')
    # plt.legend(['actual loc', 'precalc loc'])
    # plt.show()

if plot_seq:
    seq.plot()
## Prepare the sequence output for the scanner
if write_seq:
    seq.set_definition('Name', 'se_ptb')
    seq.write('./sequences/' + seq_file)

# %%
