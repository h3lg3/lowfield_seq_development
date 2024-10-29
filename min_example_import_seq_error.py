import math
import warnings
import numpy as np
from matplotlib import pyplot as plt
import pypulseq as pp


#@title 3D SE - sequence
experiment_id = 's10_3d_se'
# %% SETUP system
# choose the scanner limits
system = pp.Opts(
        max_grad=25,       # MaRCoS's limits are 25, 40 and 35 mT/m for X, Y and Z axes
        grad_unit="mT/m",
        max_slew=50,       # MaRCoS's limits are 50, 80 and 70 mT/m/ms for X, Y and Z, respectively
        slew_unit="T/m/s",
        rf_ringdown_time=15e-6,
        rf_dead_time=15e-6,
        adc_dead_time=0e-6 )
# Create a new sequence object
seq = pp.Sequence(system)

# %% DEFINE the sequence, FOV, resolution, and other parameters
fov_mm = (150, 150, 150) # Define FOV in [mm]
nRD, nPH, n3D = [40, 40, 10]     # Define resolution (matrix sizes)
TE = 15e-3 # Echo time
TR = 300e-3  # Repetition time

dG = 500e-6 # ramping time for all gradients
sampling_time = 4.2e-3
t_ex  = 60e-6 # [s], duration of the excitation pulse.
t_ref = 100e-6 # [s], duration of the refocusing pulse
fsp_r = 1 # spoiler area in the read direction in parts of the read gradient area
rf_ex_phase = np.pi / 2 # phase of the excitation pulse
rf_ref_phase = 0 # phase of the refocusing pulse

# derived and modifed parameters
fov = np.array(fov_mm) * 1e-3 # FOV in meters
TE = round(TE/system.grad_raster_time/2) * system.grad_raster_time * 2 # TE (=ESP) should be divisible to a double gradient raster, which simplifies calcuations
rf_add = math.ceil(max(system.rf_dead_time,system.rf_ringdown_time)/system.grad_raster_time)*system.grad_raster_time # round up dead times to the gradient raster time to enable correct TE & ESP calculation
t_sp = round((0.5 * (TE - sampling_time - t_ref) - rf_add)/system.grad_raster_time)*system.grad_raster_time # the duration of gradient spoiler after the refocusing pulse
t_spex = round((0.5 * (TE - t_ex - t_ref) - 2*rf_add)/system.grad_raster_time)*system.grad_raster_time # the duration of readout prephaser after the excitation pulse. note: exclude the RF ringdown time of excitation pulse and rf dead time of refocusing pulse

# ======
# CREATE EVENTS
# ======
# excitation and refocusing pulses
flip_ex = 90 * np.pi / 180
rf_ex = pp.make_block_pulse(
    flip_angle=flip_ex,
    system=system,
    duration=t_ex,
    delay=rf_add,
    phase_offset=rf_ex_phase )

d_ex=pp.make_delay(t_ex+rf_add*2)

flip_ref = 180 * np.pi / 180
rf_ref = pp.make_block_pulse(
    flip_angle=flip_ref,
    system=system,
    duration=t_ref,
    delay=rf_add,
    phase_offset=rf_ref_phase,
    use="refocusing" )

d_ref=pp.make_delay(t_ref+rf_add*2)

delta_kx = 1 / fov[0]
rd_amp = nRD * delta_kx / sampling_time

gr_acq = pp.make_trapezoid(
    channel="x",
    system=system,
    amplitude = rd_amp,
    flat_time=sampling_time,
    delay=0,
    rise_time=dG )

adc = pp.make_adc(
    num_samples=nRD,
    dwell=sampling_time/nRD,
    delay=gr_acq.rise_time )

# Readout spoiler gradient
gr_spr = pp.make_trapezoid(
    channel="x",
    system=system,
    area=gr_acq.area * fsp_r,
    duration=t_sp,
    rise_time=dG )

# hint: Phase-encoding
delta_ky = 1 / fov[1]
gp_max = pp.make_trapezoid(
                channel="y",
                system=system,
                area=delta_ky*nPH/2,
                duration=t_sp,
                rise_time=dG )

# hint: Partition encoding
delta_kz = 1 / fov[2]
gs_max = pp.make_trapezoid(
                channel="z",
                system=system,
                area=delta_kz*n3D/2,
                duration=t_sp,
                rise_time=dG )

# readout prephaser: account for readout pre-spoiler and readout gradient
agr_preph = gr_acq.area / 2 + gr_spr.area
gr_preph = pp.make_trapezoid(
    channel="x", system=system, area=agr_preph, duration=t_spex, rise_time=dG )

# Fill-times
t_ex = pp.calc_duration(d_ex) + pp.calc_duration(gr_preph)
t_ref = pp.calc_duration(d_ref) + pp.calc_duration(gr_acq) + 2 * pp.calc_duration(gr_spr)

TR_fill = TR - t_ex - t_ref
# Round to gradient raster
TR_fill = system.grad_raster_time * np.round(TR_fill / system.grad_raster_time)
if TR_fill < 0:
    TR_fill = 1e-3
    warnings.warn(
        f"TR too short, adapted to: {1000 * (t_ex + t_ref + TR_fill)} ms"
    )
else:
    print(f"TR fill: {1000 * TR_fill} ms")
delay_TR = pp.make_delay(TR_fill)

# ======
# CONSTRUCT SEQUENCE
# ======
for Cz in range(-1,n3D): # partition encoding loop, -1 for dummy scan
    if Cz >= 0:
        sl_scale = (Cz-n3D/2)/n3D*2; # from -1 to +1
        nPH_range=range(nPH)
    else:
        sl_scale = 0.0
        nPH_range=range(1) # skip the nPH loop for dummy scan(s)

    gs=pp.scale_grad(gs_max, sl_scale)

    for Cy in nPH_range: # phase encoding loop
        seq.add_block(rf_ex, d_ex)
        seq.add_block(gr_preph)

        if Cz >= 0:
            pe_scale = (Cy-nPH/2)/nPH*2; # from -1 to 1
        else:
            pe_scale = 0.0

        gp=pp.scale_grad(gp_max, pe_scale)

        seq.add_block(rf_ref, d_ref)
        if Cz >= 0:
            seq.add_block(gs, gp, gr_spr)
            seq.add_block(gr_acq, adc)
            seq.add_block(pp.scale_grad(gs, -1), pp.scale_grad(gp, -1), gr_spr)
        else:
            seq.add_block(gs, gp, gr_spr)
            seq.add_block(gr_acq)
            seq.add_block(pp.scale_grad(gs, -1), pp.scale_grad(gp, -1), gr_spr)

        seq.add_block(delay_TR)

(ok, error_report) = seq.check_timing()  # Check whether the timing of the sequence is correct
if ok:
    print("Timing check passed successfully")
else:
    print("Timing check failed. Error listing follows:")
    [print(e) for e in error_report]

# ======
# VISUALIZATION
# ======
seq.plot(time_range=(TR, TR+0.04))

# Plot k-space trajectory
k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

plt.figure()
plt.plot(k_traj[0],k_traj[1], 'b-')
plt.plot(k_traj_adc[0],k_traj_adc[1],'r.')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.title('k-space trajectory')
plt.xlabel('kx')
plt.ylabel('ky')
plt.show()

# Prepare the sequence output for the scanner
seq.set_definition('Name', experiment_id)
seq.write(experiment_id+'.seq')


# Create a Sequence object and read the sequence file
seqX = pp.Sequence()
seqX.read(experiment_id+'.seq')
# %% Analyze sequence
# Plot k-space trajectory
k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seqX.calculate_kspace()

plt.figure()
plt.title('full k-space trajectory ($k_{x}$ x $k_{y}$)')
plt.plot(k_traj[0],k_traj[1])
plt.plot(k_traj_adc[0],k_traj_adc[1],'.')
plt.show()