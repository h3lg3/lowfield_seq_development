from math import pi
import math
import numpy as np

import pypulseq as pp
from pypulseq.opts import Opts

from packages import seq_utils
from packages.seq_utils import Dimensions
from packages.seq_utils import Channels
from packages.seq_utils import raster
from packages.mr_systems import low_field as default_system

import warnings


def constructor(
    echo_time: float = 15e-3,   # should be named echo spacing (esp), sequence should calculate effective TE (sampling of k-space center)
    repetition_time: float = 600e-3,
    etl: int = 7,               # etl*esp gives total sampling duration for 1 excitation pulse, should be in the order of 2*T2? 
    dummies: int = 0,
    rf_duration: float = 400e-6,
    ramp_duration: float = 200e-6,
    ro_bandwidth: float = 20e3,
    ro_oversampling: int = 5,
    input_fov: Dimensions = Dimensions(x=220e-3, y=220e-3, z=225e-3),
    input_enc: Dimensions = Dimensions(x=70, y=70, z=49),
    trajectory: seq_utils.Trajectory = seq_utils.Trajectory.OUTIN,
    excitation_angle: float = pi / 2,
    excitation_phase: float = 0.,
    refocussing_angle: float = pi,
    refocussing_phase: float = pi / 2,
    channels: Channels = Channels(ro="y", pe1="z", pe2="x"),
    system:Opts = default_system,
) -> pp.Sequence:
    # -> tuple[pp.Sequence, list, list]:

    # Create a new sequence object
    seq = pp.Sequence(system)

    # DEFINE the sequence, FOV, resolution, and other parameters

    nRD, nPH, n3D = input_enc.x, input_enc.y, input_enc.z     # Define resolution (matrix sizes)
    fov = input_fov
    TE = echo_time # Echo time
    TR = repetition_time  # Repetition time

    dG = ramp_duration # ramping time for all gradients
    
    sampling_time = nRD / ro_bandwidth
    t_ex  = rf_duration # [s], duration of the excitation pulse.
    t_ref = rf_duration# [s], duration of the refocusing pulse
    fsp_r = 1 # spoiler area in the read direction in parts of the read gradient area
    rf_ex_phase = excitation_phase # phase of the excitation pulse
    rf_ref_phase = refocussing_phase # phase of the refocusing pulse

    # derived and modifed parameters
    TE = round(TE/system.grad_raster_time/2) * system.grad_raster_time * 2 # TE (=ESP) should be divisible to a double gradient raster, which simplifies calcuations
    rd_flattop_time = sampling_time ; # duration of the flat top of the read gradient
    rf_add = math.ceil(max(system.rf_dead_time,system.rf_ringdown_time)/system.grad_raster_time)*system.grad_raster_time # round up dead times to the gradient raster time to enable correct TE & ESP calculation
    t_sp = round((0.5 * (TE - rd_flattop_time - t_ref) - rf_add)/system.grad_raster_time)*system.grad_raster_time # the duration of gradient spoiler after the refocusing pulse
    t_spex = round((0.5 * (TE - t_ex - t_ref) - 2*rf_add)/system.grad_raster_time)*system.grad_raster_time # the duration of readout prephaser after the excitation pulse. note: exclude the RF ringdown time of excitation pulse and rf dead time of refocusing pulse

    # ======
    # CREATE EVENTS
    # ======
    # excitation and refocusing pulses
    flip_ex = excitation_angle
    rf_ex = pp.make_block_pulse(
        flip_angle=flip_ex,
        system=system,
        duration=t_ex,
        delay=rf_add,
        phase_offset=rf_ex_phase )

    d_ex=pp.make_delay(t_ex+rf_add*2)

    rf_ref = pp.make_block_pulse(
        flip_angle=refocussing_angle,
        system=system,
        duration=t_ref,
        delay=rf_add,
        phase_offset=rf_ref_phase,
        use="refocusing" )

    d_ref=pp.make_delay(t_ref+rf_add*2)

    delta_kx = 1 / fov.x
    rd_amp = nRD * delta_kx / sampling_time

    gr_acq = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        amplitude = rd_amp,
        flat_time=rd_flattop_time,
        delay=t_sp,
        rise_time=dG )

    adc = pp.make_adc(
        num_samples=nRD*ro_oversampling,
        dwell=sampling_time/nRD/ro_oversampling,
        delay=t_sp )

    gr_spr = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        area=gr_acq.area * fsp_r,
        duration=t_sp,
        rise_time=dG )

    # Phase-encoding
    delta_ky = 1 / fov.y
    gp_max = pp.make_trapezoid(
                    channel=channels.pe1,
                    system=system,
                    area=delta_ky*nPH/2,
                    duration=t_sp,
                    rise_time=dG )

    # Partition encoding
    delta_kz = 1 / fov.z
    gs_max = pp.make_trapezoid(
                    channel=channels.pe2,
                    system=system,
                    area=delta_kz*n3D/2,
                    duration=t_sp,
                    rise_time=dG )

    # Gradient surgery for readout gradients
    gc_times = np.array(
        [
            0,
            gr_spr.rise_time,
            gr_spr.flat_time,
            gr_spr.fall_time,
            gr_acq.flat_time,
            gr_spr.fall_time,
            gr_spr.flat_time,
            gr_spr.rise_time ] )
    gc_times = np.cumsum(gc_times)

    gr_amp = np.array([0, gr_spr.amplitude, gr_spr.amplitude, gr_acq.amplitude, gr_acq.amplitude, gr_spr.amplitude, gr_spr.amplitude, 0])
    gr = pp.make_extended_trapezoid(channel=channels.ro, times=gc_times, amplitudes=gr_amp)

    agr_preph = gr_acq.area / 2 + delta_kx / 2 + gr_spr.area # readout prephaser: account for readout pre-spoiler, readout gradient, and that sample are taken at the center of raster
    gr_preph = pp.make_trapezoid(
        channel="x", system=system, area=agr_preph, duration=t_spex, rise_time=dG )

    # Gradient surgery for phase encoding
    gp_amp = np.array([0, gp_max.amplitude, gp_max.amplitude, 0, 0, -gp_max.amplitude, -gp_max.amplitude, 0])
    gp_max = pp.make_extended_trapezoid(channel=channels.pe1, times=gc_times, amplitudes=gp_amp)

    # Gradient surgery for partition encoding
    gs_amp = np.array([0, gs_max.amplitude, gs_max.amplitude, 0, 0, -gs_max.amplitude, -gs_max.amplitude, 0])
    gs_max = pp.make_extended_trapezoid(channel=channels.pe2, times=gc_times, amplitudes=gs_amp)

    # Fill-times
    t_ex = pp.calc_duration(d_ex) + pp.calc_duration(gr_preph)
    t_ref = pp.calc_duration(d_ref) + pp.calc_duration(gr)

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
    for Cz in range(0-dummies,n3D):
        if Cz >= 0:
            if n3D == 1:
                sl_scale = 0.0
            else:
                sl_scale = (Cz-n3D/2)/n3D*2; # from -1 to +1
            nPH_range=range(nPH)
        else:
            sl_scale = 0.0
            nPH_range=range(1) # skip the nPH loop for dummy scan(s)

        gs=pp.scale_grad(gs_max, sl_scale)
        
        for Cy in nPH_range:
            seq.add_block(rf_ex, d_ex)
            seq.add_block(gr_preph)

            if Cz >= 0:
                pe_scale = (Cy-nPH/2)/nPH*2; # from -1 to 1
            else:
                pe_scale = 0.0
            if Cz==0:
                gp=pp.scale_grad(gp_max, pe_scale)

            seq.add_block(rf_ref, d_ref)
            if Cz >= 0:
                seq.add_block(gs, gp, gr, adc)
            else:
                seq.add_block(gs, gp, gr)

            seq.add_block(delay_TR)


    return (seq)
    # return (seq, acq_pos, [n_enc_ro, n_enc_pe1, n_enc_pe2])