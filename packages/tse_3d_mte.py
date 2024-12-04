"""Constructor for 3D TSE Imaging sequence.
QUESTIONS: Book Tofts p.159: Alternate echoes of 180 echo train with have error proportional to 180° imperfection and should not be used in fitting. 
TODO: add sampling patterns (elliptical masks, partial fourier, CS)
TODO: add optional variable refocussing pulses (pass list rather than float)
"""
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
from packages.write_seq_definitions import write_seq_definitions

import warnings

def constructor(
    echo_time: float = 15e-3,   # should be named echo spacing (esp), sequence should calculate effective TE (sampling of k-space center)
    repetition_time: float = 600e-3,
    etl: int = 7,               # number of echoes in echo train, etl*esp gives total sampling duration for 1 excitation pulse, should be in the order of 2*T2? 
    dummies: int = 0,
    rf_duration: float = 400e-6,
    gradient_correction: float = 0.,
    ro_bandwidth: float = 20e3,
    ro_oversampling: int = 5,
    input_fov: Dimensions = Dimensions(x=220e-3, y=220e-3, z=225e-3),
    input_enc: Dimensions = Dimensions(x=70, y=70, z=49),
    excitation_angle: float = pi / 2,
    excitation_phase: float = 0.,
    refocussing_angle: float = pi,
    refocussing_phase: float = pi / 2,
    channels: Channels = Channels(ro="y", pe1="z", pe2="x"),
    system:Opts = default_system,
) -> tuple[pp.Sequence, list, list]:
    """Construct 3D turbo spin echo sequence.

    Parameters
    ----------
    echo_time, optional
        Time constant between center of 90 degree pulse and center of ADC, by default 15e-3
    repetition_time, optional
        Time constant between two subsequent 90 degree pulses (echo trains), by default 600e-3
    etl, optional
        Echo train length, by default 7
    dummies, optional
        Number of dummy shots to acquire, default is 0
    rf_duration, optional
        Duration of the RF pulses (90 and 180 degree), by default 400e-6
    gradient_correction, optional
        Time constant to center ADC event, by default 510e-6
    adc_correction, optional
        Time constant which is added at the end of the ADC and readout gradient.
        This value is not taken into account for the prephaser calculation.
    ro_bandwidth, optional
        Readout bandwidth in Hz, by default 20e3
    fov, optional
        Field of view per dimension, by default default_fov
    n_enc, optional
        Number of encoding steps per dimension, by default default_encoding = Dimensions(x=70, y=70, z=49).
        If an encoding dimension is set to 1, the TSE sequence becomes a 2D sequence.
    excitation_angle, excitation_phase, optional
        set the flip angle and phase of the excitation pulse in radians, defaults to 90 degree flip angle, 0 phase
    refocussing_angle, refocussing_phase, optional
        Set the flip angle and phase of the refocussing pulse in radians,
        defaults to 180 degree flip angle, 90 degree phase
        TODO: allow this to be a list/array to vary flip angle along echo train.
    channel_ro, channel_pe1, channel_pe2, optional
        set the readout, phase1 and phase2 encoding directions, default to y, z and x.

    Returns
    -------
        Pulseq sequence and a list which describes the trajectory
    """
    #  additional ADC samples to initialize decimation filters prior to the 'sampling_time'
    nRD_pre = 0
    nRD_post = 0

    # Create a new sequence object
    seq = pp.Sequence(system)

    # map fov and n_enc with channels  
    n_enc, fov = seq_utils.map_fov_enc(channels, input_fov, input_enc)
                                           
    # derived and modifed parameters
    delta_k_ro = 1/fov['ro']
    adc_dwell_time = raster(1 / ro_bandwidth, precision=system.grad_raster_time) # sample everything on grad_raster_time
    adc_duration = n_enc['ro'] * adc_dwell_time    
    #adc_duration = raster(n_enc['ro'] / ro_bandwidth, precision=system.grad_raster_time)   # seems more robust to use dwell time instead of duration due to rounding errors when 1/ro_bandwidth is calculated
    ro_amplitude = n_enc['ro'] * delta_k_ro / adc_duration
    gradient_correction = raster(gradient_correction, precision=system.grad_raster_time)

    echo_time = round(echo_time/system.grad_raster_time/2) * system.grad_raster_time * 2 # TE (=ESP) should be divisible to a double gradient raster, which simplifies calcuations
    ro_flattop_time = adc_duration + 2*gradient_correction ; # duration of the flat top of the read gradient
    rf_add = math.ceil(max(system.rf_dead_time,system.rf_ringdown_time)/system.grad_raster_time)*system.grad_raster_time # round up dead times to the gradient raster time to enable correct TE & ESP calculation
    t_sp = round((0.5 * (echo_time - ro_flattop_time - rf_duration) - rf_add)/system.grad_raster_time)*system.grad_raster_time # the duration of gradient spoiler after the refocusing pulse
    t_spex = round((0.5 * (echo_time - rf_duration - rf_duration) - 2*rf_add)/system.grad_raster_time)*system.grad_raster_time # the duration of readout prephaser after the excitation pulse. note: exclude the RF ringdown time of excitation pulse and rf dead time of refocusing pulse

    # Definition of RF pulses
    rf_ex = pp.make_block_pulse(
        system=system,
        flip_angle=excitation_angle,
        phase_offset=excitation_phase,
        duration=rf_duration,
        delay=rf_add,
        use="excitation"
    )
    d_ex=pp.make_delay(rf_duration+rf_add*2)

    rf_ref = pp.make_block_pulse(
        system=system,
        flip_angle=refocussing_angle,
        phase_offset=refocussing_phase,
        duration=rf_duration,
        delay=rf_add,
        use="refocusing"
    )
    d_ref=pp.make_delay(rf_duration+rf_add*2)

    # Define readout gradient, spoiler and prewinder
    grad_ro = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        amplitude=ro_amplitude,
        flat_time=ro_flattop_time, # Add gradient correction time and ADC correction time
        delay=t_sp
    )
    
    # Readout spoiler gradient
    grad_ro_spr = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        area=grad_ro.area,  # grad_ro_spr.area = 0 why not same as set 0 here
        duration=t_sp,
        )
    
    # # Calculate readout prephaser without correction timess
    grad_ro_pre = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        area=grad_ro.area / 2 + grad_ro_spr.area,
        duration=t_spex,
    )
     
    adc = pp.make_adc(
        system=system,
        num_samples=(n_enc['ro']+nRD_pre+nRD_post) * ro_oversampling,
        duration=adc_duration,
        delay=t_sp-nRD_pre*adc_duration/n_enc['ro']# nRD_pre*sampling_time/nRD: delay for additional ADC samples to initialize decimation filters prior to the 'sampling_time'
    )
    
    # Phase-encoding
    delta_k_pe1 = 1 / fov['pe1']
    grad_pe_1_max = pp.make_trapezoid(
                    channel= channels.pe1,
                    system=system,
                    area=delta_k_pe1*n_enc['pe1']/2,
                    duration=t_sp,
                    )

    # Partition encoding
    delta_k_pe2 = 1 / fov['pe2']
    grad_pe_2_max = pp.make_trapezoid(
                    channel=channels.pe2,
                    system=system,
                    area=delta_k_pe2*n_enc['pe2']/2,
                    duration=t_sp)
    
    # combine parts of the read gradient
    gc_times = np.array(
        [
            0,
            grad_ro_spr.rise_time,
            grad_ro_spr.flat_time,
            grad_ro_spr.fall_time,
            grad_ro.flat_time,
            grad_ro_spr.fall_time,
            grad_ro_spr.flat_time,
            grad_ro_spr.rise_time ] )
    gc_times = np.cumsum(gc_times)

    gr_amp = np.array([0, 
                       grad_ro_spr.amplitude, 
                       grad_ro_spr.amplitude,
                       grad_ro.amplitude,
                       grad_ro.amplitude,
                       grad_ro_spr.amplitude,
                       grad_ro_spr.amplitude,
                       0])
    gr = pp.make_extended_trapezoid(
        channel=channels.ro, 
        times=gc_times, 
        amplitudes=gr_amp)

    gp_amp = np.array([0,
                       grad_pe_1_max.amplitude,
                       grad_pe_1_max.amplitude, 
                       0, 
                       0, 
                       -grad_pe_1_max.amplitude, 
                       -grad_pe_1_max.amplitude, 
                       0])
    grad_pe_1_max = pp.make_extended_trapezoid(
        channel=channels.pe1, 
        times=gc_times, 
        amplitudes=gp_amp)

    gs_amp = np.array([0, 
                       grad_pe_2_max.amplitude, 
                       grad_pe_2_max.amplitude, 
                       0, 
                       0, 
                       -grad_pe_2_max.amplitude, 
                       -grad_pe_2_max.amplitude, 
                       0])
    grad_pe_2_max = pp.make_extended_trapezoid(
        channel=channels.pe2, 
        times=gc_times, 
        amplitudes=gs_amp)
    
    # Fill-times
    t_ex = pp.calc_duration(d_ex) + pp.calc_duration(grad_ro_pre)
    t_ref = pp.calc_duration(d_ref) + pp.calc_duration(gr)

    t_train = t_ex + etl * t_ref

    TR_fill = repetition_time - t_train

    # Round to gradient raster
    TR_fill = system.grad_raster_time * np.round(TR_fill / system.grad_raster_time)
    if TR_fill < 0:
        TR_fill = 1e-3
        warnings.warn(
            f"TR too short, adapted to: {1000 * (t_train + TR_fill)} ms"
        )
    else:
        print(f"TR fill: {1000 * TR_fill} ms")
    delay_TR = pp.make_delay(TR_fill)

    # ======
    # CONSTRUCT SEQUENCE
    # ======
    for C_pe2 in range(-dummies, n_enc['pe2']):
        if C_pe2 >= 0:
            pe2_scale = (C_pe2 - n_enc['pe2']/2)/n_enc['pe2']*2; # from -1 to +1
            n_pe1_range = range(n_enc['pe1'])
        else:
            pe2_scale = 0.0
            n_pe1_range = range(1) # skip the nPH loop for dummy scan(s)

        g_pe2 = pp.scale_grad(grad_pe_2_max, pe2_scale)

        for C_pe1 in n_pe1_range:
            seq.add_block(rf_ex, d_ex)
            seq.add_block(grad_ro_pre)

            if C_pe2 >= 0:
                pe1_scale = (C_pe1 - n_enc['pe1']/2)/n_enc['pe1']*2; # from -1 to 1
            else:
                pe1_scale = 0.0

            g_pe1 = pp.scale_grad(grad_pe_1_max, pe1_scale)

            for k_echo in range(etl):

                label_pe1 = pp.make_label(type="SET", label="LIN", value=int(C_pe1))
                label_pe2 = pp.make_label(type="SET", label="PAR", value=int(C_pe2))
                label_echo = pp.make_label(type="SET", label="PAR", value=int(k_echo))

                seq.add_block(rf_ref, d_ref)
                if C_pe2 >= 0:
                    seq.add_block(g_pe2, g_pe1, gr, adc, label_pe1, label_pe2, label_echo)
                else:
                    seq.add_block(g_pe2, g_pe1, gr)

            seq.add_block(delay_TR)
    
    header = seq_utils.create_ismrmd_header(
                n_enc = n_enc,
                fov = fov,
                system = system
                )
    
    # write all required parameters in the seq-file definitions.
    write_seq_definitions(
        seq = seq,
        fov = (fov["ro"], fov["pe1"], fov["pe2"]),
        name = "tse_3d_mte",
        alpha = excitation_angle,
        slice_thickness = fov['pe2']/n_enc['pe2'],
        Nx = n_enc['ro'],
        Ny = n_enc['ro'],
        sampling_scheme = 'cartesian',
        Nz = n_enc['pe2'],
        TE = echo_time,
        TR = repetition_time,
        etl = etl,
        ro_bandwidth = ro_bandwidth,
        ro_oversampling = ro_oversampling
    )

    return (seq, header)