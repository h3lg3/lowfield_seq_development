"""Constructor for 3D TSE Imaging sequence.

TODO: add sampling patterns (elliptical masks, partial fourier, CS)
TODO: add optional inversion pulse
TODO: add optional variable refocussing pulses (pass list rather than float)
TODO: Design goal: Allow for minimal TE/Echo spacing for maximal ETL (acceleration)? IF TE:None is set, then use min TE
TODO: Optimize gradients by gradient surgery
TODO: REcalculate trajectory such that gradient scaling can be used in loop: Get max Grad Amp and scale according to encoding step: gs=pp.scale_grad(gs_max, sl_scale)
"""
from math import pi
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
    etl: int = 7,               # etl*esp gives total sampling duration for 1 excitation pulse, should be in the order of 2*T2? 
    dummies: int = 0,
    rf_duration: float = 400e-6,
    gradient_correction: float = 0, # delete 2*gradient_correction?
    ro_bandwidth: float = 20e3,
    ro_oversampling: int = 5,
    input_fov: Dimensions = Dimensions(x=220e-3, y=220e-3, z=225e-3),
    input_enc: Dimensions = Dimensions(x=70, y=70, z=49),
    echo_shift: float = 0.0,    # delete?
    trajectory: seq_utils.Trajectory = seq_utils.Trajectory.OUTIN,
    excitation_angle: float = pi / 2,
    excitation_phase: float = 0.,
    refocussing_angle: float = pi,
    refocussing_phase: float = pi / 2,
    inversion_pulse: bool = False,
    inversion_time: float = 50e-3,
    inversion_angle: float = pi,
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
    trajectroy, optional
        The k-space trajectory, by default set to in-out, other currently implemented options are...
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
    # Sequence options
    disable_pe = False
    disable_spoiler = True
    
    # Create a new sequence object
    seq = pp.Sequence(system)

    # map fov and n_enc with channels  
    n_enc, fov = seq_utils.map_fov_enc(channels, input_fov, input_enc)
                                           
    # derived and modifed parameters
    delta_k_ro = 1/fov['ro']
    adc_dwell_time = raster(1 / ro_bandwidth, precision=system.grad_raster_time) # sample everything on grad_raster_time
    adc_duration = n_enc['ro'] * adc_dwell_time
    # adc_duration = raster(n_enc['ro'] / ro_bandwidth, precision=system.grad_raster_time)   # seems more robust to use dwell time instead of duration due to rounding errors when 1/ro_bandwidth is calculated
    grad_amplitude = n_enc['ro'] * delta_k_ro / adc_duration
    gradient_correction = raster(gradient_correction, precision=system.grad_raster_time)

    pe_traj = seq_utils.get_traj(
                n_enc=n_enc,
                etl = etl,
                trajectory = trajectory,
                )
    
    trains, trains_pos = seq_utils.get_trains(
                pe_traj = pe_traj,
                n_enc=n_enc,
                fov=fov,
                etl = etl,
                trajectory = trajectory,
                )

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
    
    if inversion_pulse:
        rf_inversion = pp.make_block_pulse(
            system=system,
            flip_angle=inversion_angle,
            phase_offset=refocussing_phase,
            duration=rf_duration,
            use="refocusing"
        )

    # Define readout gradient, spoiler and prewinder
    grad_ro = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        amplitude=grad_amplitude,
        # Add gradient correction time and ADC correction time
        flat_time=adc_duration + 2 * gradient_correction # HH: why 2*gradient_correction?
    )
    ro_duration = pp.calc_duration(grad_ro)
    
    # Readout spoiler gradient
    grad_ro_spr = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        area=grad_ro.area,  # grad_ro_spr.area = 0 why not same as set 0 here
        duration=ro_duration,
        )
    ro_spr_duration = pp.calc_duration(grad_ro_spr)

    if disable_spoiler:
        warnings.warn('Spoiler are disabled')
        grad_ro_spr.amplitude = 0   # necessary for event creation
        grad_ro_spr.area = 0        # necessary for creating ro prewinder
    
    # # Calculate readout prephaser without correction timess
    grad_ro_pre = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        area=grad_ro.area / 2 + grad_ro_spr.area,
        duration=ro_duration,
    )
    ro_pre_duration = pp.calc_duration(grad_ro_pre) 
     
    adc = pp.make_adc(
        system=system,
        num_samples=n_enc['ro']*ro_oversampling,
        duration=adc_duration,
        # Add gradient correction time and ADC correction time
        delay=2 * gradient_correction + grad_ro.rise_time
    )
    
    # Calculate delays
    # Note: RF dead-time is contained in RF delay
    # Delay duration center excitation (90 degree) and center refocussing (180 degree) RF pulse
    tau_1 = raster(val=echo_time/2 - rf_duration/2 - rf_90.ringdown_time - rf_duration/2 - rf_180.delay - ro_pre_duration, precision=system.grad_raster_time)
    
    # Delay duration between center refocussing (180 degree) RF pulse and center readout
    # check if gradient_correction should be added to this delay, maybe only 1*gradient_correction? and tau_3 should get also 1*gradient_correction?
    tau_2 = raster(echo_time/2 - adc_duration/2 - 2*gradient_correction - rf_duration/2 - rf_180.ringdown_time  \
        - grad_ro.rise_time - ro_spr_duration + echo_shift, precision=system.grad_raster_time)

    # Delay duration between center readout and next center refocussing (180 degree) RF pulse 
    tau_3 = raster(echo_time/2 - adc_duration/2 - rf_duration/2 - rf_180.delay  - grad_ro.rise_time - ro_spr_duration - echo_shift, precision=system.grad_raster_time)

    # recommended_timing = seq_utils.get_esp_etl(tau_1=tau_1, tau_2=tau_2, tau_3=tau_3, echo_time=echo_time, T2=100, n_enc_pe1=n_enc['pe1'])
    # print(recommended_timing)
    
    for _ in range(dummies):
        if inversion_pulse:
            seq.add_block(rf_inversion)
            seq.add_block(pp.make_delay(raster(val=inversion_time - rf_duration, precision=system.grad_raster_time)))
        seq.add_block(rf_90)
        seq.add_block(pp.make_delay(raster(val=echo_time / 2 - rf_duration, precision=system.grad_raster_time)))
        for _ in range(etl):
            seq.add_block(rf_180)
            seq.add_block(pp.make_delay(raster(
                val=echo_time - rf_duration,
                precision=system.grad_raster_time
            )))
        if inversion_pulse:
            seq.add_block(pp.make_delay(raster(
                val=repetition_time - (etl + 0.5) * echo_time - rf_duration - inversion_time,
                precision=system.grad_raster_time
            )))
        else:
            seq.add_block(pp.make_delay(raster(
                val=repetition_time - (etl + 0.5) * echo_time - rf_duration,
                precision=system.grad_raster_time
            )))

    for train, position in zip(trains, trains_pos):     
        if inversion_pulse:
            seq.add_block(rf_inversion)
            seq.add_block(pp.make_delay(raster(
                val=inversion_time - rf_duration,
                precision=system.grad_raster_time
            )))
        seq.add_block(rf_90)
        seq.add_block(grad_ro_pre)
        seq.add_block(pp.make_delay(tau_1))
        
        for echo, pe_indices in zip(train, position):
            pe_1, pe_2 = echo
            if disable_pe:
                warnings.warn('Phase encoding is disabled')
                pe_1 = 0
                pe_2 = 0
            
            seq.add_block(rf_180)
            seq.add_block(
                pp.make_trapezoid(
                    channel=channels.pe1,
                    area=pe_1,
                    duration=ro_spr_duration,
                    system=system,
                ),
                pp.make_trapezoid(
                    channel=channels.pe2,
                    area=pe_2,
                    duration=ro_spr_duration,
                    system=system,
                ),
                grad_ro_spr                
            )          

            seq.add_block(pp.make_delay(tau_2))
            
            # Cast index values from int32 to int, otherwise make_label function complains
            label_pe1 = pp.make_label(type="SET", label="LIN", value=int(pe_indices[0]))
            label_pe2 = pp.make_label(type="SET", label="PAR", value=int(pe_indices[1]))
            seq.add_block(grad_ro, adc, label_pe1, label_pe2)

            seq.add_block(
                pp.make_trapezoid(
                    channel=channels.pe1,
                    area=-pe_1,
                    duration=ro_spr_duration,
                    system=system,
                ),
                pp.make_trapezoid(
                    channel=channels.pe2,
                    area=-pe_2,
                    duration=ro_spr_duration,
                    system=system,
                ),
                grad_ro_spr
            )

            seq.add_block(pp.make_delay(tau_3))

        # recalculate TR each train because train length is not guaranteed to be constant
        tr_delay = raster(repetition_time - echo_time * len(train) - adc_duration / 2 - ro_spr_duration \
            - tau_3 - rf_90.delay - rf_duration / 2 - grad_ro.rise_time, precision=system.grad_raster_time)

        if inversion_pulse:
            tr_delay -= inversion_time
        
        seq.add_block(pp.make_delay(tr_delay))

    print(f"TR fill: {1000 * tr_delay} ms")
    
    # Check labels
    labels = seq.evaluate_labels(evolution="adc")
    acq_pos = np.concatenate(trains_pos).T
    
    if not np.array_equal(labels["LIN"], acq_pos[0, :]):
        raise ValueError("LIN labels don't match actual acquisition positions.")
    if not np.array_equal(labels["PAR"], acq_pos[1, :]):
        raise ValueError("PAR labels don't match actual acquisition positions.")
    
    # # Calculate some sequence measures
    train_duration_tr = (seq.duration()[0]) / len(trains)
    train_duration = train_duration_tr - tr_delay
    
    header = seq_utils.create_ismrmd_header(
                n_enc = n_enc,
                fov = fov,
                system = system
                )
        
        # write all required parameters in the seq-file definitions.
    write_seq_definitions(
        seq = seq,
        fov = (fov["ro"], fov["pe1"], fov["pe2"]),
        encoding_dim = (n_enc["ro"], n_enc["pe1"], n_enc["pe2"]),
        name = "tse_3d",
        alpha = excitation_angle,
        slice_thickness = fov['pe2']/n_enc['pe2'],
        Nx = n_enc['ro'],
        Ny = n_enc['ro'],
        sampling_scheme = 'cartesian',
        Nz = n_enc['pe2'],
        TE = echo_time,
        TR = repetition_time,
        train_duration = train_duration,
        n_total_trains = len(trains),
        tr_delay = tr_delay,
        channel_order = (channels.ro, channels.pe1, channels.pe2),
        etl = etl,
        ro_bandwidth = ro_bandwidth,
        ro_oversampling = ro_oversampling
    )    

    return (seq, header)

    # # check if channel labels are valid
    # channel_ro, channel_pe1, channel_pe2 = seq_utils.validate_inputs(channel_ro, channel_pe1, channel_pe2)

    # # map fov and n_enc according to channels    
    # n_enc_ro, fov_ro, n_enc_pe1, fov_pe1, n_enc_pe2, fov_pe2 = seq_utils.calculate_enc_fov_order(
    #                                                             channel_ro, channel_pe1, n_enc, fov
    #                                                             )

    # # Map trajectory to PE1 + PE2 samples in format required by sort_kspace()
    # # Commented because calculates the same as: np.concatenate(trains_pos).T
    # k_traj_adc = seq.calculate_kspacePP()[0]
    # acq_pos = seq_utils.calculate_acq_pos(k_traj_adc,n_enc_ro,channel_pe1,fov_pe1,channel_pe2,fov_pe2)