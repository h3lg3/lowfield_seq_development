"""Constructor for 3D TSE Imaging sequence.

TODO: add sampling patterns (elliptical masks, partial fourier, CS)
TODO: add optional inversion pulse
TODO: add optional variable refocussing pulses (pass list rather than float)
TODO: Design goal: Allow for minimal TE/Echo spacing for maximal ETL (acceleration)?
TODO: Alternate phase of 180 refocussing pulse by 180° each time to compensate imperfect excitaion (!=180 degree)
TODO: Instead of mapping fov to fov_ro, fov_pe1, fov_pe2,,, just get indices (e.g.,rO=1,pe1=2,pe2=0), and index fov[id.ro],n_enc[id.ro] etc throughout the script
"""
from math import pi
import numpy as np

import pypulseq as pp
from pypulseq.opts import Opts

from packages import seq_utils
from packages.seq_utils import Dimensions
from packages.seq_utils import raster
from packages.mr_systems import low_field as default_system

import warnings

default_fov = Dimensions(x=220e-3, y=220e-3, z=225e-3)
default_encoding = Dimensions(x=70, y=70, z=49)

def constructor(
    echo_time: float = 15e-3,   # should be named echo spacing (esp), sequence should calculate effective TE (sampling of k-space center)
    repetition_time: float = 600e-3,
    etl: int = 7,               # etl*esp gives total sampling duration for 1 excitation pulse, should be in the order of 2*T2? 
    dummies: int = 0,
    rf_duration: float = 400e-6,
    ramp_duration: float = 200e-6,
    gradient_correction: float = 0.,
    ro_bandwidth: float = 20e3,
    ro_oversampling: int = 5,
    fov: Dimensions = default_fov,
    n_enc: Dimensions = default_encoding,
    echo_shift: float = 0.0,
    trajectory: seq_utils.Trajectory = seq_utils.Trajectory.OUTIN,
    excitation_angle: float = pi / 2,
    excitation_phase: float = 0.,
    refocussing_angle: float = pi,
    refocussing_phase: float = pi / 2,
    inversion_pulse: bool = False,
    inversion_time: float = 50e-3,
    inversion_angle: float = pi,
    channel_ro: str = "y",
    channel_pe1: str = "z",
    channel_pe2: str = "x",
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
    alternate_refocussing_phase = False
    disable_spoiler = False
    
    # Create a new sequence object
    seq = pp.Sequence(system)
    
    # check if channel labels are valid
    channel_ro, channel_pe1, channel_pe2 = seq_utils.validate_inputs(channel_ro, channel_pe1, channel_pe2)

    # map fov and n_enc according to channels    
    n_enc_ro, fov_ro, n_enc_pe1, fov_pe1, n_enc_pe2, fov_pe2 = seq_utils.calculate_enc_fov_order(
                                                                channel_ro, channel_pe1, n_enc, fov
                                                                )
    
    # derived and modifed parameters
    delta_k_ro = 1/fov_ro
    delta_k_pe1 = 1/fov_pe1
    delta_k_pe2 = 1/fov_pe2
    
    adc_duration = raster(n_enc_ro / ro_bandwidth, precision=system.grad_raster_time)   # sample everything on grad_raster_time
    adc_dwell = 1 / ro_bandwidth
    grad_amplitude = n_enc_ro*delta_k_ro/ adc_duration
    
    gradient_correction = raster(gradient_correction, precision=system.grad_raster_time)

    pe_traj = seq_utils.get_traj(
                n_enc_pe1 = n_enc_pe1,
                n_enc_pe2 = n_enc_pe2,
                etl = etl,
                trajectory = trajectory,
                )
    
    trains = seq_utils.get_trains(
                pe_traj = pe_traj,
                n_enc_pe1 = n_enc_pe1,
                n_enc_pe2 = n_enc_pe2,
                fov_pe1 = fov_pe1,
                fov_pe2 = fov_pe2,
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

    # Define readout gradient and prewinder
    grad_ro = pp.make_trapezoid(
        channel=channel_ro,
        system=system,
        amplitude=grad_amplitude,
        rise_time=ramp_duration,
        fall_time=ramp_duration,
        # Add gradient correction time and ADC correction time
        flat_time=adc_duration + 2 * gradient_correction # HH: why 2*gradient_correction?
    )

    # # Calculate readout prephaser without correction timess
    grad_ro_pre = pp.make_trapezoid(
        channel=channel_ro,
        system=system,
        area=grad_ro.area / 2 + delta_k_ro / 2 ,
        rise_time=ramp_duration,
        fall_time=ramp_duration,
        duration=pp.calc_duration(grad_ro) / 2,
    )
    ro_pre_duration = pp.calc_duration(grad_ro_pre) 
    
    adc = pp.make_adc(
        system=system,
        num_samples=n_enc_ro*ro_oversampling,
        duration=adc_duration,
        # Add gradient correction time and ADC correction time
        delay=2 * gradient_correction + grad_ro.rise_time
    )
    
    # ## Spoiler gradient on PE2 (used three times: before excitation (or after ADC), before refocusing, after refocusing) 
    area_pe2_sp = 4*pi/(2*pi*42.57*fov.z/n_enc.z) # unit area: mt/m*ms
    area_pe2_sp = area_pe2_sp*1e-6*system.gamma # unit area: 1/m
    area_pe2_sp = np.round(area_pe2_sp*1e3)/1e3
    if disable_spoiler:
        area_pe2_sp = 0 # Turn off spoiler
    
    grad_pe2_sp = pp.make_trapezoid(
        channel=channel_pe2, 
        area=area_pe2_sp,
        rise_time=200.e-6,
        fall_time=200.e-6,
        duration=.5e-3, 
        system=system
        )
    
    # Calculate delays
    # Note: RF dead-time is contained in RF delay
    # Delay duration center excitation (90 degree) and center refocussing (180 degree) RF pulse
    tau_1 = raster(val=echo_time/2 - rf_90.shape_dur/2 - rf_90.ringdown_time - rf_180.shape_dur/2 - rf_180.delay - ro_pre_duration - pp.calc_duration(grad_pe2_sp), precision=system.grad_raster_time)
    
    # Delay duration between center refocussing (180 degree) RF pulse and center readout
    tau_2 = raster(echo_time/2 - adc_duration/2 - rf_180.shape_dur/2 - rf_180.ringdown_time - 2*gradient_correction \
        - ramp_duration - ro_pre_duration - pp.calc_duration(grad_pe2_sp) + echo_shift, precision=system.grad_raster_time)

    # Delay duration between center readout and next center refocussing (180 degree) RF pulse 
    tau_3 = raster(echo_time/2 - adc_duration/2 - rf_180.shape_dur/2 - rf_180.delay  - ramp_duration - ro_pre_duration - pp.calc_duration(grad_pe2_sp) - echo_shift, precision=system.grad_raster_time)

    recommended_timing = seq_utils.get_esp_etl(tau_1=tau_1, tau_2=tau_2, tau_3=tau_3, echo_time=echo_time, T2=100, n_enc_pe1=n_enc_pe1)
    print(recommended_timing)
    
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

    for train in trains:
        # Reset phase of 180° refocussing pulse if alternate refocussing phase is enabled
        if alternate_refocussing_phase:
            rf_180.phase_offset = refocussing_phase
        
        if inversion_pulse:
            seq.add_block(rf_inversion)
            seq.add_block(pp.make_delay(raster(
                val=inversion_time - rf_duration,
                precision=system.grad_raster_time
            )))
        seq.add_block(rf_90)
        seq.add_block(grad_ro_pre)
        seq.add_block(pp.make_delay(tau_1))
        
        for echo in train:
            pe_1, pe_2 = echo
            if disable_pe:
                pe_1 = 0
                pe_2 = 0
            
            seq.add_block(grad_pe2_sp) 
            seq.add_block(rf_180)
            seq.add_block(grad_pe2_sp)
            seq.add_block(
                pp.make_trapezoid(
                    channel=channel_pe1,
                    area=-pe_1,
                    duration=ro_pre_duration,
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
                )
            )          

            seq.add_block(pp.make_delay(tau_2))

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

            seq.add_block(pp.make_delay(tau_3))
            
            # Alternate refocussing phase by 180° each time
            if alternate_refocussing_phase:
                rf_180.phase_offset = (rf_180.phase_offset + pi) % (2 * pi)

        seq.add_block(grad_pe2_sp) # add spoiler after last 180 pulse in echo train

        # recalculate TR each train because train length is not guaranteed to be constant
        tr_delay = raster(repetition_time - echo_time * len(train) - adc_duration / 2 - ro_pre_duration \
            - tau_3 - rf_90.delay - rf_duration / 2 - ramp_duration - pp.calc_duration(grad_pe2_sp), precision=system.grad_raster_time)

        if inversion_pulse:
            tr_delay -= inversion_time
        
        seq.add_block(pp.make_delay(tr_delay))

    # Calculate some sequence measures
    train_duration_tr = (seq.duration()[0]) / len(trains)
    train_duration = train_duration_tr - tr_delay
    
    # Map trajectory to PE1 + PE2 samples in format required by sort_kspace()
    k_traj_adc = seq.calculate_kspacePP()[0]
    acq_pos = seq_utils.calculate_acq_pos(k_traj_adc,n_enc_ro,channel_pe1,fov_pe1,channel_pe2,fov_pe2)

    # # Add measures to sequence definition
    # seq.set_definition("n_total_trains", len(trains))
    # seq.set_definition("train_duration", train_duration)
    # seq.set_definition("train_duration_tr", train_duration_tr)
    # seq.set_definition("tr_delay", tr_delay)

    return (seq, acq_pos, [n_enc_ro, n_enc_pe1, n_enc_pe2])


def sort_kspace(raw_data: np.ndarray, trajectory: np.ndarray, kdims: list) -> np.ndarray:
    """
    Sort acquired k-space lines.

    Parameters
    ----------
    kspace
        Acquired k-space data in the format (averages, coils, pe, ro)
    trajectory
        k-Space trajectory returned by TSE constructor with dimension (pe, 2)
    dim
        dimensions of kspace
    """
    n_avg, n_coil, _, _ = raw_data.shape
    ksp = np.zeros((n_avg, n_coil, kdims[2], kdims[1], kdims[0]), dtype=complex)
    for idx, kpt in enumerate(trajectory):
        ksp[..., kpt[1], kpt[0], :] = raw_data[:, :, idx, :]
    return ksp
