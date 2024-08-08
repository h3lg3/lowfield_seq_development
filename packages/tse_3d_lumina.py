"""Constructor for 3D TSE Imaging sequence.

TODO: add sampling patterns (elliptical masks, partial fourier, CS)
TODO: add optional inversion pulse
TODO: add optional variable refocussing pulses (pass list rather than float)
TODO: move trajectory calculation to seperate file to sharew with other imaging experiments (needed?)
TODO: Design goal: Allow for minimal TE/Echo spacing for maximal ETL (acceleration)?

"""
from math import pi
import math
import numpy as np
import pypulseq as pp

from console.interfaces.interface_acquisition_parameter import Dimensions
from console.utilities.sequences.system_settings import raster
from packages import def_trajectory

# low field system
# from console.utilities.sequences.system_settings import system
# Siemens system
from packages.siemens_system import system

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
    trajectory: def_trajectory.Trajectory = def_trajectory.Trajectory.OUTIN,
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
    # system.rf_ringdown_time = 0 # used for lowfield system?
    
    system.max_grad=24*1e-3*system.gamma
    system.max_slew=100*system.gamma
    
    seq = pp.Sequence(system)
    seq.set_definition("Name", "tse_3d_lumina")

    # check if channel labels are valid
    channel_valid = True
    if len(channel_ro) > 1 or len(channel_ro) == 0:
        channel_valid = False
        print("Invalid readout channel: %s" % (channel_ro))
    if len(channel_pe1) > 1 or len(channel_pe1) == 0:
        channel_valid = False
        print("Invalid pe1 channel: %s" % (channel_pe1))
    if len(channel_pe2) > 1 or len(channel_pe2) == 0:
        channel_valid = False
        print("Invalid pe2 channel: %s" % (channel_pe2))

    channel_ro = channel_ro.lower()
    channel_pe1 = channel_pe1.lower()
    channel_pe2 = channel_pe2.lower()  # set all channels to lower case

    if channel_ro not in ("x", "y", "z") or channel_pe1 not in ("x", "y", "z") or channel_pe2 not in ("x", "y", "z"):
        channel_valid = False
        print("Invalid axis orientation")
    if channel_ro == channel_pe1 or channel_ro == channel_pe2 or channel_pe1 == channel_pe2:
        channel_valid = False
        print("Error, multiple channels have the same gradient")
        print("Readout channel: %s, pe1 channel: %s, pe2 channel: %s" % (channel_ro, channel_pe1, channel_pe2))
    if not channel_valid:
        print("Defaulting to readout in y, pe1 in z, pe2 in x")
        channel_ro = "y"
        channel_pe1 = "z"
        channel_pe2 = "x"

    if (channel_ro == "x"):
        n_enc_ro = n_enc.x
        fov_ro = fov.x
        if channel_pe1 == "y":
            n_enc_pe1 = n_enc.y
            fov_pe1 = fov.y
            n_enc_pe2 = n_enc.z
            fov_pe2 = fov.z
        else:
            n_enc_pe1 = n_enc.z
            fov_pe1 = fov.z
            n_enc_pe2 = n_enc.y
            fov_pe2 = fov.y
    elif (channel_ro == "y"):
        n_enc_ro = n_enc.y
        fov_ro = fov.y
        if channel_pe1 == "x":
            n_enc_pe1 = n_enc.x
            fov_pe1 = fov.x
            n_enc_pe2 = n_enc.z
            fov_pe2 = fov.z
        else:
            n_enc_pe1 = n_enc.z
            fov_pe1 = fov.z
            n_enc_pe2 = n_enc.x
            fov_pe2 = fov.x
    else:
        n_enc_ro = n_enc.z
        fov_ro = fov.z
        if channel_pe1 == "y":
            n_enc_pe1 = n_enc.y
            fov_pe1 = fov.y
            n_enc_pe2 = n_enc.x
            fov_pe2 = fov.x
        else:
            n_enc_pe1 = n_enc.x
            fov_pe1 = fov.x
            n_enc_pe2 = n_enc.y
            fov_pe2 = fov.y
    
    pe_traj = def_trajectory.constructor(
                            n_enc_pe1 = n_enc_pe1,
                            n_enc_pe2 = n_enc_pe2,
                            etl = etl,
                            trajectory = trajectory,
                            )[0]

    # Divide all PE steps into echo trains
    if trajectory.name == 'SYMMETRIC':
        num_trains = int(np.ceil(n_enc_pe1 / etl))
        temp = [pe_traj[k::num_trains] for k in range(num_trains)]
        trains = []
        for k in np.arange(n_enc_pe2):
            for i in np.arange(num_trains):
                for j in np.arange(etl):
                    trains.append([temp[i][j]/fov_pe1, k/fov_pe2])
                    
        trains = [trains[k*etl:(k+1)*etl] for k in range(num_trains*n_enc_pe2)]      
    else:
        # # calculate the required gradient area for each k-point
        pe_traj[:, 0] /= fov_pe1
        pe_traj[:, 1] /= fov_pe2
        num_trains = int(np.ceil(pe_traj.shape[0] / etl))
        trains = [pe_traj[k::num_trains, :] for k in range(num_trains)]

    # Create a list with the kspace location of every line of kspace acquired, in the order it is acquired
    # trains_pos = [pe_order[k::num_trains, :] for k in range(num_trains)]
    # acq_pos = []
    # for train_pos in trains_pos:
    #     acq_pos.extend(train_pos)
    acq_pos = []

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
    
    adc = pp.make_adc(
        system=system,
        num_samples=int(n_enc_ro*ro_oversampling),
        duration=raster(val=adc_duration, precision=system.adc_raster_time),
        # Add gradient correction time and ADC correction time
        delay=raster(val=2 * gradient_correction + grad_ro.rise_time, precision=system.adc_raster_time)
    )

    # ## Spoiler gradient on PE2 (used three times: before excitation (or after ADC), before refocusing, after refocusing) 
    area_pe2_sp = 4*pi/(2*pi*42.57*fov.z/n_enc.z) # unit area: mt/m*ms
    area_pe2_sp = area_pe2_sp*1e-6*system.gamma # unit area: 1/m
    grad_pe2_sp = pp.make_trapezoid(
        channel=channel_pe2, area=area_pe2_sp, system=system
        )
    
    # Calculate delays
    # Note: RF dead-time is contained in RF delay
    # Delay duration center excitation (90 degree) and center refocussing (180 degree) RF pulse
    tau_1 = echo_time/2 - rf_90.shape_dur/2 - rf_90.ringdown_time - rf_180.shape_dur/2 - rf_180.delay - ro_pre_duration - pp.calc_duration(grad_pe2_sp)
    
    # Delay duration between center refocussing (180 degree) RF pulse and center readout
    tau_2 = echo_time/2 - adc_duration/2 - rf_180.shape_dur/2 - rf_180.ringdown_time - 2*gradient_correction \
        - ramp_duration - pp.calc_duration(grad_ro_pre) - pp.calc_duration(grad_pe2_sp) + echo_shift

    # Delay duration between center readout and next center refocussing (180 degree) RF pulse 
    tau_3 = echo_time/2 - adc_duration/2 - rf_180.shape_dur/2 - rf_180.delay  - ramp_duration - pp.calc_duration(grad_ro_pre) - pp.calc_duration(grad_pe2_sp) - echo_shift 

    min_esp = -(min(tau_1, tau_2, tau_3) - echo_time/2)*2   # minimum echo spacing to accomodate gradients
    min_esp = np.ceil(min_esp*1e3)*1e-3                     # round to 1 ms -> new echo time
    T2 = 80
    max_sampling = -math.log(0.1)*T2*1e-3                   # sampling duration [ms] till signal drops to 20%
    max_etl = np.floor(max_sampling/min_esp)                # maximum numbers of 180° echoes
    
    # ETL that is multiple of n_enc_pe1 and closest to max_etl
    cc = [(i, n_enc_pe1%i, np.abs(max_etl-i)) for i in np.arange(50)]
    # Filter the list to find tuples where the second element is zero
    filtered_cc = [item for item in cc if item[1] == 0]
    # Find the tuple with the minimal value in the third element
    result = min(filtered_cc, key=lambda x: x[2])           # closest ETL value that is multiple of n_enc_pe1 -> new ETL

    for dummy in range(dummies):
        if inversion_pulse:
            seq.add_block(rf_inversion)
            seq.add_block(pp.make_delay(raster(val=inversion_time - rf_duration, precision=system.grad_raster_time)))
        seq.add_block(rf_90)
        seq.add_block(pp.make_delay(raster(val=echo_time / 2 - rf_duration, precision=system.grad_raster_time)))
        for idx in range(etl):
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
        if inversion_pulse:
            seq.add_block(rf_inversion)
            seq.add_block(pp.make_delay(raster(
                val=inversion_time - rf_duration,
                precision=system.grad_raster_time
            )))
        seq.add_block(rf_90)
        seq.add_block(grad_ro_pre)
        seq.add_block(pp.make_delay(raster(val=tau_1, precision=system.grad_raster_time)))

        for echo in train:
            pe_1, pe_2 = echo
            
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

        seq.add_block(grad_pe2_sp) # add spoiler after last 180 pulse in echo train

        # recalculate TR each train because train length is not guaranteed to be constant
        tr_delay = repetition_time - echo_time * len(train) - adc_duration / 2 - ro_pre_duration \
            - tau_3 - rf_90.delay - rf_duration / 2 - ramp_duration - pp.calc_duration(grad_pe2_sp)

        if inversion_pulse:
            tr_delay -= inversion_time

        seq.add_block(pp.make_delay(raster(
            val=tr_delay,
            precision=system.block_duration_raster
        )))

    # Calculate some sequence measures
    train_duration_tr = (seq.duration()[0]) / len(trains)
    train_duration = train_duration_tr - tr_delay

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