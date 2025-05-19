"""Constructor for 3D TSE Imaging sequence."""

import warnings
from math import pi

import numpy as np
import pypulseq as pp
from pypulseq.opts import Opts

from packages import seq_utils
from packages.seq_utils import Dimensions
from packages.seq_utils import Channels
from packages.seq_utils import raster
from packages.mr_systems import low_field as system



def constructor(
    te: float = 15e-3,
    tr: float = 600e-3,
    etl: int = 7,
    n_dummies: int = 0,
    rf_ex_duration: float = 400e-6,
    rf_ref_duration: float = 400e-6,
    g_ro_correction: float = 0.0,
    adc_bandwidth: float = 20e3,
    adc_oversampling: int = 1,
    input_fov: Dimensions | None = None,
    input_enc: Dimensions | None = None,
    trajectory: seq_utils.Trajectory = seq_utils.Trajectory.OUTIN,
    rf_ex_flip_angle: float = pi / 2,
    rf_ex_phase: float = 0.0,
    rf_ref_flip_angle: float = pi,
    rf_ref_phase: float = pi / 2,
    channels: Channels | None = None,
    system: Opts = system,
) -> tuple[pp.Sequence, list]:
    """Construct 3D turbo spin echo sequence.

    Parameters
    ----------
    te
        Time between center of 90 degree pulse and center of ADC in s.
    tr
        Time between two subsequent 90 degree pulses (echo trains) in s.
    etl
        Echo train length.
    n_dummies
        Number of dummy shots to acquire.
    rf_ex_duration
        Duration of the excitation RF pulse in s.
    rf_ref_duration
        Duration of the refocussin RF pulse in s.
    g_ro_correction
        Time to center ADC event in s.
    adc_bandwidth
        Readout bandwidth in Hz.
    adc_oversampling
        Factor to oversample the ADC.
    input_fov
        Field of view per dimension in m. Set None for default values.
    input_enc
        Number of encoding steps per dimension. Set None for default values.
    trajectory
        The k-space trajectory.
    rf_ex_flip_angle
        Flip angle of the excitation pulse in radians.
    rf_ex_flip_phase
        Phase of the excitation pulse in radians.
    rf_ref_flip_angle
        Flip angle of the refocussing pulse in radians.
    rf_ref_flip_phase
        Phase of the refocussing pulse in radians.
    channels
        Channel configuration. Sets the readout, phase encoding 1 and phase encoding 2 channels.
        Set None for default values.
    system
        System object containing hardware parameters.

    Returns
    -------
    tuple[pp.Sequence, list]: A tuple containing a PyPulseq sequence and a list:
        - seq: PyPulseq sequence.
        - header: ISMRMD header.
    """
    # Sequence options
    # Toggle the use of constant PE amplitudes (avoid artifact due to GPA nonlinearity) or
    # constant duration (classical implementation) for all phase encoding gradients
    DISABLE_PE = False  # Debugging: disable phase encoding in both directions
    DISABLE_SPOILER = False  # Debugging: disable spoiler

    # Set default values if not provided
    if input_fov is None:
        input_fov = Dimensions(x=192e-3, y=192e-3, z=192e-3)
    if input_enc is None:
        input_enc = Dimensions(x=64, y=64, z=64)
    if channels is None:
        channels = Channels(ro='x', pe1='y', pe2='z')

    # Create a new sequence object
    seq = pp.Sequence(system)

    # Map fov and n_enc with channels
    n_enc, fov = seq_utils.map_fov_enc(channels, input_fov, input_enc)

    # Derive ADC timing and RO gradient amplitude
    adc_dwell = 1 / adc_bandwidth
    adc_dwell = raster(adc_dwell, precision=system.grad_raster_time)
    adc_duration = adc_dwell * n_enc['ro']
    adc_bandwidth_raster = 1 / adc_dwell

    g_ro_amplitude = n_enc['ro'] / fov['ro'] / adc_duration
    g_ro_correction = raster(g_ro_correction, precision=system.grad_raster_time)

    # Create readout gradient to use timing for other gradients
    g_ro = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        amplitude=g_ro_amplitude,
        flat_time=adc_duration + g_ro_correction,  # Add gradient correction time
    )

    # Define phase encoding gradients
    _k1_max = 1 / fov['pe1'] * n_enc['pe1'] / 2  # maximum kspace value in pe1
    _k2_max = 1 / fov['pe2'] * n_enc['pe2'] / 2  # maximum kspace value in pe2
    k_max = max(_k1_max, _k2_max)  # maximum kspace value in both directions

    # Get duration for largest phase encode gradient
    g_pe_max = pp.make_trapezoid(channel=channels.pe1, area=k_max, system=system)
    # Use this for all phase encode gradients
    g_pe_duration = pp.calc_duration(g_pe_max)

    g_pe_delay = pp.make_delay(g_pe_duration)

    # Generate k-space trajectory
    pe_traj = seq_utils.get_traj(
        n_enc=n_enc,
        etl=etl,
        trajectory=trajectory,
    )

    # Sort k-space points into echo trains
    trains, trains_pos = seq_utils.get_trains(
        pe_traj=pe_traj,
        n_enc=n_enc,
        fov=fov,
        etl=etl,
        trajectory=trajectory,
    )

    # Create excitation pulse
    rf_ex = pp.make_block_pulse(
        system=system,
        flip_angle=rf_ex_flip_angle,
        phase_offset=rf_ex_phase,
        duration=rf_ex_duration,
        delay=system.rf_dead_time,
        use='excitation',
    )

    # Create refocussing pulse
    rf_ref = pp.make_block_pulse(
        system=system,
        flip_angle=rf_ref_flip_angle,
        phase_offset=rf_ref_phase,
        duration=rf_ref_duration,
        delay=system.rf_dead_time,
        use='refocusing',
    )

    # Create spoiler gradient on readout axis before and after 180°-SE-refocusing pulse
    g_ro_spoil = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        area=g_ro.area,
    )

    # Disable spoiler if needed
    if DISABLE_SPOILER:
        warnings.warn('!!!Spoiler are disabled!!!', stacklevel=2)
        g_ro_spoil.amplitude = 0  # avoid spoiler creation during sequence building
        g_ro_spoil.area = 0  # necessary for prewinder definition

    # Create readout prewinder accounting for readout and spoiler area
    g_ro_prew = pp.make_trapezoid(
        channel=channels.ro,
        system=system,
        area=g_ro.area / 2 + g_ro_spoil.area,
    )
    g_ro_prew_duration = pp.calc_duration(g_ro_prew)

    # Create ADC event
    adc = pp.make_adc(
        system=system,
        num_samples=n_enc['ro'] * adc_oversampling,
        duration=adc_duration,
        delay=g_ro_correction + g_ro.rise_time,  # Add gradient correction time
    )

    # Calculate delays
    # Delay duration between center excitation (90 degree) and center refocussing (180 degree) RF pulse
    delay_rf_ex_ref = raster(
        val=te / 2 - rf_ex_duration / 2 - rf_ex.ringdown_time - rf_ref_duration / 2 - rf_ref.delay - g_ro_prew_duration,
        precision=system.grad_raster_time,
    )

    # Delay duration between center refocussing (180 degree) RF pulse and center readout
    delay_rf_ref_adc = raster(
        te / 2
        - adc_duration / 2
        - g_ro_correction
        - rf_ref_duration / 2
        - rf_ref.ringdown_time
        - g_ro.rise_time
        - pp.calc_duration(g_ro_spoil, g_pe_delay),
        precision=system.grad_raster_time,
    )

    # Delay duration between center ADC and next center refocussing (180 degree) RF pulse
    delay_adc_rf_ref = raster(
        te / 2
        - adc_duration / 2
        - rf_ref_duration / 2
        - rf_ref.delay
        - g_ro.rise_time
        - pp.calc_duration(g_ro_spoil, g_pe_delay),
        precision=system.grad_raster_time,
    )

    # Add dummys by adding NaN entries at beginning of trains_pos
    for _ in range(n_dummies):
        trains_pos.insert(0, np.full((etl, 2), np.nan).tolist())
        trains.insert(0, np.zeros((etl, 2), dtype=int))

    # Add imaging blocks to sequence
    for train, position in zip(trains, trains_pos, strict=True):
        _start_time = sum(seq.block_durations.values())

        # add 90° excitation rf pulse
        seq.add_block(rf_ex)
        # add readout prewinder
        seq.add_block(g_ro_prew)
        # add delay to refocussing rf pulse
        seq.add_block(pp.make_delay(delay_rf_ex_ref))

        for echo, pe_indices in zip(train, position, strict=True):
            # no labels for dummy scans
            if not np.isnan(pe_indices).any():
                # Cast index values from int32 to int, otherwise make_label function complains
                label_pe1 = pp.make_label(type='SET', label='LIN', value=int(pe_indices[0]))
                label_pe2 = pp.make_label(type='SET', label='PAR', value=int(pe_indices[1]))

            # get phase encoding values
            pe_1, pe_2 = echo

            # Disable phase encoding if needed
            if DISABLE_PE:
                warnings.warn('!!!Phase encoding is disabled!!!', stacklevel=2)
                pe_1 = 0
                pe_2 = 0

            # create phase encoding gradients and rewinder for PE1
            g_pe1_prew = pp.make_trapezoid(
                channel=channels.pe1,
                area=pe_1,
                duration=g_pe_duration,
                system=system,
            )
            g_pe1_rew = pp.scale_grad(g_pe1_prew, -1)

            # create phase encoding gradients and rewinder for PE2
            g_pe2_prew = pp.make_trapezoid(
                channel=channels.pe2,
                area=pe_2,
                duration=g_pe_duration,
                system=system,
            )
            g_pe2_rew = pp.scale_grad(g_pe2_prew, -1)

            # add refocussing rf pulse
            seq.add_block(rf_ref)
            # add phase encoding gradients and first spoiler gradient on readout axis
            seq.add_block(
                g_pe1_prew, g_pe2_prew, g_ro_spoil, g_pe_delay
            )  # ensures that event is at least g_pe_duration long
            # add delay to readout
            seq.add_block(pp.make_delay(delay_rf_ref_adc))

            # handle dummy scans
            if np.isnan(pe_indices).any():
                # only readout gradient for dummy scans
                seq.add_block(g_ro)
            else:
                # add readout gradient and ADC
                seq.add_block(g_ro, adc, label_pe1, label_pe2)

            # add phase encode rewinder and second spoiler gradient on readout axis
            seq.add_block(
                g_pe1_rew, g_pe2_rew, g_ro_spoil, g_pe_delay
            )  # ensures that event is at least g_pe_duration long
            # add delay to next refocussing rf pulse
            seq.add_block(pp.make_delay(delay_adc_rf_ref))

        # calculate TR delay
        _duration_pe_step = sum(seq.block_durations.values()) - _start_time
        delay_tr = tr - _duration_pe_step
        delay_tr = raster(delay_tr, precision=system.grad_raster_time)

        if delay_tr < 0:
            raise ValueError('TR too short')
        # add delay to next TR
        seq.add_block(pp.make_delay(delay_tr))

    print(f'TR fill: {1000 * delay_tr} ms')

    # Check labels
    labels = seq.evaluate_labels(evolution='adc')
    # handle dummy scans: remove nan from trains_pos
    trains_pos = [x for x in trains_pos if not np.isnan(x).any()]
    acq_pos = np.concatenate(trains_pos).T

    if not np.array_equal(labels['LIN'], acq_pos[0, :]):
        raise ValueError("LIN labels don't match actual acquisition positions.")
    if not np.array_equal(labels['PAR'], acq_pos[1, :]):
        raise ValueError("PAR labels don't match actual acquisition positions.")

    # calculate train duration
    train_duration = seq.duration()[0] / len(trains) - delay_tr

    # Create ISMRMD header
    header = []

    # write all required parameters in the seq-file header/definitions dictionary
    seq.set_definition('name', 'tse_3d')
    seq.set_definition('TE', te)
    seq.set_definition('TR', tr)
    seq.set_definition('FOV', (fov['ro'], fov['pe1'], fov['pe2']))
    seq.set_definition('encoding_dim', (n_enc['ro'], n_enc['pe1'], n_enc['pe2']))
    seq.set_definition('channel_order', (channels.ro, channels.pe1, channels.pe2))
    seq.set_definition('echo_train_length', int(etl))
    seq.set_definition('train_duration', int(train_duration))
    seq.set_definition('n_total_trains', len(trains))
    seq.set_definition('ro_bandwidth', adc_bandwidth_raster)
    seq.set_definition('ro_oversampling', int(adc_oversampling))
    seq.set_definition('k_space_encoding1', int(n_enc['pe1']))
    seq.set_definition('k_space_encoding2', int(n_enc['pe2']))    

    return (seq, header)
