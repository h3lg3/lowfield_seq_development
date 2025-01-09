"""Constructor for 2D SE T1 Mapping Imaging sequence.
Source: https://github.com/mritogether/ESMRMB2024_Hardware_to_Map/blob/main/02_sequence_design_for_mapping/notebooks/s30_2D_IR_SE_T1mapping.ipynb
"""
# TODO: Make spoiler longer if there is enough time between excitation pulse and 180 pulse and after readout, similar to what was done for gz_s, gx_s, gy_s
# TODO: Add multi-slice capability
# TODO: Add multi-echo capability: Does it affect T1 mapping if multiple echoes are used?
# TODO: Add multi-TI readouts: Is it possible to acquire the same k-space line for all inversion times (TI) after one inversion pulse? 
# TODO: RF gradients are automatically set to z-channel. Add option to set the channel for RF gradients.

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


def constructor(
    echo_time: float = 15e-3,   
    repetition_time: float = 5000e-3,
    TI: list = [50e-3, 100e-3, 500e-3, 1500e-3, 4500e-3],
    etl: int = 1,               
    slice_thickness: int = 8e-3,  
    ro_bandwidth: float = 20e3,
    ro_oversampling: int = 5,
    input_fov: Dimensions = Dimensions(x=220e-3, y=220e-3, z=225e-3),
    input_enc: Dimensions = Dimensions(x=70, y=70, z=49),
    excitation_angle: float = pi / 2,
    channels: Channels = Channels(ro="x", pe1="y", pe2="z"),    # channel pe2 is used for slice selection
    system:Opts = default_system,
) -> tuple[pp.Sequence, list, list]:
    """Construct 3D turbo spin echo sequence.

    Parameters
    ----------
    echo_time, optional
        Time constant between center of 90 degree pulse and center of ADC, by default 15e-3
    repetition_time, optional
        Time constant between two subsequent 90 degree pulses (echo trains), by default 600e-3
    ro_bandwidth, optional
        Readout bandwidth in Hz, by default 20e3
    fov, optional
        Field of view per dimension, by default default_fov
    n_enc, optional
        Number of encoding steps per dimension, by default default_encoding = Dimensions(x=70, y=70, z=49).
        If an encoding dimension is set to 1, the TSE sequence becomes a 2D sequence.
    excitation_angle, excitation_phase, optional
        set the flip angle and phase of the excitation pulse in radians, defaults to 90 degree flip angle, 0 phase
    channel_ro, channel_pe1, channel_pe2, optional
        set the readout, phase1 and phase2 encoding directions, default to y, z and x.

    Returns
    -------
        Pulseq sequence and a list which describes the trajectory
    """

    # increase gradient duration to relax performance
    SILENT_MODE = False

    # Create a new sequence object
    seq = pp.Sequence(system)

    # map fov and n_enc with channels  
    n_enc, fov = seq_utils.map_fov_enc(channels, input_fov, input_enc)

    # Define Geometry
    Nx, Ny, Nz = [n_enc["ro"], n_enc["pe1"], n_enc["pe2"]]

    TE  = echo_time   # Echo time [s]
    TR  = repetition_time    # Repetition time [s]

    adc_dwell_time = raster(1 / ro_bandwidth, precision=system.grad_raster_time) # sample everything on grad_raster_time
    adc_duration = n_enc['ro'] * adc_dwell_time    

    # Inversion Pulse
    rf180inv, gz180inv, _ = pp.make_adiabatic_pulse(
        pulse_type = 'hypsec',
        duration = 8e-3, 
        beta = 800,
        slice_thickness = slice_thickness, 
        mu = 4.9,
        return_gz = True, 
        system = seq.system, 
        use = "inversion"
        )    
    gz180inv.channel = channels.pe2 # change channel to pe2 (slice selection)
    
    rf_inv_dur = pp.calc_duration(gz180inv)

    # Generate sequence
    seq = pp.Sequence(system)

    # Sequence to calculate delay
    seqa = pp.Sequence(system)

    # Export delays for TR and TI
    _, _, _, time_remove_TI, remove_delayTR = SE_module(seqa, fov = fov, Nx = Nx, Nz = Nz, Ny = Ny, TE = TE,
                                TR = TR, TI = TI[0], rf_inv_dur = rf_inv_dur, ky_i = 0, 
                                ETL = etl, adc_duration = adc_duration, system = system,
                                channels = channels, silent_mode = SILENT_MODE)

    time_remove_TI = np.ceil(time_remove_TI/system.grad_raster_time)*system.grad_raster_time 
    remove_delayTR = np.ceil(remove_delayTR/system.grad_raster_time)*system.grad_raster_time 

    delay_TR = TR - remove_delayTR
    TId = TI - time_remove_TI

    for idx, nTI in enumerate(TId):
        for ky_i in range(Ny):
            # Cast index values from int32 to int, otherwise make_label function complains
            label_PE = pp.make_label(type = "SET", label = "LIN", value = idx)
            label_TI = pp.make_label(type = "SET", label = "REP", value = Ny-1 - ky_i) # PE goes from max (31) to min (-32), therefore order needs to be reversed here, such that the min line is at 0 index of array

            seq.add_block(pp.make_delay(round(delay_TR - nTI, 9)))
            seq.add_block(rf180inv, gz180inv)
            seq.add_block(pp.make_delay(nTI), label_PE, label_TI)
            ##########################################################################################
            # SE module
            seq, TR, _,_,_ = SE_module(seq, fov = fov, Nx = Nx, Nz = Nz, Ny = Ny, TE = TE,
                                TR = TR, TI = TI[0], rf_inv_dur = rf_inv_dur, ky_i = ky_i, 
                                ETL = etl, adc_duration = adc_duration, system = system,
                                channels = channels, silent_mode = SILENT_MODE)
            ##########################################################################################


    header = seq_utils.create_ismrmd_header(
                n_enc = n_enc,
                fov = fov,
                system = system
                )

    # write all required parameters in the seq-file definitions.
    write_seq_definitions(
        seq = seq,
        fov = (fov["ro"], fov["pe1"], fov["pe2"]),
        name = "se_t1_mapping",
        alpha = excitation_angle,
        slice_thickness = slice_thickness,
        Nx = n_enc['ro'],
        Ny = n_enc['ro'],
        sampling_scheme = 'cartesian',
        Nz = n_enc['pe2'],
        TE = echo_time,
        TR = repetition_time,
        ro_bandwidth = ro_bandwidth,
        ro_oversampling = ro_oversampling
    )

    return (seq, header)

def SE_module(seq, fov = 200e-3, Nx = 128, Nz = 1, Ny = 128, TE = 15e-3,
              TR = 5000e-3, TI = 50e-3, rf_inv_dur = 1e-3, ky_i = 0, slice_thickness = 5e-3, ETL = 1, 
              adc_duration = 6.4e-3, system:Opts = default_system,
              channels: Channels = Channels(ro="x", pe1="y", pe2="z"), silent_mode = True):
    """
    Creates a Spin-echo module (SE) sequence and adds it
    to the seq object

    Parameters
    ----------
    seq: object
        Sequence object
    fov : float, optional
        Field-of-view [m]. Default is 0.2 m.
    Nx, Ny, Nz : float, optional
        Ny - Number of phase encoding lines. Default is 128.
    TE: Echo time [s]
    TR: Repetition Time [s]
    ky_i: line in k-space
   Returns
    -------
    seq : SimpleNamespace
        Seq object with bSSFP readout
    TR : float
        Repetition time
    Ny : float
    time_remove_TI: float
                    time for TI delay calculation
    remove_delayTR: float
                    time for TR delay calculation
    """
    readout_time =  adc_duration  # ADC duration
    TE = TE * np.arange(1, ETL+1) # [s]
    
    # Readout gradients
    delta_kx = 1 / fov["ro"]
    delta_ky = 1 / fov["pe1"]
    kx_width = (Nx) * delta_kx

    gx = pp.make_trapezoid(channel = channels.ro, system = system, flat_area = kx_width,
                    flat_time = readout_time)
    adc = pp.make_adc(num_samples = Nx, duration = readout_time, delay = gx.rise_time)

    gx_pre = pp.make_trapezoid(channel = channels.ro, system = system, area = (gx.area) / 2)
    gx_post = pp.make_trapezoid(channel = channels.ro, system = system, area = 3 * gx.area / 2)

    # Prephase and rephase
    phase_areas = np.linspace(-0.5*Ny * delta_ky, (0.5*Ny-1) * delta_ky, Ny) # scale from -Ny/2 to Ny/2-1
    gy_pe_max = pp.make_trapezoid(channel = channels.pe1, system = system, area = np.max(np.abs(phase_areas)))
    pe_duration = pp.calc_duration(gy_pe_max)    

    rf_flip = 90
    rf_offset = 0

    flip90 = round(rf_flip * pi / 180, 3)
    flip180 = 180 * pi / 180
    rf90, gz90, gz_reph = pp.make_sinc_pulse(flip_angle = flip90, 
                                             system = system, 
                                             duration = 2.5e-3,
                                             slice_thickness = slice_thickness, 
                                             apodization = 0.5,
                                             phase_offset = pi / 2,
                                             time_bw_product = 4, 
                                             return_gz = True, 
                                             use = "excitation")
    gz90.channel = channels.pe2 # change channel to pe2 (slice selection)
    gz_reph.channel = channels.pe2 # change channel to pe2 (slice selection)

    # increase duration to relax the gradient
    if silent_mode:
        gz_reph = pp.make_trapezoid(channel = gz_reph.channel, system = system, area = gz_reph.area, duration = 1.5*pp.calc_duration(gz_reph))

    rf180, gz180, _ = pp.make_sinc_pulse(flip_angle = flip180, 
                                         system = system,
                                         duration = 2.5e-3,
                                         slice_thickness = 1.25 * slice_thickness,
                                         apodization = 0.5,
                                         time_bw_product = 4, 
                                         phase_offset = 0,
                                         return_gz = True, 
                                         use = "refocusing")
    gz180.channel = channels.pe2 # change channel to pe2 (slice selection)
    
    # Spoiler
    gss_spoil_A = 4 / slice_thickness
    gss_spoil = pp.make_trapezoid(channel = channels.pe2, system = system, area = gss_spoil_A)
    # increase duration to relax the gradient
    if silent_mode:
        gss_spoil = pp.make_trapezoid(channel = channels.pe2, system = system, area = gss_spoil_A, duration = 2*pp.calc_duration(gss_spoil))

    gss_times = np.cumsum([0, 
                           gss_spoil.rise_time, 
                           gss_spoil.flat_time,
                           gss_spoil.fall_time, 
                           gz180.flat_time,
                           gss_spoil.rise_time,
                           gss_spoil.flat_time, 
                           gss_spoil.fall_time])

    gss_amps = np.array([0, 
                         gss_spoil.amplitude, 
                         gss_spoil.amplitude, 
                         gz180.amplitude, 
                         gz180.amplitude, 
                         gss_spoil.amplitude, 
                         gss_spoil.amplitude, 
                         0])

    gss_spoil_add = pp.make_extended_trapezoid(channel = channels.pe2, amplitudes = gss_amps, times = gss_times, system = system)
    gss_spoil_duration = pp.calc_duration(gss_spoil_add)

    # End of TR spoiler
    temp = pp.make_trapezoid(channel = channels.pe2, system=system, area=2 / slice_thickness)
    gz_spoil = pp.make_trapezoid(channel = channels.pe2, system=system, area=2 / slice_thickness, duration = 3*pp.calc_duration(temp))

    # gz_spoil = pp.make_trapezoid(channel = channels.pe2, system=system, area=2 / slice_thickness)
    gz_spoil.id = seq.register_grad_event(gz_spoil)

    # Spoiler on all axes
    g_spoil_A = - 4/slice_thickness
    gz_s = pp.make_trapezoid(channel = channels.pe2, system = seq.system, area = g_spoil_A)
    gx_s = pp.make_trapezoid(channel = channels.ro, system = seq.system, area = g_spoil_A)
    gy_s = pp.make_trapezoid(channel = channels.pe1, system = seq.system, area = g_spoil_A)

    spoiler_duration = pp.calc_duration(gz_s)

    # TI delay
    time_remove_TI = pp.calc_duration(gx_s, gy_s, gz_s) + pp.calc_duration(gz90)/2 + rf_inv_dur/2
    dTI = TI - time_remove_TI

    # Make spoiler longer if there is enough time between inversion pulse and excitation pulse
    if silent_mode:
        if 2*spoiler_duration < dTI:
                gz_s = pp.make_trapezoid(channel = channels.pe2, system = seq.system, area = g_spoil_A, duration = 2*spoiler_duration)
                gx_s = pp.make_trapezoid(channel = channels.ro, system = seq.system, area = g_spoil_A, duration = 2*spoiler_duration)
                gy_s = pp.make_trapezoid(channel = channels.pe1, system = seq.system, area = g_spoil_A, duration = 2*spoiler_duration)
                time_remove_TI = pp.calc_duration(gx_s, gy_s, gz_s) + pp.calc_duration(gz90)/2 + rf_inv_dur/2

    # Phase encoding
    gy_pre = pp.make_trapezoid(channel = channels.pe1, system = system,
                                 area = -phase_areas[-ky_i - 1], duration = pe_duration)
    gy_post = pp.make_trapezoid(channel = channels.pe1, system = system,
                                  area = -phase_areas[-ky_i - 1], duration = pe_duration)

    # Calculate timing for delays
    rf180.delay = gss_spoil.rise_time + gss_spoil.flat_time + gss_spoil.fall_time

    # Echo time (TE) and repetition time (TR)
    pre180d = (TE[0] / 2
               - pp.calc_duration(gz90) / 2
               - pp.calc_duration(gx_pre, gy_pre, gz_reph)
               - gss_spoil_duration / 2
                )
    post180d = (TE[0] / 2
                - gss_spoil_duration / 2
                - (gx.rise_time + readout_time / 2)
                )
    pre180delay = pp.make_delay(pre180d)
    post180delay= pp.make_delay(post180d)

    # TI delay and TR delay
    remove_delayTR = time_remove_TI + TE[0] + pp.calc_duration(gx)/2 + pp.calc_duration(gx_post, gy_post, gz_spoil)

    # Add inversion pulse duration to delays
    remove_delayTR = remove_delayTR + rf_inv_dur/2

    # Build sequence block SE
    seq.add_block(gz_s, gx_s, gy_s)
    seq.add_block(rf90, gz90)
    seq.add_block(gx_pre, gy_pre, gz_reph)

    seq.add_block(pre180delay)
    seq.add_block(rf180, gss_spoil_add)
    seq.add_block(post180delay)

    #seq.add_block(gy_pre)

    seq.add_block(gx, adc)
    seq.add_block(gx_post, gy_post, gz_spoil)

    return seq, TR, Ny, time_remove_TI, remove_delayTR