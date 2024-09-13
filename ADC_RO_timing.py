## Timing for ADC + RO + ROpre to account for Siemens ADC samples taken at ADC raster center instead of ADC raster edges

## How to implement that for tse_3d_ptb
# Option 1)
# Shift ADC by half the dwell time to the start of RO gradient
# Problem? ADC center is not at center of readout, center sample is not at peak of RF signal
adc = pp.make_adc(
    system=system,
    num_samples=int(n_enc_ro*ro_oversampling),
    duration=raster(val=adc_duration, precision=system.adc_raster_time),
    # Add gradient correction time and ADC correction time, HH: why 2*gradient_correction?
    # delay=raster(val=2 * gradient_correction + grad_ro.rise_time, precision=system.adc_raster_time)
    # Account for Siemens  
    delay=raster(val=2 * gradient_correction + grad_ro.rise_time - 0.5 * adc_duration/int(n_enc_ro*ro_oversampling), precision=system.adc_raster_time)
)
# Option 2)
# Use RO prephasor to dephase k-space a but further such that first kspace sample is taken at the right location.
# Problem? k-space is not symmetric anymore and after each 180° pulse, additional phase flips on other side of k-space
# therefore additional rephasor after ADC is needed to compensate that

# PREPHASOR
# # Calculate readout prephaser with correction times -> add small rephasor event after adc 
# grad_ro_pre = pp.make_trapezoid(
#     channel=channel_ro,
#     system=system,
#     area=grad_ro.area / 2 + grad_ro.amplitude*adc.dwell*0.5,
#     rise_time=ramp_duration,
#     fall_time=ramp_duration,
#     duration=raster(pp.calc_duration(grad_ro) / 2, precision=system.grad_raster_time),
# )
# ro_pre_duration = pp.calc_duration(grad_ro_pre) 

# REPHASOR
# # Calculate readout prephaser with correction times -> add small rephasor event after adc 
# ,
#     pp.make_trapezoid(
#         channel=channel_ro,
#         area=grad_ro.amplitude*adc.dwell,
#         duration=ro_pre_duration/2,
#         system=system,
#         rise_time=ramp_duration,
#         fall_time=ramp_duration
#     )


import pypulseq as pp
import numpy as np

system = pp.Opts(
    max_grad=30,  
    grad_unit="mT/m",
    max_slew=75,
    slew_unit="T/m/s",
    rf_dead_time=100e-6,
    rf_ringdown_time=100e-6,
    rf_raster_time=1e-6,
    grad_raster_time=1e-5,
    block_duration_raster=1e-5,
    adc_raster_time=1e-7,
)

system.adc_dead_time = 10e-6    # time delay at the beginning of ADC event
dG = 250e-6                     # rise time readout
sampling_time= 6.4e-3

fov = 256e-3
Nx = 64

delta_k=1/fov
k_width = Nx*delta_k

## Implementation for GRE sequence 
# gr_acq.area = 259.765625
# gr_acq.flat_time = 6.4e-3
# adc.duration = gr_acq.flat_time = 6.4e-3
# adc.delay = gr_acq.rise_time = 250e-6
# gr_preph.area = gr_acq.area / 2 = -129.8828125
gr_acq = pp.make_trapezoid(
    channel="x", flat_area=k_width, flat_time=sampling_time, rise_time=dG, system=system
)
adc = pp.make_adc(
    num_samples=Nx, duration=sampling_time, delay=gr_acq.rise_time, system=system
)
gr_preph = pp.make_trapezoid(
    channel="x", area=-gr_acq.area / 2, rise_time=dG, system=system
)
    
## Implementation for 2D TSE sequence
# gr_acq.area = 259.73520249221184
# gr_acq.flat_time = 6.42e-3
# adc.duration = gr_acq.flat_time - 2*system.adc_dead_time = 6.4e-3
# adc.delay = system.adc_dead_time = 1e-5
# gr_preph.area = gr_acq.area / 2 = 129.86760124610592

readout_time = sampling_time + 2*system.adc_dead_time
gr_acq = pp.make_trapezoid(
    channel="x",
    system=system,
    flat_area=k_width,
    flat_time=readout_time,
    rise_time=dG,
)
adc = pp.make_adc(
    num_samples=Nx, duration=sampling_time, delay=system.adc_dead_time
) # doesnt the ADC start before GRacq reached full amplitude? 
    
# Readout prephasor
gr_preph = pp.make_trapezoid(
        channel="x", system=system, area=gr_acq.area/2, rise_time=dG
)

## Implementation for MPRAGE sequence
# gr_acq.area = 260.15625000000006
# gr_acq.flat_time = 6.41e-3
# adc.duration = gr_acq.flat_time - system.adc_dead_time = 6.4e-3
# adc.delay = gr_acq.rise_time = 250e-6 
# gr_preph.area = gr_acq.area / 2 = 129.86760124610592
readout_time = np.ceil(
    (sampling_time + system.adc_dead_time) / system.grad_raster_time
    )* system.grad_raster_time

gr_acq = pp.make_trapezoid(
    channel="x",
    amplitude=k_width / sampling_time,
    flat_time=readout_time,
    rise_time=dG,
    system=system,
)
adc = pp.make_adc(
    num_samples=Nx,
    duration=sampling_time,  # total duration? or duration after delay?
    delay=gr_acq.rise_time, # is ADC dead time covered in delay or added afterwards?
    system=system,
)
#  First 0.5 is necessary to account for the Siemens sampling in the center of the dwell periods
gro_pre = pp.make_trapezoid(
    channel="x",
    area=-gr_acq.amplitude
    * (adc.dwell * (adc.num_samples / 2 + 0.5) + 0.5 * gr_acq.rise_time),
    system=system,
)
