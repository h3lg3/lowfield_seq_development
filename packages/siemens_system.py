"""Global definition of siemens system settings to be imported by sequence constructors."""

from pypulseq.opts import Opts

system = Opts(
    # standard max grad is 40 mT/m
    # max_grad=24,
    # grad_unit="mT/m",
    
    # time delay at the end of RF event, SETS RF DELAY!
    rf_dead_time=100e-6,

    # Set raster times to spectrum card frequency (timing checks)
    grad_raster_time=1e-5,
    rf_raster_time=1e-6,
    block_duration_raster=1e-5,
    adc_raster_time=1e-7,

    # Time delay at the beginning of an RF event
    rf_ringdown_time=100e-6,
    
    # Time delay at the beginning of ADC event
    adc_dead_time=10e-6,

    # Set maximum slew rate
    max_slew=50,
    slew_unit="T/m/s",
)