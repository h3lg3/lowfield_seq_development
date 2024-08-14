"""Global definition of MR system settings to be imported by sequence constructors."""

from pypulseq.opts import Opts

lumina = Opts(
    max_grad=24,    # 24*system.gamma*1e-3
    grad_unit="mT/m",
    max_slew=75,
    slew_unit="T/m/s",
    
    rf_dead_time=100e-6, # time delay at the end of RF event, SETS RF DELAY!
    rf_ringdown_time=100e-6, # Time delay at the beginning of an RF event
    rf_raster_time=1e-6,
    
    grad_raster_time=1e-5,
    block_duration_raster=1e-5,
    
    adc_raster_time=1e-7,
    adc_dead_time=10e-6, # Time delay at the beginning of ADC event

    B0=3.0
)

low_field = Opts(
    max_grad=250.e3,
    max_slew=500.e3,
    rise_time=100.e-6,   # redundant with max_slew

    # time delay at the end of RF event, SETS RF DELAY!
    rf_dead_time=20.e-6,

    # Set raster times to spectrum card frequency (timing checks)
    grad_raster_time=5.e-8,
    rf_raster_time=5.e-8,
    block_duration_raster=5.e-8,
    adc_raster_time=5.e-8,

    # Time delay at the beginning of an RF event
    rf_ringdown_time=0.,
    
    # Time delay at the beginning of ADC event
    adc_dead_time=0.,
    B0=50.e-3
)