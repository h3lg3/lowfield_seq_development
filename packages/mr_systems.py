"""Global definition of MR system settings to be imported by sequence constructors."""

from pypulseq.opts import Opts
GAMMA = 42576000
lumina = Opts(
    max_grad=20,    # 24*system.gamma*1e-3
    grad_unit="mT/m",
    max_slew=40,
    slew_unit="T/m/s",
    
    rf_dead_time=100e-6, # time delay at the beginning of RF event, SETS RF DELAY!
    rf_ringdown_time=100e-6, # Time delay at the end of an RF event
    rf_raster_time=1e-6,
    
    grad_raster_time=10e-6,
    block_duration_raster=10e-6,
    
    adc_raster_time=100e-9,
    adc_dead_time=10e-6, # Time delay at the BEGINNING/END of ADC event

    B0=3.0
)

low_field = Opts(
    max_grad = 300e3/GAMMA*1e3,
    grad_unit = "mT/m",

    # Set maximum slew rate
    max_slew = 50,
    slew_unit = "T/m/s",

    # time delay at the beginning of RF event, SETS RF DELAY!
    rf_dead_time=20.e-6,

    # Set raster times to spectrum card frequency (timing checks)
    grad_raster_time=1e-6,
    rf_raster_time=1e-6,
    block_duration_raster=1e-6,
    adc_raster_time=1e-6,

    # Time delay at the end of an RF event
    rf_ringdown_time=0, # 2e-3
    
    # Time delay at the BEGINNING/END of ADC event, at beginning it is usually covered by delay
    #adc_dead_time=200e-6,
    
    B0=50.e-3
)
