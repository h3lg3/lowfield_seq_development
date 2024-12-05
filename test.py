import pypulseq as pp
import numpy as np
system = pp.Opts()

# Set max slew to 40 T/m/s
system.max_slew = 40 * system.gamma  # Hz/m/s

gs_rise_time = 500e-6
gs_fall_time = 500e-6
gs_flat_time = 440e-6
gs_amplitude = 851063.829787234

grf_rise_time = 200e-6
grf_fall_time = 200e-6
grf_flat_time = 2.5e-3
grf_amplitude = 256000
gss_times = np.cumsum([0, 
                        gs_rise_time, 
                        gs_flat_time,
                        gs_fall_time - grf_rise_time, 
                        grf_flat_time,
                        gs_rise_time - grf_fall_time,
                        gs_flat_time, 
                        gs_fall_time])

gss_amps = np.array([0, 
                        gs_amplitude, 
                        gs_amplitude, 
                        grf_amplitude, 
                        grf_amplitude, 
                        gs_amplitude, 
                        gs_amplitude, 
                        0])

gss_spoil_add = pp.make_extended_trapezoid(channel = 'z', amplitudes = gss_amps, times = gss_times, system = system)

g = pp.points_to_waveform(times=gss_times, amplitudes=gss_amps, grad_raster_time=system.grad_raster_time)

slew = np.squeeze(np.subtract(g[1:], g[:-1]) / system.grad_raster_time)
print(f"Max slew rate of gradient: {round(max(slew)/system.gamma, 2)} T/m/s")

