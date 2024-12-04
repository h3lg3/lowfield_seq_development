import pypulseq as pp
system = pp.Opts()

# Set max slew to 50 T/m/s
system.max_slew = 50 * system.gamma  # Hz/m/s

grad_spoil = pp.make_trapezoid(channel='z', system=system, area = 800, duration = 1e-3, rise_time = 200e-6)

grad_slew_rate = grad_spoil.amplitude/grad_spoil.rise_time/system.gamma
print(f"Slew rate of gradient: {round(grad_slew_rate)} T/m/s")

# Slew rate of gradient: 117 T/m/s