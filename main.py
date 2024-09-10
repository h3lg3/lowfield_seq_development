from packages.simulate_seq import simulate_seq
from packages.plot_sim import plot_sim
from packages.analyze_seq import analyze_seq

# ======
# SEQUENCE
# ======
# from write_tse_pypulseq import main as write_seq
# seq_name = "tse_2d"

from write_tse_3d_ptb import main as write_seq
seq_name = "tse_3d"

# from write_MPRAGE import main as write_seq
# seq_name = "mprage"

# ======
# SCANNER
# ======
from packages.mr_systems import lumina as system
seq_name = seq_name + '_lumina'

# from packages.mr_systems import low_field as system
# seq_name = seq_name + '_ptb'

# ======
# FOV
# ======
fov = (220e-3, 220e-3, 220e-3)
nk = (120, 120, 1)

write_sequence = False
analyze_sequence = False
simulate_sequence = True
plot_simulation = True

if write_sequence:
    if plot_simulation:
        write_seq(plot=False, write_seq=True, seq_filename=seq_name, system=system, fov=fov, nk=nk)
    else:
        write_seq(plot=True, write_seq=True, seq_filename=seq_name, system=system, fov=fov, nk=nk)
    
if analyze_sequence:
    analyze_seq(seq_filename=seq_name)
            
if simulate_sequence:
    simulate_seq(save=True, seq_filename=seq_name, fov=fov, nk=nk)
    
if plot_simulation:
    plot_sim(plot={
        "phantom": False,
        "seq": False,
        "kspace": False,
        "reco": True
        }, seq_filename=seq_name) 
    

# # Lumina test setting
# fov = (220e-3, 220e-3, 220e-3)
# nk = (120, 120, 1)

#     seq = tse_3d.constructor(
#                             echo_time=16e-3,
#                             repetition_time=2000e-3,   # 600
#                             etl=10, # define max sampling period (tmax = 200ms?), etl_max = round(tmax/esp), nr. of pe1 steps should be multiple of etl
#                             dummies=5,     # 5
#                             ro_bandwidth=20e3,
#                             ro_oversampling = 2, 
#                             rf_duration = 100e-6,
#                             fov=Dimensions(x=fov[0], y=fov[1], z=fov[2]),  
#                             n_enc=Dimensions(x=nk[0], y=nk[1], z=nk[2]),           
#                             trajectory=Trajectory.SYMMETRIC,
#                             refocussing_angle=120/180 * pi,    # pi
#                             excitation_phase=pi/2,
#                             refocussing_phase=0,
#                             channel_ro="x", 
#                             channel_pe1="y",
#                             channel_pe2="z",
#                             system=system
#                             )[0]