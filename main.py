from packages.simulate_seq import simulate_seq
from packages.plot_sim import plot_sim
from packages.analyze_seq import analyze_seq

# ======
# SEQUENCE
# ======
# from write_tse_pypulseq import main as write_seq
# seq_name = "tse_2d"

from write_tse_3d_demo import main as write_seq
seq_name = "tse_3d-demo"

from write_tse_3d_ptb import main as write_seq
seq_name = "tse_3d"

# from write_tse_3d_ptb_untouched import main as write_seq
# seq_name = "tse_3d_untouched"

# from write_MPRAGE import main as write_seq
# seq_name = "mprage"

# ======
# SCANNER
# ======
from packages.mr_systems import lumina as system
seq_name = seq_name + '_lumina'

# ======
# FOV
# ======
fov = (256e-3, 256e-3, 256e-3)
nk = (16, 16, 1)
# # Lumina test setting
# fov = (220e-3, 220e-3, 220e-3)
# nk = (120, 120, 1)

write_sequence = True
analyze_sequence = False
simulate_sequence = False
plot_simulation = False

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
        "seq": True,
        "kspace": True,
        "reco": True
        }, seq_filename=seq_name) 
