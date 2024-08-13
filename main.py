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
from packages.siemens_system import system
seq_name = seq_name + '_lumina'

# from packages.lf_system import system
# seq_name = seq_name + '_ptb'

# ======
# FOV
# ======
fov = (220e-3, 220e-3, 220e-3)
nk = (64, 64, 1)

write_sequence = True
analyze_sequence = True
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
        "phantom": True,
        "seq": True,
        "kspace": True,
        "reco": True
        }, seq_filename=seq_name) 