from packages.simulate_seq import simulate_seq
from packages.plot_sim import plot_sim
from packages.analyze_seq import analyze_seq

# ======
# SEQUENCE
# # ======
# from write_tse_3d import main as write_seq
# seq_name = "tse_3d"

from write_tse_3d_qmri import main as write_seq
seq_name = "tse_3d_qmri"

# from write_tse_3d_ptb_untouched import main as write_seq
# seq_name = "tse_3d_ptb_untouched"

# from write_tse_3d_mte import main as write_seq
# seq_name = "tse_3d_mte"

# from write_se_t1_mapping import main as write_seq
# seq_name = "se_t1_mapping"
# fov = (160e-3, 160e-3, 8e-3)
# n = 32 
# n_enc = (n, n, 1) 

# from write_tse_3d_demo import main as write_seq 
# seq_name = "tse_3d_demo"

# from write_MPRAGE import main as write_seq
# seq_name = "mprage"

# ======
# SCANNER
# ======
from packages.mr_systems import low_field as system
seq_name = seq_name + '_lowfield'


# ======
# FOV
# ======

#channels = Channels(ro = "x", pe1 = "y", pe2 = "z")
fov = (256e-3, 256e-3, 256e-3)
n_enc = (64, 64, 1) 


write_sequence = True
analyze_sequence = False
simulate_sequence = True
plot_simulation = True

if write_sequence:
    if plot_simulation or analyze_sequence:
        write_seq(plot=False, write_seq=True, seq_filename=seq_name, system=system, fov=fov, n_enc=n_enc)
    else:
        write_seq(plot=True, write_seq=True, seq_filename=seq_name, system=system, fov=fov, n_enc=n_enc)
    
if analyze_sequence:
    analyze_seq(seq_filename=seq_name, system=system)
            
if simulate_sequence:
    simulate_seq(save=True, seq_filename=seq_name, system=system)
    
if plot_simulation:
    if analyze_sequence:
        plot_sim(plot={
            "phantom": False,
            "seq": False,
            "kspace": False,
            "reco": True,
            }, seq_filename=seq_name, system=system) 
    else:
        plot_sim(plot={
            "phantom": False,
            "seq": False,
            "kspace": False,
            "reco": True,
            }, seq_filename=seq_name, system=system) 
