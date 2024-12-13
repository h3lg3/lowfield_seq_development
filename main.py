from packages.simulate_seq import simulate_seq
from packages.plot_sim import plot_sim
from packages.analyze_seq import analyze_seq

# ======
# SEQUENCE
# # ======
from write_tse_3d import main as write_seq
seq_name = "tse_3d"

# from write_tse_3d_mte import main as write_seq
# seq_name = "tse_3d_mte"

# from write_se_t1_mapping import main as write_seq
# seq_name = "se_t1_mapping"

# from write_tse_3d_demo import main as write_seq 
# seq_name = "tse_3d_demo"

# from write_tse_3d_ptb_untouched import main as write_seq
# seq_name = "tse_3d_ptb_untouched"

# from write_MPRAGE import main as write_seq
# seq_name = "mprage"

# ======
# SCANNER
# ======
from packages.mr_systems import lumina as system
seq_name = seq_name + '_lumina'

#seq_name = 'tse_3d_lumina_64_64_32_TR1500'
# ======
# FOV
# ======
fov = (160e-3, 160e-3, 160e-3)
n = 64
n_enc = (n, n, 32) 

write_sequence = True
analyze_sequence = True
simulate_sequence = False
plot_simulation = False

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
