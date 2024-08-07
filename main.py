from packages.simulate_seq import simulate_seq
from packages.plot_sim import plot_sim
from packages.analyze_seq import analyze_seq

# from write_tse_pypulseq import main as write_seq
# seq_name = "tse_pypulseq"

from write_tse_3d_console import main as write_seq
seq_name = "tse_3D_console"

# from write_MPRAGE import main as write_seq
# seq_name = "mprage_pypulseq"

# from write_3Dt1_mprage import main as write_seq
# seq_name = "3Dt1_mprage_pypulseq"

write_sequence = False
analyze_sequence = False
simulate_sequence = False
plot_simulation = True

if write_sequence:
    if plot_simulation:
        write_seq(plot=False, write_seq=True, seq_filename=seq_name)
    else:
        write_seq(plot=True, write_seq=True, seq_filename=seq_name)
    
if analyze_sequence:
    analyze_seq(seq_filename=seq_name)
            
if simulate_sequence:
    simulate_seq(save=True, seq_filename=seq_name)
    
if plot_simulation:
    plot_sim(plot={
        "phantom": False,
        "seq": True,
        "kspace": True,
        "reco": True
        }, seq_filename=seq_name) 