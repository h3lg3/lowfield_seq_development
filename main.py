from packages.simulate_seq import simulate_seq
from packages.plot_sim import plot_sim
from packages.analyze_seq import analyze_seq

# from write_tse_pypulseq import main as write_seq
# seq_name = "tse_pypulseq"

from write_tse_3d_console import main as write_seq
seq_name = "tse_3D_console"


write_sequence = True
analyze_sequence = False
simulate_sequence = True
plot_simulation = True

if write_sequence:
    write_seq(plot=True, write_seq=True, seq_filename=seq_name)
    
if analyze_sequence:
    analyze_seq(seq_filename=seq_name)
            
if simulate_sequence:
    simulate_seq(save=True, seq_filename=seq_name)
    
if plot_simulation:
    plot_sim(plot={
        "phantom": True,
        "seq": True,
        "kspace": True,
        "reco": True
        }, seq_filename=seq_name) 