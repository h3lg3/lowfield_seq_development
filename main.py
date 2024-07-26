from write_tse_pypulseq import main as write_seq
from packages.simulate_seq import simulate_seq
from packages.plot_sim import plot_sim
from packages.analyze_seq import analyze_seq

seq_name = "tse_pypulseq"
write_sequence = False
analyze_sequence = True
simulate_sequence = False
plot_simulation = False

if write_sequence:
    write_seq(plot=False, write_seq=True, seq_filename=seq_name)
    
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