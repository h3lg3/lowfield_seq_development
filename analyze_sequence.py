# %% 
import pypulseq as pp
import matplotlib.pyplot as plt
import numpy as np
import pypulseq.utils.siemens as siemens

# import specific function from package within pypulseq
#from pypulseq.utils.safe_pns_prediction import safe_example_hw
#from pypulseq.Sequence.calc_pns import calc_pns 

# Define the file name of the sequence
seq_file = './Examples/epi_se_pypulseq.seq'

# Create a Sequence object and read the sequence file
seq = pp.Sequence()
seq.read(seq_file)
# %% 
## Analyze sequence

## Plot k-space trajectory
k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

plt.figure()
plt.plot(k_traj[0],k_traj[1])
plt.plot(k_traj_adc[0],k_traj_adc[1],'.')
plt.show()

## waveforms
# waveforms = seq.waveforms_and_times()

## Very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within
## slew-rate limits
# rep = seq.test_report()
# print(rep)

## Simulate slice profile
# https://github.com/pulseq/MR-Physics-with-Pulseq/blob/main/tutorials/02_rf_pulses/notebooks/se2d_sliceprofile_exercise.ipynb

## Calculate SAR only possible if sequence has certain length (>2s?)
# if pp.Sequence.duration(seq)[0] > 2:
#     pp.SAR.SAR_calc.calc_SAR(seq)

## Calculate PNS
# use example specs
# seq.calculate_pns(pp.utils.safe_pns_prediction.safe_example_hw(), do_plots=True) 
# use PRISMA specs
# seq.calculate_pns('E:\Python\MP_GPA_K2309_2250V_951A_AS82.asc', do_plots=True)  

## Calculate mechanical resonances
# asc_dict = siemens.readasc.readasc('E:\Python\MP_GPA_K2309_2250V_951A_AS82.asc')
# resonances = siemens.asc_to_hw.asc_to_acoustic_resonances(asc_dict[0])
# seq.calculate_gradient_spectrum(plot=True, acoustic_resonances=resonances)
# plt.show()

# ## Plot the entire Sequence
# plt.figure(figsize=(12, 8))
# seq.plot()
# plt.show()

# # Define the time range for plotting
# time_start = 5e-3  # Start time in seconds (5 ms)
# time_stop = 15e-3  # Stop time in seconds (15 ms)

# # Plot the Sequence within the specific time range
# plt.figure(figsize=(12, 8))
# seq.plot(time_range=[time_start, time_stop])
# plt.show()

