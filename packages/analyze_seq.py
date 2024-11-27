# %% Setup script
# import packages and .seq file
import pypulseq as pp
import matplotlib.pyplot as plt
import numpy as np
from packages.write_seq_definitions import read_definitions

def analyze_seq(seq_filename: str, system:pp.Opts, seq_path: str = "./sequences/"):

    # Create a Sequence object and read the sequence file
    seq = pp.Sequence(system=system)
    seq.read(seq_path + seq_filename + ".seq", detect_rf_use = True)
# %% Analyze sequence
    # Plot sequence
    if seq.definitions['k_space_encoding1']*seq.definitions['k_space_encoding2'] > 900:
        print("Analyzing only first 20% of sequence")
        seq.plot(time_range = (0, round(0.2*seq.duration()[0])))
    else:
        seq.plot()
    
    # Plot k-space trajectory
    k_traj_adc, k_traj = seq.calculate_kspace()[0:2]

    plt.figure()
    plt.title('full k-space trajectory ($k_{x}$ x $k_{y}$)')
    plt.plot(k_traj[0],k_traj[1])
    plt.plot(k_traj_adc[0],k_traj_adc[1],'.')
    plt.xlabel(r'k_x [1/m]')
    plt.ylabel(r'k_y [1/m]')
    plt.show()
    
    # # %% Analyze waveforms: Gradient amplitude and slew rates:: 
    # # Combined amplitudes and slew rate for case of oblique scanning interesting
    # # Values should stay below system specs even for oblique scanning
    wf = {}
    if seq.definitions['k_space_encoding1']*seq.definitions['k_space_encoding2'] > 900:
        wf_grad = seq.waveforms(time_range = (0, round(0.2*seq.duration()[0])))
    else:
        wf_grad = seq.waveforms()

    wf['t_gx'] = wf_grad[0][0]
    wf['t_gy'] = wf_grad[1][0]
    wf['t_gz'] = wf_grad[2][0]

    wf['gx'] = wf_grad[0][1]
    wf['gy'] = wf_grad[1][1]
    wf['gz'] = wf_grad[2][1]
    min_t = np.min(np.concatenate((wf['t_gx'], wf['t_gy'], wf['t_gz'])))
    max_t = np.max(np.concatenate((wf['t_gx'], wf['t_gy'], wf['t_gz'])))
    dt = 10e-6
    wf_interp = {}
    wf_interp['t'] = np.arange(min_t, max_t, dt)
    # Gradient amplitude, unit mT/m
    wf_interp['gx'] = np.interp(wf_interp['t'], wf['t_gx'], wf['gx'])/(seq.system.gamma*1e-3)
    wf_interp['gy'] = np.interp(wf_interp['t'], wf['t_gy'], wf['gy'])/(seq.system.gamma*1e-3)
    wf_interp['gz'] = np.interp(wf_interp['t'], wf['t_gz'], wf['gz'])/(seq.system.gamma*1e-3)
    wf_interp['g_norm'] = np.sqrt(pow(wf_interp['gx'], 2.0) + pow(wf_interp['gy'], 2.0) + pow(wf_interp['gz'], 2.0))
    # Slew rate, unit T/m/s
    wf_interp['slew_x'] = np.diff(wf_interp['gx'])/dt*1e-3
    wf_interp['slew_y'] = np.diff(wf_interp['gy'])/dt*1e-3
    wf_interp['slew_z'] = np.diff(wf_interp['gz'])/dt*1e-3
    wf_interp['slew_norm'] = np.sqrt(pow(wf_interp['slew_x'], 2.0) + pow(wf_interp['slew_y'], 2.0) + pow(wf_interp['slew_z'], 2.0))
    # # Gradient defined by vertices
    # plt.figure()
    # plt.title('gradients')
    # plt.plot(wf['t_gx'],wf['gx'], '.-', c='r')
    # plt.plot(wf['t_gy'],wf['gy'], '.-', c='g')
    # plt.plot(wf['t_gz'],wf['gz'], '.-', c='b')
    # plt.xlabel("time in s")
    # plt.ylabel("gradient amplitude in Hz/m")
    # plt.legend(["x", "y", "z"])
    # plt.show()
    # Gradient amplitudes
    fig, axs = plt.subplots(2, 1, figsize=(10, 8))

    # Gradient amplitudes
    axs[0].plot(wf_interp['t'], wf_interp['gx'], c='r')
    axs[0].plot(wf_interp['t'], wf_interp['gy'], c='g')
    axs[0].plot(wf_interp['t'], wf_interp['gz'], c='b')
    axs[0].plot(wf_interp['t'], wf_interp['g_norm'], c='k')
    axs[0].set_xlabel("time [s]")
    axs[0].set_ylabel("gradient amplitude [mT/m]")
    axs[0].legend(["x", "y", "z", "norm"], loc='upper right')
    axs[0].set_title('Gradient Amplitudes')

    # Slew rates
    axs[1].plot(wf_interp['t'][0:-1], wf_interp['slew_x'], c='r')
    axs[1].plot(wf_interp['t'][0:-1], wf_interp['slew_y'], c='g')
    axs[1].plot(wf_interp['t'][0:-1], wf_interp['slew_z'], c='b')
    axs[1].plot(wf_interp['t'][0:-1], wf_interp['slew_norm'], c='k')
    axs[1].set_xlabel("time [s]")
    axs[1].set_ylabel("slew rate [T/m/s]")
    axs[1].legend(["x", "y", "z", "norm"], loc='upper right')
    axs[1].set_title('Slew Rates')

    plt.tight_layout()
    plt.show()
    # %% Sequence test report
    # For the real TE, TR or for staying within slew-rate limits
    rep = seq.test_report()
    print(rep)
    # %% Simulate slice profile
    # https://github.com/pulseq/MR-Physics-with-Pulseq/blob/main/tutorials/02_rf_pulses/notebooks/se2d_sliceprofile_exercise.ipynb
    
    # # REQUIRES SEQUENCES COMPILED WITH PYPULSEQ Version > 1.4.0
    # # %% Calculate SAR only possible if sequence has certain length (>2s?)
    
    # if pp.Sequence.duration(seq)[0] > 2:
    #     pp.SAR.SAR_calc.calc_SAR(seq)
    # %% Calculate PNS
    # use example specs
    # seq.calculate_pns(pp.utils.safe_pns_prediction.safe_example_hw(), do_plots=True) 
    # use PRISMA specs
    #seq.calculate_pns('.\data\MP_GPA_K2368_2250V_950A_GC25_Lumina.asc', do_plots=True)  
    # %% Calculate mechanical resonances
    asc_dict = pp.utils.siemens.readasc.readasc('.\data\MP_GPA_K2368_2250V_950A_GC25_Lumina.asc')
    resonances = pp.utils.siemens.asc_to_hw.asc_to_acoustic_resonances(asc_dict[0])
    seq.calculate_gradient_spectrum(plot=True, acoustic_resonances=resonances)
    plt.show()
    # %%
if __name__ == "__main__":
    analyze_seq(seq_filename = "tse_pypulseq")    