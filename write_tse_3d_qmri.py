from matplotlib import pyplot as plt
from math import pi

from packages import tse_3d_qmri
from packages.seq_utils import Trajectory, Dimensions, Channels
from packages.mr_systems import low_field as default_system

from pypulseq.opts import Opts
import os

def main(plot:bool, write_seq:bool, seq_filename:str = "tse_3d_qmri",
         system:Opts = default_system, 
         fov:tuple = (256e-3, 256e-3, 256e-3), 
         n_enc:tuple = (64, 64, 64)
         ):
    seq = tse_3d_qmri.constructor(
                            input_fov=Dimensions(x = fov[0], y = fov[1], z = fov[2]),
                            input_enc=Dimensions(x = n_enc[0], y = n_enc[1], z = n_enc[2]),
                            te=80e-3,
                            tr=2000e-3,
                            etl=10,
                            n_dummies=10,
                            adc_bandwidth=20e3,
                            adc_oversampling=1,
                            g_ro_correction=0,
                            rf_ex_duration=200e-6,
                            rf_ref_duration=200e-6,
                            rf_ex_phase=0,
                            rf_ref_phase=pi / 2,
                            trajectory=Trajectory.INOUT,
                            channels=Channels(ro = "x", pe1 = "y", pe2 = "z"),
                            system=system,
                            )[0]

    if plot:
        plot_kspace = True
        plot_seq = True
    else:
        plot_kspace = False
        plot_seq = False
        
        
    ## Check whether the timing of the sequence is compatible with the scanner
    (ok,error_report,) = seq.check_timing()  # Check whether the timing of the sequence is correct
    if ok:
        print("Timing check passed successfully")
        print("Sequence duration is: ", round(seq.duration()[0]), "s")
    else:
        print("Timing check failed. Error listing follows:")
        [print(e) for e in error_report]        
        
    if plot_kspace:
        k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

        plt.figure()
        plt.title('full k-space trajectory ($k_{x}$ x $k_{y}$)')
        plt.plot(k_traj[0],k_traj[1])
        plt.plot(k_traj_adc[0],k_traj_adc[1],'.')
        plt.xlabel('kx [1/m]')
        plt.ylabel('ky [1/m]')
        plt.show()
                        
    if plot_seq:
        if n_enc[1]*n_enc[2] > 900  or seq.duration()[0] > 600:
            print("Plotting only the first 20% of the sequence")
            seq.plot(time_range = (0, round(0.2*seq.duration()[0])))
        else:
            seq.plot()
    # rep = seq.test_report()
    # print(rep)
    # =========
    # WRITE .SEQ
    # =========    
    if write_seq:
        seq.set_definition('Name', seq_filename)
        seq.write('./sequences/' + seq_filename)
        if os.path.exists(r"C:\Users\hhert\VirtualMachines\SharedFolder\pulseq"):
            seq.write(r"C:\Users\hhert\VirtualMachines\SharedFolder\pulseq\external.seq")
        else:
            print("Shared folder not found. Sequence not written to shared folder.")
if __name__ == "__main__":
    main(plot=True, write_seq=True)        
    