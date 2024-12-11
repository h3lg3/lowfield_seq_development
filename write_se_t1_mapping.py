from matplotlib import pyplot as plt
from math import pi

from packages import se_t1_mapping
from packages.seq_utils import Dimensions, Channels
from packages.mr_systems import low_field as default_system

from pypulseq.opts import Opts

import os

def main(plot:bool, write_seq:bool, seq_filename:str = "tse_3d",
         system:Opts = default_system, 
         fov:tuple = (256e-3, 256e-3, 256e-3), 
         n_enc:tuple = (64, 64, 64)
         ):
    seq = se_t1_mapping.constructor(
                            echo_time = 15e-3, 
                            repetition_time = 5000e-3,
                            TI = [50e-3, 100e-3, 500e-3, 1500e-3],
                            slice_thickness = 8e-3,   
                            ro_bandwidth = 10e3,
                            ro_oversampling = 1, 
                            input_fov = Dimensions(x = fov[0], y = fov[1], z = fov[2]),  
                            input_enc = Dimensions(x = n_enc[0], y = n_enc[1], z = n_enc[2]),           
                            channels = Channels(ro = "x", pe1 = "y", pe2 = "z"),
                            system = system
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
        if n_enc[1]*n_enc[2] > 900 or seq.duration()[0] > 600:
            print("Plotting only the first 20% of the sequence")
            seq.plot(time_range = (0, round(0.2*seq.duration()[0])))
        else:
            seq.plot()
        
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
    