from matplotlib import pyplot as plt
from packages import tse_3d_untouched
from packages.seq_utils import Trajectory
from packages.seq_utils import Dimensions
#from console.interfaces.dimensions import Dimensions
from math import pi
from pypulseq.opts import Opts
from packages.mr_systems import low_field as default_system

def main(plot:bool, write_seq:bool, seq_filename:str = "tse_3d_ptb_untouched",
         system:Opts = default_system, 
         fov:tuple = (256e-3, 256e-3, 256e-3), 
         n_enc:tuple =(64, 64, 64)
         ):
    seq = tse_3d_untouched.constructor(
                            echo_time = 20e-3,
                            repetition_time = 500e-3,  
                            etl = 4, # define max sampling period (tmax = 200ms?), etl_max = round(tmax/esp), nr. of pe1 steps should be multiple of etl
                            dummies = 2,    
                            ro_bandwidth = 20e3,
                            rf_duration = 160e-6,
                            fov=Dimensions(x=fov[0], y=fov[1], z=fov[2]),  
                            n_enc=Dimensions(x=n_enc[0], y=n_enc[1], z=n_enc[2]),           
                            trajectory=Trajectory.LINEAR,
                            gradient_correction = 80e-6,
                            refocussing_angle = pi,  
                            excitation_phase = pi/2,
                            refocussing_phase = 0,
                            channel_ro = "x", 
                            channel_pe1 = "y",
                            channel_pe2 = "z",
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
    else:
        print("Timing check failed. Error listing follows:")
        [print(e) for e in error_report]        
        
    if plot_kspace:
        k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

        plt.figure()
        plt.plot(k_traj[0],k_traj[1])
        plt.plot(k_traj_adc[0],k_traj_adc[1],'.')
        plt.show()
                
    if plot_seq:
        seq.plot()
        
    # =========
    # WRITE .SEQ
    # =========    
    if write_seq:
        seq.set_definition('Name', seq_filename)
        seq.write('./sequences/' + seq_filename)
        seq.write(r"C:\Users\hhert\VirtualMachines\SharedFolder\pulseq\external.seq")

if __name__ == "__main__":
    main(plot=True, write_seq=True)        
    