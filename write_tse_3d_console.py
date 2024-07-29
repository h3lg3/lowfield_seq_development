from matplotlib import pyplot as plt
from packages import tse_3d
from console.interfaces.interface_acquisition_parameter import Dimensions
from math import pi

def main(plot: bool, write_seq: bool, seq_filename: str = "tse_3D_ptb_console.seq"):
    select_fov = Dimensions(x=220e-3, y=220e-3, z=220e-3)
    select_encoding = Dimensions(x=64, y=64, z=1)

    seq = tse_3d.constructor(
                             echo_time=28e-3,
                             repetition_time=2000e-3, 
                             etl=8, 
                             dummies=5, 
                             ro_bandwidth=10e3, 
                             fov=select_fov, 
                             n_enc=select_encoding,
                             trajectory=tse_3d.Trajectory.LINEAR,
                             excitation_phase=pi/2,
                             refocussing_phase=0
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
        seq.set_definition('Name', 'se_ptb_console')
        seq.write('./sequences/' + seq_filename)
        # seq.write(r"C:\Users\hhert\VirtualMachines\SharedFolder\pulseq\external.seq")

if __name__ == "__main__":
    main(plot=True, write_seq=True)        