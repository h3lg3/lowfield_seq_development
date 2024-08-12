from matplotlib import pyplot as plt
from packages import tse_3d
from packages.tse_trajectory import Trajectory
from console.interfaces.interface_acquisition_parameter import Dimensions
from math import pi

def main(plot: bool, write_seq: bool, seq_filename: str = "tse_3d_ptb.seq"):
    seq = tse_3d.constructor(
                             echo_time=14e-3,
                             repetition_time=1000e-3, 
                             etl=10, # define max sampling period (tmax = 200ms?), etl_max = round(tmax/esp), nr. of pe1 steps should be multiple of etl
                             dummies=10, 
                             ro_bandwidth=20e3,
                             ro_oversampling = 2, 
                             rf_duration = 100e-6,
                             fov=Dimensions(x=220e-3, y=220e-3, z=220e-3),  
                             n_enc=Dimensions(x=120, y=120, z=1),           
                             trajectory=Trajectory.SYMMETRIC,
                             refocussing_angle=pi/3,
                             excitation_phase=pi/2,
                             refocussing_phase=0,
                             channel_ro="x", 
                             channel_pe1="y",
                             channel_pe2="z" 
                             )[0]
    # LF setting
    # fov=Dimensions(x=140e-3, y=140e-3, z=140e-3), 
    # n_enc=Dimensions(x=1, y=120, z=120) # x along bore (HF, PE2), y bottom to top(AP, PE1), z along B0 (LR, readout),
    # channel_ro="z", 
    # channel_pe1="y",
    # channel_pe2="x"   

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
        seq.set_definition('Name', 'se_3d_ptb')
        seq.write('./sequences/' + seq_filename)
        # seq.write(r"C:\Users\hhert\VirtualMachines\SharedFolder\pulseq\external.seq")

if __name__ == "__main__":
    main(plot=True, write_seq=True)        
    
    
# k-space order for 32 phase enc steps
# # PYPULSEQ
# array([[-14., -13., -12., -11.],
#        [-10.,  -9.,  -8.,  -7.],
#        [ -6.,  -5.,  -4.,  -3.],
#        [ -2.,  -1.,   0.,   1.],
#        [  2.,   3.,   4.,   5.],
#        [  6.,   7.,   8.,   9.],
#        [ 10.,  11.,  12.,  13.],
#        [ 14.,  15., -16., -15.]])

# # LINEAR ORDER
# [array([ 15.5,  11.5,   7.5,   3.5,  -0.5,   4.5,  -8.5, -12.5]), 
# array([-14.5, -10.5,   6.5,  -2.5,   1.5,   5.5,  -9.5, -13.5]), 
# array([13.5,  9.5, -5.5, -1.5,  2.5, -6.5, 10.5, 14.5]), 
# array([ 12.5,   8.5,  -4.5,   0.5,  -3.5,  -7.5, -11.5, -15.5])]

# # OUTIN
# [array([15.5, 13.5, 11.5,  9.5,  7.5, -5.5,  3.5, -1.5]), 
# array([-15.5, -13.5, -11.5,  -9.5,  -7.5,   5.5,  -3.5,   1.5]),
#  array([-14.5,  12.5, -10.5,   8.5,   6.5,  -4.5,  -2.5,   0.5]), 
# array([ 14.5, -12.5,  10.5,  -8.5,  -6.5,   4.5,   2.5,  -0.5])]

# #INOUT
# [array([ -0.5,   2.5,   4.5,  -6.5,  -8.5,  10.5, -12.5,  14.5]), 
# array([  0.5,  -2.5,  -4.5,   6.5,   8.5, -10.5,  12.5, -14.5]), 
# array([  1.5,  -3.5,   5.5,  -7.5,  -9.5, -11.5, -13.5, -15.5]), 
# array([-1.5,  3.5, -5.5,  7.5,  9.5, 11.5, 13.5, 15.5])]