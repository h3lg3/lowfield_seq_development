from matplotlib import pyplot as plt
from packages import tse_3d
from packages.seq_utils import Trajectory
from packages.seq_utils import Dimensions
#from console.interfaces.dimensions import Dimensions
from math import pi
from pypulseq.opts import Opts
from packages.mr_systems import low_field as default_system

def main(plot:bool, write_seq:bool, seq_filename:str = "tse_2d_lumina.seq",
         system:Opts = default_system, 
         fov:tuple = (256e-3, 256e-3, 256e-3), 
         nk:tuple =(64, 64, 64)
         ):
    seq = tse_3d.constructor(
                            echo_time = 16e-3,
                            repetition_time = 2000e-3,  
                            etl = 1, # define max sampling period (tmax = 200ms?), etl_max = round(tmax/esp), nr. of pe1 steps should be multiple of etl
                            dummies = 5,    
                            ro_bandwidth = 20e3,
                            ro_oversampling = 1, 
                            rf_duration = 100e-6,
                            fov=Dimensions(x=fov[0], y=fov[1], z=fov[2]),  
                            n_enc=Dimensions(x=nk[0], y=nk[1], z=nk[2]),           
                            trajectory=Trajectory.SYMMETRIC,
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
        seq.set_definition('Name', 'se_3d_ptb')
        seq.write('./sequences/' + seq_filename)
        seq.write(r"C:\Users\hhert\VirtualMachines\SharedFolder\pulseq\external.seq")

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