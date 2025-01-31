# %% Import packages
import numpy as np
import pypulseq as pp
import MRzeroCore as mr0
import pickle
import os 

from packages.write_seq_definitions import custom_seq_definitons

# %% Create phantom, simulate sequence, reconstruct image
def simulate_seq(save: bool,
                seq_filename: str,
                system:pp.Opts,
                seq_path: str = "./sequences/",
                sim_path: str = "./simulation/",
                add_noise: bool = False
                ):
    sim_name = "sim_" + seq_filename
    seq_file = seq_path + seq_filename + ".seq"

    seq = pp.Sequence(system=system)
    seq.read(seq_file, detect_rf_use = True)

    if 'k_space_encoding2' not in seq.definitions:
        seq.definitions['k_space_encoding2'] = 1
    if 'ro_oversampling' not in seq.definitions:
        seq.definitions['ro_oversampling'] = 1
    if 'name' not in seq.definitions:
        if 'Name' in seq.definitions:
            seq.definitions['name'] = seq.definitions['Name']
        else:
            seq.definitions['name'] = ''

    try:
        n_enc = (seq.definitions['number_of_readouts'], seq.definitions['k_space_encoding1'], seq.definitions['k_space_encoding2'])
        n_enc = tuple(map(int, n_enc))
    except KeyError:
        n_enc = tuple(map(int,seq.definitions['encoding_dim']))
    
    if seq.definitions['name'] == 'tse_3d_mte':
        n_echo = int(seq.definitions['etl'])

    if 'multi_echo' in seq.definitions['name']:
        n_echo = int(seq.definitions['repetition'])
        
    fov = tuple(seq.definitions['FOV'])
    ro_oversampling = int(seq.definitions['ro_oversampling'])

    # Remove definitions from the sequence because they cause import error in mr0.Sequence.import_file
    temp_seq_file = seq_path + seq_filename + '_temp.seq'
    for definition in custom_seq_definitons:
        seq.definitions.pop(definition, None)   
    seq.write(temp_seq_file)
    # import the sequence with MRzero
    seq0 = mr0.Sequence.import_file(temp_seq_file)
    # Delete the temporary sequence file
    os.remove(temp_seq_file)

    # Reload seq file
    seq.read(seq_file, detect_rf_use = True)

    # Setup spin system/object on which we can run the MR sequence
    sel_phantom = "simbrain" # select phantom type: invivo, simbrain, pixel
    if sel_phantom == 'simbrain' and n_enc[2] > 1:
        print('3D simulation is not supported for simbrain phantom - change to invivo phantom')
        sel_phantom = "invivo"

    print('load phantom')
    if sel_phantom == 'pixel':
        obj_p = mr0.CustomVoxelPhantom(
            pos=[[0., 0., 0]],
            PD=[1.0],
            T1=[3.0],
            T2=[0.5],
            T2dash=[30e-3],
            D=[0.0],
            B0=0,
            voxel_size=0.1,
            voxel_shape="box"
        )
    elif sel_phantom == 'simbrain':
        obj_p = mr0.VoxelGridPhantom.load_mat('./data/numerical_brain_cropped.mat') # has only 1 slice
        obj_p = obj_p.interpolate(n_enc[0], n_enc[1], 1)  
        obj_p.B0[:] = 0
        obj_p.D[:] = 0
    elif sel_phantom == 'invivo':
        obj_p = mr0.VoxelGridPhantom.brainweb("./data/subject05.npz")
        obj_p.B0[:] = 0 # Remove B0 inhomogeneity
        obj_p.D[:] = 0 # Remove diffusion
        if n_enc[2] == 1:
            center_slice = obj_p.PD.shape[2]//2
            obj_p = obj_p.slices([center_slice])    # select center slice
            obj_p = obj_p.interpolate(n_enc[0], n_enc[1], 1)  # interpolate     
        else:
            range_slices = tuple(np.linspace(40, 90, n_enc[2], dtype=int))
            obj_p = obj_p.slices(range_slices)      # select slices within the range
            obj_p = obj_p.interpolate(n_enc[0], n_enc[1], n_enc[2])      # interpolate
            obj_p.size[0] = 0.22
            obj_p.size[1] = 0.22
            obj_p.size[2] = 0.032
    else:
        print('Select proper phantom')
    obj_sim = obj_p.build()
    
    # SIMULATE the external.seq file and add acquired signal to ADC plot
    graph = mr0.compute_graph(seq0, obj_sim, 200, 1e-3)
    signal = mr0.execute_graph(graph, seq0, obj_sim)
    if add_noise:         # additional noise as simulation is ideal
        signal = signal + 1e-6 * np.random.randn(signal.shape[0], 2).view(np.complex128)

# check if field seq.definitions['name'] contains 'multi_echo' or 'tse_3d_mte'
    if 'name' not in seq.definitions:
        if 'Name' in seq.definitions:
            seq.definitions['name'] = seq.definitions['Name']
        else:
            seq.definitions['name'] = ''
    if 'tse_3d_mte' in seq.definitions['name']:
        # 3D FFT
        def fft_3d(x):
            return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x, axes=(0, 1, 3)), axes=(0, 1, 3)), axes=(0, 1, 3))

        kspace = np.reshape(signal, (n_enc[2], n_enc[1], n_echo, ro_oversampling*(n_enc[0])))
        kspace = kspace[:, :, :, 0:ro_oversampling*(n_enc[0])]

        reco = fft_3d(kspace)
        reco = reco[:, :, :, round(ro_oversampling*n_enc[0]/2 - n_enc[0]/2) : round(ro_oversampling*n_enc[0]/2 + n_enc[0]/2)]
        reco = np.transpose(reco, (3, 1, 0, 2))
    if 'gre_cartesian' in seq.definitions['name']:        
        # 2D FFT
        def fft_2d(x):
            return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x, axes=(-2, -1)), axes=(-2, -1)), axes=(-2, -1))

        # kspace = np.reshape(signal, (n_enc[2], n_enc[1], ro_oversampling*(n_enc[0])))
        # kspace = kspace[:, :, 0:ro_oversampling*(n_enc[0])]

        # reco = np.zeros((n_enc[2], n_enc[1], n_enc[0]), dtype=np.complex128)
        # for i in range(n_enc[2]):
        #     reco[i, :, :] = fft_2d(np.squeeze(kspace[i, :, :]))
        # reco = np.transpose(reco, (1, 2, 0))        

        kspace = np.reshape(signal, (n_enc[2], n_enc[1], n_enc[0]))

        reco = np.zeros((n_enc[0], n_enc[1], n_enc[2]), dtype=np.complex128)
        for i in range(n_enc[2]):
            reco[:, :, i] = fft_2d(np.squeeze(kspace[i, :, :]))


    else:
        reco = mr0.reco_adjoint(signal, seq0.get_kspace(), resolution=n_enc, FOV=fov) # Recommended: RECO has same Reso and FOV as sequence


    # kk = signal[0:4096]
    # kk_r = np.reshape(kk, (64, 64))
    # rr = fft_2d(kk_r)

    # import matplotlib.pyplot as plt

    # # # Plot the first slice of the input
    # # plt.imshow(np.abs(obj_p.PD[:, :, 3]), cmap='gray')
    # # plt.title('Reconstructed Image - First Slice')
    # # plt.colorbar()
    # # plt.show()
    # # plt.pause(0.1)

    # # Plot the first slice of the reconstructed image
    # plt.imshow(np.abs(rr), cmap='gray')
    # plt.title('Reconstructed Image - First Slice')
    # plt.colorbar()
    # plt.show()
    # plt.pause(0.1)

    # # # Plot the absolute value of the signal
    # # plt.figure()
    # # plt.plot(np.abs(signal))
    # # plt.title('Absolute Value of Signal')
    # # plt.xlabel('Time')
    # # plt.ylabel('Signal Amplitude')
    # # plt.grid(True)
    # # plt.show()

    # # Plot the first slice of the reconstructed image
    # plt.imshow(np.abs(reco[:, :, 0]), cmap='gray')
    # plt.title('Reconstructed Image - First Slice')
    # plt.colorbar()
    # plt.show()
    # plt.pause(0.1)

    # %% save results
    if save:
        # Check if directory exists
        if not os.path.exists(sim_path):
            # Create the directory
            os.makedirs(sim_path)
            print(f"Directory '{sim_path}' created.")

        with open(sim_path+ sim_name + '_obj_p.pkl', 'wb') as file:
            pickle.dump(obj_p, file)

        np.save(sim_path + sim_name + '_signal.npy', signal)

        with open(sim_path + sim_name + '_reco.pkl', 'wb') as file:
            pickle.dump(reco, file)
            
        seq_file = sim_name + '_seq.seq'
        seq.write(sim_path + seq_file)

        with open(sim_path + sim_name + '_seq0.pkl', 'wb') as file:
            pickle.dump(seq0, file)

if __name__ == "__main__":
    simulate_seq(save=True, seq_filename='tse_pypulseq')