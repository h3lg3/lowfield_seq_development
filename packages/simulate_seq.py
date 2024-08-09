# %% Import packages
import numpy as np
import pypulseq as pp
import MRzeroCore as mr0
import pickle

# %% Create phantom, simulate sequence, reconstruct image
def simulate_seq(save: bool, seq_filename: str, seq_path: str = "./sequences/", sim_path: str = "./simulation/"):
    sim_name = "sim_" + seq_filename
    seq_file = seq_path + seq_filename + ".seq"

    seq = pp.Sequence()
    seq.read(seq_file)
    seq0 = mr0.Sequence.import_file(seq_file)

    # Setup spin system/object on which we can run the MR sequence
    sel_phantom = "invivo" # select phantom type: invivo, simbrain, pixel
    reso = (64, 64, 1)  # Recommended: Same as Sequence

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
        obj_p = obj_p.interpolate(reso[0], reso[1], 1)  
        obj_p.B0[:] = 0
        obj_p.D[:] = 0
    elif sel_phantom == 'invivo':
        obj_p = mr0.VoxelGridPhantom.brainweb("./data/subject05.npz")
        obj_p = obj_p.interpolate(reso[0], reso[1], 32).slices([16])    # select center slice
        # obj_p = obj_p.interpolate(reso[0], reso[1], reso[2])          # can be used for 3D simulation
    else:
        print('Select proper phantom')
    obj_sim = obj_p.build()

    # SIMULATE the external.seq file and add acquired signal to ADC plot
    graph=mr0.compute_graph(seq0, obj_sim, 200, 1e-3)
    signal=mr0.execute_graph(graph, seq0, obj_sim)
    reco = mr0.reco_adjoint(signal, seq0.get_kspace(), resolution=reso, FOV=(0.22, 0.22, 0.22)) # Recommended: RECO has same Reso and FOV as sequence
    # %% save results
    if save:
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