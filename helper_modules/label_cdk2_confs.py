import pytraj 
import numpy as np

def get_distance_ij(traj_obj: pytraj.Trajectory, 
                    atom_mask_i: str, 
                    atom_mask_j: str):
    '''
    Get the euclidean distance between two i and j atoms along a trajectory.
    Example: (traj, ':33@NZ :51@CD')
    '''
    distances = pytraj.distance(traj_obj, f'{atom_mask_i} {atom_mask_j}')
    return distances

def get_phi_angle(traj_obj: pytraj.Trajectory, 
                  residue: int):
    ''' Evaluates the angle of a dihedral given a 
        Trajectory object and a residue number'''
    phi = pytraj.calc_phi(traj_obj, 
                       resrange = str(residue))
    # phi.values is an array of n times n_frames
    return phi.to_ndarray()

def get_geom_center_distance(traj_obj: pytraj.Trajectory, 
                             resmask_1: str, 
                             resmask_2: str):
    ''' Evaluates the distance between the geometry center of two sets of atoms 
    given a Trajectory object and two atom groups defined by two Amber atom masks'''
    center_1 = pytraj.calc_center_of_geometry(traj_obj, resmask_1)
    center_2 = pytraj.calc_center_of_geometry(traj_obj, resmask_2)
    center_distance = np.linalg.norm( center_1 - center_2, axis = 1)
    return center_distance

def secondary_struc(traj_obj: pytraj.Trajectory, 
                    resmask: str):
    ''' Returns the secondary structure string of a given subsecuence 
    and a Trajectory object'''
    dssp_traj = pytraj.dssp_analysis.dssp(traj_obj,
                                       resmask)
    # dssp_traj[1] has n_frames, n_residues shape
    return dssp_traj[1]

def label_cdk2_conformations(traj_obj, 
                             saltbridge_cutoff   = 7.0,
                             dfg_angle_cutoff    = 110.0,
                             aC_b4b5_dist_cutoff = 14.5,
                             dist_PHE146_GLY11   = 12.0
                             ):

    # Evaluates if the salt bridge between LYS33-GLU51 exists
    salt_bridge_33_51 = get_distance_ij(traj_obj, ':33@NZ', ':51@CD')
    is_salt_bridge = salt_bridge_33_51 < saltbridge_cutoff
    
    # Evaluates if the structure is DFG-OUT
    # The phi angle between ASP145 and PHE146 is above abs(135.0)
    d145_angle        = np.absolute(get_phi_angle(traj_obj, residue = 145))
    phe146_gly11_dist = get_distance_ij(traj_obj, ':146@CG', ':11@CA')
    is_dfg_out = (d145_angle > dfg_angle_cutoff) & \
                 (phe146_gly11_dist < dist_PHE146_GLY11)
    
    # Evaluates if the structure is Open
    # Measures the distance between geometric centers of
    # the αC-helix (residues 46 to 57) and the β4-β5 (66-72, 75-81) sheets 
    aC_b4b5_distances = get_geom_center_distance(
                            traj_obj, 
                            ':46-57', ':66-72,75-81')
    is_open = aC_b4b5_distances > aC_b4b5_dist_cutoff
    
    # Evaluates if the A-loop has a secondary structure (helix)
    # COMMENTED: Initially used to differentiate between inact_a and inact_b confs
    # sec_struc = secondary_struc(traj_obj, ':144-157')
    # is_helix = np.char.count(sec_struc, '0').sum(axis=1) < min_num_res_helix
    
    # descriptors, the first three descriptors
    descriptors = np.vstack((is_salt_bridge, is_dfg_out, is_open)).T
    
    # Labels
    active     = (np.array([ 1, 0, 0 ]), "active")
    dfg_out    = (np.array([ 1, 1, 0 ]), "dfg_out")
    inact_open = (np.array([ 0, 0, 1 ]), "inact_open")
    inact_src  = (np.array([ 0, 0, 0 ]), "inact_src")
    
    # Labels
    conf_labels = np.repeat("undefined", traj_obj.n_frames)
    
    for conformation in [active, dfg_out, inact_open, inact_src]:
        index_conf = np.where(( descriptors == conformation[0] ).all(axis=1))[0]
        conf_labels[index_conf] = conformation[1]
    # Now that inactives has bee identified, we've to diferentiate between inactives A and B
    # thus, if the label at index j is inactive_b and the index i is true for is_helix, 
    # change the label to inact_a
    # COMMENTED: Initially used to differentiate between inact_a and inact_b confs
    # conf_labels = ["inact_a" if i and j == "inact_b" else j for i, j in zip(is_helix, conf_labels) ]
    
    return conf_labels