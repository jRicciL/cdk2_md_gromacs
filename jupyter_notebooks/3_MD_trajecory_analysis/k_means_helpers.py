import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

def get_system_chuncks(df_all_systems, 
                       used_stride):
    n_systems = len(df_all_systems.groupby(['Conf.', 'MD-Protocol']))
    reps = df_all_systems.groupby(['Conf.', 'MD-Protocol', 'Rep.']).size()
    n_reps = len(reps)
    n_confs_per_rep = reps[0]
     
    # Replace the 'Conf. type' column 
    # to the name of the reference conf.
    df_ = df_all_systems.copy()

    # Add a `frame_numbers` column to identify
    # the # of frame inside each trajectory
    tot_confs_per_rep = n_confs_per_rep * used_stride
    frame_numbers = list(range(0, tot_confs_per_rep, used_stride)) \
                            * (n_reps)
    assert len(frame_numbers) == df_.shape[0]
    # Add `frame_number` as a new column
    df_['frame_number'] = frame_numbers
    group_columns = ['Conf.', 'MD-Protocol']
    idx_columns   = ['Conf.', 'MD-Protocol', 
                     'Rep.', 'frame_number']
    # Now, split the dataframe into chuncks, each chunk belongs
    # to each of the 20 systems
    dfs_sytems = [chunck.set_index(idx_columns)
                  for sys, chunck in 
                  df_.groupby(group_columns)]
    return dfs_sytems

# For each system, perform k Means
def run_kMeans(X, n_clusters, random_state = 42):
    
    km = KMeans(n_clusters = n_clusters, 
                random_state = random_state)\
                .fit(X)
    X['k_label'] = km.labels_
    return X, km

def get_medoid_confs_per_chunk(X, km):
    medoid_confs = []
    for label in X['k_label'].unique():
        X_clust = X.query(f'k_label == {label}')
        dims    = X_clust.drop('k_label', axis = 1).values
        # Get the centroid
        centroid = km.cluster_centers_[label, :]
        # Get the closest conf to the centroid, ie., medoid
        med_idx  = np.argmin(
                        ((dims - centroid)**2).sum(axis = 1)
                       )
        med_conf_row = X_clust.iloc[med_idx]
        medoid_confs.append(med_conf_row)
    df_medoids = pd.concat(medoid_confs, axis = 1).T
    return df_medoids

def get_medoid_confs(df_all_systems, used_stride, 
                     n_clusters, random_state = 42):
    sys_chunks = get_system_chuncks(df_all_systems, 
                                    used_stride = used_stride)
    all_medoids = []
    for chunk in sys_chunks:
        X, km = run_kMeans(X = chunk, 
                   n_clusters = n_clusters,
                   random_state = random_state)
        df_medoids = get_medoid_confs_per_chunk(X, km)
        all_medoids.append(df_medoids)
        
    df_all_medoids = pd.concat(all_medoids, 
                               axis = 0).sort_index()
    df_all_medoids.index.names = ['Conf.', 'MD-Protocol', 
                                  'Rep.', 'frame_number']
    return df_all_medoids.reset_index()