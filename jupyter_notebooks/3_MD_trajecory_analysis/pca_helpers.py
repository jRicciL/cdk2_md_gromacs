import pytraj as pyt
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

def get_coords(traj, mask, ref_conf, stride = None):
    # traj alingment
    traj = pyt.align(traj, mask = mask, 
                     ref = ref_conf)[mask]
    if stride:
        traj = traj[::stride]
    # get coordinates
    traj_coords = traj.xyz
    m, n, c     = traj_coords.shape
    mtx_coords  = traj_coords.reshape((m, n*c)) 
    return mtx_coords

def get_pca_projs(mtx_coords, df_labels, 
                  stride = None, 
                  n_components = 2):
    # Perform PCA
    if stride:
        df_labels = df_labels.iloc[::stride]
    pca = PCA(n_components = n_components)
    pca_proj = pca.fit_transform(mtx_coords)
    exp_vars = pca.explained_variance_ratio_ * 100
    df = df_labels.copy()
    df = pd.concat((df, 
                    pd.DataFrame(pca_proj,
                    columns = range(1, n_components + 1),
                                )), 
                    axis = 1)
    return df, exp_vars, pca

def plot_projection(df, color_hue, plot_title = '',
                    exp_vars = (None, None), 
                    pcs_to_plot = (1, 2),
                    alpha = 0.5,
                    linewidth = 0, ax = None, **kwargs):
    if ax == None:
        fig, ax = plt.subplots(figsize = (8,8))
    sns.scatterplot(data = df, 
                     x = pcs_to_plot[0], y = pcs_to_plot[1], 
                     alpha = alpha,
                     linewidth = linewidth,
                     ax = ax,
                     hue = color_hue, **kwargs)
    ax.set_title(plot_title)
    ax.set(xlabel = f'PC{pcs_to_plot[0]} ({exp_vars[0]:.2f}%)', 
           ylabel = f'PC{pcs_to_plot[1]} ({exp_vars[1]:.2f}%)')
    ax.axvline(0, ls = ':', color = 'gray', lw = 2.5, zorder = -10)
    ax.axhline(0, ls = ':', color = 'gray', lw = 2.5, zorder = -10)
    return ax
   

def plot_explained_variance(pca_obj, title = '', ax = None):
    PC_values = np.arange(pca_obj.n_components_)
    if not ax:
        fig, ax = plt.subplots(1, 1, figsize = (10, 6))
    sns.barplot(x = PC_values + 1, 
                y = pca_obj.explained_variance_ratio_,
                palette = 'Blues_r', ax = ax)
    sns.lineplot(x = PC_values, 
                 y = pca_obj.explained_variance_ratio_.cumsum(),
                color = 'red', marker = 'o', ax = ax)
    ax.set_title(f'Scree Plot:\n{title}')
    ax.set_xlabel('Principal Components')
    ax.set_ylim((0, 1))
    ax.grid()
    ax.set_ylabel('Explained Variance')
    return ax