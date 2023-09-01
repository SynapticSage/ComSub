import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

def perform_umap(df_matrix:pd.DataFrame, n_neighbors=15, min_dist=0.1,
                 n_components=2, metric='euclidean'):
    """
    Perform UMAP on the input matrix.
    Parameters:
    - matrix (DataFrame or numpy array): The input data matrix.
    - n_neighbors (int): UMAP parameter. Size of local neighborhood.
    - min_dist (float): UMAP parameter. Minimum distance between points in the low-dimensional representation.
    - n_components (int): Number of dimensions in the low-dimensional space.
    - metric (str): Metric used to compute distances.
    Returns:
    - DataFrame: UMAP-transformed data.
    """
    import umap
    if any([x in df_matrix.columns for x in index]):
        df_matrix = df_matrix.set_index(index)
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, metric=metric)
    print(f"Performing UMAP with {n_neighbors} neighbors, {min_dist} min_dist, {n_components} components, and {metric} metric...")
    embedding = reducer.fit_transform(df_matrix.values)
    print(f"...done. UMAP embedding shape: {embedding.shape}")
    return pd.DataFrame(embedding, columns=[f"dim_{i}" for i in range(1, n_components + 1)],
                        index=df_matrix.index)

print("Preping data...")
df_matrix = prep_uv_melt(df_clean)
df_matrix_projected = prep_uv_melt(df_clean, project=True)
df_matrix.dropna(inplace=True)
um =  perform_umap(df_matrix)
um = um.reset_index()
# um_anim = perform_umap(df_matrix_projected)
# um_anim = um_anim.reset_index()

def plot_umap(matrix:pd.DataFrame, row=None, col=None, hue=None, sample=None,
              col_wrap=None, **kws):
    """
    Plot the UMAP-transformed data.
    Parameters:
    - matrix (DataFrame): UMAP-transformed data.
    - row (str, optional): Column in matrix to use for row-wise subplotting.
    - col (str, optional): Column in matrix to use for column-wise subplotting.
    - hue (str, optional): Column in matrix to use for color encoding.
    Returns:
    - FacetGrid: Seaborn FacetGrid object with scatter plots.
    """
    if sample:
        matrix = matrix.sample(sample)
    g = sns.FacetGrid(matrix, row=row, col=col, hue=hue, height=5, aspect=1,
                      col_wrap=col_wrap)
    g.map(plt.scatter, "dim_1", "dim_2", **kws)
    g.add_legend()
    return g

plot_umap(um, hue='animal', sample=10_000, alpha=0.1)
# plot_umap(um_anim, hue='animal', sample=10_000, alpha=0.1)

plt.close('all')
plot_umap(um, hue='genH', col="animal", sample=10_000, alpha=0.1, col_wrap=3)
# plot_umap(um_anim, hue='genH', col="animal", sample=10_000, alpha=0.1, col_wrap=3)

plt.close('all')
plot_umap(um, hue='genH_highlow', col="animal", sample=10_000, alpha=0.1, col_wrap=3)
# plot_umap(um_anim, hue='genH_highlow', col="animal", sample=10_000, alpha=0.1, col_wrap=3)

plt.close('all')
plot_umap(um, hue='genH_highlow', row="patterns", col="animal", sample=10_000, alpha=0.1)

um =  perform_umap(df_matrix, n_components=3)
um = um.reset_index()
# um_anim = perform_umap(df_matrix_projected, n_components=3)
# um_anim = um_anim.reset_index()

def plot_umap_3d(matrix, row=None, col=None, hue=None, sample=None,
                 col_wrap=None, **kws):
    """
    Plot the UMAP-transformed data in 3D.
    Parameters:
    - matrix (DataFrame): UMAP-transformed data.
    - row (str, optional): Column in matrix to use for row-wise subplotting.
    - col (str, optional): Column in matrix to use for column-wise subplotting.
    - hue (str, optional): Column in matrix to use for color encoding.
    Returns:
    - FacetGrid: Seaborn FacetGrid object with scatter plots.
    """
    from mpl_toolkits.mplot3d import Axes3D
    if 'magnitude_u' in matrix.columns:
        matrix = matrix.drop(columns=['magnitude_u'])
    if 'magnitude_v' in matrix.columns:
        matrix = matrix.drop(columns=['magnitude_v'])
    if sample:
        matrix = matrix.sample(sample)
    # Create a grid of plots
    g = sns.FacetGrid(matrix, row=row, col=col, hue=hue, height=5, aspect=1,
                      col_wrap=col_wrap, subplot_kws={'projection': '3d'})
    hue_mapping = {unique_value: index for index, unique_value in
                   enumerate(matrix[hue].unique())}
    hue_order = list(hue_mapping.keys())
    for ax, data in g.facet_data():
        palette = sns.utils.get_color_cycle()
        c = g._get_palette(data, hue, hue_order, palette)[ax[-1]]
        print("ax", ax)
        if g.axes.ndim == 1:
            AX = g.axes[ax[1]]
        else:
            AX = g.axes[ax[0], ax[1]]
        AX.scatter(data["dim_1"], data["dim_2"], data["dim_3"], c=c, alpha=0.1)
    # # For each axis object, set it to 3D and plot the 3D scatter plot
    # for ax in g.axes.flat:
    #     ax.remove()
    #     ax = g.fig.add_subplot(ax, projection='3d')
    #     ax.scatter(matrix["dim_1"], matrix["dim_2"], matrix["dim_3"], c=matrix[hue])
    g.add_legend()
    return g

plt.close('all')
plot_umap_3d(um, hue='animal', sample=10_000, alpha=0.1)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=animal.pdf'), dpi=300)

plt.close('all')
plot_umap_3d(um, hue='genH', col="animal", sample=10_000, alpha=0.4,
             col_wrap=3)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_col=animal.pdf'), dpi=300)

plt.close('all')
plot_umap_3d(um, hue='genH_highlow', col="animal", sample=5_000, alpha=1,
             col_wrap=3)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_col=animal.pdf'), dpi=300)

plt.close('all')
plot_umap_3d(um, hue='genH_highlow', row="patterns", col="animal",
             sample=10_000, alpha=0.1)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_row=patterns_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_hue=genH_highlow_row=patterns_col=animal.pdf'), dpi=300)

# plot_umap_3d(um_anim, hue='genH_highlow', col="animal", sample=10_000,
alpha=0.9, col_wrap=3)
plt.savefig(os.path.join(figfolder,'umap_3d_rotanim_hue=genH_highlow_col=animal.png'), dpi=300)
plt.savefig(os.path.join(figfolder,'umap_3d_rotanim_hue=genH_highlow_col=animal.pdf'), dpi=300)

# plot_umap_3d(um_anim, hue='animal', sample=10_000, alpha=0.1)
# plot_umap_3d(um_anim, hue='genH', col="animal", sample=10_000, alpha=0.1,
#              col_wrap=3)
