import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import pandas as pd
import random
import plotly.io as pio
from matplotlib.lines import Line2D
import seaborn as sns
from sklearn.metrics import adjusted_rand_score
import scanpy as sc
import warnings
from natsort import natsort_keygen



###------------setting some local variables

qc_metrics = ['total_counts', 
    'n_genes_by_counts', 
    'pct_counts_in_top_50_genes', 
    'pct_counts_in_top_500_genes',
    'ratio_spliced', 
    'pct_counts_mito', 
    'pct_counts_ribo',
    'pct_counts_hemo',
    'doublet_score']

marker_guttube = [
"gene_a"
]

marker_ve = [ 
"gene_b"]
    
marker_ve_liver = [     
"gene_c"]

marker_mesoderm = [
"gene_d"]

marker_ve_query = [
"gene_e"
]



marker_de = [
"gene_f"
]

marker_dict = {
    'mesoderm' : marker_mesoderm,
    'guttube' : marker_guttube, 
    've' : marker_ve, 
    've_liver' : marker_ve_liver, 
    've_query' : marker_ve_query, 
    'de' : marker_de
}



#####--------------------------------------
def visualize_library_size(adata, colormap, color_by=None, default_color = 'grey', threshold=2000, width=8, height=8):
    data = adata.copy()
    plt.figure(figsize=(width, height))

    # Determine groups if color_by is given
    groups = data.obs[color_by].unique() if color_by else None

    def plot_hist(ax, metric, threshold=None):
        ax.set_xlabel(metric)
        ax.set_ylabel('Count')

        group_iter = [None] if groups is None else groups

        for group in group_iter:
            if group is None:
                mask = slice(None)  # all data
                label = None
                color = default_color
            else:
                mask = data.obs[color_by] == group
                label = str(group)
                color = colormap.get(group, default_color)

            ax.hist(
                data.obs.loc[mask, metric],
                bins=100,
                alpha=0.4,
                label=label,
                color=color)

        if groups is not None:
            ax.legend(title=color_by)

        if threshold is not None:
            ax.axvline(x=threshold, color='red', linestyle="--", linewidth=1)

    # Metrics and thresholds
    metrics = [
        ('log1p_total_counts', np.log(threshold + 1)),
        ('total_counts', threshold),
        ('log1p_n_genes_by_counts', None),
        ('n_genes_by_counts', None),
        ('pct_counts_in_top_50_genes', None)
    ]

    # Plot histograms
    for i, (metric, metric_threshold) in enumerate(metrics):
        ax = plt.subplot(3, 2, i + 1)
        plot_hist(ax, metric, metric_threshold)

    # Scatter plot
    ax6 = plt.subplot(3, 2, 6)

    if color_by:
        # Precompute colors once
        point_colors = data.obs[color_by].map(lambda g: colormap.get(g, default_color))

        ax6.scatter(
            x=data.obs['log1p_total_counts'],
            y=data.obs['log1p_n_genes_by_counts'],
            c=point_colors,
            s=3,
        )

        # Build legend manually
        legend_elements = [
            Line2D([0], [0], marker='o', color='w',
                label=str(group),
                markerfacecolor=colormap.get(group, default_color),
                markersize=6)
            for group in groups
        ]

        ax6.legend(handles=legend_elements, title=color_by)
    else:
        ax6.scatter(
            x=data.obs['log1p_total_counts'],
            y=data.obs['log1p_n_genes_by_counts'],
            s=3,
            color=default_color,
        )

    ax6.set_xlabel('log1p_total_counts')
    ax6.set_ylabel('log1p_n_genes_by_counts')
    ax6.axvline(x=np.log(threshold + 1), color='red', linestyle="--", linewidth=1)

    plt.tight_layout()


##-----------------------------------------
def visualize_mito(
    adata,
    colormap,
    color_by=None,
    default_color='gray',
    obs_counts=False,
    save_plot="no",
    file_path='./mito_content.png',
    threshold=25,
    min_total_counts_per_cell=2000,
    width=8,
    height=8
):
    data = adata.copy()
    plt.figure(figsize=(width, height))

    # Determine groups if color_by is given
    groups = data.obs[color_by].unique() if color_by else None

    def plot_hist(ax, metric, threshold=None):
        ax.set_xlabel(metric)
        ax.set_ylabel('Count')

        group_iter = [None] if groups is None else groups

        for group in group_iter:
            if group is None:
                mask  = slice(None)
                label = None
                color = default_color
            else:
                mask  = data.obs[color_by] == group
                label = str(group)
                color = colormap.get(group, default_color)

            counts, bins, _ = ax.hist(
                data.obs.loc[mask, metric],
                bins=100,
                alpha=0.4,
                label=label,
                color=color
            )

            # add counts on top of bars if requested
            if obs_counts:
                for cnt, left, right in zip(counts, bins[:-1], bins[1:]):
                    if cnt > 0:
                        ax.text(
                            (left + right) / 2, cnt, str(int(cnt)),
                            ha='center', va='bottom', fontsize=7, rotation=45
                        )

        if groups is not None:
            ax.legend(title=color_by)

        if threshold is not None:
            ax.axvline(x=threshold, color='red', linestyle='--', linewidth=1)

    # 1) Histogram of mitochondrial %
    ax1 = plt.subplot(3, 1, 1)
    plot_hist(ax1, 'pct_counts_mito', threshold)

    # 2) Total counts vs mito%
    ax2 = plt.subplot(3, 1, 2)
    if groups is None:
        ax2.scatter(
            data.obs['log1p_total_counts'],
            data.obs['pct_counts_mito'],
            s=3,
            color=default_color,
            alpha=0.5
        )
    else:
        for group in groups:
            mask = data.obs[color_by] == group
            ax2.scatter(
                data.obs.loc[mask, 'log1p_total_counts'],
                data.obs.loc[mask, 'pct_counts_mito'],
                s=3,
                alpha=0.5,
                label=str(group),
                color=colormap.get(group, default_color)
            )
        ax2.legend(title=color_by)

    ax2.set_xlabel('log1p_total_counts')
    ax2.set_ylabel('pct_counts_mito')
    ax2.axhline(y=threshold, color='red', linestyle='--', linewidth=1)
    ax2.axvline(x=np.log(min_total_counts_per_cell + 1), color='red', linestyle='--', linewidth=1)

    # 3) Genes detected vs mito%
    ax3 = plt.subplot(3, 1, 3)
    if groups is None:
        ax3.scatter(
            data.obs['log1p_n_genes_by_counts'],
            data.obs['pct_counts_mito'],
            s=3,
            color=default_color,
            alpha=0.5
        )
    else:
        for group in groups:
            mask = data.obs[color_by] == group
            ax3.scatter(
                data.obs.loc[mask, 'log1p_n_genes_by_counts'],
                data.obs.loc[mask, 'pct_counts_mito'],
                s=3,
                alpha=0.5,
                label=str(group),
                color=colormap.get(group, default_color)
            )
        ax3.legend(title=color_by)

    ax3.set_xlabel('log1p_n_genes_by_counts')
    ax3.set_ylabel('pct_counts_mito')
    ax3.axhline(y=threshold, color='red', linestyle='--', linewidth=1)

    plt.tight_layout()

    if save_plot.lower() == 'yes':
        plt.savefig(file_path, dpi=300)
        plt.close()

## -------------------------------------------------------------------------

def visualize_ribo(
    adata,
    colormap,
    color_by=None,
    default_color='gray',
    width=8,
    height=8
):
    data = adata.copy()
    plt.figure(figsize=(width, height))

    # Determine groups if color_by is given
    groups = data.obs[color_by].unique() if color_by else None

    def plot_hist(ax, metric):
        ax.set_xlabel(metric)
        ax.set_ylabel('Count')

        group_iter = [None] if groups is None else groups
        for group in group_iter:
            if group is None:
                mask  = slice(None)
                label = None
                color = default_color
            else:
                mask  = data.obs[color_by] == group
                label = str(group)
                color = colormap.get(group, default_color)

            ax.hist(
                data.obs.loc[mask, metric],
                bins=100,
                alpha=0.4,
                label=label,
                color=color
            )

        if groups is not None:
            ax.legend(title=color_by)

    def plot_scatter(ax, x_metric, y_metric):
        ax.set_xlabel(x_metric)
        ax.set_ylabel(y_metric)

        if groups is None:
            ax.scatter(
                data.obs[x_metric],
                data.obs[y_metric],
                s=3,
                color=default_color,
                alpha=0.5
            )
        else:
            for group in groups:
                mask = data.obs[color_by] == group
                ax.scatter(
                    data.obs.loc[mask, x_metric],
                    data.obs.loc[mask, y_metric],
                    s=3,
                    label=str(group),
                    color=colormap.get(group, default_color),
                    alpha=0.5
                )
            ax.legend(title=color_by)

    # 1) Histogram of ribosomal %
    ax1 = plt.subplot(3, 1, 1)
    plot_hist(ax1, 'pct_counts_ribo')

    # 2) Total counts vs ribo%
    ax2 = plt.subplot(3, 1, 2)
    plot_scatter(ax2, 'log1p_total_counts', 'pct_counts_ribo')

    # 3) Genes detected vs ribo%
    ax3 = plt.subplot(3, 1, 3)
    plot_scatter(ax3, 'log1p_n_genes_by_counts', 'pct_counts_ribo')

    plt.tight_layout()


## -----------------------------------------

def visualize_splicing(
    adata,
    colormap,
    color_by=None,
    default_color='gray',
    obs_counts=False,
    width=6,
    height=3,
):
    data = adata.copy()

    # Determine groups if color_by is given
    groups = data.obs[color_by].unique() if color_by else None

    def plot_hist(ax, metric, annotate=False):
        ax.set_xlabel(metric)
        ax.set_ylabel('Count')

        group_iter = [None] if groups is None else groups
        for group in group_iter:
            if group is None:
                mask  = slice(None)
                label = None
                color = default_color
            else:
                mask  = data.obs[color_by] == group
                label = str(group)
                color = colormap.get(group, default_color)

            counts, bins, _ = ax.hist(
                data.obs.loc[mask, metric],
                bins=100,
                alpha=0.4,
                label=label,
                color=color
            )

            if annotate and obs_counts:
                for cnt, left, right in zip(counts, bins[:-1], bins[1:]):
                    if cnt > 0:
                        ax.text(
                            (left + right) / 2, cnt, str(int(cnt)),
                            ha='center', va='bottom', fontsize=7, rotation=45
                        )

        if groups is not None:
            ax.legend(title=color_by)

    # 1) ratio_spliced histogram (with optional counts)
    plt.figure(figsize=(width, height))
    ax1 = plt.gca()
    plot_hist(ax1, 'ratio_spliced', annotate=True)
    plt.tight_layout()

    # 2) spliced / unspliced / ambiguous in one row
    metrics = ['spliced', 'unspliced', 'splicing_ambiguous']
    plt.figure(figsize=(width, height))
    for i, metric in enumerate(metrics, start=1):
        ax = plt.subplot(1, 3, i)
        plot_hist(ax, metric, annotate=False)
    plt.tight_layout()

    
        
######---------------------------------

def visualize_hemo(
    adata,
    colormap,
    color_by=None,
    default_color='gray',
    obs_counts=False,
    width=6,
    height=3
):
    data = adata.copy()
    plt.figure(figsize=(width, height))

    # Determine groups if color_by is given
    groups = data.obs[color_by].unique() if color_by else None

    def plot_hist(ax, metric):
        ax.set_xlabel(metric)
        ax.set_ylabel('Count')

        group_iter = [None] if groups is None else groups
        for group in group_iter:
            if group is None:
                mask  = slice(None)
                label = None
                color = default_color
            else:
                mask  = data.obs[color_by] == group
                label = str(group)
                color = colormap.get(group, default_color)

            counts, bins, _ = ax.hist(
                data.obs.loc[mask, metric],
                bins=100,
                alpha=0.4,
                label=label,
                color=color
            )

            if obs_counts:
                for cnt, left, right in zip(counts, bins[:-1], bins[1:]):
                    if cnt > 0:
                        ax.text(
                            (left + right) / 2, cnt, str(int(cnt)),
                            ha='center', va='bottom', fontsize=7, rotation=45
                        )

        if groups is not None:
            ax.legend(title=color_by)

    # Single-panel histogram of hemoglobin counts %
    ax = plt.gca()
    plot_hist(ax, 'pct_counts_hemo')

    plt.tight_layout()


### -----------------------------------------------------------

def visualize_doublets(
    adata,
    colormap,
    color_by=None,
    default_color='gray',
    obs_counts=False,
    width=6,
    height=3
):
   
    data = adata.copy()
    plt.figure(figsize=(width, height))

    # Determine grouping categories
    groups = data.obs[color_by].unique() if color_by else None

    ax = plt.gca()
    ax.set_xlabel('doublet_score')
    ax.set_ylabel('Count')

    # Loop over each group (or all cells if no grouping)
    for group in ([None] if groups is None else groups):
        if group is None:
            mask  = slice(None)
            label = None
            color = default_color
        else:
            mask  = data.obs[color_by] == group
            label = str(group)
            color = colormap.get(group, default_color)

        counts, bins, _ = ax.hist(
            data.obs.loc[mask, 'doublet_score'],
            bins=100,
            alpha=0.4,
            label=label,
            color=color
        )

        # Annotate counts on top of bars if requested
        if obs_counts:
            for cnt, left, right in zip(counts, bins[:-1], bins[1:]):
                if cnt > 0:
                    ax.text(
                        (left + right) / 2, cnt, str(int(cnt)),
                        ha='center', va='bottom', fontsize=7, rotation=45
                    )

    # Add legend only if grouping
    if groups is not None:
        ax.legend(title=color_by)

    plt.tight_layout()


def visualize_genes(adata, threshold = 50, width = 8, height = 6):
    adata = adata.copy()

    plt.figure(figsize=(width,height))
    
    plt.subplot(3, 1, 1)
    plt.hist(adata.var['n_cells_by_counts'], bins=100, alpha=0.4)
    plt.xlabel('n_cells_by_counts') 
    plt.ylabel('Count')
    plt.legend(title='Sample')
    plt.ylim([0, 1000])
    plt.axvline(x=threshold, color = 'red', linestyle = "--", linewidth = 1)
    
    plt.subplot(3, 1, 2)
    plt.hist(adata.var['log1p_mean_counts'], bins=100, alpha=0.4)
    plt.xlabel('log1p_mean_counts') 
    plt.ylabel('Count')
    plt.legend(title='Sample')
    plt.ylim([0, 1000])
    
    plt.subplot(3, 1, 3)
    plt.hist(adata.var['log1p_total_counts'], bins=100, alpha=0.4)
    plt.xlabel('log1p_total_counts') 
    plt.ylabel('Count')
    plt.legend(title='Sample')
    plt.ylim([0, 1000])
    
    plt.tight_layout()
    plt.show()


###-------------------------------------------

def latent_space_interactive(
    data,
    title,
    voi,
    reduction_method = 'UMAP',
    reduction_column = 'X_umap',
    colormap=None,  # Make colormap optional
    meta_fields=None,  # Allow empty or no meta fields
    height=600,
    width=800,
    background_color='gray',
    background_opacity=0.2,
    background_size=3,
    foreground_size=5,
    missing_color = 'gold',
    file_format = 'ipynb'
):
    """
    Function to create a UMAP plot with foreground and background scatter traces.
    
    Args:
        data: Anndata object containing the dataset.
        title: Title for the plot.
        reduction_method: Name of the dimensionality reduction method (e.g., 'UMAP').
        reduction: Key for the reduction embedding in `data.obsm`.
        voi: Variable of interest (column in `data.obs`).
        colormap: Optional dictionary mapping unique values of `voi` to colors.
        meta_fields: List of metadata fields to include in hover information.
        height: Height of the plot.
        width: Width of the plot.
        background_color: Color for background points.
        background_opacity: Opacity for background points.
        background_size: Size of background points.
        foreground_size: Size of foreground points.
        
    Returns:
        Plotly Figure object.
    """
    ## setting the rendering right
    if file_format == 'ipynb': 
        pio.renderers.default = "iframe_connected"
    elif file_format == 'html': 
        pio.renderers.default = "notebook"
    else: 
        raise ValueError('Format type not supported')

    meta_fields = meta_fields or []  # Default to empty list if not provided

    # Function to generate a random hex color
    random_color = lambda: f'#{random.randint(0, 0xFFFFFF):06x}'

    # Set up plotting dataframe
    umap_df = pd.DataFrame(data.obsm[reduction_column], columns=[f'{reduction_method}1', f'{reduction_method}2'])
    for meta in [voi] + meta_fields:
        umap_df[meta] = data.obs[meta].values

    # Sort umap_df alphanumerically by voi
    umap_df = umap_df.sort_values(by=voi, key=lambda col: col.str.lower().map(natsort_keygen()))

    # Initialize a new Figure
    fig = go.Figure()

    # Create the background trace with go.Scattergl for efficiency with large datasets
    traces_back = go.Scattergl(
        x=umap_df[f'{reduction_method}1'],
        y=umap_df[f'{reduction_method}2'],
        mode='markers',
        marker=dict(
            color=background_color,   # Set color to gray
            size=background_size,     # Smaller point size
            opacity=background_opacity  # Semi-transparent points
        ),
        showlegend=False,  # Hide legend for background trace
        hoverinfo='skip'
    )

    # Add background trace to the figure
    fig.add_trace(traces_back)

    # Define helper function for constructing hover text
    def construct_hover_text(row, fields):
        return "<br>".join(f"{field}: {row[field]}" for field in fields)

    # Create and add foreground traces for each unique value in 'voi'
    for value in umap_df[voi].unique():
        # Filter data for the current unique value
        value_data = umap_df[umap_df[voi] == value]

        # Determine the color for the current value
        if colormap:
            value_color = colormap.get(value, missing_color)
        else:
            value_color = random_color()  # Generate random color if no colormap or key missing

        # Construct hover text for each point
        hover_text = value_data.apply(lambda row: construct_hover_text(row, [voi] + meta_fields), axis=1)

        # Create a trace for this value
        trace = go.Scattergl(
            x=value_data[f'{reduction_method}1'],
            y=value_data[f'{reduction_method}2'],
            mode='markers',
            marker=dict(
                size=foreground_size,
                color=value_color
            ),
            name=str(value),  # Legend name for this category
            hoverinfo='text',
            hovertext=hover_text
        )

        # Add the trace to the figure
        fig.add_trace(trace)

    # Update layout parameters
    fig.update_layout(
        title=title,
        xaxis_title=f'{reduction_method}1',
        yaxis_title=f'{reduction_method}2',
        legend_title='',
        xaxis=dict(
            showgrid=False,
            showticklabels=False,
        ),
        yaxis=dict(
            showgrid=False,
            showticklabels=False,
        ),
        plot_bgcolor='white',
        height=height,
        width=width
    )

    return fig

#### pairplot ----------------

def qc_metrics_pairplot(adata, colormap = None, color_by = None, metrics = qc_metrics, height = 2): 
    metrics_tmp = metrics.copy()
    
    if color_by is not None:
        metrics_tmp.append(color_by)

    # Extract relevant metrics_tmp from adata.obs
    data = adata.obs[metrics_tmp].copy()
        
    # Drop rows with missing values to avoid plotting issues
    data.dropna(inplace=True)
    

    if color_by is not None:
        g = sns.pairplot(data, corner=False, hue=color_by, palette=colormap,
                         plot_kws={'alpha': 0.5, 's': 10}, diag_kind='kde',
                         height=height, aspect=1)
    else:
        g = sns.pairplot(data, corner=False,
                         plot_kws={'alpha': 0.5, 's': 10}, diag_kind='kde',
                         height=height, aspect=1)

    # <--- Only these lines are added ---
    for ax in g.axes.flatten():
        if ax is not None:
            ax.set_xlabel(ax.get_xlabel(), fontweight='bold')
            ax.set_ylabel(ax.get_ylabel(), fontweight='bold')
    
    plt.suptitle("QC metrics correlation", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()



#####----------------------------------------------

def qc_metrics_corr_cov(adata, metrics = qc_metrics, method='spearman'):
    metrics_tmp = metrics.copy()
    
    # Extract and clean the data
    data = adata.obs[metrics_tmp].copy()
    data.dropna(inplace=True)

    # Compute correlation and covariance matrices
    corr_matrix = data.corr(method=method, numeric_only = True)
    cov_matrix = data.cov(numeric_only = True)
    cov_matrix = np.sign(cov_matrix) * np.log2(np.abs(cov_matrix) + 1)

    # Plot side-by-side
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Correlation heatmap
    sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", center=0, ax=axes[0])
    axes[0].set_title(f"{method.capitalize()} Correlation Matrix")

    # Covariance heatmap
    sns.heatmap(cov_matrix, annot=True, fmt=".2f", cmap="coolwarm", ax=axes[1])
    axes[1].set_title("Log2p1 Covariance Matrix")

    plt.tight_layout()
    plt.show()

###### -------------- leiden cluster resolution

from sklearn.metrics import adjusted_rand_score

def clustering_resolutions(
    adata, 
    flavour = 'leiden',
    resolution_end = 1, 
    resolution_start = 0, 
    resolution_stepsize = 0.1, 
    n_neighbors = 30, 
    neighbors_key = 'neighbors'):

    adata = adata.copy()
    #resolution_stepsize = np.negative(resolution_stepsize)
    resolutions = [round(r * 0.01, 2) for r in range(int(resolution_start*100), int(resolution_end*100), int(resolution_stepsize*100))] 

    cluster_label_list = []
    for resolution in resolutions:
        if flavour == 'leiden':
            sc.tl.leiden(adata, resolution=resolution, neighbors_key=neighbors_key, key_added=f'leiden_res_{resolution}', flavor='leidenalg', use_weights = True)
            cluster_label_list.append(adata.obs[f'leiden_res_{resolution}'].values)
        elif flavour == 'phenograph_leiden': 
            sc.external.tl.phenograph(adata, clustering_algo='leiden', k=n_neighbors, jaccard=True, primary_metric='euclidean', resolution_parameter=resolution, n_jobs = -1, copy = False)
            cluster_label_list.append(adata.obs[f'pheno_leiden'].values)
        else: 
            print('Choose a valid flavour - either leiden or phenograph_leiden')


        # Loop through both the resolution and cluster labels at the same time
    for resolution, cluster in zip(resolutions, cluster_label_list):
        # Get unique values and their counts
        unique, counts = np.unique(cluster, return_counts=True)
        # Print the resolution and corresponding cluster counts
        print(f"Resolution: {resolution}\t\tNumber of clusters {len(unique)}")

            
    # Now compute the ARI matrix for all pairs of resolutions
    rand_index_mat = np.zeros((len(resolutions), len(resolutions)))
    
    # Compute ARI for each pair of resolutions
    for i in range(len(resolutions)):
        for j in range(len(resolutions)):
            rand_index_mat[i, j] = adjusted_rand_score(cluster_label_list[i], cluster_label_list[j])
    
    # Plotting the ARI heatmap
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    im = ax.imshow(rand_index_mat, cmap='Reds', vmin=0, vmax=1)
    fig.colorbar(im, ax=ax)
    
    # Set ticks and labels for the axes
    ax.set_xticks(range(len(resolutions)))
    ax.set_xticklabels(resolutions, rotation=90)
    ax.set_yticks(range(len(resolutions)))
    ax.set_yticklabels(resolutions)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('Resolution')
    ax.set_title('Adjusted Rand Index between Leiden resolutions')
    
    # Add text annotations inside the heatmap rectangles
    for i in range(len(resolutions)):
        for j in range(len(resolutions)):
            ax.text(j, i, f"{rand_index_mat[i, j]:.2f}", ha='center', va='center', color='black', fontsize = 'xx-small')
    
    plt.show()

def pca_variance_cumsum(adata):
    # Compute cumulative sum
    cumulative_variance = np.cumsum(adata.uns['pca']['variance_ratio'])
    # Plot as a line
    plt.plot(cumulative_variance, marker='o')
    plt.xlabel('Principal Component')
    plt.ylabel('Cumulative Explained Variance')
    plt.title('Cumulative Variance Explained by PCA')
    plt.grid(True)
    plt.show()



def dotplot_markers(adata, marker_dict = marker_dict, groupby = "leiden", out_path = ".", save = False):
    adata = adata.copy()
    celltypes = {}
    for key in marker_dict:
        celltypes[key] = [gene for gene in marker_dict[key] if gene in adata.var_names]
        if len(celltypes[key]) == 0: 
            print(f'No marker genes of {key} recorded in adata')
            del celltypes[key]

    categories_order = sorted(adata.obs['leiden'].unique().tolist(), key=lambda x: int(x))

    if save == True: 
        os.chdir(out_path)
        sc.pl.dotplot(adata, celltypes, groupby=groupby,use_raw = False, dendrogram = False, swap_axes = False, categories_order = categories_order, save = f'marker_genes.png')
    else: 
        sc.pl.dotplot(adata, celltypes, groupby=groupby,use_raw = False, dendrogram = False, swap_axes = False, categories_order = categories_order)


def boxplot_cluster_metrics(adata, qc_metrics = qc_metrics, cluster_column = 'leiden', height = 20, width = 10, save = False, file_path = './boxplot_cluster_metrics.png'):
    columns = [cluster_column] + qc_metrics
    data = adata.obs[columns].copy()
    # Suppress warnings
    warnings.filterwarnings("ignore")


    fig = plt.figure(figsize=(width, height))  
    for j, item in enumerate(data.columns[1:]):
        ax = fig.add_subplot(len(qc_metrics), 1, j + 1)
        sns.boxplot(
            x="leiden", y=item, hue='leiden', 
            data=data, ax=ax, 
            palette=sc.pl.palettes.zeileis_28
        )
        ax.grid(True)
        ax.legend().set_visible(False)

    plt.tight_layout()
    if save == True:
        fig.savefig(file_path, dpi=300)

    plt.show()  # Ensure all plots render


def visualize_contingency_table(adata, x_column, y_column, height = 7, width = 9): 
    data = adata.obs.copy()
    contingency_table = pd.crosstab(data[y_column], data[x_column])
    contingency_log2 = np.log2(contingency_table + 1)
    
    # Step 3: Set the figure size and plot the heatmap
    plt.figure(figsize=(width, height))  

    # Create heatmap and use 'annot' to display sampling frequencies
    heatmap = sns.heatmap(contingency_log2,         # Use log-transformed values for color mapping
                        annot=contingency_table,    # Annotate the heatmap with sampling frequencies
                        fmt="",           # No specific formatting
                        cmap="coolwarm",    # Color map for better contrast
                        linewidths=.5,      # Add gridlines for clarity
                        cbar_kws={'shrink': .8},  # Adjust the color bar size
                        annot_kws={"size": 6}) 
    colorbar = heatmap.collections[0].colorbar
    colorbar.ax.set_title('count', fontsize=12) 

    plt.xlabel(x_column, fontsize=12)
    plt.ylabel(y_column, fontsize=14)

    plt.xticks(rotation=90, fontsize=8)  
    plt.yticks(rotation=0, fontsize=12)

    plt.tight_layout()

    plt.show()

def interactive_contingency_table(
    adata, x_column, y_column,
    file_format='ipynb',
    height=900,
    width=700, 
    colorscale = 'bluered'
):
    # Set renderer
    if file_format == 'ipynb': 
        pio.renderers.default = "iframe_connected"
    elif file_format == 'html': 
        pio.renderers.default = "notebook"
    else: 
        raise ValueError('Format type not supported')

    data = adata.obs[[x_column, y_column]].copy()

    # Cross-tab
    contingency = pd.crosstab(data[y_column], data[x_column])
    contingency = contingency.loc[
        sorted(contingency.index, key=natsort_keygen()),
        sorted(contingency.columns, key=natsort_keygen())
    ]

    # Compute log2 values for color
    contingency_log = np.log2(contingency + 1)

    # Prepare hover text with raw counts
    hover_text = np.array([
        [f"X: {x}<br>Y: {y}<br>Count: {contingency.loc[y, x]}"
         for x in contingency.columns]
        for y in contingency.index
    ])

    fig = go.Figure(data=go.Heatmap(
        z=contingency_log.values,
        x=contingency.columns.tolist(),
        y=contingency.index.tolist(),
        text=hover_text,
        hoverinfo='text',
        colorscale=colorscale,
        colorbar=dict(title='log2(count+1)', len=0.75),
    ))

    fig.update_layout(
        title=f"Contingency Table: {y_column} vs {x_column}",
        xaxis=dict(title=x_column, tickangle=45),
        yaxis=dict(title=y_column, autorange='reversed'),
        height=height,
        width=width,
        plot_bgcolor='white'
    )

    return fig



        