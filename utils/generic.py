import scanpy as sc 
import numpy as np
import anndata as ad
import sys
import os
sys.path.append('/data/hadjantalab/lucas/sonja_project/utils')
import visuals
import scipy

os.environ['R_HOME'] = '/home/schroel1/miniconda3/envs/gpu2/lib/R'
os.environ['R_USER'] = '/home/schroel1/miniconda3/envs/gpu2/lib/R'
import rpy2.robjects as ro
from rpy2.robjects import r
r('library(Seurat)')


def quality_control(adata):
    adata = adata.copy()
    if 'raw' in adata.layers: 
        adata.X = adata.layers['raw'].copy()
        print("adata.X set to adata.layers['raw']")
    var_names_lower = adata.var_names.str.lower()
    mask_ribo = var_names_lower.str.startswith('rps') | var_names_lower.str.startswith('rpl')
    adata.var['ribo'] = mask_ribo
    mask_mito = var_names_lower.str.startswith('mt-')
    adata.var['mito'] = mask_mito
    mask_hemo = var_names_lower.str.startswith('hb') & ~var_names_lower.isin(["hbp1", "hbegf"])
    adata.var['hemo'] = mask_hemo
    sc.pp.calculate_qc_metrics(adata, inplace = True, qc_vars=['mito', 'ribo', 'hemo'])
    if 'normalized' in adata.layers: 
        adata.X = adata.layers['normalized'].copy()
        print('Adata.X set back to adata.layers["normalized"]')
    return adata
        
def normalize(adata, target_sum = 1e4, base = 2):
    adata = adata.copy()
    print(f'Parameters set to: target_sum = {target_sum}, log_base = {base}')
    if 'raw' in adata.layers:
        print('adata.layers[raw] exists')
        all_integers = np.all(np.equal(np.mod(adata.layers['raw'].data, 1), 0))
        print(f'adata.layers[raw] contains only integers: {all_integers}')
        adata.X = adata.layers['raw'].copy()
        print("adata.X set to adata.layers['raw'] before normalization")
    else: 
        print("adata.layers[raw] does not exist - we therefore assume adata.X contains raw data")
        all_integers = np.all(np.equal(np.mod(adata.X.data, 1), 0))
        print(f'adata.X contains only integers: {all_integers}')
        adata.layers['raw'] = adata.X.copy()
        print("adata.layers[raw] is set to adata.X - and normalization will be performed on adata.X")
    sc.pp.normalize_total(adata, inplace = True, target_sum = target_sum)
    sc.pp.log1p(adata, base = base, copy = False)
    adata.layers['normalized'] = adata.X.copy()
    print("Normalization was successful - adata.X is now normalized and adata.layers[normalized] is a copy of adata.X")
    return adata
    
def reduce_dimensions(adata, n_top_genes = 4000, n_comps = 50, n_neighbors = 30, min_dist = 0.1, mask_var = 'highly_variable'): 
    """
    mask_var must be either highly_variable or None
    """
    print(f'Parameters set to: n_top_genes = {n_top_genes}, n_comps = {n_comps}, n_neighbors = {n_neighbors}, min_dist = {min_dist}')
    adata = adata.copy()
    if ('normalized' in adata.layers): 
        print("adata.layers[normalized] will be used for the pca")
        if mask_var == 'highly_variable':
            if 'raw' in adata.layers:
                print("adata.layers[raw] will be used for highly variable gene selection (seurat_ve3)")
                sc.pp.highly_variable_genes(adata, layer = 'raw', n_top_genes = n_top_genes, flavor = 'seurat_v3', inplace = True) 
                sc.pp.pca(adata, n_comps = n_comps, layer = 'normalized', mask_var = mask_var, copy = False)
 
            else: 
                raise ValueError("adata.layers['raw'] not found - hvg selection with seurat_v3 can not be performed")
        elif mask_var is None: 
            sc.pp.pca(adata, n_comps = n_comps, layer = 'normalized', mask_var = mask_var, copy = False)
        else: 
            raise ValueError("Incompatible value for mask_var received")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca', metric='euclidean', copy = False)
        sc.tl.umap(adata, min_dist=min_dist, copy = False)
        return adata
    else: 
        print("adata.layers[normalized] not found - normalize your data")

def clustering(adata, resolution, flavour = 'leiden',n_neighbors = 30, neighbors_key = 'neighbors'):
    adata = adata.copy()
    if flavour == 'leiden':
        sc.tl.leiden(adata, resolution=resolution, neighbors_key=neighbors_key, key_added='leiden', flavor='leidenalg', random_state = 1, copy = False)
    if flavour == 'phenograph_leiden': 
        sc.external.tl.phenograph(adata, clustering_algo='leiden', k=n_neighbors, jaccard=True, primary_metric='euclidean', resolution_parameter=resolution, n_jobs = -1, copy = False)
    
    return adata

def reset_to_var_in_raw(adata): 
    """
    You gonna lose all slots except .var, .obs and X
    I am not using adata.raw.obs_names because adata.raw.obs_names gets modified when adata.obs_names gets modified
    this is true for subsetting, sorting the index, renaming etc etc
    """
    adata = adata.copy()
    all_var = adata.raw.var.copy()
    complete_X = adata.raw.X.copy()
    obs = adata.obs.copy()
    adata = ad.AnnData(var = all_var, obs = obs, X = complete_X)
    return adata

def filtering_var(adata, threshold_proportion = 0.005):
    """
    Function to remove genes that are expressed in less than a percentage of cells
    while keeping a fixed set of genes of interest (manually selected genes from the visuals script)
    """

    adata = adata.copy()
    ## expressed in number of cells threshold
    threshold = adata.shape[0]*threshold_proportion
    print(f'For shape {adata.shape}, the threshold is set to {threshold} cells')
    print(f'Genes expressend in more than {threshold_proportion*100}%: {sum(adata.var["n_cells_by_counts"] >= threshold)}')
    print(f'Genes expressend in less than {threshold_proportion*100}%: {sum(adata.var["n_cells_by_counts"] < threshold)}')
    print(f"A list of selected genes (Sonja's genes of interest) will be retained regardless of expression")

    ## list of marker genes of interest
    genes_of_interest = [
        gene.lower()
        for gene in (
    visuals.marker_guttube+
    visuals.marker_mesoderm+
    visuals.marker_ve+
    visuals.marker_ve_liver+
    visuals.marker_ve_query+
    visuals.marker_de)]
    var_names_lower = adata.var_names.str.lower()
    mask_genes_of_interest = var_names_lower.isin(genes_of_interest)
    not_found = [gene for gene in genes_of_interest if gene not in var_names_lower]
    if not_found != []: 
        raise ValueError(f'Exiting loop because these genes of interest are not found: {not_found}')

    ## setting the mask to keep var
    adata.var['keep'] = adata.var["n_cells_by_counts"] >= threshold
    adata.var.loc[adata.var['ribo'], 'keep'] = False
    adata.var.loc[mask_genes_of_interest, 'keep'] = True


    adata = adata[:, adata.var.pop('keep')].copy()
     
    return adata
    
def get_size_of_adata(adata):
    total_size = sys.getsizeof(adata)  # Base size of the AnnData object
    
    # Add sizes of the attributes
    total_size += sys.getsizeof(adata.X)         # Expression matrix
    total_size += sys.getsizeof(adata.obs)       # Observations
    total_size += sys.getsizeof(adata.var)       # Variables (genes)
    total_size += sys.getsizeof(adata.obsm)      # Multi-dimensional observations
    total_size += sys.getsizeof(adata.varm)      # Multi-dimensional variables
    total_size += sys.getsizeof(adata.uns)       # Unstructured annotations
    total_size += sys.getsizeof(adata.layers)    # Layers of data
    
    # Include the size of adata.raw, if it exists
    if adata.raw is not None:
        total_size += sys.getsizeof(adata.raw)
        total_size += sys.getsizeof(adata.raw.X)
        total_size += sys.getsizeof(adata.raw.var)
    
    return total_size / (1024 ** 3)  # Convert bytes to GB



def get_gene_count_statistics(adata, population, layer_normalized = 'normalized', layer_raw = 'raw'):
    statistics_dict = {}

    raw_counts = adata.layers[layer_raw].copy()
    norm_counts = adata.layers[layer_normalized].copy()
    
    if scipy.sparse.issparse(norm_counts):
        # Convert sparse matrix to dense
        cells_per_gene = np.array((raw_counts >= 1).sum(axis=0)).flatten()
        signal_in_percent = (cells_per_gene / adata.shape[0]) * 100
        sum_norm = np.array(norm_counts.sum(axis=0)).flatten()
        mean_norm = np.array(norm_counts.mean(axis=0)).flatten()
        median_norm = np.median(norm_counts.toarray(), axis=0)

        statistics_dict['gene_symbol'] = adata.var_names.copy()
        statistics_dict['highly_variable'] = adata.var['highly_variable'].copy()
        statistics_dict[f'#_cells_{population}'] = np.repeat(adata.shape[0], adata.shape[1])
        statistics_dict[f'>=1_in_%_{population}'] = signal_in_percent
        statistics_dict[f'sum_{population}'] = sum_norm
        statistics_dict[f'mean_{population}'] = mean_norm
        statistics_dict[f'mean_if_expressed_{population}'] = np.divide(sum_norm, cells_per_gene, out=np.zeros_like(sum_norm), where=cells_per_gene != 0)
        statistics_dict[f'median_{population}'] = median_norm
    else:
        print(f'The adata for {population} does not have a sparse matrix in norm_counts')

    return statistics_dict
     

def counts_to_obs(adata, gene_symbols, layer = 'normalized'): 
    adata = adata.copy()
    gene_indices = adata.var.index.get_indexer(gene_symbols)
    gene_expression = adata.layers[layer][:, gene_indices]
    gene_expression = gene_expression.toarray()

    for i, gene in enumerate(gene_symbols): 
        adata.obs[gene] = gene_expression[:,i].copy()
    return adata


def get_cell_cycle_score(adata): 
    adata = adata.copy()
    r('s_genes_satija <- cc.genes$s.genes')
    r('g2m_genes_satija <- cc.genes$g2m.genes')
    s_genes_satija = ro.globalenv['s_genes_satija']
    g2m_genes_satija = ro.globalenv['g2m_genes_satija']
    s_genes_satija = [gene.capitalize() for gene in s_genes_satija]
    g2m_genes_satija = [gene.capitalize() for gene in g2m_genes_satija]
    if adata.raw is not None:
        print(f'{sum(adata.raw.var_names.isin(s_genes_satija))} out of {len(s_genes_satija)} s-phase genes (annotation from satija lab) found in adata.raw.var_names')
        print(f'{sum(adata.raw.var_names.isin(g2m_genes_satija))} out of {len(g2m_genes_satija)} g2m-phase genes (annotation from satija lab) found in adata.raw.var_names')
        s_genes = adata.raw.var_names[adata.raw.var_names.isin(s_genes_satija)]
        g2m_genes = adata.raw.var_names[adata.raw.var_names.isin(g2m_genes_satija)]
    else: 
        print(f'{sum(adata.var_names.isin(s_genes_satija))} out of {len(s_genes_satija)} s-phase genes (annotation from satija lab) found in adata.var_names')
        print(f'{sum(adata.var_names.isin(g2m_genes_satija))} out of {len(g2m_genes_satija)} g2m-phase genes (annotation from satija lab) found in adata.var_names')
        s_genes = adata.var_names[adata.var_names.isin(s_genes_satija)]
        g2m_genes = adata.var_names[adata.var_names.isin(g2m_genes_satija)]
    sc.tl.score_genes_cell_cycle(adata, s_genes = s_genes, g2m_genes = g2m_genes, copy = False)

    return adata


def run_magic(adata, knn = 5, t = 3, name_list = 'all_genes', solver = 'exact', n_pca = None, n_jobs = 7): 
    adata = adata.copy()
    adata.X = adata.layers['normalized'].copy()
    sc.external.pp.magic(adata, 
                         name_list = name_list, 
                         knn = knn, 
                         t = t, 
                         solver = solver, 
                         n_pca = n_pca, 
                         n_jobs = n_jobs, 
                         random_state = 3, 
                         copy = False)
    adata.layers['magic'] = adata.X.copy()
    adata.X = adata.layers['normalized'].copy()
    return adata