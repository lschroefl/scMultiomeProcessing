import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import anndata as ad
from sklearn.decomposition import PCA
import random
import os
import scanpy as sc
os.environ['R_HOME'] = '/home/schroel1/miniconda3/envs/gpu2/lib/R'
os.environ['R_USER'] = '/home/schroel1/miniconda3/envs/gpu2/lib/R'
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects.conversion as conversion
import rpy2.robjects.numpy2ri as numpy2ri
import sys 
sys.path.append('/data/hadjantalab/lucas/sonja_project/utils')
import visuals
import anndata2ri


# Activate automatic conversion for pandas <-> R
pandas2ri.activate()
numpy2ri.activate()

# Import R packages
limma = importr("limma")
base = importr("base")
stats = importr("stats")
#anndata2ri.activate() ## apparently leads to issues
#pandas2ri.activate() ## apparentl



def hist_pseudosamples(adata, sample_column = 'library', cluster_column = 'leiden', threshold = 30, logscale = True, height = 4, width = 6):
    adata = adata.copy()

    # Step 1: Group the data by 'leiden_scVI' and 'covariate_composite'
    grouped = adata.obs.groupby([cluster_column, sample_column], observed=True)

    # Step 2: Calculate the size of each group
    group_sizes = grouped.size()

    # Step 3: Plot the histogram of the group sizes with log2 x-axis scaling
    plt.figure(figsize=(width, height))
    sns.histplot(group_sizes, bins=1000, kde=False)  # Plot the raw group sizes

    # Set log2 scale for the x-axis
    if logscale == True: 
        plt.xscale('log', base=2)

    # Add a red faint vertical line at x = 30
    plt.axvline(x=threshold, color='red', linestyle='--', alpha=0.6)

    # Set plot titles and labels
    plt.title('Histogram of Group Sizes (Log2 X-axis Scale)', fontsize=16)
    plt.xlabel('Group Size (Number of Cells)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)

    # Show the plot
    plt.tight_layout()
    plt.show()



def contingency_table_pseudosamples(adata, sample_column = 'library', cluster_column = 'leiden', min_cells = 1, logscale = True, height = 4, width = 8):
    adata = adata.copy()

    # Step 1: Create the contingency table
    contingency_table = pd.crosstab(adata.obs[sample_column], adata.obs[cluster_column])
    contingency_table[contingency_table < min_cells] = 0
    annot = contingency_table.copy()
    if logscale == True:
        contingency_table = np.log2(contingency_table + 1)

    # Step 3: Set the figure size and plot the heatmap
    plt.figure(figsize=(width, height))  

    # Create heatmap and use 'annot' to display sampling frequencies
    heatmap = sns.heatmap(contingency_table,         # Use log-transformed values for color mapping
                        annot=annot,    # Annotate the heatmap with sampling frequencies
                        fmt="",           # No specific formatting
                        cmap="coolwarm",    # Color map for better contrast
                        linewidths=.5,      # Add gridlines for clarity
                        cbar_kws={'shrink': .8},  # Adjust the color bar size
                        annot_kws={"size": 6})  # Set font size for annotations

    # Step 4: Add title to the colorbar
    colorbar = heatmap.collections[0].colorbar
    if logscale == True:
        colorbar.ax.set_title('Log2 cells', fontsize=12) 
    else:
        colorbar.ax.set_title('cell count', fontsize=12) 
    # Step 5: Adjust axis labels
    plt.title('Cell number per pseudosample', fontsize=16)
    plt.xlabel(sample_column, fontsize=12)
    plt.ylabel(cluster_column, fontsize=14)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90, fontsize=8)  
    plt.yticks(rotation=0, fontsize=12)

    # Step 6: Adjust layout to ensure everything fits
    plt.tight_layout()

    # Show the plot
    plt.show()


### NOT USED -------------------------

def create_pseudobulk(adata, min_cells = 30, sample_column = 'library', cluster_column = 'leiden'):
    adata = adata.copy()
    
    grouped = adata.obs.groupby([cluster_column, sample_column], observed=True).filter(lambda x: len(x) > min_cells)
    grouped = grouped.groupby([cluster_column, sample_column], observed=True)
    
    aggregated_counts = []
    
    for group_name, group_df in grouped:
        obs_names = group_df.index
        adata_obs_indices = np.where(adata.obs.index.isin(obs_names))[0]
        
        group_data = adata.X[adata_obs_indices, :]
        group_data = group_data.toarray() if hasattr(group_data, "toarray") else group_data
        aggregated = np.sum(group_data, axis = 0)
        aggregated_counts.append((group_name, aggregated))
        
    X_pseudobulk = pd.DataFrame({group: aggregated for group, aggregated in aggregated_counts}).T
    X_pseudobulk.index = [group_name for group_name, _ in aggregated_counts]
    
    X_pseudobulk.columns = adata.var_names

    # Compute cell counts per group
    cell_counts = adata.obs.groupby([cluster_column, sample_column], observed=True).size()
    
    # Build obs_pseudobulk DataFrame
    obs_pseudobulk = pd.DataFrame({
        'cluster': [cluster for cluster, sc_sample in X_pseudobulk.index],
        'library': [sc_sample for cluster, sc_sample in X_pseudobulk.index],
        'cell_count': [cell_counts.loc[(cluster, sc_sample)] for cluster, sc_sample in X_pseudobulk.index],
        'lib_size': np.sum(X_pseudobulk, axis=1)
    }, index=X_pseudobulk.index)

    ## setting index to string otherwise I cannot create an adata object later
    obs_pseudobulk.index = obs_pseudobulk.index.map(lambda x: str(x))
    X_pseudobulk.index = X_pseudobulk.index.map(lambda x: str(x))


    return X_pseudobulk, obs_pseudobulk

def X_scaling(X_pseudobulk, target_sum = 1e6):
    X_pseudobulk = X_pseudobulk.apply(lambda x: x * ( target_sum / x.sum()), axis=1)
    return X_pseudobulk

def hist_gene_counts(X_pseudobulk, threshold = 25, logscale = True, width = 6, height = 4):
    plt.figure(figsize=(width, height))
    if logscale == False: 
        sns.histplot(X_pseudobulk.sum(axis=0), bins=1000, kde=False)  # Plot the raw group sizes
        plt.axvline(x=threshold, color='red', linestyle='--', alpha=0.6)
    else:
        sns.histplot(np.log2(X_pseudobulk.sum(axis=0)+1), bins=1000, kde=False)  # Plot the raw group sizes
        plt.axvline(x=np.log2(threshold+1), color='red', linestyle='--', alpha=0.6)
        
    # Set plot titles and labels
    plt.xlabel('log2(counts+1)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    
    # Show the plot
    plt.tight_layout()
    plt.show()

def X_filtering_genes(X_pseudobulk, threshold = 30):
    X_pseudobulk = X_pseudobulk.copy()
    print(f'X_pseudobulk shape before filtering: {X_pseudobulk.shape}')
    ## Filtering genes that have less than 30 counts across all samples
    X_pseudobulk = X_pseudobulk.loc[:, X_pseudobulk.sum(axis=0) >= threshold]
    print(f'X_pseudobulk shape after filtering: {X_pseudobulk.shape}')
    return X_pseudobulk

def X_log2(X_pseudobulk):
    X_pseudobulk = X_pseudobulk.copy()
    X_pseudobulk = np.log2(X_pseudobulk+1)
    return X_pseudobulk

def pseudobulk_to_adata_and_pca(X_pseudobulk, obs_pseudobulk, color_by = ['leiden', 'library', 'lib_size']):
    adata_pb = ad.AnnData(X = X_pseudobulk, obs = obs_pseudobulk)
    ## normalize and pca
    sc.pp.pca(adata_pb)
    adata_pb.obs["lib_size"] = np.sum(adata_pb.X, axis=1)
    # Extract variance ratios
    variance_ratios = adata_pb.uns['pca']['variance_ratio']
    pc0_variance = variance_ratios[0] * 100  # Percentage for PC1
    pc1_variance = variance_ratios[1] * 100  # Percentage for PC2
    pc2_variance = variance_ratios[2] * 100  # Percentage for PC1
    pc3_variance = variance_ratios[3] * 100  # Percentage for PC1

    sc.pl.pca(
        adata_pb,
        color=color_by,
        ncols=1,
        size=300,
        components=['0,1'],  # We are interested in PC1 and PC2
        title=f'PCA: PC0 ({pc0_variance:.2f}%) vs PC2 ({pc1_variance:.2f}%)'
    )
    sc.pl.pca(
    adata_pb,
    color=color_by,
    ncols=1,
    size=300,
    components=['2,3'],  # We are interested in PC1 and PC2
    title=f'PCA: PC2 ({pc2_variance:.2f}%) vs PC3 ({pc3_variance:.2f}%)'
    )


def pseudobulk_pca_plot(X_pseudobulk, obs_pseudobulk, color_by=['cluster', 'library', 'cell_count', 'lib_size'], width = 6, height = 4):
    # Compute library size if not already in obs
    
    # Perform PCA
    pca = PCA(n_components=10)  # you can increase this if you want more PCs
    X_pca = pca.fit_transform(X_pseudobulk)
    
    # Variance explained per PC
    variance_ratios = pca.explained_variance_ratio_
    
    # Add PC coordinates to obs DataFrame for plotting
    obs_pseudobulk = obs_pseudobulk.copy()
    obs_pseudobulk['PC1'] = X_pca[:, 0]
    obs_pseudobulk['PC2'] = X_pca[:, 1]
    obs_pseudobulk['PC3'] = X_pca[:, 2]
    obs_pseudobulk['PC4'] = X_pca[:, 3]
    
    # Plotting
    for pair, title in [ (('PC1', 'PC2'), f'PC1 ({variance_ratios[0]*100:.2f}%) vs PC2 ({variance_ratios[1]*100:.2f}%)'),
                         (('PC3', 'PC4'), f'PC3 ({variance_ratios[2]*100:.2f}%) vs PC4 ({variance_ratios[3]*100:.2f}%)') ]:
        for color in color_by:
            plt.figure(figsize=(width,height))
            sns.scatterplot(
                x=pair[0],
                y=pair[1],
                data=obs_pseudobulk,
                hue=color,
                s=100,
                alpha=0.8
            )
            plt.title(f'PCA: {title} colored by {color}')
            plt.xlabel(pair[0])
            plt.ylabel(pair[1])
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.show()


#### NOT USED OVER ------------------------------


def run_dge(X_pseudobulk, obs_pseudobulk, gene_symbols, block_variable = False, weights_variable = False, model_formula = '~0 + leiden_of_interest'):
    obs_pseudobulk = obs_pseudobulk.copy()
    X_pseudobulk = X_pseudobulk.copy()
    X_pseudobulk = X_pseudobulk.T
    gene_symbols = ro.StrVector(gene_symbols)


    # Convert pandas dataframe and numpy array to R objects
    obs_pseudobulk_r = pandas2ri.py2rpy(obs_pseudobulk)
    X_pseudobulk_r = numpy2ri.py2rpy(np.asarray(X_pseudobulk))
    
    # Assign to R globalenv so we can access from R code
    ro.globalenv['obs_pseudobulk'] = obs_pseudobulk_r
    ro.globalenv['X_pseudobulk'] = X_pseudobulk_r
    ro.globalenv['model_formula'] = model_formula
    ro.globalenv['block_variable'] = block_variable
    ro.globalenv['weights_variable'] = weights_variable
    ro.globalenv['gene_symbols'] = gene_symbols

    
    # Run R code as a multiline string
    r_script = """
    library(limma)
    leidens <- unique(obs_pseudobulk$leiden)
    dge_results <- list()
    rownames(X_pseudobulk) <- gene_symbols

    for (leiden in leidens) {
        leiden_of_interest <- ifelse(obs_pseudobulk$leiden == leiden, 1, 0)
        design <- model.matrix(as.formula(model_formula), obs_pseudobulk)
        colnames(design) <- make.names(colnames(design))
        if (block_variable == FALSE){
            if (weights_variable == FALSE){
                fit <- lmFit(X_pseudobulk, design = design)
            }
            else{
                weights <- obs_pseudobulk[[weights_variable]]
                fit <- lmFit(X_pseudobulk, design = design, weights = weights)
            }
        }
        else{
            block <- obs_pseudobulk[[block_variable]]
            cor <- duplicateCorrelation(X_pseudobulk, design, block=block)
            print(paste('Consensus Correlation', toString(cor$consensus.correlation)))
            if (weights_variable == FALSE) {
                fit <- lmFit(X_pseudobulk, design=design, block=block, correlation=cor$consensus.correlation)
            }
            else{
                weights <- obs_pseudobulk[[weights_variable]]
                fit <- lmFit(X_pseudobulk, design=design, block=block, correlation=cor$consensus.correlation, weights = weights)
            }
            }
        contrast_matrix <- makeContrasts(leiden_of_interest, levels=colnames(design))
        fit <- contrasts.fit(fit, contrast_matrix)
        fit <- eBayes(fit)
        
        results_leiden_vs_rest <- topTable(fit, adjust="BH", number=Inf)
        dge_results[[paste0("leiden_", leiden)]] <- results_leiden_vs_rest
    }

    dge_results
    """

    # 5) run & unpack
    out_r = ro.r(r_script)

    # pull back each component
    out_py   = {}
    for key, table in out_r.items():
        out_py[key] = pandas2ri.rpy2py(table)
        #out_py[key] = table 


    return out_py


def aggregate_and_filter(
    adata,
    cell_identity,
    donor_key="sample",
    cell_identity_key="cell_type",
    obs_to_keep=None,  # which additional metadata to keep, e.g. gender, age, etc.
    replicates_per_patient=1,
    NUM_OF_CELL_PER_DONOR = 3,
    layer = 'raw',
    ):
    """
    code from https://www.sc-best-practices.org/conditions/differential_gene_expression.html 
    modified by LS
    """
    # subset adata to the given cell identity
    if obs_to_keep is None:
        obs_to_keep = []
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)
    df = pd.DataFrame(
        columns=[*adata_cell_pop.var_names, *obs_to_keep, "cell_counts", "log2_cell_counts", "pseudosample"])

    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            # create replicates for each donor
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                agg_dict["cell_counts"] = "first"
                agg_dict["log2_cell_counts"] = "first"
                agg_dict["pseudosample"] = "first"
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"
                # create a df with all genes, donor and condition info
                df_donor = pd.DataFrame(adata_replicate.layers[layer].A)
                df_donor.index = adata_replicate.obs_names
                df_donor.columns = adata_replicate.var_names
                df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
                #print(f"{donor}, len: {len(df_donor.index)}")
                df_donor['cell_counts'] = len(df_donor.index)
                df_donor['log2_cell_counts'] = np.log2(df_donor['cell_counts'])
                df_donor['pseudosample'] = "donor_" + donor + "_population" + df_donor[cell_identity_key][0] + "_sample" + str(i)
                #print(df_donor)
                # aggregate
                df_donor = df_donor.groupby(donor_key).agg(agg_dict)
                df_donor[donor_key] = donor
                df.loc[f"donor_{donor}_population{df_donor[cell_identity_key][0]}_sample{i}"] = df_donor.loc[donor]


    print("\n")
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names].astype(np.float32), obs=df.drop(columns=adata_cell_pop.var_names)
    )
    return adata_cell_pop




def filtering_var(adata, threshold = 200):
    """
    Function to remove genes where total_counts < threshold
    while keeping a fixed set of genes of interest (manually selected genes from the visuals script)

    """
    adata = adata.copy()
    ## expressed in number of cells threshold
    
    print(f'Threshold for shape {adata.shape}: {threshold}')
    print(f'Genes with total counts more than {threshold}: {sum(adata.var["total_counts"] >= threshold)}')
    print(f'Genes with total counts less than {threshold}: {sum(adata.var["total_counts"] < threshold)}')

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


def ma_plot(data):
    data = data.copy()
    
    plots_per_row = 3
    num_plots = len(data)
    rows = (num_plots // plots_per_row) + (num_plots % plots_per_row > 0)
    
    # Create subplots grid
    fig, axes = plt.subplots(rows, plots_per_row, figsize=(15, 5 * rows))  # Adjust size as needed
    
    # Flatten axes array if there are multiple rows
    axes = axes.flatten()
    
    # Plot each dataframe
    for idx, (name, df) in enumerate(data.items()):
        # Create a condition where the color is assigned based on the adjusted p-value
        significant = np.where(df['adj.P.Val'] < 0.1, 'red', 'black')
    
        # Create scatter plot in the corresponding subplot
        axes[idx].scatter(df['AveExpr'], df['logFC'], color=significant, s=4)
    
        # Add labels and title for clarity
        axes[idx].set_xlabel('log2(cpm+1) expression')
        axes[idx].set_ylabel('log fold change')
        axes[idx].set_title(f'Scatter plot for {name}')
    
    # Hide any empty subplots
    for ax in axes[num_plots:]:
        ax.axis('off')
    
    # Adjust layout for better spacing
    plt.tight_layout()
    
    # Show the entire figure
    plt.show()


def fit_model_edgeR(adata): 
    """
    function taken from https://www.sc-best-practices.org/conditions/differential_gene_expression.html 

    """
    adata = adata.copy()
    adata.X = adata.layers['unscaled'].copy()
    del adata.var, adata.uns, adata.obsm, adata.varm, adata.layers, adata.obsm
    adata.obs = adata.obs[['library', 'leiden']].copy()
    print(adata)

    ro.globalenv['adata_'] = anndata2ri.py2rpy(adata)


    r_script = """
    library(edgeR)
    # create an edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "X"), group = colData(adata_)$leiden)
    # filter out genes with low counts
    print("Dimensions before subsetting:")
    print(dim(y))
    print("")
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE]
    print("Dimensions after subsetting:")
    print(dim(y))
    print("")
    # normalize
    y <- calcNormFactors(y)
    # create a vector that is a concatenation of condition and cell type that we will later use with contrasts
    group <- paste0(colData(adata_)$leiden)
    replicate <- colData(adata_)$library
    # create a design matrix: here we have multiple donors so also consider that in the design matrix
    design <- model.matrix(~ 0 + group + replicate)
    # estimate dispersion
    y <- estimateDisp(y, design = design)
    # fit the model
    fit <- glmQLFit(y, design)

    print(colnames(y$design))


    return(list("fit"=fit, "design"=design, "y"=y)) """

    # Evaluate the R code and get the results
    out = ro.r(r_script)

    fit = out.rx2('fit')
    design = out.rx2('design')
    y = out.rx2('y')

    return fit