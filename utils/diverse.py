"""
All kind of functions that do not fit in any of the other scripts
"""
import os
os.environ['R_HOME'] = '/home/schroel1/miniconda3/envs/gpu2/lib/R'
os.environ['R_USER'] = '/home/schroel1/miniconda3/envs/gpu2/lib/R'

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import anndata2ri
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import r, globalenv
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from pybiomart import Server

anndata2ri.activate()


def run_fastmnn(adata, batch_column = 'library'): 
    if 'normalized' in adata.layers: 
        adata = ad.AnnData(X = adata.layers['normalized'], 
            obs = adata.obs[[batch_column]], 
            var = pd.DataFrame(index = adata.var.index))
    else: 
        raise InputError("No adata.layers[normalized] detected - adata must be normalized before running fastmnn")

    sce_in = anndata2ri.py2rpy(adata)

    ro.globalenv['sce_in'] = sce_in

    r_script = """
    library(SingleCellExperiment)
    library(batchelor, quietly=TRUE)
    logcounts(sce_in) <- assay(sce_in, "X")
    print("X set as logcounts, now starting to compute the fastmnn corrected counts")
    sce_tmp <- fastMNN(sce_in, 
                            batch = sce_in$library, 
                            correct.all = TRUE, 
                            d = 50, 
                            k = 30)
    mat <- assay(sce_tmp, "reconstructed")
    mat <- as(mat, "dgCMatrix")
    sce_out <- SingleCellExperiment(list(counts = mat))
    rownames(sce_out) <- rownames(sce_tmp)
    colnames(sce_out) <- colnames(sce_tmp)
    reducedDim(sce_out) <- reducedDim(sce_tmp)
    reducedDimNames(sce_out) <- 'corrected_pca'

    merge_info = metadata(sce_tmp)$merge.info

    print(merge_info)

    sce_out
    """
    sce_out = ro.r(r_script)
    # Convert only sce_out to AnnData
    adata_out = anndata2ri.rpy2py(sce_out)
    return adata_out


def get_transcription_factors(go_list_path = '/data/hadjantalab/lucas/sonja_project/processing/data/transcription_factors.csv'):
    ## extracting all the transcription factors from my results
    ## list from https://amigo.geneontology.org/amigo/term/GO:0003700     
    # Load transcription factors data
    file_path = go_list_path
    transcription_factors_table = pd.read_table(file_path, sep=',')
    
    # Function to extract the correct MGI ID, handling cases like 'MGI:MGI:95388'
    def extract_mgi_id(mgi_string):
        # Find the last occurrence of 'MGI:' and extract it with the following ID
        return 'MGI:' + mgi_string.split('MGI:')[-1]
    
    # Apply the extraction of the MGI ID
    transcription_factors_table['MGI ID'] = transcription_factors_table['bioentity'].apply(extract_mgi_id)
    
    # Connect to the Ensembl BioMart server
    server = Server(host='http://www.ensembl.org')
    
    # Access the Mouse (Mus musculus) dataset
    mart = server.marts['ENSEMBL_MART_ENSEMBL']
    dataset = mart.datasets['mmusculus_gene_ensembl']
    
    # Query BioMart to get all MGI IDs and corresponding HGNC symbols
    # Get all entries with MGI ID and HGNC symbol (for human genes)
    biomart_results = dataset.query(attributes=['mgi_id', 'mgi_symbol'])
    
    # Filter BioMart results to include only MGI IDs present in the 'transcription_factors_table' DataFrame
    biomart_results = biomart_results[biomart_results['MGI ID'].isin(transcription_factors_table['MGI ID'])]
    
    # Merge the filtered BioMart results with the transcription factors DataFrame
    transcription_factors_table = transcription_factors_table.merge(biomart_results, on='MGI ID', how='left')
    
    ## extract all the unique mgi symbols from my table
    transcription_factors = set(transcription_factors_table['MGI symbol'])
    return transcription_factors