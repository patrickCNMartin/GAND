#!/usr/bin/env python
import scanpy as sc
import numpy as np
import pandas as pd
import scanorama as scan
import os
import re
import argparse
import anndata as ad
import matplotlib.pyplot as plt
#-----------------------------------------------------------------------------#
# Import Data sets and prepare for integration
#-----------------------------------------------------------------------------#

def import_data(path: str, meta: pd.DataFrame):
    data_sets = []
    for index, row in meta.iterrows():
        item = row['ID']
        cond = row['cond']
        # Get matching files
        matching_files = [f for f in os.listdir(path) if item in f]
        if not matching_files:
            raise ValueError(f"No matching files for ID {item} in {path}")
        # Use regex to find the prefix from one of the files
        for fname in matching_files:
            match = re.match(r"^(.*)_(barcodes|matrix|features)\.(tsv|mtx)\.gz$", fname)
            if match:
                prefix = match.group(1) + '_'
                break
        else:
            raise ValueError(f"No valid 10X-style files found for ID {item}")
        # Read data
        counts = sc.read_10x_mtx(path, prefix=prefix)
        counts.obs['ID'] = item
        counts.obs['cond'] = cond
        cond_merged = re.sub(r"\d+", "", cond)
        counts.obs['cond_merged'] = cond_merged
        data_sets.append(counts)
    return data_sets



#-----------------------------------------------------------------------------#
# integrate with scanorama
#-----------------------------------------------------------------------------#
def scano_integrate(data_sets:list):
    adata_concat = ad.concat(data_sets)
    sc.pp.log1p(adata_concat)
    sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'ID')
    #var_genes_batch = adata_concat.var.highly_variable_nbatches > 0
    var_select = adata_concat.var.highly_variable_nbatches > 2
    var_genes = var_select.index[var_select]
    adata_concat.raw = adata_concat
    adata_concat = adata_concat[:,var_genes]
    sc.pp.scale(adata_concat)
    sc.pp.pca(adata_concat)
    batches = adata_concat.obs['ID'].cat.categories.tolist()
    alldata = {}
    for batch in batches:
        alldata[batch] = adata_concat[adata_concat.obs['ID'] == batch,]
        
    alldata2 = dict()
    for ds in alldata.keys():
        alldata2[ds] = alldata[ds][:,var_genes]
        
    adatas = list(alldata2.values())
    scan.integrate_scanpy(adatas, dimred = 50)
    scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]
    # make into one matrix.
    all_s = np.concatenate(scanorama_int)
    print(all_s.shape)
    # add to the AnnData object, create a new object first
    adata_concat.obsm["Scanorama"] = all_s
    sc.pp.neighbors(adata_concat, n_pcs =30, use_rep = "Scanorama")
    sc.tl.leiden(adata_concat, resolution = 1)
    sc.tl.umap(adata_concat)
    sc.tl.tsne(adata_concat, n_pcs = 30, use_rep = "Scanorama")
    return adata_concat


#-----------------------------------------------------------------------------#
# Run main
#-----------------------------------------------------------------------------#
def main():
    parser = argparse.ArgumentParser(description='GAND Integration')
    parser.add_argument('--path', help='Path to scRNA data')
    parser.add_argument('--manifest_name', help = 'Name of Manifest file')
    args = parser.parse_args()
    path = args.path
    meta = args.manifest_name
    meta = os.path.join(path, meta)
    meta = pd.read_table(meta, sep = ' ', header = None)
    meta = meta.rename(columns= {0:"ID", 1:"cond", 2:"type"})
    data_sets = import_data(path, meta)
    data_sets = scano_integrate(data_sets)
    fig, axs = plt.subplots(2, 2, figsize=(10,8), constrained_layout=True)
    sc.pl.umap(data_sets, color="ID", title="Scanorama UMAP - Sample ID", ax=axs[0,0], show=False)
    sc.pl.umap(data_sets, color="leiden", title="Scanorama UMAP - Clusters", ax=axs[0,1], show=False)
    sc.pl.umap(data_sets, color="cond", title="Scanorama UMAP - Conditions", ax=axs[1,0], show=False)
    sc.pl.umap(data_sets, color="cond_merged", title="Scanorama UMAP - Conditions Merged", ax=axs[1,1], show=False)

    # Save in current directory
    file_name = 'GAND_integrated.png'
    plt.savefig(file_name)
    return 0

if __name__ == "__main__":
    main()