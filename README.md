# GAND
GAND project repository for reproducible analysis.

# Reproducible Workflow
## Nextflow
Building a nextflow pipeline that will handle everything from data download (when relevant or possible), to report building. 

## Data Download 
The data should be stored in the `data` directory under the appropriate category.
To reduce the memory footprint of this repository, we do not provide the data here.

### scRNA Data

Nextflow will handle data download automatically using a data download workflow.

We download two data sets:

1. ScRNA-seq data taken from this [publication](https://www.nature.com/articles/s41398-023-02678-x) under GEO accession number [GSE244477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244477). To facilitate, formatting and naming of data sets, we provide a data `manifest.txt` which contains the names and conditions of each individual data set. 

2. A reference data sets for cell type annotation ([GEO: GSE123335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123335))

The `nextflow.config` contains a `dwl` section which allows you to select which data sets should be downloaded. Simply toggle the booleans to `true` to selectively downloaded data sets. The download URLs are alreaddy added.

```
dwl {
        dwl_scrna = false
        scrna_url = ["https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE244477&format=file"]
        dwl_ref = false
        ref_url = ["https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123335/suppl/GSE123335%5FE14%5Fcombined%5Fmatrix.txt.gz",
                    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123335/suppl/GSE123335%5FE14%5Fcombined%5Fmatrix%5FClusterAnnotations.txt.gz"]
    }

```

### BigWig
TBD

## Containers & Environments

### Conda yml
Conda `.yml` files are available for each sub-workflow in the `envs` directory. The current workflow will automatically make use of the appropriate environment. 

NOTE: The config file contains `beforeScript` clause that will be remove in final iterations. 

### Containers
Singularity/Docker containers will be made available. 

## Data upload
[Nextflow pipeline for data upload?](https://github.com/nf-core/proposals/issues/79)





