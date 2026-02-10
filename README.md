# GAND
GAND project repository for reproducible analysis.

# Reproducible Workflow

The following analysis was wrapped in a `Nextflow` pipeline which handles data download to report building. A more extensive pipeline breakdown is available below. All path to files are all relative to the root of this directory.

IMPORTANT NOTICE: While we provide containers for this analysis, the target architecture was Linux (Amd64). We do not provide Mac images. Mac images can build using the `Dockerfile` and the `Nix` flake if required.  

## Main Workflow

### Data Download 

Data download of publically avalible data sets can be handled directly through the `workflows/dwl_data.nf` workflow. 

To modify data download parameters, open the `nextflow.config` file and:


1. Toggle data download ` run_download = true`
2. Modify data URLs and final locations (NOTE: Nextflow expects data to be placed in these directories - modify at your own risk)

```
dwl = [
        scrna_url    : [
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE244477&format=file"
        ],
        ref_url      : [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123335/suppl/GSE123335%5FE14%5Fcombined%5Fmatrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123335/suppl/GSE123335%5FE14%5Fcombined%5Fmatrix%5FClusterAnnotations.txt.gz"
        ],
        output_scrna : "${baseDir}/data/scRNA/",
        output_ref   : "${baseDir}/data/ref/"
    ]

```

NOTE: For clarity, we provide a `manifest.txt` in `data/scRNA` to facilitate the automatic naming of files. 

### scRNA Data Integration

Integration parameters can be found in the `nextflow.config` file under the following section:

```
// Integration workflow params
    integration = [
        input              : "${baseDir}/data/scRNA/",
        ref                : "${baseDir}/data/ref",
        tmp                : "${baseDir}/data/tmp_scrna/",
        manifest           : "${baseDir}/data/scRNA/manifest.txt",
        npcs               : 30, // Number of Princpal Components
        min_features       : 100, // Minimum number of feature per cell
        max_features       : 10000, // Maximum number of features per cell
        percent_mt         : 10, // Percentage of Mitochondrial RNA allowed per cell
        n_var_features     : 2000, // Number of Variable features used for PCA
        cluster_resolution : 0.4, // Louvain clustering resolution
        integration_tag    : "integrated", // Tag to name the integration objects and meta data
        integration_method : "RPCAIntegration", // Seurat method used to integrate data
    ]
```

### Cell Type Modelling

To model the expression of certain gene sets by cell type and by condition (accounting for samples), we used a *Linear Mixed Effect Model*. We also check if expression of certain gene were mutually exclusive in terms of expression patterns.

To modify, the gene sets update the parameters in the following section:

```
 // Modelling Params
    modelling = [
        tmp       : "${baseDir}/data/tmp_scrna/",
        annotated : "${params.integration.tmp}/GAND_seurat_annotated.rds",
        gene_sets : [["Chd3", "Foxp1","Foxp2","Satb2"],
                    ["Chd3", "Foxp1","Satb2"],
                    ["Chd3","Foxp2","Satb2"],
                    ["Arx"]], // Gene sets to check for enrichment by cell type
        min_cells  : 1, // Minimum number of cells used by cell type for modelling
        mut_genes : ["Foxp1","Foxp2"], // Check if two genes are mutually exclusively expressed in cells
        score_type: ["module","counts"], // mode = Seurat module score || counts = log counts => which to use for modelling
    ]

```

### Report

A final pdf report is built from all the intermediate `csv` files that are generated. These files are produced to mimic the Source Data formatting often required by journal for publication. While we provide a report building, we encourage you to use the source data to make your own plots should you wish to change the aesthetics.

The source data location is shown and handled by the `nextflow.config` file.

```
 report = [
        annotated : "${baseDir}/data/tmp_scrna/GAND_seurat_annotated.csv",
        mut_genes : "${baseDir}/data/tmp_scrna/mutually_exclusive_genes.csv",
        gene_sets  : "${baseDir}/data/tmp_scrna/*_geneset_list.csv",
        template  : "${baseDir}/bin/scRNA_report_template.Rmd",
        output    : "${baseDir}/results/scRNA",
    ]
```

## Containers

### Docker & Apptainer
A Docker image was built using the definition file found in `containers`. Currently, only a Linux (amd64) target was built since the analysis was run exclusively run on a Linux HPC. Specifially, we used the following commands from the root directory of this project. 

```
cd containers
docker buildx build --platform linux/amd64 \
  -t gand:v0.0.1 \
  --output type=tar,dest=gand_image.tar .
```

Converstion to HPC safe apptainer `.sif` file was achieved through the following command:

```
# If using a SLURM grid engine 
# module load apptainer

# Build the SIF from the docker-archive
apptainer build gand_image.sif docker-archive://gand_image.tar
```

### Nix OCI

For maximum reproducibility, we wanted provide a `flake.nix` to build `Nix` style bit-for-bit reproducible environments. Due to the challenge of building with `Nix`, we show the procedure in `nix_build.md` for Linux. Mac builds are complex and show a partial build procedure. We would recommend avoiding cross architecture builds with `Nix` unless you are ready to fight the system.





