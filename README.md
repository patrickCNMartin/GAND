# GAND
GAND project repository for reproducible analysis.

# Reproducible Workflow

The following analysis was wrapped in a `Nextflow` pipeline which handles data download to report building. A more extensive pipeline breakdown is available below. All path to files are all relative to the root of this directory.

IMPORTANT NOTICE: While we provide containers for this analysis, the target architecture was Linux (Amd64). We do not provide Mac images. Mac images can build using the `Dockerfile` and the `Nix` flake if required.  

## Main Workflow

### Data Download 
Data download of publically avalible data sets is handled directly through the `workflows/dwl_data.nf` workflow. 

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

### scRNA Data

Nextflow will handle data download automatically using a data download workflow.

We download two data sets:

1. ScRNA-seq data taken from this [publication](https://www.nature.com/articles/s41398-023-02678-x) under GEO accession number [GSE244477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244477). To facilitate, formatting and naming of data sets, we provide a data `manifest.txt` which contains the names and conditions of each individual data set. 

2. A reference data sets for cell type annotation ([GEO: GSE123335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123335))

The `nextflow.config` contains a `dwl` section which allows you to select which data sets should be downloaded. Simply toggle the booleans to `true` to selectively downloaded data sets. The download URLs are alreaddy added.

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
For highly reproducible environment builds, we also provide a `Nix` flake which will build a OCI/Docker image. Using Nix to build the docker image ensures that all package versions and underlying libraries will be as close to bit-for-bit reproducibl across machines and time. First, make sure that you have `Nix` installed. You can find it [here](https://nixos.org/download/). Note that `Nix` is not compatible with Windows but can work through `WSL`.

To build the image with `Nix` using linux as target:

```
nix run nixpkgs#darwin.linux-builder-x86_64
nix build .#default --builders 'ssh://builder x86_64-linux'
```

The resulting image with also use the `v0.0.1` nomenclature but the image will be called `gand_nix_image.tar` to clarify that this is strictly built image less prone to config drift. 






