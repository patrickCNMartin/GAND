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





