#!/bin/bash
#SBATCH --job-name=GAND
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
module load apptainer
module load nextflow
nextflow run main.nf -profile container_hpc