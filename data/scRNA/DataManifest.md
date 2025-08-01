# Single Cell RNA-seq

## Data Accession

Data taken from this [publication](https://www.nature.com/articles/s41398-023-02678-x) under GEO accession number:

[GSE244477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244477)

### On Linux
From the root of this repository, run the following:

```
cd data/scRNA/
wget -O GSE244477_RAW.tar https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE244477&format=file
tar -xvf GSE244477_RAW.tar
```

### On Mac
From the root of this repository, run the following:

```
cd data/scRNA/
curl -o GSE244477_RAW.tar https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE244477&format=file
tar -xvf GSE244477_RAW.tar
```


## Manifest file

Formatted file containing condition information for each data set.

GSM7817739	WT1 scRNAseq
GSM7817740	WT2 scRNAseq
GSM7817741	WT3 scRNAseq
GSM7817742	Het 1 scRNAseq
GSM7817743	Het 2 scRNAseq
GSM7817744	Het 3 scRNAseq
GSM7817745	Het 4 scRNAseq
GSM7817746	KO 1 scRNAseq
GSM7817747	KO 2 scRNAseq
GSM7817748	KO 3 scRNAseq

