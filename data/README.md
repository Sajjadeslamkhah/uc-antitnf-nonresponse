# Input data

Nothing in this folder is tracked by git except this file. All inputs derive from public GEO series.

## Notebook

Part 0 downloads the series-matrix files of GSE87466, GSE92415, GSE75214, GSE16879, GSE12251, GSE73661 and GSE23597,
and Part 1 the platform annotation tables of GPL13158, GPL6244 and GPL570, into `raw_geo/` of the project folder.
Part 1 writes the gene-level matrices, sample metadata and sample sheets into `01_preprocessing_DEG/expression/`.

## Revision scripts

`revision/reverse_discovery.py --data data` expects:

```
data/expression/GSE92415_gene_level_all_samples.csv.gz
data/expression/GSE16879_gene_level_all_samples.csv.gz
data/expression/GSE23597_gene_level_all_samples.csv.gz
data/samplesheets/samplesheet_L2_disc_GOL_GSE92415.csv
data/samplesheets/samplesheet_L2_placebo_GSE92415.csv
data/samplesheets/samplesheet_L2_val_IFX_GSE16879.csv
data/samplesheets/samplesheet_L2_val_IFX_GSE23597.csv
data/Supplementary_Tables_S1-S12.xlsx        (sheet S4_core_UC_genes_687)
```

Copy them from `01_preprocessing_DEG/expression/` of the notebook's project folder, and the workbook from the supplementary material.
