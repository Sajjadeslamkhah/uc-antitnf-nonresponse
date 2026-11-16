# Baseline mucosal signatures of anti-TNF non-response in ulcerative colitis

Analysis code for:

> Erçağ E, Eslamkhah S, Gülmez H, Şahin H. *Baseline mucosal signatures of anti-TNF non-response in ulcerative colitis
> largely reflect inflammatory burden: a placebo-controlled multi-cohort transcriptomic analysis.*
> *Genes* (MDPI), manuscript genes-4682835, under revision. Submitted under the title *Predictive but not independent: a
> placebo-controlled, non-redundant multi-cohort analysis of baseline mucosal signatures of anti-TNF non-response in
> ulcerative colitis*.

The study re-analyses four non-overlapping public microarray series of colonic biopsies from the Gene Expression Omnibus:

| GEO series | Platform | Cohort | Role |
|---|---|---|---|
| GSE92415 | GPL13158 | PURSUIT-SC golimumab induction (phase 2/3, placebo-controlled) | Discovery: UC vs control (Layer 1); non-response vs response under golimumab, placebo arm as negative control (Layer 2) |
| GSE75214 | GPL6244 | IBD mucosal cohort with non-IBD controls | Validation: UC vs control (Layer 1) |
| GSE16879 | GPL570 | Biopsies before a first infliximab infusion | Validation: non-response vs response under infliximab (Layer 2) |
| GSE23597 | GPL570 | Active Ulcerative Colitis Trial 1 (ACT 1), infliximab | Validation: non-response vs response under infliximab (Layer 2) |

Three further series (GSE87466, GSE12251, GSE73661) were screened and excluded because their samples duplicate, or cannot be shown to be independent of, the retained series (see Part 0 of the notebook and Supplementary Methods 2.2).

## Repository layout

```
notebooks/UC_AntiTNF_pipeline.ipynb   complete analysis, run top to bottom in Google Colab (outputs retained)
revision/reverse_discovery.py         alternative-discovery analysis added at revision (Reviewer 2, comment 7)
revision/plot_reverse_discovery.py    Supplementary Figure S13 for that analysis (Table S13 = revision/results/T_reverse_discovery_*.csv)
revision/results/                     output of the two revision scripts
data/README.md                        which input files are needed and where they come from
requirements.txt                      Python packages
```

### Notebook sections

| Part | Content | Supplementary Methods |
|---|---|---|
| 0 | Download of the seven candidate series; three-level overlap audit (GSM accessions, subject identifiers, cross-series expression similarity) | 2.2 |
| 1 | Probe annotation, gene-level collapsing, sample sheets, QC, limma differential expression for the six contrasts, Layer 1 and Layer 2 replication, signature scoring | 2.3–2.7 |
| 1b | Within-group outlier detection and sensitivity analysis, adjustment for Mayo score and transcriptomic inflammatory burden, signature decomposition, sex inference | 2.5, 2.7, 2.8 |
| 2 | Contribution of the response-specific component beyond inflammatory burden (empirical null), MCP-counter deconvolution, ORA and pre-ranked GSEA | 2.8–2.10 |
| 3 | WGCNA, module–trait association, module preservation, pre-specified candidate module selection, hub genes | 2.11 |
| 4b | STRING network centrality, druggability annotation (UniProt, ChEMBL, DGIdb, RCSB PDB, AlphaFold DB), target tiers | 2.12–2.13 |
| R.0, F1 | Collection of results tables and Figure 1 | — |

Each part writes a machine-readable run log (date, parameters, sample counts, package versions) next to its results.
Thresholds and decision rules are set in the first cell of each part, before any result is computed.
Cell `0.FIX` restores the Part 0 verification tables after an accidental overwrite on Drive; it embeds the outputs of the original run and recomputes nothing.

## Reproducing the analysis

The notebook was written for Google Colab with the project folder on Google Drive (`/content/drive/MyDrive/UC_AntiTNF`).
Part 0 downloads the series-matrix and platform files from GEO, so an internet connection is required on the first run.
Part 4b additionally needs the STRING network export for the 56-gene list (instructions in cell 4b.1) and queries public web APIs; responses are cached on first use.

Software used for the published results: R 4.6.1 with limma 3.68.5, WGCNA 1.74 and MCPcounter 1.2.0; Python 3 with gseapy 1.3.1, statsmodels 0.15.0, scikit-learn 1.6.1, SciPy 1.16.3 and NetworkX 3.7. Random seed 2026 throughout.

### Revision analysis: alternative discovery cohorts

`revision/reverse_discovery.py` repeats the Layer 2 pipeline with each infliximab cohort, and with both infliximab cohorts pooled (cohort as a blocking covariate), as the discovery set. Every resulting signature is scored in all remaining cohorts, including the placebo arm. Thresholds and scoring are unchanged from the primary analysis.

```bash
python revision/reverse_discovery.py      --data data --out revision/results
python revision/plot_reverse_discovery.py --out revision/results
```

Inputs are the gene-level matrices and sample sheets written by Part 1 of the notebook (see `data/README.md`). Requires R with limma.
With the golimumab arm as discovery the script reproduces the published differential-expression table (Table S3), the 43-gene signature and the AUCs of the primary analysis exactly; the results in `revision/results` were produced with R 4.3.3 and limma 3.58.1.

## Data availability

All expression data are public and are not redistributed here. They can be obtained from GEO under the accessions above; the notebook downloads them automatically.

## License

Code: MIT License (see `LICENSE`).
