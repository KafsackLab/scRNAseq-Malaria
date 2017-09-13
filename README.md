# scRNAseq-Malaria
This repository contains all the scripts and external data files necessary to reproduce the analysis presented in our manuscript "Single-cell RNA-seq reveals a transcriptional signature of sexual commitment in malaria parasites", Nature, 2017 DOI:####
To perform the analysis, download all files to a local folder. To avoid changing paths in the scripts, keep all files in the same folder, without subfolders.
A working version of R is required, including the packages loaded in the beginning of each R script. The analysis was performed using Rstudio and R version 3.4.0. 

The files are divided to three categories:
a) All_DGE_Matrices.tar.gz - A tar file containing all single-cell gene expression matrices. Simply untar it into the main working folder.
b) 7 R scripts
- RawMatrices_to_dataFrame.R (creates a combined data matrix from all samples, and a metadata object including collection time, strain, etc.)

- Seurat_setup_code.R (a seurat setup and time assignment code, specifically for ap2gdd cells, but can easily be modified to other subsets as well. Make sure to use Seurat v1.4. For more information, see http://satijalab.org/seurat/)

- figure1_seurat.R (cell selection and clustering for figures 1b,c)

- DifferentialExP_ByStrain.R (calculation of differentially expressed genes for either NF54 or ap2gdd by cluster)

- monocle_analysis.R (one example for monocle ordering of our data, specifically for treated ap2g-dd cells. For more information see https://bioconductor.org/packages/release/bioc/html/monocle.html)

- Main_figures_code.R (code for most of the figures. Where external objects which are calculated elsewhere, referrals are provided)

- InternalFunctions.R (a set of function written to streamline the main figures code, including various functions for visualization, differential gene expression, etc.)

c) supplementary files needed for the codes to work:
- GeneModules.tab (The sets of genes appearing in Table S2)
- Bartfai_RNASeqTC.pct (bulk RNA seq sequential data for time point assignment),
- notPolyA_PlasmoDB32.txt (a list of non-polyadenylated genes to be excluded from analysis)
- GeneInfoTable_expanded_May2017_PlasmoDB32.csv (Information regarding genes including names, description, etc.)
