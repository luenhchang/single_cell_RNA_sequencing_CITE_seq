#--------------------------------------------------------------------------------------------------------------
# File path: E:/Lab_MarkS/lunC_work/single-cell-RNA-seq/scripts/Bald-scRNAseq-HNSCC.R
# Modified from: E:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/scripts/replicate-scRNAseq_GSE139324.R
# Date created: 2021416
# Author: Chang
# Purposes:
## (1) Uncompress .gz files
## (2) read 10X files, including barcodes.tsv, genes.tsv, matrix.mtx, as Seurat objects
## (3) 

# Reference
## [RenameIdent can only rename one cluster #771](https://github.com/satijalab/seurat/issues/771)
## [Changing R plot options in Jupyter](https://notebook.community/andrie/jupyter-notebook-samples/Changing%20R%20plot%20options%20in%20Jupyter)
## [Using a new windows version of R in Jupyter notebooks](https://stackoverflow.com/questions/51647561/using-a-new-windows-version-of-r-in-jupyter-notebooks)
## [scRNAseq_Braun_etal_Immunity_Script.Rmd](https://github.com/BaldLab/2020_Braun_et_al_CD226_Immunity/blob/master/scRNAseq_Mouse/R/scRNAseq_Braun_etal_Immunity_Script.Rmd)
## [genename to ensembl id #1867](https://github.com/satijalab/seurat/issues/1867)
## [Introduction to scRNA-seq integration](https://satijalab.org/seurat/articles/integration_introduction.html)
## [Using dittoSeq to visualize (sc)RNAseq data](https://bioconductor.org/packages/devel/bioc/vignettes/dittoSeq/inst/doc/dittoSeq.html)
## [Setup a Seurat object, add the RNA and protein data](https://satijalab.org/seurat/articles/multimodal_vignette.html)
## [scRNA-seq/lessons](https://github.com/hbctraining/scRNA-seq/tree/master/lessons)
## [How to solve Error: cannot allocate vector of size 1.2 Gb in R?](https://www.researchgate.net/post/How_to_solve_Error_cannot_allocate_vector_of_size_12_Gb_in_R)
#----------------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------------------------
# 20210520    This script file saved as E:/Lab_MarkS/lunC_work/single-cell-RNA-seq/scripts/Bald-scRNAseq-HNSCC.R
# 20210416    Jason said to see DNAM expression in CD8-T cells, not globally
#----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# Type      location
# Reference https://github.com/hbctraining/scRNA-seq/blob/master/lessons/02_SC_generation_of_count_matrix.md
# Reference https://github.com/hbctraining/scRNA-seq/blob/master/lessons/03_SC_quality_control-setup.md
# Reference https://github.com/hbctraining/scRNA-seq/blob/master/lessons/04_5SC_quality_control.md
# Reference https://github.com/hbctraining/scRNA-seq/blob/master/lessons/05_normalization_and_PCA.md
# Reference https://github.com/hbctraining/scRNA-seq/blob/master/lessons/06_SC_SCT_and_integration.md
# Reference https://github.com/hbctraining/scRNA-seq/blob/master/lessons/07_SC_clustering_cells_SCT.md
# Reference https://github.com/hbctraining/scRNA-seq/blob/master/lessons/08_SC_clustering_quality_control.md

# Reference https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD226

# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX1_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX1_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/features.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX1_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/matrix.mtx.gz

# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX2_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX2_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/features.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX2_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/matrix.mtx.gz

# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX3_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX3_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/features.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX3_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/matrix.mtx.gz

# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX4_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX4_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/features.tsv.gz
# Input     /mnt/lustre/working/genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX4_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/matrix.mtx.gz

# Input     [human cell cycle markers](https://www.dropbox.com/s/hus4mrkueh1tfpr/cycle.rda?dl=1)

# Output    C:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CountsOnly_no_ADT/analysis-data-sets/.*.png
# Output    C:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CountsOnly_no_ADT/analysis-data-sets/merged_seurat.RData
# Output    C:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CountsOnly_no_ADT/analysis-data-sets/integrated_seurat.rds


# Input     C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/GSE139324_RAW/.*\\.tsv.gz
# Input     C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/GSE139324_RAW/.*\\.mtx.gz

# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/Seurat_Read10x_format/GSM.*/barcodes.tsv
# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/Seurat_Read10x_format/GSM.*/genes.tsv
# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/Seurat_Read10x_format/GSM.*/natrix.mtx

# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/analysis-data-sets/cell-counts-by-conditions.png
# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/analysis-data-sets/number-unique-molecular-identifiers-by-conditions.png
# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/analysis-data-sets/distribution-genes-detected-per-cell-by-conditions.png

# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/analysis-data-sets/merged_seurat.RData
# Output    C:/Lab_MarkS/lunC_work/NCBI-NLM-NIH-GOV-GEO-GSE/GSE139324/analysis-data-sets/filtered_seurat.RData
#--------------------------------------------------------------------------------------------------------------

#---------------------------------
# Set up directory in system drive
#---------------------------------
dir.C <- "C:"

dir.R.packages <- file.path(dir.C,"Program Files","R","R-4.1.0","library") #"C:/Program Files/R/R-4.0.2/library"

#---------------------------------
# Set up directory in storage drive
#---------------------------------
dir.E <- "E:"

dir.Bald.scRNA.202002 <- file.path(dir.E,"backup-genomeinfo","share","analysis","Bald_-_RNASeq","scRNA_Feb2020_CITE_VDJ","GEX_CiteSeq") 

dir.Bald.scRNA.202002.analysis <- file.path(dir.Bald.scRNA.202002,"analysis-results")
dir.Bald.scRNA.202002.out <- file.path(dir.Bald.scRNA.202002,"output")
dir.Bald.scRNA.202002.log <- file.path(dir.Bald.scRNA.202002,"error-log")

dir.data <- file.path(dir.E,"backup-genomeinfo","share","analysis","data")

dir.Bald.scRNA.202002.filtered.sample.1 <- file.path(dir.Bald.scRNA.202002,"GEX1_20200121_BaldNextSeq","outs","filtered_feature_bc_matrix")
dir.Bald.scRNA.202002.filtered.sample.2 <- file.path(dir.Bald.scRNA.202002,"GEX2_20200121_BaldNextSeq","outs","filtered_feature_bc_matrix")
dir.Bald.scRNA.202002.filtered.sample.3 <- file.path(dir.Bald.scRNA.202002,"GEX3_20200121_BaldNextSeq","outs","filtered_feature_bc_matrix")
dir.Bald.scRNA.202002.filtered.sample.4 <- file.path(dir.Bald.scRNA.202002,"GEX4_20200121_BaldNextSeq","outs","filtered_feature_bc_matrix")

dir.Bald.scRNA.202002.raw.sample.1 <- file.path(dir.Bald.scRNA.202002,"GEX1_20200121_BaldNextSeq","outs","raw_feature_bc_matrix")
dir.Bald.scRNA.202002.raw.sample.2 <- file.path(dir.Bald.scRNA.202002,"GEX2_20200121_BaldNextSeq","outs","raw_feature_bc_matrix")
dir.Bald.scRNA.202002.raw.sample.3 <- file.path(dir.Bald.scRNA.202002,"GEX3_20200121_BaldNextSeq","outs","raw_feature_bc_matrix")
dir.Bald.scRNA.202002.raw.sample.4 <- file.path(dir.Bald.scRNA.202002,"GEX4_20200121_BaldNextSeq","outs","raw_feature_bc_matrix")

dir.reference.Zhang2018 <- file.path(dir.E, "Lab_MarkS","lunC_work","cancer-immunology-references","Zhang et al 2018 Lineage tracking reveals dynamic relationships of T cells in colorectal cancer")

#dir.create(dir.GSE139324.10X)
#dir.create(dir.Bald.scRNA.202002.analysis)
#dir.create(dir.Bald.scRNA.202002.log)
#dir.create(dir.Bald.scRNA.202002.out, recursive = TRUE)

dir.exists(c(  dir.Bald.scRNA.202002
              ,dir.data
              ,dir.Bald.scRNA.202002.filtered.sample.1
              ,dir.Bald.scRNA.202002.filtered.sample.2
              ,dir.Bald.scRNA.202002.filtered.sample.3
              ,dir.Bald.scRNA.202002.filtered.sample.4
              ,dir.Bald.scRNA.202002.raw.sample.1
              ,dir.Bald.scRNA.202002.raw.sample.2
              ,dir.Bald.scRNA.202002.raw.sample.3
              ,dir.Bald.scRNA.202002.raw.sample.4
              ,dir.Bald.scRNA.202002.analysis
              ,dir.reference.Zhang2018)) # All TRUE

lapply(c( dir.Bald.scRNA.202002.raw.sample.1
         ,dir.Bald.scRNA.202002.raw.sample.2
         ,dir.Bald.scRNA.202002.raw.sample.3
         ,dir.Bald.scRNA.202002.raw.sample.4)
       ,function(x) list.files(path = x))

#--------------------------------------------------------------------------------------
# Update or install packages, dependent packages in 
## E:/Lab_MarkS/lunC_work/single-cell-RNA-seq/scripts/packages.R
#--------------------------------------------------------------------------------------

#--------------
# Load packages
#--------------
library(vctrs, lib.loc = dir.R.packages)
library(dplyr, lib.loc = dir.R.packages)
library(tidyverse, lib.loc = dir.R.packages)
library(R.utils, lib.loc = dir.R.packages)
library(htmltools, lib.loc = dir.R.packages)
library(Seurat, lib.loc = dir.R.packages)
library(foreach, lib.loc = dir.R.packages)
library(doParallel, lib.loc = dir.R.packages)
library(parallel, lib.loc = dir.R.packages)
library(ggplot2, lib.loc = dir.R.packages)
library(RCurl, lib.loc = dir.R.packages)
library(cowplot, lib.loc = dir.R.packages)

library(dittoSeq, lib.loc = dir.R.packages)
library(scRNAseq, lib.loc = dir.R.packages)
library(SingleCellExperiment, lib.loc = dir.R.packages)

library(mygene, lib.loc = dir.R.packages)
library(org.Hs.eg.db, lib.loc = dir.R.packages)
library(biomaRt, lib.loc = dir.R.packages)
library(ensembldb,lib.loc = dir.R.packages)
library(EnsDb.Hsapiens.v79, lib.loc = dir.R.packages)

library(EnhancedVolcano, lib.loc = dir.R.packages)

#---------------------------------------
# Decompress all gz files in a folder
#---------------------------------------

# Get file paths of barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
file.paths.10X.filtered.sample.1 <- list.files(path = dir.Bald.scRNA.202002.filtered.sample.1, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.filtered.sample.1) 3
file.paths.10X.filtered.sample.2 <- list.files(path = dir.Bald.scRNA.202002.filtered.sample.2, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.filtered.sample.2) 3
file.paths.10X.filtered.sample.3 <- list.files(path = dir.Bald.scRNA.202002.filtered.sample.3, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.filtered.sample.3) 3
file.paths.10X.filtered.sample.4 <- list.files(path = dir.Bald.scRNA.202002.filtered.sample.4, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.filtered.sample.4) 3

file.paths.10X.raw.sample.1 <- list.files(path = dir.Bald.scRNA.202002.raw.sample.1, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.raw.sample.1) 3
file.paths.10X.raw.sample.2 <- list.files(path = dir.Bald.scRNA.202002.raw.sample.2, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.raw.sample.2) 3
file.paths.10X.raw.sample.3 <- list.files(path = dir.Bald.scRNA.202002.raw.sample.3, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.raw.sample.3) 3
file.paths.10X.raw.sample.4 <- list.files(path = dir.Bald.scRNA.202002.raw.sample.4, pattern = "\\.gz", full.names = TRUE) # length(file.paths.10X.raw.sample.4) 3

# Create sample-level meta data
metaData <- data.frame(folder.path.filtered=c(  dir.Bald.scRNA.202002.filtered.sample.1
                                               ,dir.Bald.scRNA.202002.filtered.sample.2
                                               ,dir.Bald.scRNA.202002.filtered.sample.3
                                               ,dir.Bald.scRNA.202002.filtered.sample.4)
                       ,folder.path.raw=c( dir.Bald.scRNA.202002.raw.sample.1
                                           ,dir.Bald.scRNA.202002.raw.sample.2
                                           ,dir.Bald.scRNA.202002.raw.sample.3
                                           ,dir.Bald.scRNA.202002.raw.sample.4)
                       ,condition=c("unstimulated","unstimulated","stimulated","stimulated")
                       ,object.name.10X=c("GEX1","GEX2","GEX3","GEX4")
                       ,sample=c("unstimulated.1","unstimulated.2","stimulated.1","stimulated.2")
                       ,stringsAsFactors = F) %>%
  dplyr::mutate(folder.name.filtered=  basename(dirname(dirname(folder.path.filtered)))
                ,folder.name.raw=  basename(dirname(dirname(folder.path.raw)))
                ) # dim(metaData) 4 7

#--------------------------------------------------------------------------------
# Decompress .gz files
# Read 10X files (barcodes.tsv, genes.tsv, matrix.mtx) files as dgCMatrix objects
# Read dgCMatrix objects as Seurat objects
#--------------------------------------------------------------------------------
str(metaData)

# Find out column positions of gene symbols and Ensembl ID in the features.tsv 
file.exists("E:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX1_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/features.tsv")

genes <- read.table("E:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX1_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/features.tsv", stringsAsFactors = FALSE)
str(genes)
# 'data.frame':	33544 obs. of  4 variables:
# $ V1: chr  "ENSG00000243485" "ENSG00000237613" "ENSG00000186092" "ENSG00000238009" ...
# $ V2: chr  "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
# $ V3: chr  "Gene" "Gene" "Gene" "Gene" ...
# $ V4: chr  "Expression" "Expression" "Expression" "Expression" ...
unique(genes$V3) # [1] "Gene"     "Antibody"
unique(genes$V4) # [1] "Expression" "Capture"

file.exists("E:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX1_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/barcodes.tsv")

barcodes <- read.table("E:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/GEX1_20200121_BaldNextSeq/outs/filtered_feature_bc_matrix/barcodes.tsv", stringsAsFactors = FALSE)

str(barcodes)
# 'data.frame':	3055 obs. of  1 variable:
#   $ V1: chr  "AAACCTGAGTACGATA-1" "AAACCTGCACAGTCGC-1" "AAACCTGCAGGCTCAC-1" "AAACCTGTCATGCATG-1" ...

source(file = "E:/Lab_MarkS/lunC_work/single-cell-RNA-seq/scripts/RFunction_Run-Seurat-singleCellRNAsequencing-CITEsequencing.R")

# Run Seurant function on sample 1 and export the object as a rds file
RunSeuratScRNAseqCITEseq( object.10X.name=metaData$object.name.10X[1]
                          ,input.data.dir=metaData$folder.path.filtered[1]
                          ,input.data.sample.name=metaData$object.name.10X[1]
                          ,output.folder.path=dir.Bald.scRNA.202002.analysis
                          ,vector.genes.of.interest=c("CD226", "EOMES", "RORA","GZMB", "GZMK"))

# Run Seurant function on sample 2 and export the object as a rds file
RunSeuratScRNAseqCITEseq( object.10X.name=metaData$object.name.10X[2]
                          ,input.data.dir=metaData$folder.path.filtered[2]
                          ,input.data.sample.name=metaData$object.name.10X[2]
                          ,output.folder.path=dir.Bald.scRNA.202002.analysis
                          ,vector.genes.of.interest=c("CD226", "EOMES", "RORA","GZMB", "GZMK"))

# Run Seurant function on sample 3 and export the object as a rds file
RunSeuratScRNAseqCITEseq( object.10X.name=metaData$object.name.10X[3]
                          ,input.data.dir=metaData$folder.path.filtered[3]
                          ,input.data.sample.name=metaData$object.name.10X[3]
                          ,output.folder.path=dir.Bald.scRNA.202002.analysis
                          ,vector.genes.of.interest=c("CD226", "EOMES", "RORA","GZMB", "GZMK"))

# Run Seurant function on sample 4 and export the object as a rds file
RunSeuratScRNAseqCITEseq( object.10X.name=metaData$object.name.10X[4]
                          ,input.data.dir=metaData$folder.path.filtered[4]
                          ,input.data.sample.name=metaData$object.name.10X[4]
                          ,output.folder.path=dir.Bald.scRNA.202002.analysis
                          ,vector.genes.of.interest=c("CD226", "EOMES", "RORA","GZMB", "GZMK"))

#-----------------------------------------------------------------------------------------------
# Merge all Seurat objects as a single Seurat object if they are samples of the same tissue type
#-----------------------------------------------------------------------------------------------
# Read rds files from individual sample Seurat objects

GEX1 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filteredGEX1.rds")) # class(GEX1) "SeuratObject" "Bald-HNSCC-scRNAseq-CITEseq_filteredGEX1"

GEX2 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filteredGEX2.rds")) 
GEX3 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filteredGEX3.rds")) 
GEX4 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filteredGEX4.rds")) # class(GEX4)

# Edit meta data
GEX1@meta.data$condition <- "unstimulated"
GEX2@meta.data$condition <- "unstimulated"
GEX3@meta.data$condition <- "stimulated"
GEX4@meta.data$condition <- "stimulated"

# Merge all samples of the same tissue type
merged_seurat <- merge( x=GEX1
                        ,y=c(GEX2, GEX3, GEX4)) 

unique(merged_seurat@meta.data$orig.ident) # [1] "GEX1" "GEX2" "GEX3" "GEX4"
unique(merged_seurat@meta.data$condition) # [1] "unstimulated" "stimulated"

# @assays a list of 7
# merged_seurat@assays$RNA
# merged_seurat@assays$ADT

# Create .RData object to load at any time
save(merged_seurat, file=file.path(dir.Bald.scRNA.202002.analysis,"merged_seurat.rds"))

#---------------------------------------------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
## reference [Using Seurat with multimodal data](https://satijalab.org/seurat/articles/multimodal_vignette.html)
#---------------------------------------------------------------------------------------------------------------

# See the keys for all keyed components (assays, dimensional reduction, spatial images) of a Seurat object using the Key functio
Seurat::Key(merged_seurat)
#    RNA    ADT 
# "rna_" "adt_" 

Seurat::DefaultAssay(merged_seurat) # "RNA"
merged_seurat <- Seurat::NormalizeData(object = merged_seurat
                                       ,assay="RNA"
                                       ,verbose=FALSE) %>%
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData() %>%
  Seurat::RunPCA()

# Clustering and dimension (PCs) selection
# Check dimensions
## A dimension is a PC. Choose PC23 as the break point. Take the more conservative PC.
Seurat::ElbowPlot(object = merged_seurat
                  ,ndims = 30
                  ,reduction = "pca")

merged_seurat <- Seurat::FindNeighbors(object = merged_seurat
                                       ,dims=1:23
                                       ,verbose=FALSE)

merged_seurat <- Seurat::FindClusters(object = merged_seurat, resolution=c(0.5, 0.9), verbose=FALSE)

merged_seurat <- Seurat::RunTSNE(merged_seurat, assay="RNA",dims=1:23, verbose=FALSE) 

Seurat::DimPlot(merged_seurat
                ,reduction = "tsne"
                ,label = TRUE)

merged_seurat <- Seurat::RunUMAP(merged_seurat, assay="RNA", dims=1:23)

#-------------------------
# Normalise assay ADT
#-------------------------
merged_seurat <- Seurat::NormalizeData(object = merged_seurat
                                       , normalization.method="CLR"
                                       , assay="ADT"
                                       , margin=2 # normalise data across features (margin=1) or cells (margin=2)
                                       ) %>%
  Seurat::ScaleData(verbose=FALSE) %>%
  Seurat::RunPCA(verbose=FALSE)

#-----------------------------------------------------------------------
# Categorise expression levels of proteins (ADT) of cells into meta data
#-----------------------------------------------------------------------

# Get cell barcodes (ID) from the RNA assay
.cell.ids <- rownames(merged_seurat@meta.data) # class(.cell.ids) "character" # length(.cell.ids) 10178

# Get names of ADTs
ADT.names <- rownames(merged_seurat@assays$ADT)

for(i in 1:length(ADT.names)){
  # Get name of a ADT, replacing dashes with underscores
  .ADT.name <- stringr::str_replace_all(string = ADT.names[i], pattern = "-", replacement = "_")
  cat("Processing protein", .ADT.name, "\n")
  #----------------------------------------------------------------------
  # Subset cells based on the expression levels of ADT protein.
  ## 1st Qu: bottom 25% of cells with the lowest expression
  ## 3rd Qu: top 25% of cells with the highest expression
  ## between 1st Qu and 3rd Qu: middle 50% of cells
  #----------------------------------------------------------------------
  
  # Subset expression of a protein from the ADT assay data
  .ADT.expression <- merged_seurat@assays$ADT@data[i, ] # class(.ADT.expression) [1] "numeric" # length(.ADT.expression) 1611
  
  # Add ADT expression data to metadata
  merged_seurat <- Seurat::AddMetaData(object=merged_seurat
                                       ,metadata=.ADT.expression
                                       ,col.name=.ADT.name)
  # Make a matrix
  .ADT.expression.matrix <- t(as.matrix(.ADT.expression)) # dim(.ADT.expression.matrix) 1 1611
  # Get bottom and top 25% quantile cut off values from ADT data from a given protein
  .q.vals <- quantile(.ADT.expression.matrix, probs = c(0.25, 0.75))
  print(.q.vals)
  cat("Calculated cut off values for bottom 25% and top 25% as", .q.vals[1], "and", .q.vals[2], "\n")
  
  # Check number of cells in each subset
  ## (i=1) 1*1611 = 410+794+407
  .num.cells.1st.Qu <- sum(.ADT.expression.matrix <= .q.vals[1])
  .num.cells.middle50 <-sum(.ADT.expression.matrix > .q.vals[1] & .ADT.expression.matrix < .q.vals[2])
  .num.cells.3rd.Qu <- sum(.ADT.expression.matrix >= .q.vals[2]) # 407
  
  cat(  "Number of cells that have expression of ADT",.ADT.name, "<=", .q.vals[1], "is",.num.cells.1st.Qu, "\n"
        ,"Number of cells that have expression of ADT",.ADT.name, "between", .q.vals[1],"and", .q.vals[2], "is", .num.cells.middle50, "\n"
        ,"Number of cells that have expression of ADT",.ADT.name, ">=", .q.vals[2], "is",.num.cells.3rd.Qu, "\n")
  
  # Get logical variables of which cells fit in each subset
  .logical.cell.bottom25 <- .ADT.expression.matrix <= .q.vals[1]
  .logical.cell.middle50 <- .ADT.expression.matrix > .q.vals[1] & .ADT.expression.matrix < .q.vals[2]
  .logical.cell.top25 <- .ADT.expression.matrix >= .q.vals[2]
  
  cat("Evaluated protein expression against thresholds as logicals")
  
  # Get cell IDs for each subset
  .cell.ID.bottom25 <- colnames(.ADT.expression.matrix)[.logical.cell.bottom25] # class(.cell.ID.bottom25) "character" # length(.cell.ID.bottom25) 410
  
  .cell.ID.middle50 <- colnames(.ADT.expression.matrix)[.logical.cell.middle50]
  .cell.ID.top25 <- colnames(.ADT.expression.matrix)[.logical.cell.top25]
  
  .cells.in.bottom25 <- .cell.ids %in% .cell.ID.bottom25 # class(.cells.in.bottom25) [1] "logical" # length(.cells.in.bottom25) 1611
  .cells.in.middle50 <- .cell.ids %in% .cell.ID.middle50
  .cells.in.top25 <- .cell.ids %in% .cell.ID.top25
  
  # Classify cell protein expression into bottom 25%, middle 50%, and top 25%
  .cells.ADT.clustered <- .cell.ids
  .cells.ADT.clustered[.cells.in.bottom25] <- "b25"
  .cells.ADT.clustered[.cells.in.middle50] <- "m50"
  .cells.ADT.clustered[.cells.in.top25] <- "t25"
  
  # Add ADT clusters info to metadata
  .column.name.ADT.cluster.quantiles <- paste0("expr.cluster.", .ADT.name)
  
  cat("Name the column as",.column.name.ADT.cluster.quantiles)
  
  # Add categorised ADT expression data to metadata
  merged_seurat <- Seurat::AddMetaData(object=merged_seurat
                                       ,metadata=.cells.ADT.clustered
                                       ,col.name=.column.name.ADT.cluster.quantiles)
  cat("Categorised expression levels of protein"
      ,.ADT.name
      ,"as column"
      ,.column.name.ADT.cluster.quantiles
      , "in meta data slot", "\n")   
  
  levels(merged_seurat@meta.data[,.column.name.ADT.cluster.quantiles])
  # [1] "b25"  "m50"  "t25"
  
  # Calculate ADT cluster numbers and %
  table(merged_seurat@meta.data[,.column.name.ADT.cluster.quantiles])
  # b25  m50 t25 (i=1)
  # 410  794  407
  
  prop.table(table(merged_seurat@meta.data[,.column.name.ADT.cluster.quantiles]))*100
  #      b25      m50     t25 (i=1)
  # 25.45003 49.28616 25.26381    
}

#-----------------------------------------------------------------------------
# Categorise expression levels of selected genes (RNA) of cells into meta data
#-----------------------------------------------------------------------------
# |Column name| Value|Definition|
# | --- | --- | --- |
# | rna.expr.cluster.CD226 | CD226- | zero expression of CD226 |
# | rna.expr.cluster.CD226 | CD226+ | expression of CD226 greater than zero |
# | rna.expr.cluster.PDCD1 | PDCD1- | zero expression of PDCD1|
# | rna.expr.cluster.PDCD1 | PDCD1+ | expression of PDCD1 greater than zero |

genes.to.categorise <- c("CD226","PDCD1")
genes.to.categorise %in% rownames(merged_seurat)
RNA.expression <- merged_seurat@assays$RNA@data[genes.to.categorise, ]
cat("Data type of expression data: ",class(RNA.expression),"\n"
    ,"Dimension of expression data: ", dim(RNA.expression),"\n"
    ,"Rownames of expression data: ", head(rownames(RNA.expression)),"\n"
)

for(i in 1:length(genes.to.categorise)){
  gene <- genes.to.categorise[i]
  # Subset expression of genes from the RNA assay data
  # Make a one-row matrix
  RNA.expression.matrix <- t(as.matrix(RNA.expression[i,]))
  cat( "Processing gene ",gene,"\n"
       ,"Data type of the expression: ",class(RNA.expression.matrix),"\n"
       ,"Number of items: ", length(RNA.expression.matrix),"\n"
       ,"Rownames of expression matrix: ", head(rownames(RNA.expression)),"\n")
  
  # Group cells into zero, and non-zero expression of a gene
  RNA.expression.logical.negative <- RNA.expression.matrix == 0
  RNA.expression.logical.positive <- RNA.expression.matrix > 0    
  cat("Number of FALSE and TRUE: ",table(RNA.expression.logical.negative),"\n")
  cat("Number of FALSE and TRUE: ",table(RNA.expression.logical.positive),"\n")
  
  # Get cell IDs for each subset
  cell.IDs <- colnames(RNA.expression.matrix)
  cell.IDs.negative <- cell.IDs[RNA.expression.logical.negative]
  cell.IDs.positive <- cell.IDs[RNA.expression.logical.positive]
  cat(head(cell.IDs),"\n")
  cat(head(cell.IDs.negative),"\n")
  cat(head(cell.IDs.positive),"\n")
  
  logical.cell.neg <- cell.IDs %in% cell.IDs.negative
  cat(head(logical.cell.neg),"\n")
  logical.cell.pos <- cell.IDs %in% cell.IDs.positive
  cat(head(logical.cell.pos),"\n")
  
  # Categorise expression as neg (=0) or positive (>0)
  cells.clustered <- cell.IDs
  cells.clustered[logical.cell.neg] <- paste0(gene,"-") #"CD226-"
  cells.clustered[logical.cell.pos] <- paste0(gene,"+") #"CD226+"
  cat(head(cells.clustered),"\n")
  cat(tail(cells.clustered),"\n")
  
  # Add categories to meta data
  colname.cell.clustered <- paste0("rna.expr.cluster.",gene)
  merged_seurat <- Seurat::AddMetaData(object = merged_seurat
                                       ,metadata=cells.clustered
                                       ,col.name= colname.cell.clustered
  )
} # End the for loop

#-----------------------------------------------------------------------------
# Categorise expression levels of selected genes (RNA) of cells into meta data
#-----------------------------------------------------------------------------
# |Column name| Value|Definition|
# | --- | --- | --- |
# | rna.expr.cluster.nonzero.median.CD226 | CD226low | expression < median of non-zero CD226 |
# | rna.expr.cluster.nonzero.median.CD226 | CD226high | expression >= median of non-zero CD226 |
# | rna.expr.cluster.nonzero.median.PDCD1 | PDCD1low | expression < median of non-zero PDCD1|
# | rna.expr.cluster.nonzero.median.PDCD1 | PDCD1high | expression >= median of non-zero PDCD1 |

genes.to.categorise <- c("CD226","PDCD1")
genes.to.categorise %in% rownames(merged_seurat)
RNA.expression <- merged_seurat@assays$RNA@data[genes.to.categorise, ]
cat("Data type of expression data: ",class(RNA.expression),"\n"
    ,"Dimension of expression data: ", dim(RNA.expression),"\n"
    ,"Rownames of expression data: ", head(rownames(RNA.expression)),"\n"
)

for(i in 1:length(genes.to.categorise)){
  gene <- genes.to.categorise[i]
  # Subset expression of genes from the RNA assay data
  # Make a one-row matrix
  RNA.expression.matrix <- t(as.matrix(RNA.expression[i,]))
  # Get median of non-zero expression
  median.non.zero.expression <- median(RNA.expression.matrix[RNA.expression.matrix > 0])    
  cat( "Processing gene ",gene,"\n"
       ,"Data type of the expression: ",class(RNA.expression.matrix),"\n"
       ,"Number of items: ", length(RNA.expression.matrix),"\n"
       ,"Median of non-zero expression: ", median.non.zero.expression,"\n")
  
  # Group cells into smaller than (ST) median, and greater than or equal (GE) median
  RNA.expression.logical.ST <- RNA.expression.matrix < median.non.zero.expression
  RNA.expression.logical.GE <- RNA.expression.matrix >= median.non.zero.expression
  
  cat("Number of FALSE and TRUE: ",table(RNA.expression.logical.ST),"\n")
  cat("Number of FALSE and TRUE: ",table(RNA.expression.logical.GE),"\n")
  
  # Get cell IDs for each subset
  cell.IDs <- colnames(RNA.expression.matrix)
  cell.IDs.ST <- cell.IDs[RNA.expression.logical.ST]
  cell.IDs.GE <- cell.IDs[RNA.expression.logical.GE]
  cat(head(cell.IDs),"\n")
  cat(head(cell.IDs.ST),"\n")
  cat(head(cell.IDs.GE),"\n")
  
  logical.cell.ST <- cell.IDs %in% cell.IDs.ST
  cat(head(logical.cell.ST),"\n")
  logical.cell.GE <- cell.IDs %in% cell.IDs.GE
  cat(head(logical.cell.GE),"\n")
  
  # Categorise expression as neg (=0) or positive (>0)
  cells.clustered <- cell.IDs
  cells.clustered[logical.cell.ST] <- paste0(gene,"low") 
  cells.clustered[logical.cell.GE] <- paste0(gene,"high") 
  cat(head(cells.clustered),"\n")
  cat(tail(cells.clustered),"\n")
  
  # Add categories to meta data
  colname.cell.clustered <- paste0("rna.expr.cluster.nonzero.median.",gene)
  merged_seurat <- Seurat::AddMetaData(object = merged_seurat
                                       ,metadata=cells.clustered
                                       ,col.name= colname.cell.clustered
  )
} # End the for loop

#-------------------------------------------------------
# Concatenate CD226 and PDCD1 to make 4 groups of cells
#-------------------------------------------------------
# |Column name| Value|
# | --- | --- |
# |rna.expr.cluster.nonzero.median.CD226.PDCD1|CD226low_PDCD1low|
# |rna.expr.cluster.nonzero.median.CD226.PDCD1|CD226low_PDCD1high|
# |rna.expr.cluster.nonzero.median.CD226.PDCD1|CD226high_PDCD1low|
# |rna.expr.cluster.nonzero.median.CD226.PDCD1|CD226high_PDCD1high|

# Concatenate CD226 and PDCD1
rna.expr.cluster.CD226<- merged_seurat@meta.data$rna.expr.cluster.nonzero.median.CD226
rna.expr.cluster.PDCD1<- merged_seurat@meta.data$rna.expr.cluster.nonzero.median.PDCD1
#str(rna.expr.cluster.CD226) # chr [1:10178] "CD226-" "CD226-" "CD226-" "CD226-" "CD226-" "CD226-" ...
rna.expr.cluster.CD226.PDCD1 <- paste0(rna.expr.cluster.CD226,"_",rna.expr.cluster.PDCD1)
str(rna.expr.cluster.CD226.PDCD1) # chr [1:10178] "CD226-PDCD1-" "CD226-PDCD1-" "CD226-PDCD1-" "CD226-PDCD1-" ...
# Add concatenated CD226 and PDCD1 to meta data
merged_seurat <- Seurat::AddMetaData(object=merged_seurat
                                     ,metadata=rna.expr.cluster.CD226.PDCD1
                                     ,col.name="rna.expr.cluster.nonzero.median.CD226.PDCD1")
unique(merged_seurat@meta.data$rna.expr.cluster.nonzero.median.CD226.PDCD1)

genes.selected <- c("CD226","PDCD1","CD28","ENTPD1","CD96","TIGIT") # ,"EOMES"
genes.selected %in% rownames(merged_seurat) # Should be all TRUE
genes.selected.RNA <- paste0(Seurat::Key(merged_seurat[["RNA"]]), genes.selected) # Include only genes that can be found in the rownames() #

# protein feature names
protein.selected <- rownames(merged_seurat@assays$ADT)
protein.selected.ADT <- paste0(Seurat::Key(merged_seurat[["ADT"]]), protein.selected)

#-------------------------------------------------------
# Identifying clusters of cells by different resolutions
#-------------------------------------------------------
## The resolution argument determines the number of clusters. Larger resolution values find more clusters than small resolution

merged_seurat_0.1 <- Seurat::FindClusters(merged_seurat, resolution = 0.1, verbose=FALSE)
merged_seurat_0.2 <- Seurat::FindClusters(merged_seurat, resolution = 0.2, verbose=FALSE)
merged_seurat_0.3 <- Seurat::FindClusters(merged_seurat, resolution = 0.3, verbose=FALSE)
merged_seurat_0.4 <- Seurat::FindClusters(merged_seurat, resolution = 0.4, verbose=FALSE)
merged_seurat_0.5 <- Seurat::FindClusters(merged_seurat, resolution = 0.5, verbose=FALSE)
merged_seurat_0.6 <- Seurat::FindClusters(merged_seurat, resolution = 0.6, verbose=FALSE)
merged_seurat_0.7 <- Seurat::FindClusters(merged_seurat, resolution = 0.7, verbose=FALSE)
merged_seurat_0.8 <- Seurat::FindClusters(merged_seurat, resolution = 0.8, verbose=FALSE)
merged_seurat_0.9 <- Seurat::FindClusters(merged_seurat, resolution = 0.9, verbose=FALSE)
merged_seurat_1.0 <- Seurat::FindClusters(merged_seurat, resolution = 1.0, verbose=FALSE)
merged_seurat_1.1 <- Seurat::FindClusters(merged_seurat, resolution = 1.1, verbose=FALSE)
merged_seurat_1.2 <- Seurat::FindClusters(merged_seurat, resolution = 1.2, verbose=FALSE)

p1 <- Seurat::DimPlot(merged_seurat_0.1, reduction = "tsne")+ ggplot2::ggtitle("resol_0.1")
p2 <- Seurat::DimPlot(merged_seurat_0.2, reduction = "tsne")+ ggplot2::ggtitle("resol_0.2")
p3 <- Seurat::DimPlot(merged_seurat_0.3, reduction = "tsne")+ ggplot2::ggtitle("resol_0.3")
p4 <- Seurat::DimPlot(merged_seurat_0.4, reduction = "tsne")+ ggplot2::ggtitle("resol_0.4")
p5 <- Seurat::DimPlot(merged_seurat_0.5, reduction = "tsne")+ ggplot2::ggtitle("resol_0.5")
p6 <- Seurat::DimPlot(merged_seurat_0.6, reduction = "tsne")+ ggplot2::ggtitle("resol_0.6")
p7 <- Seurat::DimPlot(merged_seurat_0.7, reduction = "tsne")+ ggplot2::ggtitle("resol_0.7")
p8 <- Seurat::DimPlot(merged_seurat_0.8, reduction = "tsne")+ ggplot2::ggtitle("resol_0.8")
p9 <- Seurat::DimPlot(merged_seurat_0.9, reduction = "tsne")+ ggplot2::ggtitle("resol_0.9")
p10 <- Seurat::DimPlot(merged_seurat_1.0, reduction = "tsne")+ ggplot2::ggtitle("resol_1.0")
p11 <- Seurat::DimPlot(merged_seurat_1.1, reduction = "tsne")+ ggplot2::ggtitle("resol_1.1")
p12 <- Seurat::DimPlot(merged_seurat_1.2, reduction = "tsne")+ ggplot2::ggtitle("resol_1.2")

Seurat::CombinePlots(plots = list(p1, p2,p3,p4,p5,p6))
Seurat::CombinePlots(plots = list(p7,p8,p9,p10,p11,p12))

#--------------------------------------------------------
# Plot cells clustered by UMAP with different resolutions
#--------------------------------------------------------
p1 <- Seurat::DimPlot(merged_seurat_0.1, reduction = "umap")+ ggplot2::ggtitle("resol_0.1")
p2 <- Seurat::DimPlot(merged_seurat_0.2, reduction = "umap")+ ggplot2::ggtitle("resol_0.2")
p3 <- Seurat::DimPlot(merged_seurat_0.3, reduction = "umap")+ ggplot2::ggtitle("resol_0.3")
p4 <- Seurat::DimPlot(merged_seurat_0.4, reduction = "umap")+ ggplot2::ggtitle("resol_0.4")
p5 <- Seurat::DimPlot(merged_seurat_0.5, reduction = "umap")+ ggplot2::ggtitle("resol_0.5")
p6 <- Seurat::DimPlot(merged_seurat_0.6, reduction = "umap")+ ggplot2::ggtitle("resol_0.6")
p7 <- Seurat::DimPlot(merged_seurat_0.7, reduction = "umap")+ ggplot2::ggtitle("resol_0.7")
p8 <- Seurat::DimPlot(merged_seurat_0.8, reduction = "umap")+ ggplot2::ggtitle("resol_0.8")
p9 <- Seurat::DimPlot(merged_seurat_0.9, reduction = "umap")+ ggplot2::ggtitle("resol_0.9")
p10 <- Seurat::DimPlot(merged_seurat_1.0, reduction = "umap")+ ggplot2::ggtitle("resol_1.0")
p11 <- Seurat::DimPlot(merged_seurat_1.1, reduction = "umap")+ ggplot2::ggtitle("resol_1.1")
p12 <- Seurat::DimPlot(merged_seurat_1.2, reduction = "umap")+ ggplot2::ggtitle("resol_1.2")

Seurat::CombinePlots(plots = list(p1, p2,p3,p4,p5,p6))
Seurat::CombinePlots(plots = list(p7,p8,p9,p10,p11,p12))

#----------------------------------------------------------------------------------------------
# Find markers for every cluster compared to all remaining cells, report only the positive ones
#----------------------------------------------------------------------------------------------
cat("Identify differentially expressed genes for each cluster")
markers_0.5 <- Seurat::FindAllMarkers(object= merged_seurat_0.5
                                      ,min.pct=0.25
                                      ,only.pos = TRUE
                                      ,logfc.threshold = 0.25
                                      ,verbose=FALSE)

markers_0.5.top10 <- markers_0.5 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Copy genes highlighted in yellow in the file E:\backup-genomeinfo\share\analysis\Bald_-_RNASeq\scRNA_Feb2020_CITE_VDJ\GEX_CiteSeq\analysis-results\Top Ranked genes for each modlue from tobis data analysis.xlsx
# Sergio picked the genes
rep.genes.cluster.0 <- c("IL7R","SELL","ZFP36","LTB","CCR7") # T cell Na?ve 
rep.genes.cluster.1 <- c("CCL5","GNLY","SOX4") # Cytotoxic
rep.genes.cluster.2 <- c("GZMK","EOMES","SAMD3","CD27") # Cytotoxic
rep.genes.cluster.3 <- c("XCL2","XCL1","IFNG","TNF") # Cytokine markers expressed in the CD 8 cluster
rep.genes.cluster.4 <- c("CCL3","CCL4","CXCL13","TIGIT","HAVCR2","ZBTB32") # Cytokine
rep.genes.cluster.5 <- c("IFI6","IFIT3","IFIT1","STAT1") # Cytokine markers in the interferon family
rep.genes.cluster.6 <- c("FOS") # Exhausted T cells
rep.genes.cluster.7 <- c("TRBV12-4","TRAV12-3","CTSW") # Exhausted T cells
rep.genes.cluster.8 <- c("CD300A","TCF7","KLRD1","IL7R") # T cell Na?ve
rep.genes.cluster.9 <- c() # proliferative T cells
rep.genes.cluster.10 <- c("KLRB1","RGS2","TRAV1-2") # MAIT
rep.genes.cluster.11 <- c("GNLY","S100A4") # Cytotoxic
rep.genes.cluster.12 <- c("CXCL8","TXN","FCER1G") # Cytotoxic

top.genes <- unique(c( rep.genes.cluster.0, rep.genes.cluster.1, rep.genes.cluster.2, rep.genes.cluster.3
                      ,rep.genes.cluster.4, rep.genes.cluster.5, rep.genes.cluster.6, rep.genes.cluster.7
                      ,rep.genes.cluster.8, rep.genes.cluster.10,rep.genes.cluster.11,rep.genes.cluster.12))
cat("Selected top genes:",top.genes)
length(top.genes)

# Choose 0.5 as resolution
merged_seurat_0.5 <- Seurat::FindClusters(merged_seurat, resolution = 0.5, verbose=FALSE)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
cat("Identify differentially expressed genes for each cluster")
merged_seurat.markers <- Seurat::FindAllMarkers(object= merged_seurat_0.5
                                                ,min.pct=0.25
                                                ,only.pos = TRUE
                                                ,logfc.threshold = 0.25
                                                ,verbose=FALSE)
head(merged_seurat.markers)

#---------------------------------------------------------------
# Assign cell types to clusters 
#---------------------------------------------------------------

# | Cell type | Short hand | Assigned to clusters |
# | --- | --- | --- |
# | Naive CD8 T cell | TN | 0, 8 |
# | Cytotoxic CD8 T cells | Cytotoxic | 1, 2, 11, 12 |
# | Cytokine CD8 T cells | Cytokine | 3 |
# | Exhausted CD8 T cells | Tex | 4, 6, 7 |
# | CD8 T cell producing interferon gamma | CD8 IFN | 5 |
# | Proliferative CD8 T cells | Proliferative | 9 |
# | Mucosal-associated invariant T cells | MAIT_cells | 10|

# Assign cell types to clusters
celltypes <- c("0" =  'TN'
               ,"1" =  'Cytotoxic'
               ,"2" =  'Cytotoxic'
               ,"3" =  'Cytokine'
               ,"4" =  'Tex' # Exhausted CD8 T cells
               ,"5" =  'CD8 IFN' # Cd8 T cell interferon gamma
               ,"6" =  'Tex'
               ,"7" =  'Tex'
               ,"8" =  'TN' # Naive CD8 T cells
               ,"9" =  'proliferative'# proliferative CD8 T cells
               ,"10" = 'MAIT_cells'# MAIT cells
               ,"11" ='Cytotoxic' # Cytotoxic T cell
               ,"12"= 'Cytotoxic')

names(celltypes) <- levels(merged_seurat_0.5)

# Rename clusters with cell types
merged_seurat_0.5 <- Seurat::RenameIdents(merged_seurat_0.5, celltypes)

merged_seurat_0.5[["celltypes"]] <- merged_seurat_0.5@active.ident

#--------------------------------------------------------------
# Subset top N differentially expressed genes in each cell type
#--------------------------------------------------------------
merged_seurat_0.5.markers <- Seurat::FindAllMarkers(merged_seurat_0.5
                                                    , min.pct = 0.1
                                                    , logfc.threshold = 0.25
                                                    , only.pos = TRUE) #only.pos = TRUE

# Get top 10 differentially expressed genes in each cell type (cluster here) 
topN <- merged_seurat_0.5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

Seurat::DoHeatmap( object=merged_seurat_0.5
                   ,features=topN$gene
                   ,group.colors= NULL
                   ,angle = 90
                   ,size = 11) +
  ggplot2::ggtitle("Expression heatmap for top 10 genes in each cell type")+
  ggplot2::theme(text = ggplot2::element_text(size = 15)) + 
  ggplot2::scale_fill_gradientn(colors = c("steelblue1", "white", "red3"))

#--------------------------------------------------------
# Plot expression levels of selected genes  in cell types
#--------------------------------------------------------
# Genes that Jason picked
genes.of.interest <- c("CD226","ENTPD1","PDCD1","CD96","TIGIT","TOX","EOMES","CXCL13","NKG7")
genes.of.interest %in% rownames(merged_seurat)

dot.plot.1 <- Seurat::DotPlot(merged_seurat_0.5, features = genes.of.interest) + 
  Seurat::RotatedAxis() + 
  ggplot2::labs(title = "Expression of selected genes in different cell types"
                ,x="Genes"
                ,y="Cell types")

dot.plot.2 <- Seurat::DotPlot(merged_seurat_0.5, features = genes.of.interest, split.by="condition") + 
  Seurat::RotatedAxis() + 
  ggplot2::labs(title = "Expression of selected genes in different cell types and conditions"
                ,x="Genes"
                ,y="Cell types")

patchwork::wrap_plots(dot.plot.1,dot.plot.2, ncol = 2)

#-------------------------------------------------------------------------------------
# Cluster cells into cell types based on the expression of selected genes and proteins
#-------------------------------------------------------------------------------------
Seurat::FeaturePlot(object=merged_seurat_0.5
                    ,features=c(paste0("rna_",genes.of.interest)
                                ,protein.selected.ADT)
                    ,reduction ="tsne"
                    ,label=TRUE
                    ,cols=c("lightgrey","red"))

Seurat::DoHeatmap( object=merged_seurat_0.5
                   ,features=topN$gene
                   ,group.by="rna.expr.cluster.nonzero.median.CD226.PDCD1"
                   ,group.colors= NULL
                   ,angle = 45
                   ,size = 11) +
  ggplot2::theme(text = ggplot2::element_text(size = 15)) + 
  ggplot2::scale_fill_gradientn(colors = c("steelblue1", "white", "red3"))

#------------------------------------------
# Complex heatmap with multiple annotations
#------------------------------------------
# Expression matrix
mat <- as.matrix(Seurat::GetAssayData(object= merged_seurat_0.5
                                      , slot= "counts")) # class(mat) # matrix # dim(mat) # 16298 10178
type <- colnames(mat)
class(type) # [1] "character"

# ComplexHeatmap package provides very flexible supports for setting annotations and defining new annotation graphics. The annotations can be put on the four sides of the heatmap, by top_annotation, bottom_annotation, left_annotation and right_annotation arguments.
# The value for the four arguments should be in the HeatmapAnnotation class and should be constructed by HeatmapAnnotation() function, or by rowAnnotation() function if it is a row annotation.
ha <- ComplexHeatmap::HeatmapAnnotation(type = type, annotation_name_side = "left")

ht_list <- ComplexHeatmap::Heatmap( matrix= mat
                                    ,name = "expression" # Name the matrix
                                    , km = 5 # Apply k-means clustering on rows
                                    , top_annotation = ha
                                    ,show_row_names = FALSE
                                    ,show_column_names = FALSE)
draw(ht_list)





