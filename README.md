###### tags: #R #scRNA-seq #CITE-Seq #PCA #t-SNE #UMAP

# Single-Cell RNA Sequencing and CITE-Seq Analysis

This repository contains resources and scripts for analyzing single-cell RNA sequencing (scRNA-seq) data, focusing on integrating RNA and surface protein expression data using Cellular Indexing of Transcriptomes and Epitopes by Sequencing (CITE-Seq). The primary goal is to provide reproducible workflows for scRNA-seq data analysis with practical insights into handling biological datasets.

## Overview

The main analysis script, [`Bald-scRNAseq-HNSCC.R.md`](https://github.com/luenhchang/single_cell_RNA_sequencing_CITE_seq/blob/main/Bald-scRNAseq-HNSCC.R.md), describes the step-by-step approach for analyzing scRNA-seq data using R. It includes:
- Preprocessing and quality control of scRNA-seq data.
- Dimensionality reduction techniques like Principal Component Analysis (PCA), t-distributed Stochastic Neighbor Embedding (t-SNE), Uniform Manifold Approximation and Projection (UMAP).
- Clustering and identification of cell types.
- Integration of CITE-Seq protein markers.
- Visualization of results through high-quality plots.

## Repository Structure

- **`Bald-scRNAseq-HNSCC.R.md`**: Markdown version of the R script, which provides a detailed narrative of the analysis pipeline, with code snippets and outputs.
- **`Bald-scRNAseq-HNSCC.R_images/`**: A directory containing images and plots generated during the analysis, referenced in the `.md` file.

## How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/luenhchang/single_cell_RNA_sequencing_CITE_seq.git
