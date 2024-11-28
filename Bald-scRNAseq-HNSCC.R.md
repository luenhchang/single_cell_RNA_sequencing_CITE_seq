# Analyse single cell RNA seqeuncing and CITE sequencing on human head and neck cancer

#### File name: Bald-scRNAseq-HNSCC.R.ipynb

#### Date created: 20210513

#### Programmer: Chang

#### Reference

[Chapter 3 Heatmap Annotations](https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html)

[Changing R plot options in Jupyter](https://notebook.community/andrie/jupyter-notebook-samples/Changing%20R%20plot%20options%20in%20Jupyter)

[Changing R plot options in Jupyter](https://notebook.community/andrie/jupyter-notebook-samples/Changing%20R%20plot%20options%20in%20Jupyter)

[Using a new windows version of R in Jupyter notebooks](https://stackoverflow.com/questions/51647561/using-a-new-windows-version-of-r-in-jupyter-notebooks)

[scRNAseq_Braun_etal_Immunity_Script.Rmd](https://github.com/BaldLab/2020_Braun_et_al_CD226_Immunity/blob/master/scRNAseq_Mouse/R/scRNAseq_Braun_etal_Immunity_Script.Rmd)

[Introduction to scRNA-seq integration](https://satijalab.org/seurat/articles/integration_introduction.html)

[Using dittoSeq to visualize (sc)RNAseq data](https://bioconductor.org/packages/devel/bioc/vignettes/dittoSeq/inst/doc/dittoSeq.html)

[Setup a Seurat object, add the RNA and protein data](https://satijalab.org/seurat/articles/multimodal_vignette.html)

[scRNA-seq/lessons](https://github.com/hbctraining/scRNA-seq/tree/master/lessons)

[How to solve Error: cannot allocate vector of size 1.2 Gb in R?](https://www.researchgate.net/post/How_to_solve_Error_cannot_allocate_vector_of_size_12_Gb_in_R)


```R
dir.C <- "C:"

dir.R.packages <- "C:/Program Files/R/R-4.0.3/library" #"C:/Program Files/R/R-4.0.2/library"

library(dplyr, lib.loc = dir.R.packages)
```

    
    Attaching package: 'dplyr'
    
    
    The following objects are masked from 'package:stats':
    
        filter, lag
    
    
    The following objects are masked from 'package:base':
    
        intersect, setdiff, setequal, union
    
    
    


```R
dir.Bald.scRNA.202002.analysis <- "E:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/analysis-results"
```


```R
# merged_seurat <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"merged_seurat.rds"))

GEX1 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filtered_ADT-quantiledGEX1.rds"))
GEX2 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filtered_ADT-quantiledGEX2.rds")) 
GEX3 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filtered_ADT-quantiledGEX3.rds")) 
GEX4 <- readRDS(file = file.path(dir.Bald.scRNA.202002.analysis,"Bald-HNSCC-scRNAseq-CITEseq_filtered_ADT-quantiledGEX4.rds")) # class(GEX4)

# Edit meta data
GEX1@meta.data$condition <- "unstimulated"
GEX2@meta.data$condition <- "unstimulated"
GEX3@meta.data$condition <- "stimulated"
GEX4@meta.data$condition <- "stimulated"
```

    Loading required package: SeuratObject
    
    Warning message:
    "package 'SeuratObject' was built under R version 4.0.5"
    


```R
# Merge all samples of the same tissue type
merged_seurat <- merge( x=GEX1
                        ,y=c(GEX2, GEX3, GEX4))
```

    Warning message in CheckDuplicateCellNames(object.list = objects):
    "Some cell names are duplicated across objects provided. Renaming to enforce unique cell names."
    


```R
#---------------------------------------------------------------------------------------------------------------
# Cluster cells on the basis of their scRNA-seq profiles
## reference [Using Seurat with multimodal data](https://satijalab.org/seurat/articles/multimodal_vignette.html)
#---------------------------------------------------------------------------------------------------------------

# View the keys for all keyed components (assays, dimensional reduction, spatial images) of a Seurat object using the Key functio
cat("Keys to Seurat elements",Seurat::Key(merged_seurat))
#    RNA    ADT 
# "rna_" "adt_" 

# DefaultAssay(merged_seurat) # "RNA"
merged_seurat <- Seurat::NormalizeData(object = merged_seurat
                                       ,assay="RNA"
                                       ,verbose=FALSE) %>%
  Seurat::FindVariableFeatures(verbose=FALSE) %>%
  Seurat::ScaleData(verbose=FALSE) %>%
  Seurat::RunPCA(verbose=FALSE)

# Change default R plot size
options(repr.plot.width = 4*5, repr.plot.height = 3*5)

# Clustering and dimension (PCs) selection
# Check dimensions
## A dimension is a PC. Choose PC23 as the break point. Take the more conservative PC.
Seurat::ElbowPlot(object = merged_seurat
                  ,ndims = 30
                  ,reduction = "pca")
```

    Keys to Seurat elements rna_ adt_


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_5_1.png)
    



```R
merged_seurat <- Seurat::FindNeighbors(object = merged_seurat, dims=1:23, verbose=FALSE)

merged_seurat <- Seurat::FindClusters(object = merged_seurat, resolution=c(0.5, 0.9), verbose=FALSE)

merged_seurat <- Seurat::RunTSNE(merged_seurat, assay="RNA",dims=1:23, verbose=FALSE) 

merged_seurat <- Seurat::RunUMAP(merged_seurat, assay="RNA", dims=1:23, verbose=FALSE)
```

    Warning message:
    "The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session"
    


```R
expansion.factor <- 5
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 3*expansion.factor)
cat("Visualising clustering by t-SNE")
Seurat::DimPlot(merged_seurat, reduction = "tsne", label = TRUE, pt.size=2, label.size=10)
```

    Visualising clustering by t-SNE


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_7_1.png)
    



```R
cat("Visualising clustering by UMAP")
Seurat::DimPlot(merged_seurat, reduction = "umap", label = TRUE, pt.size=2, label.size=10)
```

    Visualising clustering by UMAP


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_8_1.png)
    



```R
#-------------------------
# Normalise assay ADT
#-------------------------
merged_seurat <- Seurat::NormalizeData(object = merged_seurat
                                       , normalization.method="CLR"
                                       , assay="ADT"
                                       , margin=2 # normalise data across features (margin=1) or cells (margin=2)
                                       , verbose=FALSE
                                       ) %>%
  Seurat::ScaleData(verbose=FALSE) %>%
  Seurat::RunPCA(reduction.name='adtpca',verbose=FALSE)
```

    Warning message:
    "Cannot add objects with duplicate keys (offending key: PC_), setting key to 'adtpca_'"
    


```R
#-------------------------------------------
# Visualize multiple modalities side-by-side
#-------------------------------------------
# Find out the prefixes of different assays
## Features (i.e., genes, protein) can be referred as assayKey_featureName
Seurat::Key(merged_seurat[["RNA"]]) # [1] "rna_"
Seurat::Key(merged_seurat[["ADT"]]) # [1] "adt_"

# Compare gene expression between conditions and clusters
genes.selected <- c("CD226","PDCD1","CD28","ENTPD1","CD96","TIGIT") # ,"EOMES"
genes.selected %in% rownames(merged_seurat) # Should be all TRUE
genes.selected.RNA <- paste0(Seurat::Key(merged_seurat[["RNA"]]), genes.selected) # Include only genes that can be found in the rownames() # ,"CD279" ,"CD39"

# protein feature names
protein.selected <- rownames(merged_seurat@assays$ADT)
protein.selected.ADT <- paste0(Seurat::Key(merged_seurat[["ADT"]]), protein.selected)
```


'rna_'



'adt_'



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>




```R
# Compare gene expression of all genes between conditions
cat("Compare gene expression of all genes between conditions stimulated and unstimulated")
dittoSeq::dittoDimPlot(object=merged_seurat
                       ,var = "condition"
                       ,reduction.use = "umap"
                       ,assay = "RNA"
                      ,size=2)
```

    Compare gene expression of all genes between conditions stimulated and unstimulated


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_11_1.png)
    



```R
# Plot expression of CD226 at protein and RNA levels
## gene name is specified as assayName_featureName (e.g., rna_CD226, adtcd226_CD226)
featurePlot.1.rnaCD226 <- Seurat::FeaturePlot(object = merged_seurat
                                              ,features = genes.selected.RNA[1]
                                              , cols = c("lightgrey", "darkgreen")) + ggplot2::ggtitle("CD226 RNA")

featurePlot.1.adtCD226 <- Seurat::FeaturePlot(object = merged_seurat
                                              ,features = protein.selected.ADT[1]
                                              , cols = c("lightgrey", "darkgreen")) + ggplot2::ggtitle("CD226 protein")
# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 2*expansion.factor)

featurePlot.1.adtCD226 | featurePlot.1.rnaCD226
```


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_12_0.png)
    



```R
#------------------------------------------------
# Expression of single gene between conditions in assay RNA
#------------------------------------------------
## Arrange multiple violin plots into a grid using cowplot::plot_grid()
cowplot::plot_grid(
     dittoSeq::dittoPlot(object=merged_seurat
                    ,var=genes.selected[1] # var argument cannot take assayName_featureName
                    ,assay = "RNA"
                    ,group.by = "condition"
                    ,plots = c("vlnplot", "jitter"))
    ,dittoSeq::dittoPlot(object=merged_seurat
                    ,var=protein.selected[1]
                    ,assay = "ADT"
                    ,group.by = "condition"
                    ,plots = c("vlnplot", "jitter"))
    ,ncol=2)
```


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_13_0.png)
    



```R
# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 6*expansion.factor)

# Plot RNA and protein side by side
cat("Violin plots for expression levels of RNA and protein side by side")
## Arrange multiple violin plots into a grid using cowplot::plot_grid()
cowplot::plot_grid(
  # CD226 RNA and protein
   Seurat::VlnPlot(object=merged_seurat, features = genes.selected.RNA[1])
  ,Seurat::VlnPlot(object=merged_seurat, features = protein.selected.ADT[1])
  # PDCD1 RNA and protein CD279-PD1-TotalSeqC
  ,Seurat::VlnPlot(object=merged_seurat, features = genes.selected.RNA[2])
  ,Seurat::VlnPlot(object=merged_seurat, features = protein.selected.ADT[2])
  # CD28 RNA and protein
  ,Seurat::VlnPlot(object=merged_seurat, features = genes.selected.RNA[3])
  ,Seurat::VlnPlot(object=merged_seurat, features = protein.selected.ADT[3])
  ,ncol=2)
```

    Violin plots for expression levels of RNA and protein side by side


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_14_1.png)
    



```R
# Plot RNA and protein side by side
cat("Violin plots for expression levels of RNA and protein side by side")
## Arrange multiple violin plots into a grid using cowplot::plot_grid()
cowplot::plot_grid(
  # ENTPD1 RNA and protein CD39-TotalSeqC
   Seurat::VlnPlot(object=merged_seurat, features = genes.selected.RNA[4])
  ,Seurat::VlnPlot(object=merged_seurat, features = protein.selected.ADT[4])
   # CD 96 RNA and CD96TACTILE-TotalSeqC protein 
  ,Seurat::VlnPlot(object=merged_seurat, features = genes.selected.RNA[5])
  ,Seurat::VlnPlot(object=merged_seurat, features = protein.selected.ADT[5])
   # TIGIT RNA and protein TIGITVSTM3-TotalSeqC
  ,Seurat::VlnPlot(object=merged_seurat, features = genes.selected.RNA[6])
  ,Seurat::VlnPlot(object=merged_seurat, features = protein.selected.ADT[6]) 
  ,ncol=2)
```

    Violin plots for expression levels of RNA and protein side by side


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_15_1.png)
    



```R
# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 3*expansion.factor)
cat("Violin plots for RNA gene expression between conditions stimulated and unstimulated")
dittoSeq::multi_dittoPlot(object= merged_seurat
                          , assay="RNA"
                          , var= genes.selected # c(genes.selected.RNA, genes.selected.ADT) 
                          , group.by = "condition"
                          , vlnplot.lineweight = 0.3
                          , jitter.size = 1
                          ,legend.show = FALSE
                          ,xlab = "Conditions"
                          ,ylab = "Expression") # CD226, TIGIT have higher expression in stim than unsti; EOMES
```

    Violin plots for RNA gene expression between conditions stimulated and unstimulated


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_16_1.png)
    



```R
# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 6*expansion.factor)

# CD226, CDPD1 gene and protein expression between conditions
cat("tSNE plots for expression of RNA and protien between conditions stimulated and unstimulated")
Seurat::FeaturePlot(object=merged_seurat
                    , features = c( genes.selected.RNA[1], protein.selected.ADT[1]
                                   ,genes.selected.RNA[2], protein.selected.ADT[2])
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    , split.by = "condition"
                    ) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```

    tSNE plots for expression of RNA and protien between conditions stimulated and unstimulated


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_17_1.png)
    



```R
# CD28, ENTPD1 gene and protein expression between conditions
cat("tSNE plots for expression of RNA and protien between conditions stimulated and unstimulated")
Seurat::FeaturePlot(object=merged_seurat
                    , features = c( genes.selected.RNA[3], protein.selected.ADT[3]
                                   ,genes.selected.RNA[4], protein.selected.ADT[4])
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    , split.by = "condition"
                    ) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```

    tSNE plots for expression of RNA and protien between conditions stimulated and unstimulated


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_18_1.png)
    



```R
# CD96, TIGIT gene and protein expression between conditions
cat("tSNE plots for expression of RNA and protien between conditions stimulated and unstimulated")
Seurat::FeaturePlot(object=merged_seurat
                    , features = c( genes.selected.RNA[5], protein.selected.ADT[5]
                                   ,genes.selected.RNA[6], protein.selected.ADT[6])
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    , split.by = "condition"
                    ) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```

    tSNE plots for expression of RNA and protien between conditions stimulated and unstimulated


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_19_1.png)
    


## Identify potential marker genes for each cluster

This type of analysis is typically recommended for when evaluating a single sample group/condition. With the ` FindAllMarkers()` function we are comparing each cluster against all other clusters to identify potential marker genes. The cells in each cluster are treated as replicates, and essentially a differential expression analysis is performed with some statistical test.

NOTE: The default is a Wilcoxon Rank Sum test, but there are other options available.

![marker_ident_function1.png](Bald-scRNAseq-HNSCC.R_images/marker_ident_function1.png)


```R
# Find markers for every cluster compared to all remaining cells, report only the positive ones
cat("Identify differentially expressed genes for each cluster")
markers <- Seurat::FindAllMarkers(object=merged_seurat
                                  ,min.pct=0.25
                                  ,only.pos = TRUE
                                  ,logfc.threshold = 0.25
                                  ,verbose=FALSE)
```

    Identify differentially expressed genes for each cluster

## Identify gene markers that are conserved between the groups

Since we have samples representing different conditions in our dataset, our best option is to find conserved markers. This function internally separates out cells by sample group/condition, and then performs differential gene expression testing for a single specified cluster against all other clusters (or a second cluster, if specified). Gene-level p-values are computed for each condition and then combined across groups using meta-analysis methods from the MetaDE R package.

![marker_ident_function2.png](Bald-scRNAseq-HNSCC.R_images/marker_ident_function2.png)


```R
cluster0_conserved_markers <- Seurat::FindConservedMarkers(merged_seurat
                                                           ,ident.1 = 0
                                                           ,grouping.var = "condition"
                                                           ,only.pos = TRUE
                                                           ,logfc.threshold = 0.25)
head(cluster0_conserved_markers)
```

    Testing group unstimulated: (0) vs (10, 4, 7, 8, 6, 9, 5, 12, 13, 3, 14, 15, 1, 2, 11)
    
    Testing group stimulated: (0) vs (2, 7, 3, 5, 1, 6, 10, 11, 4, 13, 9, 14, 15, 12, 8)
    
    


<table class="dataframe">
<caption>A data.frame: 6 × 12</caption>
<thead>
	<tr><th></th><th scope=col>unstimulated_p_val</th><th scope=col>unstimulated_avg_log2FC</th><th scope=col>unstimulated_pct.1</th><th scope=col>unstimulated_pct.2</th><th scope=col>unstimulated_p_val_adj</th><th scope=col>stimulated_p_val</th><th scope=col>stimulated_avg_log2FC</th><th scope=col>stimulated_pct.1</th><th scope=col>stimulated_pct.2</th><th scope=col>stimulated_p_val_adj</th><th scope=col>max_pval</th><th scope=col>minimump_p_val</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GZMK</th><td>5.659935e-265</td><td>1.741423</td><td>0.825</td><td>0.296</td><td>9.224562e-261</td><td>1.522696e-71</td><td>1.9109346</td><td>0.732</td><td>0.266</td><td>2.481690e-67</td><td>1.522696e-71</td><td>1.131987e-264</td></tr>
	<tr><th scope=row>CMC1</th><td>8.986527e-163</td><td>1.509320</td><td>0.618</td><td>0.242</td><td>1.464624e-158</td><td>3.240914e-48</td><td>1.3044352</td><td>0.614</td><td>0.243</td><td>5.282042e-44</td><td>3.240914e-48</td><td>1.797305e-162</td></tr>
	<tr><th scope=row>CST7</th><td>9.525578e-154</td><td>1.074360</td><td>0.951</td><td>0.770</td><td>1.552479e-149</td><td>1.327403e-34</td><td>0.9268426</td><td>0.950</td><td>0.837</td><td>2.163401e-30</td><td>1.327403e-34</td><td>1.905116e-153</td></tr>
	<tr><th scope=row>EOMES</th><td>1.408120e-151</td><td>1.300099</td><td>0.538</td><td>0.175</td><td>2.294955e-147</td><td>5.301662e-37</td><td>1.2112108</td><td>0.423</td><td>0.140</td><td>8.640648e-33</td><td>5.301662e-37</td><td>2.816241e-151</td></tr>
	<tr><th scope=row>RPS15A</th><td>1.118059e-122</td><td>0.462811</td><td>0.999</td><td>0.997</td><td>1.822213e-118</td><td>3.793634e-33</td><td>0.5146104</td><td>1.000</td><td>0.996</td><td>6.182865e-29</td><td>3.793634e-33</td><td>2.236118e-122</td></tr>
	<tr><th scope=row>DKK3</th><td>1.390575e-108</td><td>1.218504</td><td>0.329</td><td>0.072</td><td>2.266358e-104</td><td>1.375378e-53</td><td>1.2423714</td><td>0.309</td><td>0.055</td><td>2.241592e-49</td><td>1.375378e-53</td><td>2.781149e-108</td></tr>
</tbody>
</table>




```R
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
```


```R
# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 2.5*expansion.factor)

Seurat::CombinePlots(plots = list(p1, p2,p3,p4,p5,p6))
```

    Warning message:
    "CombinePlots is being deprecated. Plots should now be combined using the patchwork system."
    


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_25_1.png)
    



```R
Seurat::CombinePlots(plots = list(p7,p8,p9,p10,p11,p12))
```

    Warning message:
    "CombinePlots is being deprecated. Plots should now be combined using the patchwork system."
    


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_26_1.png)
    



```R
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
```


```R
Seurat::CombinePlots(plots = list(p1, p2,p3,p4,p5,p6))
```

    Warning message:
    "CombinePlots is being deprecated. Plots should now be combined using the patchwork system."
    


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_28_1.png)
    



```R
Seurat::CombinePlots(plots = list(p7,p8,p9,p10,p11,p12))
```

    Warning message:
    "CombinePlots is being deprecated. Plots should now be combined using the patchwork system."
    


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_29_1.png)
    



```R
# Find markers for every cluster compared to all remaining cells, report only the positive ones
cat("Identify differentially expressed genes for each cluster")
markers_0.5 <- Seurat::FindAllMarkers(object= merged_seurat_0.5
                                  ,min.pct=0.25
                                  ,only.pos = TRUE
                                  ,logfc.threshold = 0.25
                                  ,verbose=FALSE)
head(markers_0.5)
```

    Identify differentially expressed genes for each cluster


<table class="dataframe">
<caption>A data.frame: 6 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>IL7R</th><td>9.678977e-268</td><td>1.3660222</td><td>0.856</td><td>0.519</td><td>1.577480e-263</td><td>0</td><td>IL7R </td></tr>
	<tr><th scope=row>RPS12</th><td>1.887499e-254</td><td>0.5962669</td><td>0.998</td><td>0.998</td><td>3.076247e-250</td><td>0</td><td>RPS12</td></tr>
	<tr><th scope=row>RPL32</th><td>9.104873e-237</td><td>0.5486433</td><td>0.995</td><td>0.998</td><td>1.483912e-232</td><td>0</td><td>RPL32</td></tr>
	<tr><th scope=row>RPL10</th><td>2.111811e-234</td><td>0.5273115</td><td>0.995</td><td>0.999</td><td>3.441830e-230</td><td>0</td><td>RPL10</td></tr>
	<tr><th scope=row>RPS3A</th><td>1.218826e-225</td><td>0.5881000</td><td>0.991</td><td>0.996</td><td>1.986442e-221</td><td>0</td><td>RPS3A</td></tr>
	<tr><th scope=row>RPL34</th><td>9.950714e-223</td><td>0.5902106</td><td>0.984</td><td>0.995</td><td>1.621767e-218</td><td>0</td><td>RPL34</td></tr>
</tbody>
</table>




```R
write.table(markers_0.5, "E:/backup-genomeinfo/share/analysis/Bald_-_RNASeq/scRNA_Feb2020_CITE_VDJ/GEX_CiteSeq/analysis-results/FindAllMarkers_resolution-0.5.txt"
            , sep = "\t", row.names = F)
```


```R
top10 <- markers_0.5 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
head(top10)
```


<table class="dataframe">
<caption>A grouped_df: 6 × 7</caption>
<thead>
	<tr><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>9.678977e-268</td><td>1.3660222</td><td>0.856</td><td>0.519</td><td>1.577480e-263</td><td>0</td><td>IL7R   </td></tr>
	<tr><td>1.887499e-254</td><td>0.5962669</td><td>0.998</td><td>0.998</td><td>3.076247e-250</td><td>0</td><td>RPS12  </td></tr>
	<tr><td>6.742727e-164</td><td>1.0233487</td><td>0.942</td><td>0.918</td><td>1.098930e-159</td><td>0</td><td>PABPC1 </td></tr>
	<tr><td>5.589620e-146</td><td>1.3356441</td><td>0.982</td><td>0.985</td><td>9.109963e-142</td><td>0</td><td>FTH1   </td></tr>
	<tr><td>4.438425e-135</td><td>0.8232958</td><td>0.850</td><td>0.680</td><td>7.233745e-131</td><td>0</td><td>ZFP36L2</td></tr>
	<tr><td>2.077494e-100</td><td>0.7098626</td><td>0.818</td><td>0.665</td><td> 3.385900e-96</td><td>0</td><td>ZFP36  </td></tr>
</tbody>
</table>




```R
Seurat::DoHeatmap( subset(merged_seurat_0.5, downsample = 100)
                  ,features = top10$gene
                  , size = 9
                  ,assay = "RNA")
```

    Warning message in Seurat::DoHeatmap(subset(merged_seurat_0.5, downsample = 100), :
    "The following features were omitted as they were not found in the scale.data slot for the RNA assay: SPOCK2, SPINT2, SAMD3, CLDND1, CST7, ADK, RPS12"
    


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_33_1.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C01-LEF1 cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
genes.repre.CD8.cluster.1 <- c("CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1")
genes.repre.CD8.cluster.1 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 4*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C01-LEF1 cluster")
Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.1
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
                    ) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C01-LEF1 cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_34_2.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C02-GPR183 cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
genes.repre.CD8.cluster.2 <- c("CCR7","SELL","IL7R","CD27","CD28","PRF1","GZMA","CCL5","GPR183","S1PR1")
genes.repre.CD8.cluster.2 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 5*expansion.factor, repr.plot.height = 4*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C02-GPR183 cluster")

Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.2
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C02-GPR183 cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_35_2.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C03-CX3CR1 cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
genes.repre.CD8.cluster.3 <- c("KLRG1","CX3CR1","FCGR3A","FGFBP2","PRF1","GZMH","TBX21","EOMES","S1PR1","S1PR5")
genes.repre.CD8.cluster.3 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 5*expansion.factor, repr.plot.height = 4*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C03-CX3CR1 cluster")

Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.3
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C03-CX3CR1 cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_36_2.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C04-GZMK cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
genes.repre.CD8.cluster.4 <- c("GZMK","CXCR4","CXCR3","CD44")
genes.repre.CD8.cluster.4 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 4*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C04-GZMK cluster")

Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.4
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C04-GZMK cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_37_2.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C05-CD6 cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
# The gene NR4A1/2/3 is not found. Close genes are [1] "NR4A2" "NR4A3" "NR4A1"
genes.repre.CD8.cluster.5 <- c("CD6","XCL1","XCL2","MYADM","CAPG","RORA","NR4A1","NR4A2","NR4A3","CD69","ITGAE")
genes.repre.CD8.cluster.5 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 5*expansion.factor, repr.plot.height = 4*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C05-CD6 cluster")

Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.5
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C05-CD6 cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_38_2.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C06-CD160 cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
grep(pattern = "KLRC.*", x=rownames(merged_seurat), value = TRUE)
# The gene KLRC1/2/3 is not found. Close genes are [1] "KLRC4-KLRK1" "KLRC4" "KLRC3" "KLRC2" "KLRC1"  
# The gene NR4A1/2/3 is not found. Close genes are [1] "NR4A2" "NR4A3" "NR4A1"

genes.repre.CD8.cluster.6 <- c("CD160","KIR2DL4","TMIGD2","KLRC1","KLRC2","KLRC3","NR4A1","NR4A2","NR4A3","IKZF2","ENTPD1","CD69","ITGAE")
genes.repre.CD8.cluster.6 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 4*expansion.factor, repr.plot.height = 4*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C06-CD160 cluster")

Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.6
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'KLRC4-KLRK1'</li><li>'KLRC4'</li><li>'KLRC3'</li><li>'KLRC2'</li><li>'KLRC1'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C06-CD160 cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_39_3.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C07-LAYN cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
genes.repre.CD8.cluster.7 <- c("HAVCR2","CXCL13","PDCD1","LAYN","TOX","IFNG","GZMB","MIR155HG","TNFRSF9","ITGAE")
genes.repre.CD8.cluster.7 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 5*expansion.factor, repr.plot.height = 4*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C07-LAYN cluster")

Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.7
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C07-LAYN cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_40_2.png)
    



```R
#---------------------------------------------------------------------------------------------------
# Plot expression of representative genes in CD8_C08-SLC4A10 cluster 
# reported in Supplementary table 5, Zhang et al 2018 paper
#---------------------------------------------------------------------------------------------------
genes.repre.CD8.cluster.8 <- c("SLC4A10","KLRB1","ZBTB16","NCR3","RORC","RORA")
genes.repre.CD8.cluster.8 %in% rownames(merged_seurat) # All TRUE

# Chang R plot sizes
expansion.factor <- 4
options(repr.plot.width = 3*expansion.factor, repr.plot.height = 3.5*expansion.factor)

cat("Expression of representative genes reported in Zhang's CD8_C08-SLC4A10 cluster")

Seurat::FeaturePlot(object=merged_seurat
                    , features = genes.repre.CD8.cluster.8
                    , reduction = "tsne"
                    , label = TRUE
                    , label.size = 7.5
                    #, split.by = "condition"
) & ggplot2::theme(text = ggplot2::element_text(size = 25, face = "bold")) +
  ggplot2::theme(legend.text=ggplot2::element_text(size=25, color="black"))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li><li>TRUE</li></ol>



    Expression of representative genes reported in Zhang's CD8_C08-SLC4A10 cluster


    
![png](Bald-scRNAseq-HNSCC.R_images/Bald-scRNAseq-HNSCC.R_41_2.png)
    

