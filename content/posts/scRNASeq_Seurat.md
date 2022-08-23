---
title: "scRNASeq analysis with Seurat "
date: 2020-05-02
draft: false
tags: ["R", "scRNASeq", "tutorial"]
categories: ["Tutorial"]
author: "Yaoyu E. Wang"

autoCollapseToc: true

---


# scRNASeq analysis with Seurat and other tools

## Seurat install

The Seurat package will need to be installed. You can find further instructions for install found [here](https://satijalab.org/seurat/install.html). It is installed via a separate package manager (in R) Bioconductor and R's package manager CRAN.


```R
# Check if Bioconductor is installed. If not, install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# Install a necessary package for Seurat
BiocManager::install('multtest')
BiocManager::install('MAST')
# Install Seurat
install.packages('Seurat')
```

    Updating HTML index of packages in '.Library'
    
    Making 'packages.html' ...
     done
    
    Bioconductor version 3.10 (BiocManager 1.30.10), R 3.6.2 (2019-12-12)
    
    Installing package(s) 'BiocVersion', 'multtest'
    
    also installing the dependencies 'BiocGenerics', 'Biobase'
    
    
    Bioconductor version 3.10 (BiocManager 1.30.10), R 3.6.2 (2019-12-12)
    
    Installing package(s) 'MAST'
    
    Updating HTML index of packages in '.Library'
    
    Making 'packages.html' ...
     done

## Setup of the scRNASeq data into a Seurat object

For this first part of the tutorial, we will be analyzing a small dataset of peripheral blood mononuclear cells (PBMC) provided by 10X Genomics. This is a data set of 2,700 single cells sequenced on Illumina NextSeq 500. You can get the raw data [here](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).

When starting your RStudio session, please move to a separate folder for this project. I will assume you named it **seurat_tutorial**. Then please make a folder named **data** in the **seurat_tutorial** folder and place the downloaded 10X data in the **data** folder. It should be named **pbmc3k_filtered_gene_bc_matrices.tar.gz**. The file will need to be untarred.

* On Windows, the free and open source program 7zip can be used to decompress this file
* On Mac (and Linux), you should be able to decompress this in Finder in the GUI by default.

```R
setwd("~/")
```

*Note that* `setwd` *stands for set working directory.*

Now that the environment is set up, we can start by loading the raw data into Seurat so that it is formatted to work with the rest of the Seurat toolset.


```R
# Load the libraries

# Library for multtest; unneeded as Seurat should automatically load it
library(multtest)

# Library for Seurat, the scRNASeq analysis toolkit
library(Seurat)

# Library for MAST, the differential gene analysis algorithm for
# scRNASeq. Needs to be loaded to work in Seurat
library(MAST)

# Library for plotting in R
library(ggplot2)

# Library for making heatmaps with Seurat
library(patchwork)

# Library Seurat uses on top of ggplot
# to plot things easily onto a grid
library(cowplot)

```


```R
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
# counts       : your loaded count data
# project      : the name of this object
# min.cells    : assumed lower threshold for cells
# min.features : assumed lower threhold for genes
pbmc <- CreateSeuratObject(
    counts = pbmc.data,
    project = "pbmc3k",
    min.cells = 3,
    min.features = 200
)
```

    Warning message:
    “Feature names cannot have underscores ('_'), replacing with dashes ('-')”


The data for the PBMC scRNASeq is split between three files:
* barcodes.tsv
* genes.tsv
* matrix.mtx

The three combined work together to define a count matrix. Seurat then converts the data into a single object representing the count matrix. Now that the data has been loaded into a Seurat object, what does the data look like.


```R
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
```

       [[ suppressing 30 column names ‘AAACATACAACCAC’, ‘AAACATTGAGCTAC’, ‘AAACATTGATCAGC’ ... ]]
    



    3 x 30 sparse Matrix of class "dgCMatrix"
                                                                       
    CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
    TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
    MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .


The `.` values represent 0s (or no molecules detected). scRNASeq is naturally sparse - i.e. we expect to see a lot of 0 values. This has major effects on the type of analyses we can perform on the data, which is what scRNASeq analyses attempt to deal with.

Additionally, we should become familiar with the meta data features of a Seurat object. The meta data are annotations applied to each cell. We can view them with:


```R
head(pbmc@meta.data)
```


<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_RNA</th><th scope=col>nFeature_RNA</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACATACAACCAC</th><td>pbmc3k</td><td>2419</td><td> 779</td></tr>
	<tr><th scope=row>AAACATTGAGCTAC</th><td>pbmc3k</td><td>4903</td><td>1352</td></tr>
	<tr><th scope=row>AAACATTGATCAGC</th><td>pbmc3k</td><td>3147</td><td>1129</td></tr>
	<tr><th scope=row>AAACCGTGCTTCCG</th><td>pbmc3k</td><td>2639</td><td> 960</td></tr>
	<tr><th scope=row>AAACCGTGTATGCG</th><td>pbmc3k</td><td> 980</td><td> 521</td></tr>
	<tr><th scope=row>AAACGCACTGGTAC</th><td>pbmc3k</td><td>2163</td><td> 781</td></tr>
</tbody>
</table>



These are the default annotations in the Seurat object. We can see we've labeled every cell with the project `pbmc3k`. Additionally, the total number of RNA molecules sequenced for that cell and the total number of genes (a.k.a features) detected within the cell.

The real utility comes from adding annotations to these cells. Let us add our cell type annotations from scMatch.


```R
# First read the file in as a data frame with R's built in
# read.csv() function.
pbmc.scmatch <- read.csv(
    "./data/pbmc3k_filtered_gene_bc_matrices/pbmc1k_scmatch.csv",
    row.names = 1
)
# Now we must rename the cells, as scMatch adds a '-1' to each cell ID
# We reassign all rownames of the pbmc3k.scmatch data frame (via row.names())
# and by sending it the new list.
# The new list is composed of sub(), which replaces all instances of "-1" 
# with nothing for all rownames
# Note that row.names() is different than rownames() and works in a subtlely 
# different manner. For our purposes here, they work exactly the same.
row.names(pbmc.scmatch) <- sub("-1", "", rownames(pbmc.scmatch))
```

The makeup of the scmatch data looks as so:


```R
head(pbmc.scmatch)
```


<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>cell.type</th><th scope=col>top.sample</th><th scope=col>top.correlation.score</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACATACAACCAC</th><td>CD8+ T cell       </td><td>CD8+ T Cells (pluriselect), donor090309, donation3.CNhs12180.12196-129B9  </td><td>0.2694394</td></tr>
	<tr><th scope=row>AAACATTGAGCTAC</th><td>B cell            </td><td>CD19+ B Cells (pluriselect), donor090309, donation2.CNhs12179.12194-129B7 </td><td>0.3398503</td></tr>
	<tr><th scope=row>AAACATTGATCAGC</th><td>CD4+ T cell       </td><td>CD4+CD25-CD45RA- memory conventional T cells, donor2.CNhs13237.11798-124C7</td><td>0.3238483</td></tr>
	<tr><th scope=row>AAACCGTGCTTCCG</th><td>monocyte          </td><td>CD14-CD16+ Monocytes, donor1.CNhs13229.11790-124B8                        </td><td>0.3081723</td></tr>
	<tr><th scope=row>AAACCGTGTATGCG</th><td>gamma-delta T cell</td><td>gamma delta positive T cells, donor1.CNhs13914.11937-126A2                </td><td>0.2184953</td></tr>
	<tr><th scope=row>AAACGCACTGGTAC</th><td>CD4+ T cell       </td><td>CD4+CD25+CD45RA- memory regulatory T cells, donor2.CNhs13206.11797-124C6  </td><td>0.2730377</td></tr>
</tbody>
</table>



What is required to add any annotations to the Seurat object is that we have some reference to the cell IDs. This scmatch dataframe contains the information for the cell IDs as the rownames. When we have this information, we can add the annotations to the Seurat object as metadata.


```R
# Now we add these annotations to the Seurat object
# The subset function used here only select the column in the scMatch
# dataset that we want to add. We can add multiple columns if we so wish.
pbmc <- AddMetaData(
    object = pbmc,
    metadata = subset(pbmc.scmatch, select = c("cell.type")),
    col.name = "scmatch_celltype"
)
```

Now let us look at what it has added to the metadata feature.


```R
head(pbmc@meta.data)
```


<table>
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_RNA</th><th scope=col>nFeature_RNA</th><th scope=col>scmatch_celltype</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACATACAACCAC</th><td>pbmc3k</td><td>2419</td><td> 779</td><td>CD8+ T cell       </td></tr>
	<tr><th scope=row>AAACATTGAGCTAC</th><td>pbmc3k</td><td>4903</td><td>1352</td><td>B cell            </td></tr>
	<tr><th scope=row>AAACATTGATCAGC</th><td>pbmc3k</td><td>3147</td><td>1129</td><td>CD4+ T cell       </td></tr>
	<tr><th scope=row>AAACCGTGCTTCCG</th><td>pbmc3k</td><td>2639</td><td> 960</td><td>monocyte          </td></tr>
	<tr><th scope=row>AAACCGTGTATGCG</th><td>pbmc3k</td><td> 980</td><td> 521</td><td>gamma-delta T cell</td></tr>
	<tr><th scope=row>AAACGCACTGGTAC</th><td>pbmc3k</td><td>2163</td><td> 781</td><td>CD4+ T cell       </td></tr>
</tbody>
</table>



We have annotated cells with the cell typing predictions from scMatch. If a cell was in the scMatch data, but not in the Seurat object, the subset invoked in the AddMetaData() function would remove those scMatch cells. In the reverse scenario where Seurat has cells that the scMatch data does not have, cells would be annotated with a null value (i.e. N/A).

## Quality control visualization

We will now perform some basic visualization of QC.

First, we will look at the overall distribution of:
* number of genes detected per cell
* number of total molecules detected per cell
* the percent of mitochonrial genes detected of total genes detected per cell

We want to look at the mitochondrial genes because low-quality or dying cells tend to exhibit extensive mitchondrial contamination. So this needs to be treated as a variable to consider when normalizing the data (by a linear regression) later.


```R
# Store the mitochondrial percentage into the Seurat object as meta data
# Here, we overwrite pbmc with the returned Seurat object from PercentageFeatureSet
# object   : Seurat object
# pattern  : regular expression pattern to match genes / features to
#            In this case, we look for genes prepended with MT-
#            '^' means only consider if the following characters are at the beginning
# col.name : Name of the new meta data to add to the Seurat object
# Note: In mouse datasets, the mitochondrial genes are listed by "mt-" instead.
#       This is important because the pattern is case sensitive.
pbmc <- PercentageFeatureSet(
    object = pbmc, 
    pattern = "^MT-", 
    col.name = "percent.mt"
)

# Visualize QC metrics as a violin plot
# object : Seurat object
# features : features in Seurat object to display
#            nFeature_RNA = # of genes detected in each cell
#            nCount_RNA   = # of total molecules detected in each cell
#            percent.mt   = percentage of detected genes that are mitochondrial
# ncol : # of columns to split the plots across
VlnPlot(
    object = pbmc, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    ncol = 3)
```


![png](/images/scRNASeq_Seurat/output_20_0.png)


Furthermore, it is good to visualize the relationship between these features. *Note that the Pearson correlation (r value) is shown as the title for each plot.*


```R
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(
    pbmc, 
    feature1 = "nCount_RNA", 
    feature2 = "percent.mt"
)
plot2 <- FeatureScatter(
    pbmc, 
    feature1 = "nCount_RNA", 
    feature2 = "nFeature_RNA"
)
plot1 + plot2
```


![png](/images/scRNASeq_Seurat/output_22_0.png)


## Pre-processing workflow

By default, we are going to use sctransform. This replaces Seurat's previous normalization method.
* The transformed data is stored in the Seurat object under the SCT assay
    * This becomes the default data matrix for Seurat to use for further analysis: PCA, UMAP, etc.
* During the normalization, it is important to remove confounding sources of variation. In the typical human case, this is the mitochondrial mapping percentage
    * Low-quality or dying cells tend to exhibit extensive mitchondrial contamination. So this needs to be treated as a variable to consider when normalizing the data (by a linear regression).
    
Outside of mitochondrial mapping percentage, there are a few other variables that can be included in the normalization. These have more importance in the non-sctransform normalization method shows later.
* The number of unique genes detected in each cell
    * Low-quality cells or empty droplets will often have very few genes
    * Cell doublets or multiplets may exhibit an aberrantly high gene count
* It is possible to include the total RNA sequencing to adjust normalization so that total RNA is considered as a regressor.


```R
# run sctransform
# object          : Seurat object
# vars.to.regress : variables to regress out
# verbose         : printing messages and progress bars
pbmc <- SCTransform(
    object = pbmc, 
    vars.to.regress = "percent.mt",
    verbose = FALSE
)
```

### The alternative normalization method

It is no longer recommended to perform the older normalization method. It scales all the data to the same range, removing the effect of total cell RNA. scTransform preserves the information by the total cell RNA. However in the case that sctransform is not working or you have deliberate reason to want to use this normalization method, here is how perform the following.

** *Please do not actually run the following during this workshop.* **

#### Filtration of data

From the above visuablized QC metrics, we will use these to filter cells.
* We filter cells that have unique feature counts over 2,500 or less than 200
* We filter cells that have >5% mitochondrial counts

These are typical filters to apply for most scRNASeq purposes, but there may be other concerns with your specific library that you may want to further filter on.


```R
# Filter out the cells that fulfill any of the features
pbmc <- subset(
    pbmc, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)
```

#### Data normalization

Next, we normalize the data on the passing cells. By default, the process is to employ a global-scaling normalization method "LogNormalize." This normalized the gene expression measurements for each cell by the total expression, and multiplies this by a scale factor (10,000 by default) so that data is in a reasonably sized value. This is important to avoid floating point errors when dealing with very small numbers. Then the results are log transformed. The normalized data is stored in the Seurat object under 'data': `pbmc[["RNA"]]@data`.

The result of LogNormalize is that all genes are reported as a proportion expression of total expression. This way, samples and cells can be compared directly without concern for variance in library size, etc.

The other options for normalization available are:
* LogNormalize : Gene counts for each cell are divided by the total counts, and multiplied by scale.factor. Then finally log-transformed with *e* (natural number)
* CLR: Applied a centered log ratio transformation - i.e. the genes are divided by the geometric mean of all genes.
* RC: LogNormalize without the final log transformation.


```R
pbmc <- NormalizeData(
    pbmc, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
)
```

#### Identification of highly variable features

We want to find a subset of features that exhibit high cell-to-cell variation in the data set - i.e. highly expressed in some cells, and lowly expressed in others - as focusing on these genes better powers downstream analysis. Theoretically, these highly variable genes are probably more indicative of biological function in single cell.

We will be pulling the top 2,000 most variable genes. There is no particular reason for this number other than it generically captures most variation for most data sets. However, you can alter this as you see fit based on your data. By default, it is suggested to use `vst` as the selection method because it is a fast method to estimate your gene variance.


```R
# Pull the top 2,000 genes
# selection.method : methods to choose the variable features
#                    vst, 
# nfeatures        : number of features to select
pbmc <- FindVariableFeatures(
    pbmc, 
    selection.method = "vst", 
    nfeatures = 2000
)

# Identify the 10 most highly variable genes
# to label in the below plot
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc) + NoLegend()
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + NoLegend()
plot1 + plot2
```

    Warning message:
    “Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
    Please use `as_label()` or `as_name()` instead.
    [90mThis warning is displayed once per session.[39m”
    When using repel, set xnudge and ynudge to 0 for optimal results
    
    Warning message:
    “Transformation introduced infinite values in continuous x-axis”
    Warning message:
    “Transformation introduced infinite values in continuous x-axis”



![png](/images/scRNASeq_Seurat/output_30_1.png)


#### Scaling the data

Finally, we scale the data. This step is to make it easier and faster for dimensional reductions techniques (e.g. PCA) to do their calcuations. The function `ScaleData` works by:
* Centers the gene expression across cells to a mean of 0.
* Scales the expression so that the variance across cells per gene is 1.

This can be a very slow process.


```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

    Centering and scaling data matrix
    


*We are finished with the old normalization methods, and will be continuing from here with the results from sctransform.*

One thing to note in differences from sctransform and the previous normalization method is that sctransform preserves the information from all genes, and not just the most variable genes.

## Dimensional reduction with PCA

We want to first decorrelate the data. This results in:
1. Genes with high correlation being compressed together.
2. The ability to ignore higher dimensions as not contributing much to the data - i.e. noise.
    * This results in reducing the dimensions
We do this decorrelation with principal compoenent analysis.

One advantage of using sctransform over the previous normalization method is that more PCA components can be used in downstream analyses. In the alternative and previously used normalization methodology, you will have to try to minimize the dimensions of PCA used. Whereas, sctransform can let you use a significant number of components - e.g. 40.


```R
pbmc <- RunPCA(pbmc, verbose = FALSE)
```

Once we have the PCA information, we can start to do so preliminary analysis to gain some insights into what is globally driving variation across the cell population. One of the first things we can analyze is what genes are contributing the to each principal component - i.e. what are the most variant genes that are also correlated.


```R
# This function visualizes the 'loadings' that go in the components of PCA 
# - i.e. the gene contributions to that PC.
# dim       : number of dimensions to display
#             e.g. in PCA the components to dis
# reduction : the dimensionality reduction to display
#             while in this example we use pca, 
#             we could use UMAP, TSNE, etc. that Seurat offers
VizDimLoadings(
    pbmc,
    dims = 1:2,
    reduction = 'pca'
)
```


![png](/images/scRNASeq_Seurat/output_36_0.png)


Then we can look at a simple scatter plot of the PCA components. This will give you an idea of clustering, the extent of variation in each component, and the extent of heterogeneity along these components. This last one can give you an idea of how many components to include in further downstream analyses.


```R
DimPlot(pbmc, reduction = 'pca')
```


![png](/images/scRNASeq_Seurat/output_38_0.png)


### Identifying the dimensionality of the data (PCA)

As a means to reducing the 'curse of dimensionality' for clustering, PCA is used to reduce the dimensionality. Thereby, the PCA components are used as input for further downstream analysis: t-SNE, UMAP, clustering, etc. However, it is a good idea to get an idea of the extent of the dimensionality still present after PCA. To this end we have a visualization of the PCA data to identify the extent of variance explained by PCA components.

A simple and less computationally taxing procedure is the elbow plot. Here, the principal components are ranked by percentage of variance explained, and this variance plotted. After the inflection point, there are diminishing returns of a component's contribution to explaining variance in the total dat set.

We want to ensure that the major PCA components are relatively low (under 15) so that downstream analyses can work effectively.


```R
# We run the elbow plot for PCA
ElbowPlot(pbmc, ndims=15, reduction = 'pca')
```


![png](/images/scRNASeq_Seurat/output_40_0.png)


Furhtermore, we explore the PCA heterogeniety further with a heatmap. The cells and genes are ordered according to their PCA scores. In the case of this example heatmap, we are going to limit the cells to the most extreme variance in the component. Heatmaps in general are very computationally taxing to run, so it is not recommended to run on the full gene and cell set.

Below, we create a heatmap of the first principal component only. See how the data tends to cluster.


```R
# dims     : the PCA components to display
# cells    : default is null (i.e. all cells), otherwise top # of cells
# balanced : plot and equal number of genes with both + and - scores
DimHeatmap(
    pbmc, 
    dims = 1, 
    cells = 500, 
    balanced = TRUE
)
```


![png](/images/scRNASeq_Seurat/output_42_0.png)


Next, we may look consecutively at the components 1-8. Here we can see the degradation in the extent of variance contributed by genes in the further components values.


```R
DimHeatmap(
    pbmc, 
    dims = 1:8, 
    cells = 500, 
    balanced = TRUE
)
```


![png](/images/scRNASeq_Seurat/output_44_0.png)


### Further dimensionality reduction with UMAP

Then we want to further reduce dimensions - for visualization - with UMAP (additionally t-SNE can be performed in a similar manner). The outputs will be stored in the Seurat object under reductions. These outputs will be used for visualization, and not as input to downstream analysis with clustering.


```R
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')

# For the input dims for UMAP, you can base it off your data - i.e. from your elbow plot
# dims : number of components from PCA to use
#        in the example below, components 1 - 40
pbmc <- RunUMAP(
    pbmc, 
    dims = 1:40, 
    seed.use = 42
)
# Note, that I have fixed the random seed to 42 so that all UMAPs produced
# this session should look the same.

# t-SNE is invoked similarly to UMAP.
# For now I have commented it out so we can focus on the rest of the downstream
# analysis, as the t-SNE process is fairly slow.
# pbmc <- runTSNE(pbmc, dims = 1:40)

# Plot the UMAP embeddings
DimPlot(pbmc, reduction = 'umap')
```

    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    10:44:02 UMAP embedding parameters a = 0.9922 b = 1.112
    
    10:44:02 Read 2700 rows and found 40 numeric columns
    
    10:44:02 Using Annoy for neighbor search, n_neighbors = 30
    
    10:44:02 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    
    |
    
    10:44:04 Writing NN index file to temp file /tmp/RtmpTm68GR/file130e1f814f4b
    
    10:44:04 Searching Annoy index using 1 thread, search_k = 3000
    
    10:44:05 Annoy recall = 100%
    
    10:44:05 Commencing smooth kNN distance calibration using 1 thread
    
    10:44:06 Initializing from normalized Laplacian + noise
    
    10:44:06 Commencing optimization for 500 epochs, with 115840 positive edges
    
    10:44:16 Optimization finished
    



![png](/images/scRNASeq_Seurat/output_46_1.png)


Additionally, the default values will only give you an approximate visualization of UMAP's projected low dimensional mapping. To get the best compression of the higher dimensional PCA inputs to UMAP, you will need to optimize parameters. However, this is outside the scope of this workshop. Suffice to say, the default configuration for UMAP will give you a rough estimate of the visualization of clusters, but there is the potential that it is not entirely optimal.

## Clustering the cells

Here, we apply a graph-based clustering approach to cluster the cell from the data. This is the KNN-SNN based approach mentioned the other day. This method has two advantages:

1. It is non-parametric in it's cluster partitioning.
2. It automatically determines the number of clusters.

However, one should optimize some of the parameters for best usage. Typically, it has been found that for a approximately 3K cell data set, the resolution parameter is best set for between 0.4 - 1.2. A smaller value, at 3K cells, is tuning the clustering process to identify more clusters. With larger data sets, a larger resolution value may be required for similar results. You can use your UMAP visualization to aid in the choice of optimal resolution paramter value.


```R
# This is a two step process
# First, identify the neighbors for each cell, this is the KNN part
# FindNeighbors by default uses the PCA reduction
pbmc <- FindNeighbors(pbmc, reduction = 'pca', dims = 1:40)
# Second, partition the clusters from this KNN graph
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Then visualize these clusters against the UMAP.
DimPlot(pbmc, reduction = 'umap', label = TRUE) + NoLegend()
```

    Computing nearest neighbor graph
    
    Computing SNN
    


    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 2700
    Number of edges: 130387
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8766
    Number of communities: 10
    Elapsed time: 0 seconds


    Warning message:
    “Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
    Please use `as_label()` or `as_name()` instead.
    [90mThis warning is displayed once per session.[39m”



![png](/images/scRNASeq_Seurat/output_48_3.png)


Finally, we should perform a sanity check for the results of the clustering - if possible. If we know of some canonical biomarkers, we can check to ensure that the clustering - at least for the cells expressing the biomarker - makes sense.


```R
# Visualize canonical marker genes on the sctransform embedding.
#
# object   : Seurat object
# features : list of genes to plot
# pt.size  : size of the individual dots representing each cell
# ncol     : # of columns to display side by side
FeaturePlot(
    object = pbmc, 
    features = c("CD8A", "GZMK", "CCL5", "CCR7"),
    pt.size = 0.2, 
    ncol = 2
)
```


![png](/images/scRNASeq_Seurat/output_50_0.png)


## Fine tuning clusters and visualization

There are a number of considerations to account for when clustering the cells:
1. Number of PCA components to use in computation
2. The parameters used in UMAP
3. The parameters used in the KNN / graph-based clustering

### Choosing the appropriate number of PCA components to use

Choosing the appropriate number of PCA components to use is an important choice for the loadings into UMAP and your graph-based clustering. After a certain point, PCA provides diminishing returns for explaining differences between cells. At some point, PCA is simply returning noise. It is prudent to optimize the choice in the number of components to use.

Let us revisit the elbow plot of the PCA variance.


```R
# We are running the elbow plot for PCA with more dimensions
ElbowPlot(pbmc, ndims=40, reduction = 'pca')
```


![png](/images/scRNASeq_Seurat/output_52_0.png)


We originally used 40 dimensions into UMAP, but what happens if we were to use only the first few PCA dimensions to load into UMAP.


```R
pbmc <- RunUMAP(
    pbmc, 
    dims = 1:8, 
    seed.use = 42
)
DimPlot(pbmc, reduction = 'umap', label = FALSE)
```

![png](/images/scRNASeq_Seurat/output_54_1.png)


We see that when we use less information from reducing the number of PCA dimensions used, that we are losing a lot of fine grained subclustering. The data starts to look more like a plot of the PCA components. Looking back at the elbow plot, we can see that components 30-40 are not contributing much to the variance. The standard Seurat normalization states that the suggested maximum components to use in UMAP is 20. With scTransform, we can go up to 40. We can see from the elbow plot that 20-30 has a lot of small variance, so without scTransform we would lose that information. For these samples, we should probably aim for about 32 components, as there is little benefit to be found in the 30-40 component range and we risk further information being due to noise.


```R
pbmc <- RunUMAP(
    pbmc, 
    dims = 1:32, 
    seed.use = 42
)
DimPlot(pbmc, reduction = 'umap', label = FALSE)
```

![png](/images/scRNASeq_Seurat/output_56_1.png)


This plot now looks strikingly similar to the 40 dimensions plot, which is to be expected. It should be stated that there is no "right" answer to the number of components to use. There are guidelines, like I have shown above, but it will take some testing and educated guesswork. Additionally, the difference between 32 and 33 components to use is minimal. As you can see from the previous plots, there is minimal effect from simply using 40 components. However, every data set is different, and some may have particularly high variance that requires using more components.

### Tuning parameters for UMAP

We have two major parameters for adjusting UMAP. However, the parameter with the largest effect is `min.dist`. We will focus on adjusting `min.dist`. By default, RunUMAP() has a default `min.dist` of 0.3. Let us see the effect of increasing min.dist further. We want to try expanding some of these really tightly packed clusters so we can see more clearly some of the subclusters' densities.


```R
pbmc <- RunUMAP(
    pbmc, 
    dims = 1:32, 
    min.dist = 1.0,
    seed.use = 42,
)
DimPlot(pbmc, reduction = 'umap', label = T)
```


![png](/images/scRNASeq_Seurat/output_58_1.png)


How about if we really want to identify the contours of various subclusters? We will push min.dist lower, to increase the density of idividual clusters.


```R
pbmc <- RunUMAP(
    pbmc, 
    dims = 1:32, 
    min.dist = 0.08,
    seed.use = 42,
)
DimPlot(pbmc, reduction = 'umap', label = T)
```



![png](/images/scRNASeq_Seurat/output_60_1.png)


After loooking at the default, the increased `min.dist`, and the decreased, let us go with a value of 0.7. It will let us view the spread of cells within the clusters better, while not blurring the boundaries between clusters and subclusters too much.


```R
pbmc <- RunUMAP(
    pbmc, 
    dims = 1:32, 
    min.dist = 0.7,
    seed.use = 42,
)
DimPlot(pbmc, reduction = 'umap', label = T)
```



![png](/images/scRNASeq_Seurat/output_62_1.png)


### Tuning parameters for KNN / graph-based clustering

We can see that the default parameters for KNN / graph-based clustering could probably use refinement. Cluster 2 in particular could probably be separated in two clusters. Possibly the same for cluster 3. Additionally, the tag sticking out of cluster 4 may be able to be separated.

For KNN, we will stick with the default `k.param` of 20. It's a good balance of information and computational speed. Increasing k will only result in promoting global structure to the clustering - losing granularity. However, if the next steps do not suffice, you may consider changing `k.param`.


```R
# Running the KNN portion that builds out the graph
# We are running mostly default parameters - with 
# 32 components being used.
pbmc <- FindNeighbors(
    pbmc, 
    reduction = 'pca', 
    dims = 1:32,
    k.param = 20,
    seed.use = 42
)
```

    Warning message:
    “The following arguments are not used: seed.use”
    Warning message:
    “The following arguments are not used: seed.use”
    Computing nearest neighbor graph
    
    Computing SNN
    


The majority of the cluster tuning is best done in the FindClusters section. This adjust the granularity with which Seurat will try clustering from the graph created by `FindNeighbors()`. We know that increasing resolution will increase granularity, so let us try to increase resolution from 0.5 to 1.0.


```R
pbmc <- FindClusters(
    pbmc, 
    resolution = 1.0
)
DimPlot(pbmc, reduction = 'umap', label = T)
```

    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 2700
    Number of edges: 120190
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8022
    Number of communities: 13
    Elapsed time: 0 seconds



![png](/images/scRNASeq_Seurat/output_66_1.png)


This may be too granular. With resolution, a little goes a long way. Let us try to find a middle ground.


```R
pbmc <- FindClusters(
    pbmc, 
    resolution = 0.65
)
DimPlot(pbmc, reduction = 'umap', label = T)
```

    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 2700
    Number of edges: 120190
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8474
    Number of communities: 11
    Elapsed time: 0 seconds



![png](/images/scRNASeq_Seurat/output_68_1.png)


This looks pretty good. Maybe we believe that cluster 2 should be split, however, cluster 0 will take priority. A dense cluster of points like cluster 0 will be prioritized when resolution increases. For now, we will stick with a resolution of 0.65.

### Labeling cell clusters by cell type

Let us use our scMatch cell type annotations to label the cells. In the scMatch annotations, each individual cell is annotated. This may not align perfectly with our clusters, though in aggregate it should.


```R
# We will define our own list of colors to improve clarity
colors = c(
    "#730000", "#ff8800", "#595843", "#00ffaa", 
    "#23698c", "#acb4e6", "#8c697c", 
    "#d2d96c", "#269973", "#005fb3", 
    "#a66cd9", "#330014", "#d9a3a3", "#d9a66c", 
    "#bfff40", "#00ffee", "#262d33", "#603973", 
    "#d93677", "#ff2200", "#ffbf40"
)
DimPlot(
    pbmc,
    reduction = 'umap',
    group.by = 'scmatch_celltype',
    cols = colors
)
```


![png](/images/scRNASeq_Seurat/output_70_0.png)


However, quite a few of these cell type calls are wrong classifications. This will happen at a low rate, so we can identify these and just filter them.


```R
# Here we get the cell types and their respective count of cells
table(pbmc@meta.data$scmatch_celltype)
```


    
         acute lymphoblastic leukemia (B-ALL) cell line 
                                                      2 
                                                 B cell 
                                                    346 
                                            CD4+ T cell 
                                                    375 
                                            CD8+ T cell 
                                                   1202 
                                     gamma-delta T cell 
                                                     74 
                                hematopoietic stem cell 
                                                      5 
                                             macrophage 
                                                      2 
                                               monocyte 
                                                    683 
                 natural killer cell leukemia cell line 
                                                      2 
                            plasmacytoid dendritic cell 
                                                      4 
    T cell chronic lymphocytic leukemia (CLL) cell line 
                                                      5 


Only B cells, CD4+ T cell, CD8+ T cell, gamma-delta T cell, and monocyte were found in any appreciable quantity. We will only look at those, to clean up the UMAP plot and call clusters.


```R
celltypes_of_interest <- c(
    "B cell",
    "CD4+ T cell",
    "CD8+ T cell",
    "gamma-delta T cell",
    "monocyte"
)
pbmc.sub <- subset(
    pbmc,
    subset = scmatch_celltype %in% celltypes_of_interest
)
DimPlot(
    pbmc.sub, 
    reduction = 'umap', 
    group.by = 'scmatch_celltype',
    cols = colors
)
```


![png](/images/scRNASeq_Seurat/output_74_0.png)


Let us try to validate how accurate these cell type calls are with a few biomarkers.


```R
markers = c(
    "CD8A", "CD4", "PTPRC",
    "CD14", "CCL5", "FCGR3A"
)
FeaturePlot(
    pbmc.sub,
    features = markers,
    ncol = 2
)
```


![png](/images/scRNASeq_Seurat/output_76_0.png)


Let us assume that the scMatch cell typing is mostly correct. Let us try to change our cluster names from numbers to cell types. First let us plot the two plots side by side, so that we can reconcile the clusters.


```R
sc <- DimPlot(
    pbmc.sub, 
    reduction = 'umap', 
    group.by = 'scmatch_celltype',
    cols = colors
)
clust <- DimPlot(pbmc, reduction = 'umap', label = T)
clust + sc
```


![png](/images/scRNASeq_Seurat/output_78_0.png)


Looking over the two plots, we can presume there is some differentiation occuring along the CD8+ T cells. Let's assume that cluster 0 (and 7) are possible progenitors, en route to differentiation to CD8+ T cells. We can also presume that clusters 8, 9, and 10 are likely multiplets. We can confirm this with software, such as Scrublet. For now, we will label these clusters as "multiplet". So we will rename the cluster identifying names to something more human readable.


```R
pbmc <- RenameIdents(
    pbmc,
    '0' = "prog 1",
    '1' = 'monocyte',
    '2' = 'B cell',
    '3' = 'CD8+ T cell',
    '4' = 'gamma-delta T cell',
    '5' = 'monocyte',
    '6' = 'CD8+ T cell',
    '7' = 'prog 2',
    '8' = 'multiplet',
    '9' = 'multiplet',
    '10' = 'multiplet'
)
DimPlot(
    pbmc,
    reduction = 'umap',
    label = T
)
```


![png](/images/scRNASeq_Seurat/output_80_0.png)


## Isolating a single cluster

Now suppose we want to identify fine grained analysis of a single cluster. Will re-running PCA -> UMAP result in getting even more granular clusters? Let's find out.


```R
mono <- subset(
    pbmc,
    idents = 'monocyte'
)
mono <- RunPCA(mono, verbose = F)
mono <- RunUMAP(mono, dims = 1:40, seed.use = 42)
DimPlot(mono, reduction = 'umap')
```


![png](/images/scRNASeq_Seurat/output_82_1.png)


We used the whole data set to aid in normalization, so we did not re-run normalization. We then re-ran PCA on the monocytes only, and finally re-ran UMAP. The result is that the cluster looks the same, even the scale for UMAP is roughly the same. This means that there is little advantage to re-running the dimensionality reduction. They are robust against purely local / relative values after normalization.

So, to zoom in on a particular cluster, it is better to simple filter for that cluster, and plot.


```R
DimPlot(
    subset(
        pbmc,
        idents = 'monocyte'
    ),
    reduction = 'umap'
)
```


![png](/images/scRNASeq_Seurat/output_84_0.png)


## Identifying marker genes from a cluster

Once we have our clusters, we may be interested in identifying the marker genes for these clusters. Let us revisit cluster x from our UMAP plot. We know that it expresses CD8A fairly uniquely compared to other cell clusters. So what other genes define the expression profile of this cluster, and where does CD8A rank in that?

We are going to find the markers for cluster 1 by differential expression against the expression of all other cells not in cluster 1 - i.e. comparing cluster 1 against background. As for the statistical test, we are going to specifically invoke MAST. MAST is a differential gene expression analysis program specifically designed for scRNASeq. However, for simply finding marker genes from background, a Wilcoxon Ran Sum test may suffice if you need to run it for computational speed reasons.


```R
# Find all biomarkers for cluster 1
# ident.1  : name of vector for cell names belonging to the first group
#            in the comparison
# ident.2  : by default, all other cells - i.e. background
#            if explicitly named, the other group to compare
# test.use : Denotes which statistical test to use.
cluster1.markers <- FindMarkers(
    pbmc,
    ident.1 = "monocyte",
    test.use = "MAST"
)

# Display the top genes
head(cluster1.markers, n = 5)
```


<table>
<caption>A data.frame: 5 × 5</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>S100A9</th><td>0</td><td>3.044202</td><td>0.992</td><td>0.210</td><td>0</td></tr>
	<tr><th scope=row>LYZ</th><td>0</td><td>2.746862</td><td>1.000</td><td>0.534</td><td>0</td></tr>
	<tr><th scope=row>S100A8</th><td>0</td><td>2.716942</td><td>0.968</td><td>0.104</td><td>0</td></tr>
	<tr><th scope=row>FTL</th><td>0</td><td>1.743696</td><td>1.000</td><td>0.990</td><td>0</td></tr>
	<tr><th scope=row>CST3</th><td>0</td><td>1.721325</td><td>0.998</td><td>0.264</td><td>0</td></tr>
</tbody>
</table>



Again, it is a good idea to perform a sanity check on the results. Here we determine just how the distribution of gene expression for cluster 1's marker genes compares to the rest of the cells. So we create a violin plot to show the distribution.


```R
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
```


![png](/images/scRNASeq_Seurat/output_88_0.png)


Additionally, you can just find all markers for all clusters at once. This can then be used to get a more detailed overview of the top marker genes.


```R
# find markers for every cluster compared to all remaining cells, report only the positive ones
# only.pos : only return the positive / significant genes
# min.pct  : only test genes that are detected in a minimum fraction of cells in any population
#            default = 0.1
pbmc.markers <- FindAllMarkers(
    pbmc, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25,
    test.use = "wilcox"
)
# For speed reasons in this workshop, I am using the Wilcoxon Rank Sum test.
# This is the default DGE algorithm in Seurat.
```

    Calculating cluster 0
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Calculating cluster 8
    
    Calculating cluster 9
    



```R
# Sort all the markers by average log fold change
pbmc.markers.srt <- pbmc.markers[order(pbmc.markers$avg_logFC), ]
# Now apply the function head to all factors of 'cluster'
pbmc.markers.out <- by(pbmc.markers.srt, pbmc.markers.srt['cluster'], head, n=3)
# Let's take a look
pbmc.markers.out
#do.call("rbind", pbmc.markers.out)
```


    cluster: 0
                   p_val avg_logFC pct.1 pct.2    p_val_adj cluster    gene
    C6orf48 2.861384e-30 0.2500030 0.641 0.435 3.597332e-26       0 C6orf48
    EEF1B2  2.491012e-41 0.2513740 0.962 0.874 3.131700e-37       0  EEF1B2
    RPL19   6.479311e-74 0.2521032 1.000 0.999 8.145789e-70       0   RPL19
    ------------------------------------------------------------ 
    cluster: 1
                  p_val avg_logFC pct.1 pct.2    p_val_adj cluster   gene
    BNIP3L 1.171551e-36 0.2505443 0.341 0.117 1.472874e-32       1 BNIP3L
    OAZ2   1.690030e-28 0.2508106 0.317 0.119 2.124705e-24       1   OAZ2
    CSTB   1.394387e-23 0.2519033 0.537 0.309 1.753023e-19       1   CSTB
    ------------------------------------------------------------ 
    cluster: 2
                   p_val avg_logFC pct.1 pct.2    p_val_adj cluster    gene
    TMEM50A 3.273501e-03 0.2562152 0.308 0.245 1.000000e+00       2 TMEM50A
    CD21    1.553084e-17 0.2579718 0.481 0.279 1.952537e-13       2     CD2
    LITAF   4.709034e-08 0.2583921 0.407 0.281 5.920198e-04       2   LITAF
    ------------------------------------------------------------ 
    cluster: 3
                   p_val avg_logFC pct.1 pct.2    p_val_adj cluster    gene
    GNB2L1  3.873942e-19 0.2515910 1.000 0.994 4.870319e-15       3  GNB2L1
    RPS81   2.282538e-30 0.2538568 1.000 0.997 2.869606e-26       3    RPS8
    PLEKHF2 5.231658e-48 0.2604962 0.283 0.055 6.577240e-44       3 PLEKHF2
    ------------------------------------------------------------ 
    cluster: 4
                   p_val avg_logFC pct.1 pct.2    p_val_adj cluster   gene
    DBI     9.026568e-10 0.2507227 0.545 0.321 1.134820e-05       4    DBI
    GLIPR21 4.881219e-13 0.2560070 0.435 0.196 6.136669e-09       4 GLIPR2
    EIF3G   8.391257e-07 0.2567494 0.675 0.529 1.054949e-02       4  EIF3G
    ------------------------------------------------------------ 
    cluster: 5
                  p_val avg_logFC pct.1 pct.2    p_val_adj cluster   gene
    MSN    2.218052e-16 0.2502262 0.552 0.249 2.788535e-12       5    MSN
    RNF149 5.505231e-26 0.2509843 0.429 0.123 6.921176e-22       5 RNF149
    FCGRT1 9.467426e-24 0.2515108 0.571 0.201 1.190245e-19       5  FCGRT
    ------------------------------------------------------------ 
    cluster: 6
                   p_val avg_logFC pct.1 pct.2    p_val_adj cluster   gene
    RGCC1   3.189222e-07 0.2522290 0.424 0.216 4.009489e-03       6   RGCC
    RPL23A1 6.947722e-13 0.2556288 1.000 0.998 8.734676e-09       6 RPL23A
    RPL10   5.031896e-17 0.2569154 1.000 1.000 6.326100e-13       6  RPL10
    ------------------------------------------------------------ 
    cluster: 7
                   p_val avg_logFC pct.1 pct.2    p_val_adj cluster    gene
    NDUFV3  6.278988e-12 0.2514490 0.382 0.070 7.893943e-08       7  NDUFV3
    CCDC88A 7.512920e-24 0.2518290 0.324 0.026 9.445242e-20       7 CCDC88A
    SPATS2L 3.485685e-21 0.2529233 0.294 0.025 4.382203e-17       7 SPATS2L
    ------------------------------------------------------------ 
    cluster: 8
                  p_val avg_logFC pct.1 pct.2 p_val_adj cluster   gene
    SNRPD2 0.0053251920 0.2505309 0.810 0.536         1       8 SNRPD2
    KARS   0.0049869961 0.2529626 0.429 0.193         1       8   KARS
    RBMS1  0.0003004107 0.2545346 0.381 0.125         1       8  RBMS1
    ------------------------------------------------------------ 
    cluster: 9
                   p_val avg_logFC pct.1 pct.2    p_val_adj cluster    gene
    LDLRAP1 2.868751e-03 0.2778558 0.385 0.122 1.0000000000       9 LDLRAP1
    AGPAT1  1.408362e-08 0.2785303 0.308 0.031 0.0001770593       9  AGPAT1
    ARF3    2.823711e-03 0.2925713 0.308 0.086 1.0000000000       9    ARF3


Now that we have the markers for every cluster, we can see how these patterns present themselve in aggregate over all cells / clusters. We will look at everything in the top 10 against all the cells, grouped by their clusters.


```R
# Select from each cluster the top 10 genes of differential expression
# Now again apply the function head to all factors of 'cluster'
pbmc.markers.out <- by(pbmc.markers.srt, pbmc.markers.srt['cluster'], head, n=10)
# by() outputs a list of the dataframes for each cluster
# we are going to merge them with rbind
# then get the unique vector of genes
top10.genes <- do.call("rbind", pbmc.markers.out)$gene

# Plot a heatmap of all of these top 10 genes
# features : list of genes to filter for
DoHeatmap(
    pbmc,
    features = top10.genes
) + NoLegend()
```

    Warning message in DoHeatmap(pbmc, features = top10.genes):
    “The following features were omitted as they were not found in the scale.data slot for the SCT assay: ABCC3, ABHD16A, CYB5R3, STOM, ARF3, TRAT1, EVL, ELF1, SNRPD2, BASP1, LGALS9, FUOM, CCDC88A, NDUFV3, RPL14, RPL36A, RPL10, RPL23A, MTSS1, PCGF5, TCF7L2, SLC2A6, C4orf48, RNF149, MSN, EIF3G, GLIPR2, DBI, MARCH1, CD40, RPS8, LCK, RNASE6, BNIP3L, TOMM7, RPL37, LEF1, PRKCQ-AS1, SOD1”



![png](/images/scRNASeq_Seurat/output_93_1.png)


## Differential gene expression analysis between clusters

Now, we can use this differential analysis to identify what makes genes differentiate two clusters of cells. For example, cluster 1 and cluster 2 share CD8A expression. Considering that it is an important biomarker, we can use differential gene analysis to determine what genes differentiate the two clusters.


```R
# Find all markers of genes from cluster 1 that are differentially expressed
# from cluster 2
# ident.1  : name of vector for cell names belonging to the first group
#            in the comparison
# ident.2  : by default, all other cells - i.e. background
#            if explicitly named, the other group to compare
# logfc.threshold
#          : the minimum log fold change to consider for results of DGE
#            default is 0.25
# test.use : Denotes which statistical test to use.
cluster1.markers <- FindMarkers(
    pbmc,
    ident.1 = 1,
    ident.2 = 2,
    logfc.threshold = 0.25,
    test.use = "MAST"
)

# Display the top genes
head(cluster1.markers, n = 5)
```

    Assuming data assay in position 1, with name et is log-transformed.
    
    
  
    
    Done!
    
    Warning message in melt(llrt):
    “The melt generic in data.table has been passed a list and will attempt to redirect to the relevant reshape2 method; please note that reshape2 is deprecated, and this redirection is now deprecated as well. To continue using melt methods from reshape2 while both libraries are attached, e.g. melt.list, you can prepend the namespace like reshape2::melt(llrt). In the next version, this warning will become an error.”



<table>
<caption>A data.frame: 5 × 5</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>LYZ</th><td> 0.000000e+00</td><td>3.180927</td><td>1.000</td><td>0.470</td><td> 0.000000e+00</td></tr>
	<tr><th scope=row>FTL</th><td> 0.000000e+00</td><td>2.350735</td><td>1.000</td><td>0.992</td><td> 0.000000e+00</td></tr>
	<tr><th scope=row>FTH1</th><td>4.206178e-319</td><td>2.054131</td><td>1.000</td><td>0.992</td><td>5.288008e-315</td></tr>
	<tr><th scope=row>TYROBP</th><td>1.238976e-283</td><td>2.266187</td><td>0.998</td><td>0.165</td><td>1.557641e-279</td></tr>
	<tr><th scope=row>CST3</th><td>2.904036e-270</td><td>2.381567</td><td>0.998</td><td>0.159</td><td>3.650955e-266</td></tr>
</tbody>
</table>



Again, we can refer to the violin plots to confirm the extent of the differential expression of these genes.


```R
# Violin plot
VlnPlot(pbmc, features = c("NKG7", "PF4"))
```


![png](/images/scRNASeq_Seurat/output_97_0.png)



```R
# Scatter plot
FeaturePlot(
    pbmc, 
    features = c("NKG7", "PF4")
)
```


![png](/images/scRNASeq_Seurat/output_98_0.png)


## Integrating multiple data sets for a stimulated vs. control experiment

With this next section we aim to address two goals:
1. Provide an example of analysis comparing two conditions of otherwise similar cell cohorts.
2. How to use Seurat to combine and annotate two separate data sets.

Here, we are going to load another 10X data set of two sets of cells: a control and a stimulated set exposed for 6 hours to IFN-beta. You can get the raw data [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583).

The first steps involve loading the data like normal. Then making a list (in R) of the objects, to which we apply normalization separately.


```R
# Load the 10X data
stim.data <- Read10X(data.dir = "./data/stim")
control.data <- Read10X(data.dir = "./data/control")

# Load the 10X data into a Seurat object
stim <- CreateSeuratObject(counts = stim.data, project = "stim-control")
control <- CreateSeuratObject(counts = control.data, project = "stim-control")
```

    Warning message:
    “Feature names cannot have underscores ('_'), replacing with dashes ('-')”
    Warning message:
    “Feature names cannot have underscores ('_'), replacing with dashes ('-')”



```R
# We are going to reduce this 14K cell data set to 3K

# Identify a random assortment of cells
cells.to.sample <- sample(rownames(stim@meta.data), 3000)
# Subset the 14K dataset for the randomly selected 3K cells
stim.small <- subset(stim, cells = cells.to.sample)

# Repeat for the control
cells.to.sample <- sample(rownames(control@meta.data), 3000)
control <- subset(control, cells = cells.to.sample)
```


```R
# Annotate the objects with the type, by a user-defined group name
control$stim <- "CONTROL"
stim$stim <- "STIM"


# Merge the Seurat objects into one object
immune.combined <- merge(
    control,
    y = c(stim),
    add.cell.ids = c("CONTROL", "STIM"),
    project = "STIM-CONTROL"
)

# Split the object by stim
ifnb.list <- SplitObject(immune.combined, split.by = "stim")

# Then, for each split object, normalize each data set individually
# We are going to use a flat normalization method for 
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
```

Once we have normalized each file separately, we can move ahead with integration of the normalized points.


```R
# Find integration points
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)

# Integrate the objects based on the determined anchor points.
# This will override the previous immune.combined object.
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
```

    Computing 2000 integration features
    
    Scaling features for provided objects
    
    Finding all pairwise anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 11892 anchors
    
    Filtering anchors
    
    	Retained 3353 anchors
    
    Extracting within-dataset neighbors
    
    Merging dataset 1 into 2
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    
    Warning message:
    “Adding a command log without an assay associated with it”


Once we have an integrated data set, we can do our typical downstream analysis.


```R
# This points Seurat to use the integrated data set, and not the
# individually normalized data
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering

# We have to scale the data with a non-sctransform results first
# (not stricly necessary, but it greatly aids PCA calculation)
immune.combined <- ScaleData(immune.combined, verbose = FALSE)

# PCA
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

# UMAP and Clustering
# With this normalization method, it is important to send fewer dimensions of
# PCA to UMAP and clustering. 
# sctransform is much more lenient to using more dimensions (up to 40).
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
```

  

    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 17446
    Number of edges: 624273
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9031
    Number of communities: 13
    Elapsed time: 5 seconds



```R
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2, ncol = 1)
```


![png](/images/scRNASeq_Seurat/output_107_0.png)



```R
# To visualize the two conditions side-by-side, we can use the split.by 
# argument to show each condition colored by cluster.
DimPlot(immune.combined, reduction = "umap", split.by = "stim")
```


![png](/images/scRNASeq_Seurat/output_108_0.png)


### Identifying conserved cell type markers

We may first want to identify what genes are not changed by the stimulus condition. Seurat offers software to do this via the function `FindConservedMarkers`. Please note that this uses a different algorithmic methods from `FindMarkers` and lacks MAST as an option, so we are left with more traditional statistical approaches to dealing with single cell data. However, it is still useful for identifying genes that did not change between conditions.


```R
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(
    immune.combined, 
    ident.1 = 7, 
    grouping.var = "stim", 
    verbose = FALSE
)
head(nk.markers)
```


<table>
<caption>A data.frame: 6 × 12</caption>
<thead>
	<tr><th></th><th scope=col>CONTROL_p_val</th><th scope=col>CONTROL_avg_logFC</th><th scope=col>CONTROL_pct.1</th><th scope=col>CONTROL_pct.2</th><th scope=col>CONTROL_p_val_adj</th><th scope=col>STIM_p_val</th><th scope=col>STIM_avg_logFC</th><th scope=col>STIM_pct.1</th><th scope=col>STIM_pct.2</th><th scope=col>STIM_p_val_adj</th><th scope=col>max_pval</th><th scope=col>minimump_p_val</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>RP11-51J9.5</th><td>2.781586e-53</td><td>1.247707</td><td>0.200</td><td>0.010</td><td>9.912181e-49</td><td>3.662331e-275</td><td>1.349391</td><td>0.222</td><td>0.011</td><td>1.305072e-270</td><td>2.781586e-53</td><td>7.324661e-275</td></tr>
	<tr><th scope=row>IL23A</th><td>2.094611e-41</td><td>1.107370</td><td>0.255</td><td>0.024</td><td>7.464146e-37</td><td>1.555251e-235</td><td>1.126471</td><td>0.189</td><td>0.009</td><td>5.542137e-231</td><td>2.094611e-41</td><td>3.110502e-235</td></tr>
	<tr><th scope=row>SRSF2</th><td>1.894999e-45</td><td>1.716273</td><td>0.673</td><td>0.188</td><td>6.752829e-41</td><td>1.660106e-227</td><td>1.781830</td><td>0.668</td><td>0.168</td><td>5.915789e-223</td><td>1.894999e-45</td><td>3.320213e-227</td></tr>
	<tr><th scope=row>RSRC2</th><td>8.381521e-46</td><td>1.648335</td><td>0.545</td><td>0.113</td><td>2.986755e-41</td><td>3.452786e-220</td><td>1.598169</td><td>0.535</td><td>0.103</td><td>1.230400e-215</td><td>8.381521e-46</td><td>6.905572e-220</td></tr>
	<tr><th scope=row>EIF1</th><td>1.986893e-30</td><td>0.905377</td><td>0.982</td><td>0.920</td><td>7.080295e-26</td><td>7.110044e-188</td><td>1.030753</td><td>0.994</td><td>0.931</td><td>2.533664e-183</td><td>1.986893e-30</td><td>1.422009e-187</td></tr>
	<tr><th scope=row>HSPH1</th><td>6.784131e-78</td><td>2.112116</td><td>0.555</td><td>0.069</td><td>2.417525e-73</td><td>2.420243e-184</td><td>1.788474</td><td>0.467</td><td>0.094</td><td>8.624537e-180</td><td>6.784131e-78</td><td>4.840487e-184</td></tr>
</tbody>
</table>




```R
# Visualize some of these markers and the clusters they belong to
# min.cutoff : minimum threshold for showing features (q9 = quantile 9)
FeaturePlot(
    immune.combined, 
    features = c(
        "CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"
    ), 
    min.cutoff = "q9"
)
```


![png](/images/scRNASeq_Seurat/output_111_0.png)


Now that we "know" what these cell types belonging to each cluster are, we can rename the cluster IDs. This simplifies remembering what cluster is what.


```R
# Rename the clusters
immune.combined <- RenameIdents(
    immune.combined, 
    `0` = "CD14 Mono", 
    `1` = "CD4 Naive T", 
    `2` = "CD4 Memory T", 
    `3` = "CD16 Mono", 
    `4` = "B", 
    `5` = "CD8 T", 
    `6` = "T activated", 
    `7` = "NK", 
    `8` = "DC", 
    `9` = "B Activated", 
    `10` = "Mk", 
    `11` = "pDC", 
    `12` = "Eryth", 
    `13` = "Mono/Nk Doublets"
)

# Let's visualize the clusters with their new labels
DimPlot(immune.combined, label = TRUE)
```

    Warning message:
    “Cannot find identity 13”



![png](/images/scRNASeq_Seurat/output_113_1.png)


Additionally, a dotplot can be used with `split.by` to get a better overall sense of the gene expression across cell types than looking at a scatterplot of expression. The scatter plot can be difficult to ascertain non-binary shifts in expression.


```R
# We want to add the new cell types as factors to the Seurat object
Idents(immune.combined) <- factor(
    Idents(immune.combined), 
    levels = c(
        "pDC", "Eryth", "Mk", "DC", 
        "CD14 Mono", "CD16 Mono", "B Activated", "B", 
        "CD8 T", "NK", "T activated", "CD4 Naive T", 
        "CD4 Memory T"
    )
)

# Tell what markers for DotPlot to plot
markers.to.plot <- c(
    "CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", 
    "GNLY", "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", 
    "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", 
    "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", "HBB", 
    "TSPAN13", "IL3RA", "IGJ"
)

# Plot the dot plot
DotPlot(
    immune.combined, 
    features = rev(markers.to.plot), 
    cols = c("blue", "red"), 
    dot.scale = 8, 
    split.by = "stim"
) + RotatedAxis()
```


![png](/images/scRNASeq_Seurat/output_115_0.png)


### Identify differentially expresed genes across conditions

Now we want to perform a statistical analysis of the differences. Again, we find ourselves using Seurat's `FindMarkers` tool to handle this. We will want to merge the data labels into a single label encompassing both condition and cell type. Then we can perform differential gene expression. This is because Seurat was originally designed around a single dat set, so combining data sets must work around the construction of the code usage.


```R
# Take all the celltypes and combine them with the stimulated condition
# and add under celltype.stim
# Example: B_STIM
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")

# Create a celltype with all the identifiers in immune.combined
immune.combined$celltype <- Idents(immune.combined)

# Relabel the identifiers with celltype.stim identifiers
Idents(immune.combined) <- "celltype.stim"
```


```R
# Run the differential expression response on the selected identifiers
b.interferon.response <- FindMarkers(
    immune.combined, 
    ident.1 = "B_STIM", 
    ident.2 = "B_CONTROL", 
    test.use = "MAST",
    verbose = FALSE
)

# Let us view the top results.
head(b.interferon.response, n = 15)
```

    Assuming data assay in position 1, with name et is log-transformed.
    
    
  
    
     Completed [============================================] 100% with 0 failures
                                                                                  
    
    
    Done!
    
    Combining coefficients and standard errors
    
    Warning message in melt(coefAndCI, as.is = TRUE):
    “The melt generic in data.table has been passed a array and will attempt to redirect to the relevant reshape2 method; please note that reshape2 is deprecated, and this redirection is now deprecated as well. To continue using melt methods from reshape2 while both libraries are attached, e.g. melt.list, you can prepend the namespace like reshape2::melt(coefAndCI). In the next version, this warning will become an error.”
    Calculating log-fold changes
    
    Warning message in melt(lfc):
    “The melt generic in data.table has been passed a list and will attempt to redirect to the relevant reshape2 method; please note that reshape2 is deprecated, and this redirection is now deprecated as well. To continue using melt methods from reshape2 while both libraries are attached, e.g. melt.list, you can prepend the namespace like reshape2::melt(lfc). In the next version, this warning will become an error.”
    Calculating likelihood ratio tests
    
    Refitting on reduced model...
    
    
 
    
     Completed [============================================] 100% with 0 failures
                                                                                  
    
    
    Done!
    
    Warning message in melt(llrt):
    “The melt generic in data.table has been passed a list and will attempt to redirect to the relevant reshape2 method; please note that reshape2 is deprecated, and this redirection is now deprecated as well. To continue using melt methods from reshape2 while both libraries are attached, e.g. melt.list, you can prepend the namespace like reshape2::melt(llrt). In the next version, this warning will become an error.”



<table>
<caption>A data.frame: 15 × 5</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ISG15</th><td>3.661441e-216</td><td>2.9481852</td><td>0.989</td><td>0.274</td><td>1.304755e-211</td></tr>
	<tr><th scope=row>ISG20</th><td>2.088407e-189</td><td>2.2445984</td><td>0.985</td><td>0.439</td><td>7.442039e-185</td></tr>
	<tr><th scope=row>IFI6</th><td>2.376781e-172</td><td>3.0057743</td><td>0.954</td><td>0.071</td><td>8.469659e-168</td></tr>
	<tr><th scope=row>IFIT3</th><td>7.716120e-129</td><td>2.9012107</td><td>0.877</td><td>0.057</td><td>2.749639e-124</td></tr>
	<tr><th scope=row>IFIT1</th><td>9.268491e-124</td><td>2.7636979</td><td>0.842</td><td>0.028</td><td>3.302827e-119</td></tr>
	<tr><th scope=row>LY6E</th><td> 4.190222e-99</td><td>1.6798024</td><td>0.930</td><td>0.325</td><td> 1.493185e-94</td></tr>
	<tr><th scope=row>MX1</th><td> 1.016617e-92</td><td>2.2014597</td><td>0.811</td><td>0.099</td><td> 3.622714e-88</td></tr>
	<tr><th scope=row>PLSCR1</th><td> 2.511780e-77</td><td>2.1052458</td><td>0.738</td><td>0.075</td><td> 8.950727e-73</td></tr>
	<tr><th scope=row>OAS1</th><td> 1.147920e-62</td><td>2.0385737</td><td>0.651</td><td>0.061</td><td> 4.090612e-58</td></tr>
	<tr><th scope=row>IFIT2</th><td> 2.219838e-62</td><td>2.0480978</td><td>0.610</td><td>0.033</td><td> 7.910391e-58</td></tr>
	<tr><th scope=row>B2M</th><td> 5.625221e-60</td><td>0.3279358</td><td>1.000</td><td>1.000</td><td> 2.004547e-55</td></tr>
	<tr><th scope=row>IFITM3</th><td> 1.994747e-59</td><td>2.4469463</td><td>0.617</td><td>0.061</td><td> 7.108279e-55</td></tr>
	<tr><th scope=row>TNFSF10</th><td> 6.027435e-58</td><td>1.9220032</td><td>0.684</td><td>0.118</td><td> 2.147877e-53</td></tr>
	<tr><th scope=row>IRF7</th><td> 2.860615e-54</td><td>1.5765888</td><td>0.712</td><td>0.151</td><td> 1.019380e-49</td></tr>
	<tr><th scope=row>MT2A</th><td> 2.561631e-52</td><td>1.8159387</td><td>0.604</td><td>0.066</td><td> 9.128371e-48</td></tr>
</tbody>
</table>



We can now visually explore the underlying distribution of expression for these statistically significant genes with scatterplots or violin plots of their expression across conditions.


```R
# Scatter plot
FeaturePlot(
    immune.combined, 
    features = c("GNLY", "IFI6"), 
    split.by = "stim", 
    max.cutoff = 3, 
    cols = c("grey", "red")
)
```


![png](/images/scRNASeq_Seurat/output_121_0.png)



```R
# Violin plot
plots <- VlnPlot(
    immune.combined, 
    features = c("LYZ", "ISG15", "CXCL10"), 
    split.by = "stim", 
    group.by = "celltype", 
    pt.size = 0, 
    combine = FALSE
)
CombinePlots(plots = plots, ncol = 1)
```

    Warning message:
    “CombinePlots is being deprecated. Plots should now be combined using the patchwork system.”



![png](/images/scRNASeq_Seurat/output_122_1.png)


## Integrating multiple data sets with sctransform

Above, we used normalization assuming a fixed scale. sctransform is a more elegant and pertinent solution to normalization, as shown above in the single data set examples. However, integration with sctransform is a fairly new feature added. In our experience, it can be extremely computationally expensive, depending on the data set. In our case with the stimulated vs control experiment, it falls under the latter. Because of the computational limits of the workshop, we used the other methodology.

However, we are providing you with the code, if you want to try it on other data.


```R
# Load the 10X data
stim.data <- Read10X(data.dir = "./data/stim")
control.data <- Read10X(data.dir = "./data/control")

# Load the 10X data into a Seurat object
stim <- CreateSeuratObject(counts = stim.data, project = "stim-control")
control <- CreateSeuratObject(counts = control.data, project = "stim-control")
```


```R
# We are going to reduce this 14K cell data set to 3K

# Identify a random assortment of cells
cells.to.sample <- sample(rownames(stim@meta.data), 3000)
# Subset the 14K dataset for the randomly selected 3K cells
stim.small <- subset(stim, cells = cells.to.sample)

# Repeat for the control
cells.to.sample <- sample(rownames(control@meta.data), 3000)
control <- subset(control, cells = cells.to.sample)
```

In the above code, we loaded the data into separate Seurat object. We labeled them according to the same project. Next, we are going to combine the data sets into one object, and label each data set accordingly. Note that the new object will override the project names from before. Additionally, I used the format `y = c(stim)`, because this can be generalized to `y = c(stim, stim2, stim3)` or however many sets you want to merge. You would then add more IDs to the object: `add.cell.ids = c("CONTROL", "STIM", "STIM2", "STIM3")`.


```R
# Annotate the objects with the type, by a user-defined group name
control$cond <- "CONTROL"
stim$cond <- "STIM"


# Merge the Seurat objects into one object
combined <- merge(
    control,
    y = c(stim),
    add.cell.ids = c("CONTROL", "STIM"),
    project = "STIM-CONTROL"
)
```

Now that they are merged together in one object, we can normalize the data and start visualizations. Note that now the visualizations are split by or color labeled by the cell IDs that we annotated.


```R
# Flag the mitochondrial genes
combined <- PercentageFeatureSet(combined, pattern = "^Mt-", col.name = "percent.mt")

# Visualize QC metrics as a violin plot
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


```R
# Visualize the feature-feature relationship with scatter plots
plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
```

Now, we plan to make a new `combined` object, as to perform a comparison across the data sets, we have to integrate them together - i.e. we need to find anchor points of comparison between the data sets. To do this, the process is roughly as follows:

1. sctransform each data set individually
2. Identify the anchor features (genes) to compare dataset to dataset
3. Re-normalize the data by these anchor points

This step is of particular importance when working with different technologies being merged together - i.e. comparing 10X prepared samples with inDrop prepared samples. However, even within the same technologies it is good practice, even if not strictly as necessary.


```R
# Remake the new combined object so 'state' is clean
combined <- merge(
    control,
    y = c(stim),
    add.cell.ids = c("CONTROL", "STIM"),
    project = "STIM-CONTROL"
)

# Split the object into a list of Seurat objects
combined.list <- SplitObject(combined, split.by = 'cond')
```


```R
# Now we sctransform each data set individually
for (i in 1 : length(combined.list)) {
    combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE)
}
```


```R
# Prepare the features for anchor identification
combined.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 2000)
```


```R
# necessary to allow R to handle the size of the data
options(future.globals.maxSize=891289600) 

combined.list <- PrepSCTIntegration(
    object.list = combined.list,
    anchor.features = combined.features,
    verbose = FALSE
)
```


```R
# Identify the anchors and integrate the datasets.
# Here we must re-normalize with sctransform via 'SCT' parameter
combined.anchors <- FindIntegrationAnchors(
    object.list = combined.list,
    normalization.method = "SCT",
    anchor.features = combined.features,
    verbose = FALSE
)
combined.integrated <- IntegrateData(
    anchorset = combined.anchors,
    normalization.method = "SCT", 
    verbose = FALSE
)
```

Once we have normalized the data, we can perform our downstream analyses as before: PCA, UMAP, and clustering.


```R
# Dimensionality reduction with PCA
combined.integrated <- RunPCA(combined.integrated, verbose = FALSE)

# Run UMAP
combined.integrated <- RunUMAP(
    combined.integrated, 
    dims = 1:30, 
    verbose = FALSE
)

# Cluster the cells
combined.integrated <- FindNeighbors(
    combined.integrated, 
    dims = 1:30, 
    verbose = FALSE
)
combined.integrated <- FindClusters(
    combined.integrated, 
    verbose = FALSE
)

# Plot the UMAP of the combined data
DimPlot(
    combined.integrated, 
    reduction="umap", 
    label = TRUE
) + NoLegend()
```

### Visualization with annotations

However, these are combined data sets, and with the way we annotated the Seurat object, we can split these splits by their annotations. In the below example, we are going to display these plots by annotation in two ways:
1. an overlay of the groups
2. splitting the groups into separate plots
We are going to do this with the **cond** group annotation that we set for each data set prior to merging.


```R
# First, plot the UMAP but with the conditions split into an overlay
DimPlot(combined.integrated, reduction = "umap", group.by = "cond")
```


```R
# Then, we plot the UMAP but with the conditions split into separate side
# by side plots
DimPlot(combined.integrated, reduction = "umap", split.by = "cond")
```

Additionally, we can chain both `split.by` and `group.by` for cases where we have multiple annotations.

We can then move forward with identifying canonical biomarkers and how they may differ between conditions.


```R
FeaturePlot(combined.integrated, features = c("CD8A"), split.by = "cond")
```

### Differential analysis for conserved markers between conditions

We may first want to identify what genes are not changed by the stimulus condition. Seurat offers software to do this via the function `FindConservedMarkers`. Please note that this uses a different algorithmic methods from `FindMarkers` and lacks MAST as an option, so we are left with more traditional statistical approaches to dealing with single cell data. However, it is still useful for identifying genes that did not change between conditions.


```R
DefaultAssay(combined.integrated) <- "RNA"
nk.markers <- FindConservedMarkers(
    combined.integrated, 
    ident.1 = 7, 
    grouping.var = "cond", 
    verbose = FALSE
)
head(nk.markers)
```


```R
# We can then visualize the expression differences between the two conditions for some genes
FeaturePlot(combined.integrated, features = c("CD8A", "GNLY", "CD79A", "FCGR3A", 
    "CCL2", "PPBP"), min.cutoff = "q9", split.by = "cond")
```

Because this is a previously published data set, we can know some things about the data. This provides an example of how one can annotate your clusters with more human readable notations than "cluster 1", "cluster 2", etc.


```R
combined.integrated <- RenameIdents(
    combined.integrated, 
    `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", 
    `7` = "NK", `8` = "DC", `9` = "B Activated", `10` = "Mk", 
    `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets"
)

DimPlot(combined.integrated, label = TRUE)
```

A better overview of the expression across clusters and conditions can be done with a dot plot. This provides us a simple graphic showcasing the expressions ordered across the multiple cell types.


```R
Idents(combined.integrated) <- factor(
    Idents(combined.integrated),
    levels = c(
        "Mono/Mk Doublets", "pDC", "Eryth", "Mk", "DC", 
        "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", 
        "NK", "T activated", "CD4 Naive T", "CD4 Memory T"
    )
)
markers.to.plot <- c(
    "CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY",
    "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", 
    "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", 
    "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ"
)
DotPlot(
    combined.integrated, 
    features = rev(markers.to.plot), 
    cols = c("blue", "red"), 
    dot.scale = 8, 
    split.by = "cond"
) + RotatedAxis()
```

### Identify differentially expressed genes between conditions

Now we get to the typical central question for a single cell conditional experiment. How do individual cell types differ across some condition. First thing we will do is visually inspect how the gene expression compares beween the two conditions with a scatterplot of gene expression between the two conditions.


```R
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# First subset just the CD4 T cell expressions
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "cond"
# Get the average expression for genes across the T cells for all genes
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
avg.t.cells$gene <- rownames(avg.t.cells)

# Then subset just the CD14 Mono cell expressions
cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "cond"
# Again, gets the average expression for all genes in the cell type
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

# Here we are going to select a few genes that we already know are
# differentially expressed to highlight them
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

# Create a scatter plot of the CD4 T cells, CONTROL expression on the 
# x-axis and STIM on the y-axis
p1 <- ggplot(avg.t.cells, aes(CONTROL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
# Create a scatter plot of the CD14 Mono cells, CONTROL expression on the 
# x-axis and STIM on the y-axis
p2 <- ggplot(avg.cd14.mono, aes(CONTROL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
# Combine the two plots into a side by side grid and display
plot_grid(p1, p2)
```

Let us find the markers between B cells stimulated and B cells under control.


```R
# Create new labels of "CD14 Mono_CONTROL", "CD14 Mono_STIM", 
# "CD4 Naive T_CONTROL", etc.
combined.integrated$celltype.stim <- paste(
    Idents(combined.integrated), 
    combined.integrated$cond, 
    sep = "_"
)
# Then we merge cell typings plus their condition into a 
# new ident under celltype.stim
combined.integrated$celltype <- Idents(combined.integrated)
Idents(combined.integrated) <- "celltype.stim"

# Now we find some markers 
interferon.response <- FindMarkers(
    combined.integrated, 
    ident.1 = "B_STIM", 
    ident.2 = "B_CONTROL",
    test.use = "MAST",
    verbose = FALSE
)
head(interferon.response, n = 15)
```

Now let's compare one of the differentially expressed genes response visually compared to a conserved gene we identified previously. We will start with an overlay of the gene expression over the scatter plot of clusters.


```R
FeaturePlot(
    combined.itegrated, 
    features = c("CD3D", "IFI6"), 
    split.by = "stim", 
    max.cutoff = 3, 
    cols = c("grey", "red")
)
```

We can also visualize the distribution of expression more clearly with a violin plot.


```R
plots <- VlnPlot(
    combined.integrated, 
    features = c("CD3D", "IFI6"), 
    split.by = "cond", 
    group.by = "celltype", 
    pt.size = 0, 
    combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
```
