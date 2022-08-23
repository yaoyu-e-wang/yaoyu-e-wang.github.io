---
title: "Basic RNASeq DGE Analysis"
date: 2022-05-01
draft: false
tags: ["R", "RNASeq", "DGE", "tutorial"]
categories: ["Tutorial"]
author: "Yaoyu E. Wang"

autoCollapseToc: true

---

Begin RNA-Seq DGE Analysis
==========================

We first load the RNA-Seq count matrix from ‘RNASeqData.Rdata’. This is
a simple 3 vs 3 RNA-Seq experiment with Control vs Experiment in count
matrix by gene.

    # One can set up the path
    # Set working directory to PROJECT_DIR
    setwd("~/Dropbox (Harvard University)/qBRC-Share/Presentations/QBRC-Workshops/Intro to RNASeq Analysis/2022-0424-IID/RNASeq_Exploratory")
    # setwd("~/modules/rnaseq_dge/")
    load('RNASeq_exploratory.RData') 


    ls()                       # list all current variables

    ## [1] "count_data"          "exploratory_results" "groups"

    #head(count_data)           # The head(count_data) provides the first 6 row of the 'count_data'
    #head(groups)           # The head(count_data) provides the first 6 row of the 'count_data'

Principle component analysis
----------------------------

Once the matrix is loaded, we can perform principle component analysis
to visualize sample distribution.  
We first group samples into two vectors/variables: CONTROL and TREATED.

Calculate PCA use ‘prcomp’ command. Since prcomp compute PCA by the
rows, we will need to tranpose the count matrix such that the samples
are represented by rows and genes by the columns.

    CONTROL=groups$control$sample
    TREATED=groups$treated$sample

    pca=prcomp(t(count_data))        # compute PCA on transposed count_data
    plot(pca$x[,1], pca$x[,2], xlab='pc1', ylab='pc2')
    points(pca$x[CONTROL,1], pca$x[CONTROL,2],
           col="red", pch=18)
    points(pca$x[TREATED,1], pca$x[TREATED,2],
           col="blue", pch=18)

![PCA](/images/RNASeq_DGE_files/figure-markdown_strict/unnamed-chunk-1-1.png)

Hierarchical Clustering
-----------------------

We first calculate similarity matrix using Euclidean distance. The
matrix is computed by ‘dist’ function. ‘hclust’ function computes
hierarhical cluster tree and save it into an variable object ‘hctree’
that can be ploted by using generic ‘plot’ function.

    # compute similarity matrix on the transposed count_data matrix using Euclidean distance
    dist_matrix=dist(t(count_data), method="euclidean")
    # compute hierarchical clustering tree and save it into variable called 'hctree'
    hctree=hclust(dist_matrix)
    # plot hctree
    plot(hctree)

![](/images/RNASeq_DGE_files/figure-markdown_strict/unnamed-chunk-2-1.png)

Differential Gene Expression
----------------------------

We will use DESeq2 to perform differential gene expression. Since DESeq2
is a specialized bioconductor package/library, we will need to install
it before loading in the package to use.

The following command check if DESeq2 and Bioconductor Manager packages
have been installed, and install the packages if they have not been
installed.

    # Check for Bioconductor, install it if it is not yet installed
      if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')

    # Check for DESeq2, install it if it is not yet installed
      if(!requireNamespace('DESeq2', quietly = TRUE))
        BiocManager::install('DESeq2', update=FALSE)

Once the packages are installed, we can call ‘library’ function to load
‘DESeq2’ library to use its functions.

We can then use the functions within DESeq2 to perform DGE analysis

First, we define which sampels are control group and which samples are
treated group The samples are ordered as YEW1\_Control, YEW2\_Control,
YEW3\_Control, YEW1\_Treated, YEW2\_Treated, YEW3\_Treated, so we can
use the following to define conditions of each sample

    condition <- factor(c("control", "control", "control",
                          "treated", "treated", "treated"))
    # alternatively, since the count_data is ordered in the same way as groups variable
    # we can declare condition from groups:
    # 
    # condition <- factor(c(as.character(groups$control$group), as.character(groups$treated$group)))
    #

We then use the ‘condition’ to perform contrast on the count\_data.
Noted that count\_data is still raw count matrix, D ESeq2 performs
normalization and DGE together in one function.

If you want to know how to use a specific function in R, just type
‘?function’ in Console. For example type:

> ?DESeqDataSetFromMatrix

    # Returns DGE gene list results in data frame format
    dge_results=results(dds, contrast=c("condition", "treated", "control"))

    # make sure there is no missing values
    dge_results=na.omit(dge_results)  

    # filter for results with padj<0.001 and absolute log2FoldChange>2
    filtered_results=subset(dge_results, padj<0.001 & abs(log2FoldChange)>2)

    # Obtain normalized count values
    norm_count <- counts(dds, normalized = T)

    # output filtered_results
    filtered_results

    ## log2 fold change (MLE): condition treated vs control 
    ## Wald test p-value: condition treated vs control 
    ## DataFrame with 30 rows and 6 columns
    ##                       baseMean    log2FoldChange             lfcSE
    ##                      <numeric>         <numeric>         <numeric>
    ## RP1-137D17.1  29.8839883584897  2.18080051514365 0.525441950834757
    ## AP000344.4      73.14903147842  2.12391047870749 0.333193245882354
    ## CTB-134H23.2  16.4631893142531 -4.06989298788851 0.902046395861882
    ## C12orf40      23.5935805512499 -2.69347617414478 0.622311943521575
    ## RP11-395P13.1 25.2555673433736  2.25522540567288  0.55251100155043
    ## ...                        ...               ...               ...
    ## RP11-413E6.7  62.9678012294187 -3.57669120749119 0.409461162253874
    ## CTD-2144E22.8 79.2703453044705 -2.48457438278124 0.339133969829025
    ## RP11-747D18.1  1062.4407047496  2.01064482797773  0.12042205283902
    ## RP11-864I4.4  25.0810724471142  2.31146493910142  0.55115609828326
    ## MSMP          44.7309465700437  2.16345289754786 0.400073627551603
    ##                            stat               pvalue                 padj
    ##                       <numeric>            <numeric>            <numeric>
    ## RP1-137D17.1   4.15041188028299 3.31877521851314e-05  0.00026025079086948
    ## AP000344.4      6.3744103608193 1.83668201976646e-10 3.60334193573427e-09
    ## CTB-134H23.2  -4.51184440906705 6.42663252483387e-06   5.821921826631e-05
    ## C12orf40      -4.32817689293055 1.50348695874195e-05  0.00012712159813363
    ## RP11-395P13.1  4.08177466031332 4.46931192753386e-05 0.000338444791804363
    ## ...                         ...                  ...                  ...
    ## RP11-413E6.7  -8.73511711783197 2.43403807314299e-18 1.20260893386853e-16
    ## CTD-2144E22.8 -7.32623270984574 2.36712611316642e-13 6.83284951698415e-12
    ## RP11-747D18.1  16.6966496632104 1.38633818987756e-62  7.2031449110832e-60
    ## RP11-864I4.4   4.19384807008606 2.74261710788168e-05 0.000220325854147882
    ## MSMP            5.4076368662136 6.38617383047045e-08 8.34139160148505e-07

We can compute PCA on normalized matrix again

    CONTROL=groups$control$sample
    TREATED=groups$treated$sample

    pca=prcomp(t(norm_count))        # compute PCA on transposed normalized_data
    plot(pca$x[,1], pca$x[,2], xlab='pc1', ylab='pc2')
    points(pca$x[CONTROL,1], pca$x[CONTROL,2],
           col="red", pch=18)
    points(pca$x[TREATED,1], pca$x[TREATED,2],
           col="blue", pch=18)

![](/images/RNASeq_DGE_files/figure-markdown_strict/norm_pca-1.png)

    # compute similarity matrix on the transposed count_data matrix using Euclidean distance
    dist_matrix=dist(t(norm_count), method="euclidean")
    # compute hierarchical clustering tree and save it into variable called 'hctree'
    hctree=hclust(dist_matrix)
    # plot hctree
    plot(hctree)

![](/images/RNASeq_DGE_files/figure-markdown_strict/norm_pca-2.png)

Heatmap
-------

Generate heatmap using filtered dge results ‘filtered\_results’. The
heatmap function is core function, but we want to have better coloring,
so we install and load gplots.

    # Check and install "gplots" if not installed
    if(!requireNamespace('gplots', quietly = TRUE))
      install.packages('gplots')
    library(gplots)

    ## Warning: package 'gplots' was built under R version 3.5.2

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     space

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     space

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

    DGE_matrix=norm_count[rownames(filtered_results),]
    heatmap(DGE_matrix,col=colorpanel(50, 'blue', 'white', 'red'),
            Rowv=TRUE, Colv=TRUE)

![](/images/RNASeq_DGE_files/figure-markdown_strict/heatmap-1.png)

    # call the commmand again to write it into a file.
    png("heatmap.png")
    heatmap(DGE_matrix,col=colorpanel(50, 'blue', 'white', 'red'),
            Rowv=TRUE, Colv=TRUE)
    dev.off()

    ## quartz_off_screen 
    ##                 2

Generate Volcano Plot
---------------------

Again, we first check for the package ‘EnhancedVolcano’, and install it
if it is not already installed.

Load the library and run EnhancedVolcano function to generate volcano
plot.

    library(EnhancedVolcano)

    ## Warning: package 'EnhancedVolcano' was built under R version 3.5.2

    ## Loading required package: ggplot2

    ## Loading required package: ggrepel

    ## Warning: package 'ggrepel' was built under R version 3.5.2

    volcano_plot<-EnhancedVolcano(dge_results,
        lab = rownames(dge_results),
        x = 'log2FoldChange',
        y = 'padj',
        xlim = c(-5, 8),
        ylim = c(0, 20),
        #pCutoff = 10e-16,
        FCcutoff = 1.5
    )

    volcano_plot

    ## Warning: Removed 235 rows containing missing values (geom_point).

    ## Warning: Removed 11 rows containing missing values (geom_text).

![](/images/RNASeq_DGE_files/figure-markdown_strict/volcano-1.png)

Gene Set Erichment Analysis
---------------------------

There are many different packages retreiving MSigDB data to perform GSEA
analysis. Most of the packages perform these tasks in very similar ways
with different application programming interface (API). We will use the
following packages:

-   migdbr
    (<a href="https://github.com/igordot/msigdbr" class="uri">https://github.com/igordot/msigdbr</a>)
-   fgsea
    (<a href="https://bioconductor.org/packages/release/bioc/html/fgsea.html" class="uri">https://bioconductor.org/packages/release/bioc/html/fgsea.html</a>)

First, we install **misgdbr**

    # install misgdbr from CRAN
    if(!require(msigdbr)) install.packages("msigdbr")

    ## Loading required package: msigdbr

    library(msigdbr)

The **misgdbr** retrieve data from MSigDB database hosted at Broad
Institute
(<a href="https://www.gsea-msigdb.org/gsea/msigdb/index.jsp" class="uri">https://www.gsea-msigdb.org/gsea/msigdb/index.jsp</a>).
The function **misgdbr** runs as the following:

\#\#\#\#Usage msigdbr(species = “Homo sapiens”, category = NULL,
subcategory = NULL)

\#\#\#\#Arguments |Parameters|Description| |:—:|—| |**species**|species
name, such as Homo sapiens, Mus musculus, etc.|
|**category**|collection, such as H, C1, C2, C3, C4, C5, C6, C7.|
|**subcategory**|sub-collection, such as CGP, MIR, BP, etc.|

    m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
    # reformat the data frame into list for fgsea input
    m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name) 


    # install fgsea
    if(!require(fgsea)) BiocManager::install("fgsea")

    ## Loading required package: fgsea

    ## Loading required package: Rcpp

    ## Warning: package 'Rcpp' was built under R version 3.5.2

    library(fgsea)

    # reformat
    gene_stat=dge_results$stat
    names(gene_stat)=rownames(dge_results)
    # sort genes
    gene_stat=gene_stat[order(gene_stat)]

    fgsea_results=fgsea(pathways=m_list, gene_stat, nperm=50, minSize=20)
    fgsea_results=fgsea_results[order(fgsea_results$pval),]
    head(fgsea_results)

    ##                                        pathway       pval padj         ES
    ## 1:    KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION 0.03703704    1 -0.4975169
    ## 2:                          KEGG_AXON_GUIDANCE 0.03846154    1 -0.4746919
    ## 3:  KEGG_PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM 0.07692308    1 -0.4796340
    ## 4: KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION 0.11428571    1  0.3375730
    ## 5:                             KEGG_CELL_CYCLE 0.12903226    1  0.4658802
    ## 6:                 KEGG_SMALL_CELL_LUNG_CANCER 0.15384615    1 -0.4090299
    ##          NES nMoreExtreme size                                 leadingEdge
    ## 1: -1.487950            0   24      PDIA3,KLRC4,HLA-DOA,RFX5,KIR3DL1,IFNA4
    ## 2: -1.430389            0   25        CDK5,EPHA5,SEMA3E,PPP3R2,EFNA4,UNC5D
    ## 3: -1.503249            1   27 INPP5E,CALML6,CDS1,INPP5A,PLCE1,PIK3C2G,...
    ## 4:  1.271000            3   84  IL20,CCL28,TNFRSF9,ACVR2B,CXCL6,TNFSF4,...
    ## 5:  1.349746            3   37    CCNB3,SMAD4,RAD21,ZBTB17,MCM2,CDKN1C,...
    ## 6: -1.242841            3   26      MAX,TRAF2,COL4A1,LAMA5,PIAS1,ITGB1,...

Now write out the GSEA results

    library(writexl)    # load  library to write excel file

    ## Warning: package 'writexl' was built under R version 3.5.2

    # Current directory is Data directory, change to output
    if(!dir.exists('outputs')){
      dir.create('outputs')
    }

    # using ggplots's ggsave function to save the plot
    print("Save the volcano plot in ggplot format")

    ## [1] "Save the volcano plot in ggplot format"

    ggsave("outputs/volcano.png", plot=volcano_plot)

    ## Saving 7 x 5 in image

    ## Warning: Removed 235 rows containing missing values (geom_point).

    ## Warning: Removed 11 rows containing missing values (geom_text).

    print('Write the final DGE data to Excel file.')

    ## [1] "Write the final DGE data to Excel file."

    # start a new excel file and write final results 
    write_xlsx(as.data.frame(dge_results), "outputs/dge_results.xlsx", col_names=TRUE)


    print('Write the final GSEA data to Excel file.')

    ## [1] "Write the final GSEA data to Excel file."

    # start a new excel file and write final results 
    write_xlsx(fgsea_results, "outputs/gsea_results.xlsx", col_names=TRUE)
