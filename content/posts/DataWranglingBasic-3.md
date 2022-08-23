---
title: "Data Wrangling with R Basics - 3"
date: 2020-04-26
draft: false
tags: ["R", "Data Wrangling", "tutorial"]
categories: ["Tutorial"]
author: "Yaoyu E. Wang"

autoCollapseToc: true
---
Data Wrangling with R Basics - 3 
-------------------------

Let's load the dataset


    data=list(
        dge_results=read.table('scripts and data/data/differential_results.csv', sep=',', header=T),
        expr=read.table('scripts and data/data/gene_expression.tsv', header=T, row.names=1),  # first column being the rowname
        annot=read.delim('scripts and data/data/gene_annotations.tsv')
    )

find genes with greater than 0 mean expression between a given chr coordinate
=============================================================================

Let’s find out mean first
-------------------------

You can use **rowMeans**, **rowMax**, **rowMin** to determine means for
a data matrix

    expr_mean=rowMeans(data$expr)
    gene_expr=cbind(data$expr, mean=expr_mean)

We can use the function **subset** to find the given with a given
coordinate

    # select all genes on chromosome 20
    head(subset(data$annot, chrom=='chr20'))

    ##       chrom  start    end strand       name
    ## 53155 chr20  87250  97094      +    DEFB125
    ## 53156 chr20 142369 145751      +    DEFB126
    ## 53157 chr20 157470 159163      +    DEFB127
    ## 53158 chr20 173229 173296      - AL360078.1
    ## 53159 chr20 187853 189681      -    DEFB128
    ## 53160 chr20 227258 229886      +    DEFB129

    # select all genes on chromosome 20 and on positive strand
    head(subset(data$annot, chrom == 'chr20' & strand=="+"))

    ##       chrom  start    end strand          name
    ## 53155 chr20  87250  97094      +       DEFB125
    ## 53156 chr20 142369 145751      +       DEFB126
    ## 53157 chr20 157470 159163      +       DEFB127
    ## 53160 chr20 227258 229886      +       DEFB129
    ## 53161 chr20 257736 261096      +       DEFB132
    ## 53162 chr20 267186 268857      + RP5-1103G7.10

    # select all genes on chromosome 20 and on positive strand and between 10mb to 20mb
    head(subset(data$annot, chrom == 'chr20' & strand=="+" & start>10000000 & end<20000000))

    ##       chrom    start      end strand         name
    ## 53383 chr20 10025917 10026168      +  Metazoa_SRP
    ## 53384 chr20 10173520 10196990      + RP11-416N4.4
    ## 53385 chr20 10180235 10185775      + RP11-416N4.1
    ## 53386 chr20 10218830 10307418      +       SNAP25
    ## 53387 chr20 10334419 10334698      +    HIGD1AP15
    ## 53389 chr20 10385779 10387806      +      SDAD1P2

    # select all genes on chromosome 20 and on positive strand and between 10mb to 20mb
    head(subset(data$annot, chrom == 'chr20' & strand=="+" & start>10000000 & end<20000000, select=c("chrom", "name")))

    ##       chrom         name
    ## 53383 chr20  Metazoa_SRP
    ## 53384 chr20 RP11-416N4.4
    ## 53385 chr20 RP11-416N4.1
    ## 53386 chr20       SNAP25
    ## 53387 chr20    HIGD1AP15
    ## 53389 chr20      SDAD1P2

    # store into a variable
    chr20_subset=subset(data$annot, chrom == 'chr20' & strand=="+" & start>10000000 & end<20000000, select=c("chrom", "name"))

    # merge with gene expression data with means
    chr20_expr=merge(chr20_subset, gene_expr, by.x="name", by.y=0)  #by.y=0 means using row names

    head(chr20_expr)

    ##         name chrom SW1_Control SW2_Control SW3_Control SW4_Treated
    ## 1    AIMP1P1 chr20       0.000       0.000       0.000       0.000
    ## 2 AL035045.1 chr20       0.000       0.000       0.000       0.000
    ## 3 AL136090.1 chr20       0.000       0.000       0.000       0.000
    ## 4      BANF2 chr20       0.000       0.000       0.000       0.000
    ## 5      BTBD3 chr20    2309.398    2782.401    2790.122    2900.808
    ## 6  C20orf187 chr20       0.000       0.000       0.000       0.000
    ##   SW5_Treated  SW6_Treated         mean
    ## 1       0.000    0.9892566    0.1648761
    ## 2       0.000    0.0000000    0.0000000
    ## 3       0.000    0.0000000    0.0000000
    ## 4       0.000    0.0000000    0.0000000
    ## 5    3200.091 2960.8449159 2823.9441132
    ## 6       0.000    0.0000000    0.0000000

We can put these steps into **function** to make it more reusable:
**function** takes the form of: function( arglist ) expr return(value)

Let start with a simple function

    add<-function(a,b){
      x=a+b
      return(x)
    }

    add(1,2)

    ## [1] 3

    add(3,4)

    ## [1] 7

    add(5,7)

    ## [1] 12

Let’s put everything together

    curr_coord=list(
      chr="chr20",
      strand="+",
      start=10000000,
      end=20000000
    )

    find_mean<-function(data, gene_expr, curr_coord){
      selected_genes=subset(data$annot, chrom==curr_coord$chr & 
                            strand==curr_coord$strand & 
                              start>curr_coord$start & 
                              end< curr_coord$end,
                             select=c("chrom", "name"))
      gene_mean_results=merge(selected_genes, gene_expr, by.x="name", by.y=0)
      return(gene_mean_results)
    }

    merge_results=find_mean(data, gene_expr, curr_coord)
    head(merge_results==chr20_expr)

    ##      name chrom SW1_Control SW2_Control SW3_Control SW4_Treated
    ## [1,] TRUE  TRUE        TRUE        TRUE        TRUE        TRUE
    ## [2,] TRUE  TRUE        TRUE        TRUE        TRUE        TRUE
    ## [3,] TRUE  TRUE        TRUE        TRUE        TRUE        TRUE
    ## [4,] TRUE  TRUE        TRUE        TRUE        TRUE        TRUE
    ## [5,] TRUE  TRUE        TRUE        TRUE        TRUE        TRUE
    ## [6,] TRUE  TRUE        TRUE        TRUE        TRUE        TRUE
    ##      SW5_Treated SW6_Treated mean
    ## [1,]        TRUE        TRUE TRUE
    ## [2,]        TRUE        TRUE TRUE
    ## [3,]        TRUE        TRUE TRUE
    ## [4,]        TRUE        TRUE TRUE
    ## [5,]        TRUE        TRUE TRUE
    ## [6,]        TRUE        TRUE TRUE

Lastly, we filter out genes with means of 0

    final_results=subset(merge_results, mean>0, select=c("name", "chrom", "mean"))
    head(final_results)

    ##       name chrom         mean
    ## 1  AIMP1P1 chr20 1.648761e-01
    ## 5    BTBD3 chr20 2.823944e+03
    ## 7  CSRP2BP chr20 1.096049e+03
    ## 8     DSTN chr20 3.369902e+04
    ## 9     DTD1 chr20 3.715940e+03
    ## 10 GAPDHP2 chr20 8.870973e+00

    library(writexl)
    write_xlsx(final_results, "Outputs/final_results.xlsx")
