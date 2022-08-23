---
title: "Data Wrangling with R Basics - 2"
date: 2022-04-25
draft: false
tags: ["R", "Data Wrangling", "tutorial"]
categories: ["Tutorial"]
author: "Yaoyu E. Wang"

autoCollapseToc: true

---


Working with genomic data with R
--------------------------------

We provided 5 files in the **Data** directory: -
differntial\_results.csv - gene\_expression.tsv - gene\_annotations.tsv
- my\_pathway\_genes.txt - mutations.tsv

These files are stored in either comma-delimited (csv) or tab-delimited
(tsv) formats. Files can be read with the following functions available
from default R libraries. These are built-in “functions” from R. All
read files, but with slightly different default behavior.

<table>
<colgroup>
<col style="width: 35%" />
<col style="width: 35%" />
<col style="width: 28%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: center;">Function</th>
<th style="text-align: center;">package</th>
<th style="text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: center;">read.table(<em>file</em>)</td>
<td style="text-align: center;">util</td>
<td style="text-align: left;">Reads a file in table format and creates a data frame</td>
</tr>
<tr class="even">
<td style="text-align: center;">read.delim(<em>file</em>)</td>
<td style="text-align: center;">util</td>
<td style="text-align: left;">Same as read.table with " as default delimiter</td>
</tr>
<tr class="odd">
<td style="text-align: center;">read.csv(<em>file</em>)</td>
<td style="text-align: center;">util</td>
<td style="text-align: left;">Same as read.table with “,” as default delimiter</td>
</tr>
</tbody>
</table>

    # setwd("~/modules/data_wrangling")
    # First, load all the data:
    print('Loading data...')

    ## [1] "Loading data..."

    dge_results <- read.table('scripts and data/data/differential_results.csv', sep=',', header=T)
    expressions = read.table('scripts and data/data/gene_expression.tsv', header=T)
    annotations = read.delim('scripts and data/data/gene_annotations.tsv')
    pathways = read.table('scripts and data/data/my_pathway_genes.txt', sep='\t', col.names=c('gene_name','pathway'))
    mutations = read.delim('scripts and data/data/mutations.tsv')
    print('Done loading data.')

    ## [1] "Done loading data."

Reading files with read.table (or equivalents) automatically creates a
dataframe. R’s version of a Excel spreadsheet- contains different data
(of different types!) in the columns and “observations” in the rows.

Data Inspection
---------------

Data can be inspected after loaded with these common commands:  
**class()**: data type (e.g. character, numeric, etc.) of vectors and
data structure of dataframes, matrices, and lists. **summary()**:
detailed display, including descriptive statistics, frequencies
**head()**: will display the top 5 entries for the variable **tail()**:
will display the bottom 5 entries for the variable

Specifically, for Dataframe and matrix data structure: **dim()**:
returns dimensions of the dataset **nrow()**: returns the number of rows
in the dataset **ncol()**: returns the number of columns in the dataset
**rownames()**: returns the row names in the dataset **colnames()**:
returns the column names in the dataset

    print("Show the first six rows of each variable")

    ## [1] "Show the first six rows of each variable"

    head(dge_results)

    ##         gene   baseMean        C          T foldChange log2FoldChange
    ## 1    DDX11L1  0.0000000  0.00000  0.0000000         NA             NA
    ## 2     WASH7P 35.6756001 36.87866 34.4725431  0.9347559    -0.09733839
    ## 3 MIR1302-10  0.0000000  0.00000  0.0000000         NA             NA
    ## 4    FAM138A  0.0000000  0.00000  0.0000000         NA             NA
    ## 5     OR4G4P  0.0000000  0.00000  0.0000000         NA             NA
    ## 6    OR4G11P  0.4759207  0.00000  0.9518413        Inf            Inf
    ##        pval padj
    ## 1        NA   NA
    ## 2 0.8764378    1
    ## 3        NA   NA
    ## 4        NA   NA
    ## 5        NA   NA
    ## 6 0.6676182    1

    head(expressions)

    ##         gene SW1_Control SW2_Control SW3_Control SW4_Treated SW5_Treated
    ## 1    DDX11L1     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 2     WASH7P    45.77888    22.84061    42.01649    50.16839   29.507082
    ## 3 MIR1302-10     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 4    FAM138A     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 5     OR4G4P     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 6    OR4G11P     0.00000     0.00000     0.00000     0.00000    2.855524
    ##   SW6_Treated
    ## 1     0.00000
    ## 2    23.74216
    ## 3     0.00000
    ## 4     0.00000
    ## 5     0.00000
    ## 6     0.00000

    head(annotations)

    ##   chrom start   end strand         name
    ## 1  chr1 11869 14409      +      DDX11L1
    ## 2  chr1 14404 29570      -       WASH7P
    ## 3  chr1 17369 17436      -    MIR6859-3
    ## 4  chr1 29554 31109      + RP11-34P13.3
    ## 5  chr1 30366 30503      +    MIR1302-9
    ## 6  chr1 34554 36081      -      FAM138A

    head(pathways)

    ##   gene_name    pathway
    ## 1      PCK2 glycolysis
    ## 2      PCK1 glycolysis
    ## 3      FBP2 glycolysis
    ## 4      BPGM glycolysis
    ## 5   ALDH3A2 glycolysis
    ## 6   ALDH3A1 glycolysis

    head(mutations)

    ##   chrom      pos ref alt
    ## 1     5 20578198   C   T
    ## 2     2  4642922   T   A
    ## 3     1 15947000   C   A
    ## 4    16  3573172   G   C
    ## 5    13 14306989   C   T
    ## 6    11  5028652   C   A

    print("Done showing the first row of each data")

    ## [1] "Done showing the first row of each data"

Slicing data
------------

To refer to a cell in Excel, you’re used to A1 (first row, first
column), C2 (third row, second column). Similar with a DataFrame:

    head(expressions)

    ##         gene SW1_Control SW2_Control SW3_Control SW4_Treated SW5_Treated
    ## 1    DDX11L1     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 2     WASH7P    45.77888    22.84061    42.01649    50.16839   29.507082
    ## 3 MIR1302-10     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 4    FAM138A     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 5     OR4G4P     0.00000     0.00000     0.00000     0.00000    0.000000
    ## 6    OR4G11P     0.00000     0.00000     0.00000     0.00000    2.855524
    ##   SW6_Treated
    ## 1     0.00000
    ## 2    23.74216
    ## 3     0.00000
    ## 4     0.00000
    ## 5     0.00000
    ## 6     0.00000

    # These two lines return the same number
    expressions[2,2]

    ## [1] 45.77888

    expressions[2, 'SW1_Control']

    ## [1] 45.77888

    # Select First Column:
    head(expressions[1])

    ##         gene
    ## 1    DDX11L1
    ## 2     WASH7P
    ## 3 MIR1302-10
    ## 4    FAM138A
    ## 5     OR4G4P
    ## 6    OR4G11P

    head(expressions['gene'])  

    ##         gene
    ## 1    DDX11L1
    ## 2     WASH7P
    ## 3 MIR1302-10
    ## 4    FAM138A
    ## 5     OR4G4P
    ## 6    OR4G11P

    head(expressions[,1])

    ## [1] DDX11L1    WASH7P     MIR1302-10 FAM138A    OR4G4P     OR4G11P   
    ## 56638 Levels: 5S_rRNA 7SK A1BG A1BG-AS1 A1CF A2M A2M-AS1 A2ML1 ... ZZZ3

    head(expressions$gene)

    ## [1] DDX11L1    WASH7P     MIR1302-10 FAM138A    OR4G4P     OR4G11P   
    ## 56638 Levels: 5S_rRNA 7SK A1BG A1BG-AS1 A1CF A2M A2M-AS1 A2ML1 ... ZZZ3

    # select First Row:
    head(expressions[1,])

    ##      gene SW1_Control SW2_Control SW3_Control SW4_Treated SW5_Treated
    ## 1 DDX11L1           0           0           0           0           0
    ##   SW6_Treated
    ## 1           0

    # select first three columns
    head(expressions[1:3])

    ##         gene SW1_Control SW2_Control
    ## 1    DDX11L1     0.00000     0.00000
    ## 2     WASH7P    45.77888    22.84061
    ## 3 MIR1302-10     0.00000     0.00000
    ## 4    FAM138A     0.00000     0.00000
    ## 5     OR4G4P     0.00000     0.00000
    ## 6    OR4G11P     0.00000     0.00000

    # select all columns except first:
    head(expressions[2:ncol(expressions)])

    ##   SW1_Control SW2_Control SW3_Control SW4_Treated SW5_Treated SW6_Treated
    ## 1     0.00000     0.00000     0.00000     0.00000    0.000000     0.00000
    ## 2    45.77888    22.84061    42.01649    50.16839   29.507082    23.74216
    ## 3     0.00000     0.00000     0.00000     0.00000    0.000000     0.00000
    ## 4     0.00000     0.00000     0.00000     0.00000    0.000000     0.00000
    ## 5     0.00000     0.00000     0.00000     0.00000    0.000000     0.00000
    ## 6     0.00000     0.00000     0.00000     0.00000    2.855524     0.00000

Using boolean
-------------

1.  Let’s try a little example to select genes located on chromosome 7
    and 3 from small gene annotation file ‘demo\_annotations.tsv’. This
    can be done in two ways:

#### Excel

    Open “demo_annotations.tsv” in excel
    Insert table
    Select ‘chrom’ column
    Set ‘chrom’ column equal ‘chr7’
    Set ‘chrom’ column equal ‘chr3’

#### R

In R using boolean comparison to create logical vector

    df=read.table("scripts and data/data/demo_annotations.tsv", header=T)
    chroms=df$chrom
    is_chr7 = chroms == 'chr7'
    is_chr7

    ## [1] FALSE  TRUE FALSE  TRUE

    df[is_chr7,]

    ##   chrom start end strand name
    ## 2  chr7   200 275      - TP53
    ## 4  chr7   500 600      + CD44

    is_chr3 <- chroms == 'chr3'
    df[is_chr3,]

    ## [1] chrom  start  end    strand name  
    ## <0 rows> (or 0-length row.names)

1.  Select genes from gene annotations that are in your oncogene list
    (‘KRAS’, ‘TP53’)

#### Excel

    #select gene annotations for (‘KRAS’, ‘TP53’)
    Go to ‘name’ column
    Find ‘KRAS’
    Click on ‘KRAS’
    Find ‘TP53’
    Click on ‘TP53’

#### R

    # select gene that are annotated as oncogenes use %in%
    oncogenes = c('KRAS', 'TP53')  
    is_oncogene = df$name %in% oncogenes
    df[is_oncogene,]

    ##   chrom start end strand name
    ## 1  chr1   100 150      + KRAS
    ## 2  chr7   200 275      - TP53

1.  Select genes from gene annotations that are in your oncogene list
    (‘KRAS’, ‘TP53’) and on chrom 7

#### Excel

    Go to ‘name’ column
    Find ‘KRAS’
    Click on ‘KRAS’
    Find ‘TP53’
    Click on ‘TP53’
    Set ‘chrom’ column equal ‘chr7’

#### R

    # Combine multiple criteria with AND (&) operation
    selection_criteria = is_oncogene & is_chr7
    data_subset = df[selection_criteria, ]
    data_subset

    ##   chrom start end strand name
    ## 2  chr7   200 275      - TP53

End of toy data set…

Real Size Data Set
------------------

#### 1) Use the data set loaded earlier in the script to identify ras signaling genes that are significantly up-regulated from our result files:

Steps to find these genes:

-   Get all the gene names that are in the pathway of interest
    (ras\_signaling) from ‘my\_pathway\_genes.tsv’
-   Find differential gene expression(DGE) results for those RAS genes
    from ‘differential\_results.csv’
-   Find all the differentially expressed genes an adjusted p &lt; 0.05
    threshold from ‘differential\_results.csv’
-   Find all the differentially expressed genes that are up-regulated
    with log fold change (log2FoldChange) &gt; 0 from
    ‘differential\_results.csv’
-   Keep only genes that pass all 3 “tests” (are True for all 3
    conditions)
-   Just in case, remove missing data

<!-- -->

    print('Filtering for ras signaling genes that are significantly upregulated.')

    ## [1] "Filtering for ras signaling genes that are significantly upregulated."

    # Get a vector of gene names that are in the pathway of interest
    ras_genes = pathways[pathways['pathway'] == 'ras_signaling', 'gene_name']

    #Create a boolean (True or False) vector for those RAS genes
    is_ras_gene = dge_results$gene %in% ras_genes

    # Create a boolean vector for whether the gene is significantly changed
    # at a p < 0.05 threshold
    is_significant = dge_results$padj < 0.05

    # Create a boolean vector for whether the gene is upregulated
    is_upregulated = dge_results$log2FoldChange > 0

    # Keep only rows that pass all 3 "tests" (are True for all 3)
    selected_rows = dge_results[is_ras_gene & is_significant & is_upregulated, ]

    # Just in case, remove missing data
    selected_rows = na.omit(selected_rows)

    # We now have the rows/genes-- now, keep only a subset of the columns:
    upregulated_ras_genes = selected_rows[c('gene','baseMean','log2FoldChange','padj')]
    print('Done filtering genes.')

    ## [1] "Done filtering genes."

#### 2) Find ras signaling genes that are significantly up-regulated

with gene coordinates and gene expression values from our result files

-   Merge result with ‘gene\_annotations.tsv’ on the selected genes
-   Merge result with ‘gene\_expression.tsv’ on the selected genes

<!-- -->

    # Merge with gene coordinates
    ras_up_genes_w_coords = merge(upregulated_ras_genes, annotations, by.x = 'gene', by.y='name')
    # Merge with gene expression
    ras_up_genes_w_coords_expression = merge(ras_up_genes_w_coords, expressions, by='gene')

    # examine the variable
    head(ras_up_genes_w_coords_expression)

    ##     gene   baseMean log2FoldChange         padj chrom     start       end
    ## 1 ANGPT1  1111.0617      0.5136897 2.783113e-04  chr8 107249482 107498055
    ## 2 BCL2L1 10000.0591      0.3472426 8.742192e-08 chr20  31664452  31723989
    ## 3  CALM2 73953.5081      0.2790400 5.369188e-06  chr2  47160082  47176601
    ## 4   FGF1   955.6761      0.2797321 2.553878e-02  chr5 142592178 142698070
    ## 5  FGF14   419.1095      0.9456159 7.832240e-05 chr13 101710804 102402457
    ## 6   FGF2 13101.6427      0.9214510 1.696858e-21  chr4 122826708 122898236
    ##   strand SW1_Control SW2_Control SW3_Control SW4_Treated SW5_Treated
    ## 1      -    765.5787    993.5663    986.8196   1237.1883   1274.5156
    ## 2      -   8944.4139   8877.7279   8584.9902  11402.5583  10930.9460
    ## 3      -  62166.7448  70176.7209  68128.0296  78758.1010  81115.9194
    ## 4      -    819.1497    879.3633    891.4309   1094.7459   1041.3144
    ## 5      -    216.2322    297.9661    345.2165    510.6425    502.5722
    ## 6      +   7597.3460   9885.8292   9679.6900  16386.2503  17452.9629
    ##   SW6_Treated
    ## 1   1408.7014
    ## 2  11259.7183
    ## 3  83375.5330
    ## 4   1008.0524
    ## 5    642.0275
    ## 6  17607.7777

Integrating the SNP data with expression data
---------------------------------------------

**Problem**: The chromosome notations are different and can not be
merged. chrom has chromosome as ‘chr8’ but the *mutations* has chrom as
‘8’.

Strings can be concatenated by *paste()* or *poste0()* | Function |
package |Description |—|—|—| |paste()|base|Concatenate vectors after
converting to character with " " as separator |paste0()|base| paste with
no separator, i.e. paste0(‘a’,‘b’) becomes ‘ab’

    # Select the first row from mutation and change chrom from '5' to 'chr5'
    i=1
    chrom =mutations[i,'chrom']
    chrom_w_prefix = paste('chr', chrom, sep='')  
    mutations[i,'chrom'] =chrom_w_prefix

We can change one instance. but then do we want to change each row
one-by-one for 7479 times? There is a better way to do this…

### Use flow control - the “for loop”

    print('Example using the for loop flow control')

    ## [1] "Example using the for loop flow control"

    my_vector <- c(10,11,12,13) 
    for (item in my_vector){
      print(item)
      #do actual operations here
    }

    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13

Everything between { and } is repeated. The script within the { } is
called a “block” of code. The indentation within the block is not
necessary, but helps with reading. The loop variable *item* is
arbitrary.

### Changing chromosome names in mutation dataframe

    # We need to add the 'chr' prefix to the chromosome names that are in the mutations dataframe.
    # This way, all the chromosome names are consistent.  
    # The method below is a slow way to do this, but easier to understand.
    for ( i in 2:nrow(mutations) )
    {
        chrom = mutations[i,'chrom']
        # Add ‘chr’ to each chrom value
        chrom_w_prefix = paste('chr', chrom, sep='')
        mutations[i, 'chrom'] = chrom_w_prefix
    }

Filter to keep only genes that are mutated. The idea is For each row (a
gene), see if any of the mutations are “inside” that gene e.g. Here, the
first mutation is on the KRAS gene

    # Again, this a slow, but clear way to do this.
    ras_up_mutated_genes = data.frame()                         # define an empty result variable
    for ( i in 1:nrow(ras_up_genes_w_coords) )
    {
      gene_info = ras_up_genes_w_coords[i,]               # get row i and save as gene_info
      same_chrom = gene_info$chrom == mutations$chrom       # get mutation row with same chrom
      past_start = mutations$pos >gene_info$start             # get position > then gene_info start
      before_end = mutations$pos <gene_info$end             # get position < then gene_info start   
      overlap = same_chrom & past_start & before_end        # get mutations fit all 3 conditions
      if ( any(overlap) )                                             # if there is any mutation quality
      {
        ras_up_mutated_genes = rbind(ras_up_mutated_genes, gene_info) # use rbind() to add to result variable
      }
    }

### Writing files

Now, we want to save our clean results to a file. Look at the help from
R’s built-in write.table function.

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr class="header">
<th>Function</th>
<th>package</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>write.table()</td>
<td>utils</td>
<td>write a variable into a file with provided filename with " as delimitor</td>
</tr>
<tr class="even">
<td>writexl()</td>
<td>writexl</td>
<td>Writes a data frame to an xlsx file</td>
</tr>
</tbody>
</table>

    print('Merge to add in the expression data')

    ## [1] "Merge to add in the expression data"

    # note no by.x or by.y since the 'gene' column is in each dataframe
    final_data = merge(ras_up_mutated_genes, expressions)
    head(final_data)

    ##     gene    baseMean log2FoldChange         padj chrom     start       end
    ## 1  CALM2 73953.50812      0.2790400 5.369188e-06  chr2  47160082  47176601
    ## 2 MAPK10    99.11630      0.8205749 1.241397e-02  chr4  86015123  86594131
    ## 3  RAB5A  4408.97751      0.3199285 2.044410e-05  chr3  19947079  19985175
    ## 4   RGL1  3355.69859      0.2658060 1.414847e-03  chr1 183636085 183928531
    ## 5   RRAS  2591.26593      0.3950969 3.138126e-03 chr19  49635292  49640201
    ## 6    TEK    33.28648      1.4060046 1.770961e-02  chr9  27109141  27230175
    ##   strand SW1_Control SW2_Control SW3_Control SW4_Treated SW5_Treated
    ## 1      - 62166.74477 70176.72087 68128.02960 78758.10104 81115.91944
    ## 2      -    70.12935    77.86570    66.99926   145.12998   122.78753
    ## 3      +  3531.79188  4172.56326  4061.97223  4860.06277  4877.23505
    ## 4      +  3176.27504  3054.41182  2911.62896  3743.81609  3554.17558
    ## 5      -  2598.68194  2074.34222  2042.90971  3009.20753  2851.71667
    ## 6      +    22.40243    18.68777    13.62697    51.06425    51.39943
    ##   SW6_Treated
    ## 1 83375.53303
    ## 2   111.78599
    ## 3  4950.23988
    ## 4  3693.88403
    ## 5  2970.73748
    ## 6    42.53803

We will create an ‘Outputs’ directory, if not exist, to store all of the
output files.

    # get current working directory. 
    getwd()

    ## [1] "/Users/yaoyuwang/Dropbox (Harvard University)/qBRC-Share/Presentations/QBRC-Workshops/Data Wrangling with R/2022-0424-IID"

    # Current directory is Data directory, change to output
    if(!dir.exists('Outputs')){
      dir.create('Outputs')
    }
    # write final data out as a tab-seperated-value file (.tsv)
    print('Write the final data to file.')

    ## [1] "Write the final data to file."

    write.table(final_data, 'Outputs/final_data.tsv', sep='\t', quote=F)

write out file in Excel format using **writexl**. Check if **writexl**
package is available, if not install the package using function
*install.packages()*.

    if(!require(readxl)) install.packages("readxl")
    if(!require(writexl)) install.packages("writexl")

    library(readxl)    # load  library to write excel file
    library(writexl)    # load  library to write excel file

    print('Write the final data to Excel file.')
    # start a new excel file and write final results 
    write_xlsx(final_data, "Outputs/final_results.xlsx", col_names=TRUE)
    # add a new sheet with new data onto the file we just created
    write_xlsx(ras_up_genes_w_coords, "Outputs/ras_up_genes_w_coords.xlsx", col_names=TRUE)
