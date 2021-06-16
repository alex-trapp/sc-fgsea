sc-fgsea v1.0 (6/16/2021)

Developed by Alexandre Trapp

############# 

Load packages \#\#\#\#\#\#\#\#\#\#\#\#\#

``` r
library(fgsea)
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.0.5

``` r
library(hash)
library(tictoc)
```

    ## Warning: package 'tictoc' was built under R version 4.0.5

``` r
library(xlsx)
```

    ## Warning: package 'xlsx' was built under R version 4.0.5

############ 

Load markers \#\#\#\#\#\#\#\#\#\#\#\#

``` r
# read-in input marker table
S.ALL.markers <- read.table("./data/hglBM.NOmp.markers.clust.txt",
                            header = TRUE)

# ensure clusters is of datatype double (starting with 1)
S.ALL.markers$cluster <- as.numeric(S.ALL.markers$cluster)

cluster_list <- sort(unique(S.ALL.markers$cluster))

cat("Number of clusters: ")
```

    ## Number of clusters:

``` r
cat(length((unique(S.ALL.markers$cluster))))
```

    ## 14

``` r
cat("\nList: ")
```

    ## 
    ## List:

``` r
cat(cluster_list)
```

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14

########## 

Load atlas \#\#\#\#\#\#\#\#\#\#

``` r
# load gmt file with geneset information
atlas <- gmtPathways("./data/atlasV14.gmt")
```

################## 

P-value processing \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

``` r
# initialize cluster dataframe dictionaries
clust_dict <- hash(keys = sprintf("Cluster_%s", cluster_list), values = cluster_list)
clust_dict_pval <- hash(keys = sprintf("Cluster_%s", cluster_list), values = cluster_list)

# internally rank p-values with a value of 0 based on fold change
for (x in cluster_list) {
  cat("\nCluster =", x)
  
  # filter dataframe into single celltype elements
  clusterdf <- dplyr::filter(S.ALL.markers, cluster == x)
  
  # add sign column
  clusterdf$fcSign <- sign(clusterdf$avg_log2FC)
  
  # fix p-values equal to 0
  
  # first, check if p-values of 0 are detected in given cluster
  if (0 %in% clusterdf$p_val_adj) { 
    
    # if so, count how many such values exist
    first_nonzero_p <- 0
    for (p in clusterdf$p_val_adj) {
      if (p == 0) { 
        first_nonzero_p <- first_nonzero_p + 1
      }
    }
    
    cat("\nNumber of (p_val_adj = 0) values:", first_nonzero_p, "\n")
    
    # create decreasing sequence from the last nonzero value 
    # (i.e. zero p-value closest to nonzero p-value)
    zero_values <- seq(first_nonzero_p, 1, by =-1) 
    
    # loop through all zero p-values starting with the one adjacent to nonzero value
    # for each p-value of 0, divide the previous non-zero value by the reduction
    # factor. this enable a valid nonzero p-value for each gene in each cluster
    reduction_factor <- 1.01
    for (zero_ind in zero_values) {
      clusterdf$p_val_adj[zero_ind] <- (clusterdf$p_val_adj[zero_ind+1]) / reduction_factor
    }
  
    cat("Modified!\n")  
  }
  
  else { # if there are no non-zero values, all is well
    cat("\nNo (p-value = 0) genes in this cluster \n")
  }
  
  # add logP column (-log10 transformation of the p-value)
  clusterdf$logP <- -log10((clusterdf$p_val_adj)) 
  
  # add metric column (-log10 p-value with sign information)
  clusterdf$metric <- clusterdf$logP/clusterdf$fcSign 
  
  # assign dataframe as value in hash table to corresponding key ("Cluster_x")
  clust_dict[[sprintf("Cluster_%s", x)]] <- clusterdf
  
  # assign metric/gene vectors in hash table to corresponding key ("Cluster_x")
  clust_dict_pval[[sprintf("Cluster_%s", x)]] <- setNames(clusterdf$metric, clusterdf$gene)
}
```

    ## 
    ## Cluster = 1
    ## Number of (p_val_adj = 0) values: 137 
    ## Modified!
    ## 
    ## Cluster = 2
    ## No (p-value = 0) genes in this cluster 
    ## 
    ## Cluster = 3
    ## Number of (p_val_adj = 0) values: 179 
    ## Modified!
    ## 
    ## Cluster = 4
    ## Number of (p_val_adj = 0) values: 1 
    ## Modified!
    ## 
    ## Cluster = 5
    ## No (p-value = 0) genes in this cluster 
    ## 
    ## Cluster = 6
    ## Number of (p_val_adj = 0) values: 161 
    ## Modified!
    ## 
    ## Cluster = 7
    ## Number of (p_val_adj = 0) values: 179 
    ## Modified!
    ## 
    ## Cluster = 8
    ## Number of (p_val_adj = 0) values: 574 
    ## Modified!
    ## 
    ## Cluster = 9
    ## Number of (p_val_adj = 0) values: 30 
    ## Modified!
    ## 
    ## Cluster = 10
    ## Number of (p_val_adj = 0) values: 4 
    ## Modified!
    ## 
    ## Cluster = 11
    ## Number of (p_val_adj = 0) values: 93 
    ## Modified!
    ## 
    ## Cluster = 12
    ## Number of (p_val_adj = 0) values: 135 
    ## Modified!
    ## 
    ## Cluster = 13
    ## Number of (p_val_adj = 0) values: 14 
    ## Modified!
    ## 
    ## Cluster = 14
    ## Number of (p_val_adj = 0) values: 79 
    ## Modified!

##### 

fgsea \#\#\#\#\#

``` r
# initialize hash dictionary to store fgsea results
fgsea_dict <- hash(keys = sprintf("Cluster_%s", cluster_list), values = cluster_list)

for (clust in cluster_list) {
  tic(sprintf("Cluster %s runtime", clust))
  
  # run fgsea analysis
  fgsea_dict[[sprintf("Cluster_%s", clust)]] <- fgsea(pathways = atlas, 
                                                      stats = clust_dict_pval[[sprintf("Cluster_%s", clust)]], 
                                                      nperm = 100000)
  
  # add cluster label column to the generated dataframe
  fgsea_dict[[sprintf("Cluster_%s", clust)]]$cluster <- rep(sprintf("%s", clust), 
                                                            dim(fgsea_dict[[sprintf("Cluster_%s",
                                                                                    clust)]])[1])
  toc()
}
```

    ## Cluster 1 runtime: 30.38 sec elapsed
    ## Cluster 2 runtime: 31.79 sec elapsed
    ## Cluster 3 runtime: 24.97 sec elapsed
    ## Cluster 4 runtime: 29.17 sec elapsed
    ## Cluster 5 runtime: 23.55 sec elapsed
    ## Cluster 6 runtime: 23.78 sec elapsed
    ## Cluster 7 runtime: 24.1 sec elapsed
    ## Cluster 8 runtime: 23.06 sec elapsed
    ## Cluster 9 runtime: 20.64 sec elapsed
    ## Cluster 10 runtime: 20.53 sec elapsed
    ## Cluster 11 runtime: 24.86 sec elapsed
    ## Cluster 12 runtime: 21.91 sec elapsed
    ## Cluster 13 runtime: 22.79 sec elapsed
    ## Cluster 14 runtime: 20.77 sec elapsed

################## 

Combine dataframes \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

``` r
# initialize a dataframe containing data for the first cluster/celltype
first_cluster_df = fgsea_dict[[sprintf("Cluster_%s", cluster_list[1])]]

# use rbind to progressively add clusters to this dataframe
for (clust in cluster_list[-1]) {
  first_cluster_df <- rbind(first_cluster_df, fgsea_dict[[sprintf("Cluster_%s", clust)]])
}

# dataframe with all clusters combined
head(first_cluster_df, 3)
```

    ##         pathway         pval         padj         ES       NES nMoreExtreme
    ## 1:   EBERT_hHSC 8.292888e-01 0.8845294543 -0.4376689 -0.747799        36530
    ## 2:    EBERT_hBC 1.329628e-05 0.0005647463  0.7081955  2.496312            0
    ## 3: EBERT_hCD34+ 3.449820e-01 0.5593397688  0.3999120  1.103458        23494
    ##    size                               leadingEdge cluster
    ## 1:    2                                    SPINK2       1
    ## 2:   17 MS4A1,CD79B,HLA-DQA1,CD79A,LY86,BANK1,...       1
    ## 3:    8                  SYPL1,BCL11A,RPL41,RPL31       1

########################## 

Transform leadingEdge data
\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

``` r
# copy previous dataframe
final_cluster_df_modLE <- first_cluster_df

# loop through every row of this dataframe to fix leadingEdge problem
for (row in 1:nrow(final_cluster_df_modLE)) {
  
  # get list of leading edge genes as a string
  leading_Edge_character = as.character(final_cluster_df_modLE[row, leadingEdge])
  
  # modify output string to remove unecessary components
  
  # remove "
  new_string_v1 <- gsub('"', "", leading_Edge_character, fixed = TRUE)
  
  # remove (
  new_string_v2 <- gsub("(", "", new_string_v1, fixed = TRUE)
  
  # remove )
  new_string_v3 <- gsub(")", "", new_string_v2, fixed = TRUE)
  
  # remove c
  new_string_v4 <- gsub("c", "", new_string_v3, fixed = TRUE)
  
  # set leadingEdge column to the new transformed strings
  final_cluster_df_modLE[row, "leadingEdge"] = new_string_v4
}

# Changes the final datatype to character to allow export
final_cluster_df_modLE$leadingEdge <- as.character(final_cluster_df_modLE$leadingEdge)

# final modified gsea dataframe before export
head(final_cluster_df_modLE, 3)
```

    ##         pathway         pval         padj         ES       NES nMoreExtreme
    ## 1:   EBERT_hHSC 8.292888e-01 0.8845294543 -0.4376689 -0.747799        36530
    ## 2:    EBERT_hBC 1.329628e-05 0.0005647463  0.7081955  2.496312            0
    ## 3: EBERT_hCD34+ 3.449820e-01 0.5593397688  0.3999120  1.103458        23494
    ##    size
    ## 1:    2
    ## 2:   17
    ## 3:    8
    ##                                                                                           leadingEdge
    ## 1:                                                                                             SPINK2
    ## 2: MS4A1, CD79B, HLA-DQA1, CD79A, LY86, BANK1, HLA-DMB, RALGPS2, BLNK, HLA-DOB, CD37, CD40, CR2, CD22
    ## 3:                                                                        SYPL1, BCL11A, RPL41, RPL31
    ##    cluster
    ## 1:       1
    ## 2:       1
    ## 3:       1

############### 

Export to Excel \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

``` r
sheet_name = "hglBM"
output_path = "./output/hglBM.cluster.markers.output.xlsx"

write.xlsx(final_cluster_df_modLE, 
           output_path,
           sheetName = sheet_name,
           row.names = FALSE)

cat("Output file written to Excel!")
```

    ## Output file written to Excel!
