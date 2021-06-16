# sc-fgsea

sc-fgsea is an automated tool for cluster annotation of single-cell RNA sequencing data. The tool runs on the fGSEA framework and leverages cell-type specific genesets to enable itterative, supervised cluster annotation based on computed p-values and enrichment scores.

The documented code for sc-fgsea is provided as an R markdown file, which can easily be loaded into RStudio.

## Data

Required data files (located in the `data` folder) are also provided to run the tool and reproduce the analysis.

The input to sc-fgsea consists of two components:
* a .txt markers file, outputted by the function FindAllMarkers by [Seurat](https://satijalab.org/seurat/) (example provided in `data/hglBM.NOmp.markers.clust.txt)
* a .gmt atlas file (see [here for a description of this file format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)), containing a variety of genesets specific to particular cell types

The output of sc-fgsea is an Excel file, with the following format:



containing p-values and enrichment scores for each geneset and each cluster in the input files. An example output file can be found in the `output` folder.

