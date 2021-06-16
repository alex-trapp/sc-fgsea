# sc-fgsea

sc-fgsea is an automated tool for cluster annotation of single-cell RNA sequencing data. The tool runs on the fGSEA framework and leverages cell-type specific genesets to enable itterative, supervised cluster annotation based on computed p-values and enrichment scores.

The documented code for sc-fgsea is provided as an R markdown file, along with the required data files (located in the `data` folder) to reproduce the analysis.

The output of sc-fgsea is an Excel file containing p-values and enrichment scores for each geneset and each cluster in the input files. An example output file can be found in the `output` folder.

