# sc-fgsea

sc-fgsea is an automated tool for cluster annotation of single-cell RNA sequencing data. The tool runs on the fGSEA framework and leverages cell-type specific genesets to enable itterative, supervised annotation based on computed p-values and enrichment scores.

The documented code for sc-fgsea is provided as an R markdown file, which can easily be loaded into RStudio.

## Input

Required data files (located in the `data` folder) are also provided to run the tool and reproduce the analysis.

The input to sc-fgsea consists of two components:
* a .txt markers file, outputted by the function [Seurat](https://satijalab.org/seurat/)'s FindAllMarkers function
* a .gmt atlas file (see [here for a description of this file format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)), containing a variety of genesets specific to particular cell types

Example markers and atlas files are provided in `data`:
* atlasV14.gmt
* hglBM.NOmp.markers.clust.txt

Note: some genes in the markers file may contain p-values of 0, resulting from underflow errors during processing. To address this, genes with p-values of 0 are ranked in a cluster-wise basis by their fold change, and p-values of these genes are modified to be non-zero while preserving the integrity of the rank. This enables piping our markers data into fgsea.

## Output

The output of sc-fgsea is an Excel file, containing p-values and enrichment scores for each geneset/cluster pair from the input marker file and atlas. An example output file (`hglBM.cluster.markers.output.xlsx`) can be found in the `output` folder, with the following format:

pathway |	pval | padj |	ES | NES | nMoreExtreme | size | leadingEdge | cluster |
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
ROSSI_mLTHSC | 0.009482058 | 0.039587592 | -0.854460094 | -1.738297467 | 384 | 3 | CALML4, UPP1, CLU | 1 |
