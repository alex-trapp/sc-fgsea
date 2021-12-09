# sc-fgsea

sc-fgsea is an automated tool for cluster annotation of single-cell RNA sequencing data. The tool runs on the fGSEA framework and leverages cell-type specific genesets to enable itterative, supervised annotation based on computed p-values and enrichment scores.

The documented code for sc-fgsea is provided as an R markdown file, which can easily be loaded into RStudio.

This tool was developed by Alexandre Trapp within the context of the [Emmrich et al. preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2021.08.01.454652v1):

```
Naked mole-rat hematopoietic stem and progenitors are highly quiescent with an inherent myeloid bias
Stephan Emmrich, Alexandre Trapp, Frances Tolibzoda Zakusilo, Marco Mariotti, Maggie E. Straight,
Zhihui Zhang, Spencer Gray, Michael G. Drage, Zhonghe Ke, Masaki Takasugi, Jan-Henning Klusmann,
Vadim N. Gladyshev, Andrei Seluanov, Vera Gorbunova
bioRxiv 2021.08.01.454652
```

## Input

Required data files (located in the `data` folder) are also provided to run the tool and reproduce the analysis.

The input to sc-fgsea consists of two components:
* a .txt markers file, outputted by [Seurat](https://satijalab.org/seurat/)'s FindAllMarkers function
* a .gmt atlas file (see [here for a description of this file format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)), containing a variety of genesets specific to particular cell types

Example markers and atlas files are provided in `data`:
* atlasV14.gmt
* hglBM.NOmp.markers.clust.txt

## Output

The output of sc-fgsea is an Excel file, containing p-values and enrichment scores for each geneset/cluster pair from the input marker file and atlas. An example output file (`hglBM.cluster.markers.output.xlsx`) can be found in the `output` folder, with the following format:

pathway |	pval | padj |	ES | NES | nMoreExtreme | size | leadingEdge | cluster |
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
ROSSI_mLTHSC | 0.009482058 | 0.039587592 | -0.854460094 | -1.738297467 | 384 | 3 | CALML4, UPP1, CLU | 1 |
