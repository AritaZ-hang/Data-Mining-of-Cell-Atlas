# Data-Mining-of-Cell-Atlas

## Introduction
A code repository of my bachelor's graduation project: `Data Mining of Cell Atlas`, which aims to construct a multi-omics cell atlas for Mus Musculus. \
\
*Chosen Tissues*: \
Brain(WholeBrain in scATAC-seq), Bone Marrow, Kidney, Liver, Lung, Small Intestine, Spleen, Testis, Thymus \
\
Here comes the data resources:
1. scRNA-seq data(batch effect has already been removed): https://figshare.com/s/865e694ad06d5857db4b, for the expression matrix of 9 tissues over 20 samples.
2. scATAC-seq data: http://atlas.gs.washington.edu/mouse-atac/data, for the peak count matrix and gene activity matrix. 


## Update
1. Pseudo-time analysis(see the folder`Pseudotime`, updated on 2022/05/08)
2. Doublet Removement(see the folder `RemoveDoublet`, updated on 2022/05/12)
3. Multi-omics analysis, including processing ATAC peak count matrix(filter sites and liftover), Seurat/Signac pipeline, label transfer, Spearman correlation and so on(see the folder `Multi-Omics-Analysis`, updated on 2022/05/12)
