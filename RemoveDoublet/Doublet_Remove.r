# assuming there are 2 samples ready to remove the doublet. The pipeline is adapted if there is metadata(for example, cell types).
# Finally, we will obtain a h5seurat file. Load the h5seurat file and save as RData we will obtain Tissue_unintegrated.RData occurred in other original code files.

rm(list = ls())
setwd('YOURPATH')
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(DoubletFinder)
library(data.table)
library(sctransform)
set.seed(1234)

Convert('Tissue1_matrix.h5ad', dest = 'h5seurat', overwrite = F)
Convert('Tissue2_matrix.h5ad', dest = 'h5seurat', overwrite = F)

Tissue1 = LoadH5Seurat('Tissue1_matrix.h5seurat')
Tissue2 = LoadH5Seurat('Tissue2_matrix.h5seurat')

## preprocessing
if(T)
{
  Tissue1 = NormalizeData(Tissue1)
  Tissue1 = FindVariableFeatures(Tissue1, selection.method = 'vst',nfeatures = 2000)
  Tissue1 = ScaleData(Tissue1)
  Tissue1 = RunPCA(Tissue1)
  Tissue1 = RunUMAP(Tissue1, dims = 1:20) 
}
if(T)
{
  Tissue2 = NormalizeData(Tissue2)
  Tissue2 = FindVariableFeatures(Tissue2, selection.method = 'vst',nfeatures = 2000)
  Tissue2 = ScaleData(Tissue2)
  Tissue2 = RunPCA(Tissue2)
  Tissue2 = RunUMAP(Tissue2, dims = 1:20)
}


## remove doublet
#find PK
if(T)
{
  sweep.res.list.Tissue1 = paramSweep_v3(Tissue1, PCs = 1:20, sct = FALSE)
  sweep.stats_Tissue1 = summarizeSweep(sweep.res.list.Tissue1, GT = FALSE)
  bcmvn_Tissue1 = find.pK(sweep.stats_Tissue1)
  
}
if(T)
{
  sweep.res.list.Tissue2 = paramSweep_v3(Tissue2, PCs = 1:20, sct = FALSE)
  sweep.stats_Tissue2 = summarizeSweep(sweep.res.list.Tissue2, GT = FALSE)
  bcmvn_Tissue2 = find.pK(sweep.stats_Tissue2)
}

## Tissue1 pK = 0.01 nExp.poi = 158
if(T)
{
  nExp.poi1 = round(0.075 * nrow(Tissue1@meta.data)) 
  #assuming 7.5% are doublets 
  Tissue1 = doubletFinder_v3(Tissue1, PCs = 1:20, pN = 0.25, 
                            pK = 0.01, nExp = nExp.poi1, 
                            reuse.pANN = FALSE, sct = FALSE)
}
## Tissue2 pk = 0.05 nExp.poi = 69
if(T)
{
  nExp.poi2 = round(0.075 * nrow(Tissue2@meta.data)) 
  #assuming 7.5% are doublets 
  Tissue2 = doubletFinder_v3(Tissue2, PCs = 1:20, pN = 0.25, 
                            pK = 0.05, nExp = nExp.poi2, 
                            reuse.pANN = FALSE, sct = FALSE)
}

DimPlot(Tissue1, reduction = 'umap',
        group.by = 'DF.classifications_0.25_0.01_158')
DimPlot(Tissue2, reduction = 'umap',
        group.by = 'DF.classifications_0.25_0.05_69')

Tissue1.filtered = subset(Tissue1,
                         `DF.classifications_0.25_0.01_158` == 'Singlet')
Tissue2.filtered = subset(Tissue2,
                         `DF.classifications_0.25_0.05_69` == 'Singlet')

save(Tissue1.filtered, 
     Tissue2.filtered,
     file = 'Tissue_unintegrated.RData')

Tissue.list = list(Tissue1.filtered,
                  Tissue2.filtered)

SelectIntegrationFeatures(Tissue.list)
Tissue.anchors = FindIntegrationAnchors(Tissue.list,
                                       dims = 1:20,
                                       k.anchor = 5,
                                       k.filter = 30)
Tissue.integrated = IntegrateData(Tissue.anchors,
                                 dims = 1:20)

Tissue.integrated = ScaleData(Tissue.integrated)
Tissue.integrated = RunPCA(Tissue.integrated)
Tissue.integrated = FindNeighbors(Tissue.integrated, dims = 1:20)
Tissue.integrated = FindClusters(Tissue.integrated, resolution = 0.5)
Tissue.integrated = RunUMAP(Tissue.integrated, dims = 1:20)

DimPlot(Tissue.integrated,reduction = 'umap')

## save files & convert to h5ad
setwd('YOURPATH')
SaveH5Seurat(Tissue.integrated, 
             filename = 'Tissue_unintegrated.h5seurat')