setwd('YOURPATH')
set.seed(1234)
library(data.table)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(Matrix)

atac_raw = readRDS('atac_matrix.binary.qc_filtered.rds')
metadata.raw = fread('cell_metadata.tissue_freq_filtered.txt')

## metadata:
## BoneMarrow, Cerebellum, Kidney, Liver, Lung
## PreFrontalCortex, SmallIntestine, Spleen
## Testes, Thymus, WholeBrain

tissue_list = c('BoneMarrow','Kidney',
                'Lung','SmallIntestine',
                'Spleen','Testes','Thymus','WholeBrain','Liver')
metadata.extract = metadata.raw[metadata.raw$tissue %in% tissue_list,]

cell_names = metadata.extract$cell

## extract needy cells
atac.extract = atac_raw[,cell_names]
## filter sites: threshold for site-frequency is 0.05
site_freq_threshold = 0.05
num_cells_ncounted = rowSums(atac.extract)
threshold = ncol(atac.extract) * site_freq_threshold
atac_site_filtered = atac.extract[num_cells_ncounted >= threshold,]

# liftover from mm9 genome to mm10 genome
peaks_mm9 = StringToGRanges(rownames(atac_site_filtered),sep = c('_',"_"))
mm9_mm10 = rtracklayer::import.chain('YOURPATH/mm9ToMm10.over.chain')
peaks_mm10 = rtracklayer::liftOver(peaks_mm9, chain = mm9_mm10)
names(peaks_mm10) = rownames(atac_site_filtered)

correspondence <- S4Vectors::elementNROWS(peaks_mm10)
peaks_mm10 <- peaks_mm10[correspondence == 1]
peaks_mm10 <- unlist(peaks_mm10)
atac_site_filtered.transfer <- atac_site_filtered[names(peaks_mm10), ]
rowname.atac = gsub("-","_",GRangesToString(grange = peaks_mm10))
rownames(atac_site_filtered.transfer) = rowname.atac

saveRDS(atac_site_filtered.transfer, file = 'atac_site_filtered.transfer.rds')