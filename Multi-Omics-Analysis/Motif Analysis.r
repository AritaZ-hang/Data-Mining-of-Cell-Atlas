library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(chromVAR)
library(viridis)
set.seed(1234)

DefaultAssay(atac_Tissue) = 'peaks'

pfm = getMatrixSet(x = JASPAR2020, opts = list(collection = 'CORE',
                                               tax_group = 'vertebrates', all_versions = FALSE))

atac_Tissue = AddMotifs(object = atac_Tissue,
                       genome = BSgenome.Mmusculus.UCSC.mm10,
                       pfm = pfm)

# Find Differentially Accessible Peaks
da_peaks_pertype = function(atac)
{
  types = unique(atac$refined.anno)
  da_peaks_all = list()
  for(i in types)
  {	
    da_peaks <- FindMarkers(
      object = atac,
      group.by = 'refined.anno',
      ident.1 = i,
      only.pos = TRUE,
      test.use = 'LR',
      min.pct = 0.05,
      latent.vars = 'nCount_peaks'
    )
    top.da.peak = rownames(da_peaks[da_peaks$p_val < 0.05,])
    top.da.peak = list(top.da.peak)
    da_peaks_all = c(da_peaks_all, top.da.peak)
    
  }
  names(da_peaks_all) = types
  return(da_peaks_all)
}

da_peaks_all = da_peaks_pertype(atac_Tissue)

# Motif Enrichment
motif_enrichment_analysis = function(atac, da_peaks_all)
{
  types = unique(atac$refined.anno)
  len = length(da_peaks_all)
  for (i in 1:len)
  {
    if(length(da_peaks_all[[i]]) > 0)
    {
      enriched.motifs = FindMotifs(atac,
                                   features = da_peaks_all[[i]])
      enriched.motifs = enriched.motifs[order(-enriched.motifs[,6]),]
      enriched.motifs.extract = enriched.motifs[enriched.motifs$fold.enrichment >= 1,]
      filename = paste(types[i], '_enriched_motifs_extracted_new.csv', sep = '')
      write.csv(enriched.motifs.extract, file = filename)
    }
  }
}
motif_enrichment_analysis(atac_Tissue, da_peaks_all)

atac_Tissue = RunChromVAR(object = atac_Tissue,
                         genome = BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(atac_Tissue) = 'chromvar'

fun1 = function(x){fread(x)}
filename = list.files(pattern = '_extracted_new.csv$')
Files = lapply(filename, fun1)
exact_types = gsub('_enriched_motifs_extracted_new.csv','', filename)
names(Files) = exact_types

# draw top 6 motif activity feature plot for each cell type
for(i in 1:length(Files))
{
  motif.activity.plot = FeaturePlot(atac_Tissue,
                                    features = head(Files[[i]]$motif),
                                    min.cutoff = 'q30', 
                                    max.cutoff = 'q90',
                                    pt.size = 0.2,
                                    coord.fixed = T, 
                                    order = T, 
                                    cols = c("lightgrey" ,"#DE1F1F"),label = F)
  TosaveName = paste(names(Files)[i],'_motif_activity_plot.png', sep = '')
  ggsave(motif.activity.plot, filename = TosaveName, width = 6, height = 6, dpi = 600)
}
