library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v75)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(future)
plan("multicore", workers = 12)
options(future.globals.maxSize = 50000 * 1024^2)

set.seed(1234)
setwd('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/')
project_name <- c('Young_1','Old_resistant_1','Old_sensitive_1','Old_resistant_2','Old_sensitive_2','Young_2')
pbmc_list <- list()
peak_list <- list()
meta_list <- list()
frag_list <- list()
gr_list <- list()

for(j in 1:6){
  project_name_j <- tolower(project_name[j])
  mex_dir_path <- paste("/Users/weizhao/Documents/Hep_Seq/ATAC-seq/",project_name_j, "_filtered_peak_bc_matrix/",sep = '')
  mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
  feature_path <- paste(mex_dir_path, "peaks.bed", sep = '/')
  barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')
  features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
  barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

  # mtx <- Matrix::readMM(mtx_path) %>%
  #   magrittr::set_rownames(features$feature) %>%
  #   magrittr::set_colnames(barcodes$barcode)

  peak <- read.table(
    file = feature_path,
    col.names = c("chr", "start", "end")
  )

  metadata <- read.csv(
    file = paste(project_name_j, "_filtered_peak_bc_matrix/",project_name_j, "_singlecell.csv",sep=''),
    header = TRUE,
    row.names = 1
  )[-1,]
  metadata <- metadata[metadata$passed_filters > 500, ]

  frag <- CreateFragmentObject(
    path = paste(project_name_j, "_filtered_peak_bc_matrix/",project_name_j,"_fragments.tsv.gz",sep = ''),
    cells = rownames(metadata)
  )

  peak_list[[j]] <- peak
  meta_list[[j]] <- metadata
  frag_list[[j]] <- frag
  gr_list[[j]] <- makeGRangesFromDataFrame(peak)
  # chrom_assay <- CreateChromatinAssay(
  #   counts = mtx,
  #   sep = c("_", "_"),
  #   genome = 'mm10',
  #   fragments =paste(project_name_j, "_filtered_peak_bc_matrix/",project_name_j,"_fragments.tsv.gz",sep = ''),
  #   min.cells = 10,
  #   min.features = 200
  # )
  # pbmc_list[[j]] <- CreateSeuratObject(
  # counts = chrom_assay,
  # assay = "peaks",
  # meta.data = metadata,project = project_name[j],
  # )
}

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr_list[[1]],gr_list[[2]],gr_list[[3]],gr_list[[4]],gr_list[[5]],gr_list[[6]]))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
# combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

pbmc_list <- list()
for(j in 1:6){
  count <- FeatureMatrix(
    fragments = frag_list[[j]],
    features = combined.peaks,
    cells = rownames(meta_list[[j]])
  )
  chrom_assay <- CreateChromatinAssay(counts = count, sep = c("-", "-"), fragments = frag_list[[j]])
  pbmc_list[[j]] <- CreateSeuratObject(chrom_assay, assay = "ATAC", meta.data=meta_list[[j]])
  pbmc_list[[j]]$dataset <- project_name[[j]]
}

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = pbmc_list[[1]],
  y = list(pbmc_list[[2]], pbmc_list[[3]], pbmc_list[[4]],pbmc_list[[5]],pbmc_list[[6]]),
  add.cell.ids = project_name
)
combined[["ATAC"]]

saveRDS(combined,file = 'ATAC_combined_original.rds')
# combined <- readRDS(file = 'ATAC_combined.rds')

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v75)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
# add the gene information to the object
Annotation(combined) <- annotations

# hhd <- subset(x = combined,subset= (dataset %in% 'Old_resistant_2'))
# compute nucleosome signal score per cell
combined <- NucleosomeSignal(object = combined)
# compute TSS enrichment score per cell
combined <- TSSEnrichment(object = combined, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(combined, group.by = 'high.tss') + NoLegend()

combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = combined, group.by = 'nucleosome_group',region = "chr1-1-20000000")
VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

# saveRDS(combined,file = 'ATAC_combined.rds')
combined <- readRDS(file = 'ATAC_combined.rds')
# combined$nucleosome_signal[combined$dataset=='Old_resistant_2'] <- -1
# combined$TSS.enrichment [combined$dataset=='Old_resistant_2'] <- 1e5

combined$percent.mt <- combined$mitochondrial/combined$total*100

set.seed(1234)
FeatureScatter(combined,'peak_region_fragments','nCount_ATAC')
RidgePlot(combined, 'nFeature_ATAC') + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10,sigma=1))
plot(density(log10(1+combined$peak_region_fragments)))
# subset cells
combined <- subset(
  x = combined,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 10000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.04 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
    & is__cell_barcode == 1
    #& percent.mt < 40
)
# subset peaks
# combined <- combined[rowSums(combined@assays$ATAC@counts>0) >= 20,]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
DimPlot(combined, split.by = 'dataset', pt.size = 0.1)

# ## find integration anchors
# integration.anchors <- FindIntegrationAnchors(
#   object.list = SplitObject(combined, split.by = "dataset"),
#   anchor.features = rownames(combined),
#   reduction = "rlsi",
#   dims = 2:30
# )
#
# # integrate LSI embeddings
# integrated <- IntegrateEmbeddings(
#   anchorset = integration.anchors,
#   reductions = combined[["lsi"]],
#   new.reduction.name = "integrated_lsi",
#   dims.to.integrate = 1:30
# )
#
# ndim = 30 # ndim = 18 removes the mini cluster 10
# integrated <- RunUMAP(integrated, dims = 2:ndim, reduction = 'integrated_lsi', min.dist = 0.05, spread=2, n.neighbors=30)
#
# DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)
# DimPlot(integrated, split.by = 'dataset', pt.size = 0.1)
#
# integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:ndim)
# integrated <- Seurat::FindClusters(object = integrated, algorithm = 1, resolution = 0.3, group.singletons=T, n.iter = 50)
# # integrated <- Seurat::FindSubCluster(object = hm.integrated, cluster = 1, algorithm = 1, resolution = 0.18,graph.name = 'ATAC_snn')
# p1 <- DimPlot(object = integrated, label = TRUE, pt.size = 0.5, label.size = 8) + NoLegend(); p1

# https://www.10xgenomics.com/resources/analysis-guides/batch-effect-correction-in-chromium-single-cell-atac-data
library(harmony)
hm.integrated <- RunHarmony(object = combined, group.by.vars = 'dataset', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE, theta= 2, lambda = 1, sigma= 0.3)
ElbowPlot(hm.integrated,ndims = 50,reduction = 'harmony')
ndim = 20 # ndim = 18 removes the mini cluster 10
hm.integrated <- RunUMAP(hm.integrated, dims = 2:ndim, reduction = 'harmony', min.dist = 0.05, spread=2, n.neighbors=30)
# hm.integrated <- RunTSNE(hm.integrated, dims = 2:ndim, reduction = 'harmony')

DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1)
DimPlot(hm.integrated, split.by = 'dataset', pt.size = 0.1)

hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'harmony', dims = 2:ndim)
hm.integrated <- Seurat::FindClusters(object = hm.integrated, algorithm = 1, resolution = 0.3, group.singletons=T, n.iter = 50)
hm.integrated <- Seurat::FindSubCluster(object = hm.integrated, cluster = 1, algorithm = 1, resolution = 0.18,graph.name = 'ATAC_snn')
p1 <- DimPlot(object = hm.integrated, label = TRUE, pt.size = 0.5, label.size = 8) + NoLegend(); p1
DimPlot(hm.integrated,group.by = 'sub.cluster')
# p2 <- DimPlot(object = hm.integrated, label = TRUE, reduction = 'tsne') + NoLegend()
ggsave('umap.pdf',plot =  p1, width=25,height=25,units='cm')

p2 <- VlnPlot(object = hm.integrated,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal','percent.mt'),
  pt.size = 0.1,
  ncol = 3
)

ggsave('QC_by_cluster.pdf',plot =  p2, width=50,height=25,units='cm')


DimPlot(hm.integrated, split.by = 'dataset', pt.size = 0.1)

DefaultAssay(hm.integrated) <- 'ATAC'
da_peaks <- FindAllMarkers(  object = hm.integrated, min.pct = 0.1, only.pos = T, logfc.threshold = 0.25)# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
da_peaks %>%
  group_by(cluster) %>%
  slice_max(n = 1e5, order_by = avg_log2FC) -> da_peaks
head(da_peaks)
ClosestFeature_da_peaks <- ClosestFeature(hm.integrated, regions = (da_peaks$gene))
head(ClosestFeature_da_peaks)
da_peaks$closest_gene <- ClosestFeature_da_peaks$gene_name[match(da_peaks$gene,ClosestFeature_da_peaks$query_region)]
write.csv(da_peaks,file='full_markers_ATAC_peak_closest_gene.csv')

da_peaks %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC) -> top10_peak
print(top10_peak,n=50*length(unique(hm.integrated@active.ident)))

DefaultAssay(hm.integrated) <- 'ATAC'
gene.activities <- GeneActivity(hm.integrated,extend.upstream = 0)
# gene.activities <- GeneActivity(hm.integrated,extend.upstream = 3000)

hm.integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
hm.integrated <- NormalizeData(
  object = hm.integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm.integrated$nCount_RNA)
)

hm.integrated <- readRDS('integrated_ATAC.rds')

DefaultAssay(hm.integrated) <- 'RNA'
hm.integrated <- Seurat::ScaleData(hm.integrated)
hm.integrated.markers <- FindAllMarkers(hm.integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
hm.integrated.markers <- hm.integrated.markers[grep('^mt',hm.integrated.markers$gene,invert = T),]
hm.integrated.markers <- hm.integrated.markers[grep('^Rp[sl]',hm.integrated.markers$gene,invert = T),]
hm.integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10
print(top10,n=10*length(unique(hm.integrated@active.ident)))

hm.integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1e5, order_by = avg_log2FC) -> hm.integrated.markers
write.csv(hm.integrated.markers,file='full_markers_ATAC.csv')

# saveRDS(hm.integrated,file = 'integrated_ATAC.rds')
pheatmap <- DoHeatmap(hm.integrated, top10$gene)
ggsave('markers_heatmap.pdf',plot =  pheatmap, width=25,height=25,units='cm')

FeaturePlot(
  object = hm.integrated,
  features = c('Cyp2e1', 'Cyp1a2','Cyp2c29','Slc1a2','Cyp2f2', 'Pck1','Glul','Hamp','Hamp2','Igfbp2','Pon1','Wfdc17', 'Clec4f','Ptprb', 'Dcn','Ptprc','Alb','Serpina3k','Cps1'),
  pt.size = 0.1,
  max.cutoff = 'q5',
  ncol = 4,order = F
)

VlnPlot(
  object = hm.integrated,
  features = c('Cyp2e1', 'Cyp1a2','Cyp2c29','Slc1a2','Cyp2f2', 'Pck1','Glul','Hamp','Hamp2','Igfbp2','Pon1','Wfdc17', 'Clec4f','Ptprb', 'Dcn','Ptprc','Alb','Serpina3k','Cps1'),
  pt.size = 0,
  ncol = 4
)

data_plot <- Matrix::t(hm.integrated@assays$RNA@counts[c('Glul','Gulo','Oat','Cyp2e1','Cyp2f2','Cyp1a2','Cyp2c29','Pigr','Pon1','Rgn','Arg1','Igfbp2','Sds',
                                                           'Fbp1','Pck1','Cps1','Alb','Serpina3k','Hamp','Hamp2'),])
data_plot <- as.data.frame(data_plot)
data_plot_subset <- data_plot[,c('Cyp2e1','Cyp2f2')]
data_plot_subset <- cbind(data_plot_subset,hm.integrated$seurat_clusters)
names(data_plot_subset)[3] <- 'cluster'
aggregate_by_group <-  stats::aggregate(data_plot_subset[,1:2],by=list(data_plot_subset$cluster),FUN=function(x){mean(x)})
rownames(aggregate_by_group) <- aggregate_by_group[,1]
aggregate_by_group <- aggregate_by_group[,-1]
aggregate_by_group <- as.data.frame(t(t(aggregate_by_group)/apply(aggregate_by_group,2,max)))
aggregate_by_group$diff <- (aggregate_by_group[,1]-aggregate_by_group[,2])/apply(aggregate_by_group,1,sum)
aggregate_by_group$cluster <- rownames(aggregate_by_group)
aggregate_by_group <- aggregate_by_group[(rownames(aggregate_by_group) %in% c(0,1,2,3,5)),]
ggplot(data=aggregate_by_group, aes(x=reorder(cluster,diff), y=diff)) +  geom_bar(stat="identity")


data_plot <- as.data.frame(data_plot)
x <- 'Cps1'
y <- 'Cyp2e1'
ggplot(data_plot, aes_string(x=x, y=y) ) +
  geom_bin2d() +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

## statistics & table
ident_cluster <- data.frame(ident = hm.integrated@meta.data$dataset,cluster=hm.integrated@meta.data$seurat_clusters,value=1,row.names = colnames(hm.integrated))
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(ident_cluster$ident))
ident_cluster$condition <- tmp[ident_cluster$ident]
ident_cluster$condition <- factor(ident_cluster$condition, levels = (c('Young','Old sensitive','Old resistant')))

ident_cluster_new <- ident_cluster %>%
  group_by(condition, cluster) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))
ident_cluster_new2 <- ident_cluster %>%
  group_by(ident, cluster) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

reshape2::dcast( ident_cluster_new, condition ~ cluster, value.var = 'perc')
library(ggplot2)
pdistribution <- ggplot(ident_cluster_new, aes(fill=cluster, y=perc*100, x=condition)) +
  geom_bar(stat="identity")+labs(x = "Condition", y = "Percent")+
  geom_text(aes(label = paste0(round(perc,4)*100,"%")),
            position = position_stack(vjust = 0.5), size = 4)
ggsave('props.pdf',plot =  pdistribution, width=30,height=55,units='cm')


# integration with scRNA-seq
hm.integrated <- readRDS('integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- "RNA"
# hm.integrated <- NormalizeData(hm.integrated)
# hm.integrated <- ScaleData(hm.integrated,  features = rownames(hm.integrated))
# hm.integrated <- FindVariableFeatures(hm.integrated,nfeatures = 3000)

# seurat_sp <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
immune.combined <- readRDS('/Users/weizhao/Documents/Hep_Seq/Hep_integration_NatMeta_SCT_low2000.rds')
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- subset(immune.combined,subset =  percent.mt > 0.1)
# Run the standard workflow for visualization and clustering
DefaultAssay(immune.combined) <- "integrated"
npcs = 30;
immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 2)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs,k.param=10)
immune.combined <- FindClusters(immune.combined, resolution = 0.5) # common results 10 clusters
DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 10,dims = c(1,2))
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- NormalizeData(immune.combined)
immune.combined <- ScaleData(immune.combined,  features = rownames(immune.combined))
immune.combined <- FindVariableFeatures(immune.combined,nfeatures = 3000)
seurat_sc <- immune.combined

anchors <- FindTransferAnchors(reference = seurat_sc , query = hm.integrated, reference.assay = 'RNA', query.assay = 'RNA',reduction = 'cca')
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat_sc$seurat_clusters, prediction.assay = TRUE,
                                    weight.reduction = hm.integrated[['harmony']],dims = 1:20) ## 'cca',hm.integrated[['lsi']], or hm.integrated[['harmony']]

sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
sp_to_sc$sc_cluster <- paste('sc-',sp_to_sc$sc_cluster,sep = '')
sp_to_sc$sp_cluster <- paste('ATAC-',hm.integrated$seurat_clusters[rownames(sp_to_sc)],sep='')
sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))
# sp_to_sc$condition <- sub('-.*-.*','',rownames(sp_to_sc))
# sp_to_sc$condition <- sub('_[ABC][12]','',sp_to_sc$condition)

df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
# df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  paste('sc-',c('0','1','4','6','9','2','2_1','3','5','8','7','10'),sep=''))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))
df_long$Total <- table(hm.integrated$seurat_clusters[rownames(sp_to_sc)])[sub('ATAC-','',df_long$sp_cluster)] + 1e-6
df_long$ratio <- df_long$Freq/df_long$Total
df_long <- df_long[df_long$ratio>0.1,]
library(ggalluvial)
pmapping <- ggplot(data = df_long,
                   aes(axis1 = sp_cluster, axis2 = sc_cluster, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/6) +
  geom_stratum(width = 1/6) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void()
pmapping
# ggsave('sc_ATAC_mapping.pdf',plot =  pmapping, width=25,height=30,units='cm')

############################################
# integration with scRNA-seq,  dataset-wsie
project_name <- c('Young_1','Old_resistant_1','Old_sensitive_1','Old_resistant_2','Old_sensitive_2','Young_2')
hm.integrated <- readRDS('integrated_ATAC.rds');
immune.combined <- readRDS('/Users/weizhao/Documents/Hep_Seq/Hep_integration_NatMeta_SCT_low2000.rds')
# DefaultAssay(immune.combined) <- "RNA"
immune.combined <- subset(immune.combined,subset =  percent.mt >= 0.1)
# Run the standard workflow for visualization and clustering
DefaultAssay(immune.combined) <- "integrated"
npcs = 30;
immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 2)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs, k.param=10)
immune.combined <- FindClusters(immune.combined, resolution = 0.5) # common results 10 clusters
DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 10,dims = c(1,2))

pmapping_list <- list()
df_long_list <- list()
mapping_table_list <- list()
for(j in 1:length(project_name)) {
  dataset_j <- project_name[j]
  sc_j <- subset(immune.combined, subset = orig.ident == dataset_j)
  DefaultAssay(sc_j) <- "RNA"
  sc_j <- NormalizeData(sc_j)
  sc_j <- ScaleData(sc_j,  features = rownames(sc_j))
  sc_j <- FindVariableFeatures(sc_j,nfeatures = 3000)

  atac_j <- subset(hm.integrated, subset = dataset == dataset_j)
  DefaultAssay(atac_j) <- "RNA"
  atac_j <- NormalizeData(atac_j)
  atac_j <- ScaleData(atac_j,  features = rownames(atac_j))
  atac_j <- FindVariableFeatures(atac_j,nfeatures = 3000)

  anchors <- FindTransferAnchors(reference = sc_j , query = atac_j, reference.assay = 'RNA', query.assay = 'RNA',normalization.method = 'LogNormalize',reduction = 'cca')
  predictions.assay <- TransferData(anchorset = anchors, refdata = sc_j$seurat_clusters, prediction.assay = TRUE,
                                    weight.reduction = 'cca',dims = 1:20)
                                    # weight.reduction = atac_j[['lsi']],dims = 2:50)

  sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
  sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
  sp_to_sc$sc_cluster <- paste('sc-',sp_to_sc$sc_cluster,sep = '')
  sp_to_sc$sp_cluster <- paste('ATAC-',atac_j$seurat_clusters[rownames(sp_to_sc)],sep='')
  sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))

  mapping_table_list[[j]] <- table(sp_to_sc[,c(1,2)])

  df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
  # df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  paste('sc-',c('0','1','4','6','9','2','2_1','3','5','8','7','10'),sep=''))
  df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))
  df_long$Total <- table(atac_j$seurat_clusters[rownames(sp_to_sc)])[sub('ATAC-','',df_long$sp_cluster)] + 1e-6
  df_long$ratio <- df_long$Freq/df_long$Total

  df_long_list[[j]] <- df_long

  df_long <- df_long[df_long$ratio>0.05,]
  library(ggalluvial)
  pmapping <- ggplot(data = df_long,
                     aes(axis1 = sp_cluster, axis2 = sc_cluster, y= Freq)) +
    geom_alluvium(aes(fill = sp_cluster),width = 1/6) +
    geom_stratum(width = 1/6) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum))) +
    theme_void()
  pmapping_list[[j]] <- pmapping
}
for(j in 1:6){
  pmapping_list[[j]] <- pmapping_list[[j]] + ggtitle(project_name[j])
}
p_sc <- DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 10,dims = c(1,2))
pmapping_list_plus_sc <- list()
for(j in 1:6){
  pmapping_list_plus_sc[[j]] <- pmapping_list[[j]] + p_sc
}
mapping_table_sum <- matrix(0,nrow = length(unique(immune.combined$seurat_clusters)),
                            ncol=length(unique(hm.integrated$seurat_clusters)),
                            dimnames = list(paste('sc-',c(0,1,2,5,8,4,3,6,11,9,7,10,12,13), sep=''),paste('ATAC-',c(0,1,2,3,4,5,6,7), sep=''))
                            )
mapping_table_0 <- mapping_table_sum
for(j in 1:6){
  mapping_table_j <- mapping_table_0
  mapping_table_j[rownames(mapping_table_list[[j]]),colnames(mapping_table_list[[j]])] <- mapping_table_list[[j]]
  mapping_table_sum <- mapping_table_sum + mapping_table_j
}

library(circlize)
col_fun = colorRamp2(c(0,  1), c("blue","red"))
ComplexHeatmap::Heatmap((t(mapping_table_sum)/colSums(mapping_table_sum)),col = col_fun(seq(0, 1,by=0.1)),
                         cluster_rows = F, cluster_columns =F,
                        row_names_side = 'left', column_names_rot = 45,
                        column_names_side = 'top')
df_long <- reshape2::melt(mapping_table_sum); names(df_long) <- c('sc_cluster','sp_cluster','Freq')

# integration between ATAC and spatial
hm.integrated <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
DefaultAssay(hm.integrated) <- "RNA"
hm.integrated <- NormalizeData(hm.integrated)
hm.integrated <- ScaleData(hm.integrated,  features = rownames(hm.integrated))
hm.integrated <- FindVariableFeatures(hm.integrated,nfeatures = 3000)

immune.combined <- readRDS('integrated_ATAC.rds')
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- NormalizeData(immune.combined)
immune.combined <- ScaleData(immune.combined,  features = rownames(immune.combined))
immune.combined <- FindVariableFeatures(immune.combined,nfeatures = 3000)
seurat_sc <- immune.combined

anchors <- FindTransferAnchors(reference = seurat_sc , query = hm.integrated, reference.assay = 'RNA', query.assay = 'RNA',normalization.method = 'LogNormalize',reduction = 'cca')
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat_sc$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = 'cca',dims = 1:20)

sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
sp_to_sc$sc_cluster <- paste('ATAC-',sp_to_sc$sc_cluster,sep = '')
sp_to_sc$sp_cluster <- paste('sp-',hm.integrated$seurat_clusters[rownames(sp_to_sc)],sep='')
sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('sp-',c(0,1,2,3,4,5,6,7),sep=''))
# sp_to_sc$condition <- sub('-.*-.*','',rownames(sp_to_sc))
# sp_to_sc$condition <- sub('_[ABC][12]','',sp_to_sc$condition)

df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
# df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  paste('sc-',c('0','1','4','6','9','2','2_1','3','5','8','7','10'),sep=''))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('sp-',c(0,1,2,3,4,5,6,7),sep=''))
df_long$Total <- table(hm.integrated$seurat_clusters[rownames(sp_to_sc)])[sub('sp-','',df_long$sp_cluster)]
df_long$ratio <- df_long$Freq/df_long$Total
df_long <- df_long[df_long$ratio>0.1,]
library(ggalluvial)
pmapping <- ggplot(data = df_long,
                   aes(axis1 = sp_cluster, axis2 = sc_cluster, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/6) +
  geom_stratum(width = 1/6) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void()
pmapping

# integration between ATAC and spatial, condition-wise
project_name <- c('Young_1','Old_resistant_1','Old_sensitive_1','Old_resistant_2','Old_sensitive_2','Young_2')
hm.integrated <- readRDS('integrated_ATAC.rds');
seurat_sp <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
sample_name <- c('A1','B2','A2','C1','B1','C2')
exp_condition <- structure(rep(c('Young','Old_sensitive','Old_resistant'),each=2),names=sample_name)

pmapping_list <- list()
df_long_list <- list()
mapping_table_list <- list()
for(j in 1:length(unique(exp_condition))) {
  sp_j <- subset(seurat_sp, subset = orig.ident %in% names(which(exp_condition == unique(exp_condition)[j])))
  DefaultAssay(sp_j) <- "RNA"
  sp_j <- NormalizeData(sp_j)
  sp_j <- ScaleData(sp_j,  features = rownames(sp_j))
  sp_j <- FindVariableFeatures(sp_j,nfeatures = 3000)

  dataset_j <- paste(unique(exp_condition)[j],c(1,2),sep='_')
  atac_j <- subset(hm.integrated, subset = dataset == dataset_j)
  DefaultAssay(atac_j) <- "RNA"
  atac_j <- NormalizeData(atac_j)
  atac_j <- ScaleData(atac_j,  features = rownames(atac_j))
  atac_j <- FindVariableFeatures(atac_j,nfeatures = 3000)

  anchors <- FindTransferAnchors(reference =  sp_j, query = atac_j, reference.assay = 'RNA', query.assay = 'RNA',normalization.method = 'LogNormalize',reduction = 'cca')
  predictions.assay <- TransferData(anchorset = anchors, refdata = sp_j$seurat_clusters, prediction.assay = TRUE,
                                    weight.reduction = 'cca',dims = 1:10)
  # weight.reduction = atac_j[['lsi']],dims = 2:50)

  sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
  sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
  sp_to_sc$sc_cluster <- paste('sp-',sp_to_sc$sc_cluster,sep = '')
  sp_to_sc$sp_cluster <- paste('ATAC-',atac_j$seurat_clusters[rownames(sp_to_sc)],sep='')
  sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))

  mapping_table_list[[j]] <- table(sp_to_sc[,c(1,2)])

  df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
  # df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  paste('sc-',c('0','1','4','6','9','2','2_1','3','5','8','7','10'),sep=''))
  df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('ATAC-',c(0,1,2,5,3,4,6,7),sep=''))
  df_long$Total <- table(atac_j$seurat_clusters[rownames(sp_to_sc)])[sub('ATAC-','',df_long$sp_cluster)] + 1e-6
  df_long$ratio <- df_long$Freq/df_long$Total

  df_long_list[[j]] <- df_long

  df_long <- df_long[df_long$ratio>0.05,]
  library(ggalluvial)
  pmapping <- ggplot(data = df_long,
                     aes(axis1 = sp_cluster, axis2 = sc_cluster, y= Freq)) +
    geom_alluvium(aes(fill = sp_cluster),width = 1/6) +
    geom_stratum(width = 1/6) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum))) +
    theme_void()
  pmapping_list[[j]] <- pmapping
}
for(j in 1:3){
  pmapping_list[[j]] <- pmapping_list[[j]] + ggtitle(project_name[j])
}
p_sc <- DimPlot(hm.integrated, reduction = "umap",pt.size = 0.5,label=T, label.size = 10,dims = c(1,2))
pmapping_list_plus_sc <- list()
for(j in 1:3){
  pmapping_list_plus_sc[[j]] <- pmapping_list[[j]] + p_sc
}
mapping_table_sum <- matrix(0,ncol = length(unique(seurat_sp$seurat_clusters)),
                            nrow=length(unique(hm.integrated$seurat_clusters)),
                            dimnames = list(paste('ATAC-',c(0,1,2,3,4,5,6,7), sep=''), paste('sp-',c(0,1,2,5,3,4,6,7), sep=''))
)
mapping_table_0 <- mapping_table_sum
for(j in 1:3){
  mapping_table_j <- mapping_table_0
  mapping_table_j[rownames(mapping_table_list[[j]]),colnames(mapping_table_list[[j]])] <- mapping_table_list[[j]]
  mapping_table_sum <- mapping_table_sum + mapping_table_j
}

library(circlize)
col_fun = colorRamp2(c(0,  1), c("white","blue"))
ComplexHeatmap::Heatmap((t(mapping_table_sum)/colSums(mapping_table_sum)),col = col_fun(seq(0, 1,by=0.1)),
                        cluster_rows = F, cluster_columns =F,
                        row_names_side = 'left', column_names_rot = 45, name= 'proportion',
                        column_names_side = 'top')
########## SCT results not good
# mapping between ATAC-derived activity & RNA-seq using SCT
hm.integrated <- readRDS('integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- "RNA"
hm.integrated.list <- SplitObject(hm.integrated, split.by = "dataset")
hm.integrated.list <- lapply(X = hm.integrated.list, FUN = function(x) {
  x <- SCTransform(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = hm.integrated.list,nfeatures = 3000)
hm.integrated.list <- PrepSCTIntegration(object.list = hm.integrated.list, anchor.features = features)
hm.integrated.anchors <- FindIntegrationAnchors(object.list = hm.integrated.list, normalization.method = "SCT",anchor.features = features)
hm.integrated.activity <- IntegrateData(anchorset = hm.integrated.anchors, normalization.method = "SCT")
# saveRDS(hm.integrated.activity,file='integrated_ATAC_activity.rds')

seurat_sp <- readRDS('integrated_ATAC_activity.rds')
seurat_sp <- RunPCA(seurat_sp, npcs = 30, verbose = FALSE)

immune.combined <- readRDS('/Users/weizhao/Documents/Hep_Seq/Hep_integration_NatMeta_SCT_low2000.rds')
immune.combined <- subset(immune.combined,subset =  percent.mt >= 0.1)
DefaultAssay(immune.combined) <- "integrated"
npcs = 30;
immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 2)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs,k.param=10)
immune.combined <- FindClusters(immune.combined, resolution = 0.5) # common results 10 clusters
DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 10,dims = c(1,2))
seurat_sc <- immune.combined

anchors <- FindTransferAnchors(reference = seurat_sc , query = seurat_sp, reference.assay = 'integrated', query.assay = 'integrated',normalization.method = 'SCT')
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat_sc$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = seurat_sp[["pca"]], dims = 1:20,
                                  sd.weight = 1)
seurat_sp[["predictions"]] <- predictions.assay

## Alluvial plot
sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
sp_to_sc$sc_cluster <- paste('sc-',sp_to_sc$sc_cluster,sep = '')
sp_to_sc$sp_cluster <- paste('ATAC-',seurat_sp$seurat_clusters[rownames(sp_to_sc)],sep='')
sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))

df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))
df_long$Total <- table(seurat_sp$seurat_clusters[rownames(sp_to_sc)])[sub('ATAC-','',df_long$sp_cluster)]
df_long$ratio <- df_long$Freq/df_long$Total
df_long <- df_long[df_long$ratio>0.1,]
library(ggalluvial)
ggplot(data = df_long,
       aes(axis1 = sp_cluster, axis2 = sc_cluster, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/6) +
  geom_stratum(width = 1/6) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void()


# mapping between ATAC-derived activity & sp data using SCT
seurat_sp <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
seurat_sc <- readRDS('integrated_ATAC_activity.rds')

anchors <- FindTransferAnchors(reference = seurat_sc , query = seurat_sp, reference.assay = 'integrated', query.assay = 'integrated',normalization.method = 'SCT')
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat_sc$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = seurat_sp[["pca"]], dims = 1:10,
                                  sd.weight = 1, n.trees=50,k.weight = 50)
seurat_sp[["predictions"]] <- predictions.assay

## Alluvial plot
sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
sp_to_sc$sc_cluster <- paste('ATAC-',sp_to_sc$sc_cluster,sep = '')
sp_to_sc$sp_cluster <- paste('sp-',seurat_sp$seurat_clusters[rownames(sp_to_sc)],sep='')
sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('sp-',c(0,1,2,5,3,4,6,7),sep=''))

df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('sp-',c(0,1,2,5,3,4,6,7),sep=''))
df_long$Total <- table(seurat_sp$seurat_clusters[rownames(sp_to_sc)])[sub('sp-','',df_long$sp_cluster)]
df_long$ratio <- df_long$Freq/df_long$Total
df_long <- df_long[df_long$ratio>0.1,]
library(ggalluvial)
ggplot(data = df_long,
       aes(axis1 = sp_cluster, axis2 = sc_cluster, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/6) +
  geom_stratum(width = 1/6) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void()


### motif enrichment
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)
hm.integrated <- readRDS('integrated_ATAC.rds')
p1 <- DimPlot(hm.integrated, label = TRUE, pt.size = 0.1) + NoLegend(); p1

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
DefaultAssay(hm.integrated) <- 'ATAC'
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(hm.integrated)) %in% main.chroms)
hm.integrated <- hm.integrated[keep.peaks, ]

hm.integrated <- AddMotifs(
  object = hm.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

# da_peaks <- FindAllMarkers(
#   object = hm.integrated,
#   only.pos = TRUE,
#   test.use = 'LR',
#   min.pct = 0.05,
#   #latent.vars = 'nCount_ATAC'
# )
da_peaks <- FindAllMarkers(
  object = hm.integrated,
  only.pos = TRUE,
  # test.use = 'LR',
  min.pct = 0.05,
  #latent.vars = 'nCount_ATAC'
)
# get top differentially accessible peaks
# da_peaks <- da_peaks[da_peaks$p_val < 0.005, ]
top.da.peak <- da_peaks$gene
enriched.motifs <- FindMotifs(
  object = hm.integrated,
  features = top.da.peak
)
enriched.motifs.list <- lapply(0:7,FUN = function(x){
  enriched.motifs <- FindMotifs(
    object = hm.integrated,
    features = top.da.peak[da_peaks$cluster==x]
  )
  enriched.motifs %>%
    slice_max(n = 1e5, order_by = -pvalue) -> enriched.motifs
})
names(enriched.motifs.list) <- paste('ATAC_',0:7,sep='')
for(j in names(enriched.motifs.list)){
  write.csv(enriched.motifs.list[[j]],file = paste('enriched.motifs_',j,'.csv'))
}


enriched.motifs.list2 <- lapply(0:7,FUN = function(x){
  object_tmp <- subset(hm.integrated,subset= seurat_clusters == x)
  enriched.motifs <- FindMotifs(
    object = object_tmp,
    features = names(which(rowSums(object_tmp) > quantile(rowSums(object_tmp),0.99)))
  )
  enriched.motifs %>%
    slice_max(n = 1e5, order_by = -pvalue) -> enriched.motifs
})
names(enriched.motifs.list2) <- paste('ATAC_',0:7,sep='')
for(j in names(enriched.motifs.list2)){
  write.csv(enriched.motifs.list2[[j]],file = paste('enriched.motifs_abundance_',j,'.csv'))
}

## DE peak across conditions, and then motif enrichment; correlation with RNA-seq
# hm.integrated <- readRDS('integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- 'ATAC'
unique_cluster <- 0:7
for(j in 1:length(unique_cluster)){
  hm.integrated_j <- subset(hm.integrated, subset = seurat_clusters == unique_cluster[j])
  Idents(hm.integrated_j)  <- ident_cluster[Cells(hm.integrated_j),'condition']
  hm.integrated_j@active.ident <- factor(x = hm.integrated_j@active.ident, levels = c('Young','Old sensitive','Old resistant'))

  da_peaks_j <- FindAllMarkers(  object = hm.integrated_j, min.pct = 0.1, only.pos = T, logfc.threshold = 0.25)# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
  da_peaks_j %>%
    group_by(cluster) %>%
    slice_max(n = 1e5, order_by = avg_log2FC) -> da_peaks_j

  ClosestFeature_da_peaks_j <- ClosestFeature(hm.integrated_j, regions = (da_peaks_j$gene))
  da_peaks_j$closest_gene <- ClosestFeature_da_peaks_j$gene_name[match(da_peaks_j$gene,ClosestFeature_da_peaks_j$query_region)]
  write.csv(da_peaks_j,file=paste('ATAC_peak_closest_gene_',unique_cluster[j],'_DA.csv',sep=''))

  enriched.motifs.j <- list()

  for(jj in c('Young','Old sensitive','Old resistant')){
    if(length(da_peaks_j$gene[da_peaks_j$cluster==jj])==0){
      enriched.motifs.jj <- NULL
    } else{
    enriched.motifs.jj <- FindMotifs(
      object = hm.integrated_j,
      features = da_peaks_j$gene[da_peaks_j$cluster==jj]
    )
    enriched.motifs.jj %>%
      slice_max(n = 1e5, order_by = -pvalue) -> enriched.motifs.jj
    enriched.motifs.jj$condition <- jj
    enriched.motifs.j[[jj]] <- enriched.motifs.jj
    }
  }
  enriched.motifs.j.df <- enriched.motifs.j[[1]]
  for(jjj in 2:length(enriched.motifs.j)){
    enriched.motifs.j.df <- rbind(enriched.motifs.j.df,enriched.motifs.j[[jjj]])
  }
  write.csv(enriched.motifs.j.df,file = paste('enriched.motifs_ATAC_',unique_cluster[j],'_DA.csv',sep = ''))
}

### DEG analysis based on activity assay
hm.integrated <- readRDS('integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- 'RNA'
unique_cluster <- 0:7
for(j in 1:length(unique_cluster)){
  hm.integrated_j <- subset(hm.integrated, subset = seurat_clusters == unique_cluster[j])
  Idents(hm.integrated_j)  <- ident_cluster[Cells(hm.integrated_j),'condition']
  hm.integrated_j@active.ident <- factor(x = hm.integrated_j@active.ident, levels = c('Young','Old sensitive','Old resistant'))

  deg_j <- FindAllMarkers(  object = hm.integrated_j, min.pct = 0.1, only.pos = T, logfc.threshold = 0.25)# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
  deg_j %>%
    group_by(cluster) %>%
    slice_max(n = 1e5, order_by = avg_log2FC) -> deg_j
  write.csv(deg_j,file = paste('DEG_',unique_cluster[j],'_.csv',sep = ''))
}

### DEG analysis based on activity assay, 1 vs 1
hm.integrated <- readRDS('integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- 'RNA'

ident_cluster <- data.frame(ident = hm.integrated@meta.data$dataset,cluster=hm.integrated@meta.data$seurat_clusters,value=1,row.names = colnames(hm.integrated))
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(ident_cluster$ident))
ident_cluster$condition <- tmp[ident_cluster$ident]
ident_cluster$condition <- factor(ident_cluster$condition, levels = (c('Young','Old sensitive','Old resistant')))

unique_cluster <- 0:7
for(j in 1:length(unique_cluster)){
  hm.integrated_j <- subset(hm.integrated, subset = seurat_clusters == unique_cluster[j])
  Idents(hm.integrated_j)  <- ident_cluster[Cells(hm.integrated_j),'condition']
  hm.integrated_j@active.ident <- factor(x = hm.integrated_j@active.ident, levels = c('Young','Old sensitive','Old resistant'))

  cellnumber <- table(hm.integrated_j@active.ident)

  if(cellnumber[1] >0 & cellnumber[2] >0){
  deg_j_1 <- FindMarkers(  object = hm.integrated_j, ident.1 ='Old sensitive', ident.2  = 'Young', min.pct = 0.1, only.pos = F, logfc.threshold = 0.25)# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
  # deg_j <- FindMarkers(  object = hm.integrated_j,  ident.1='Young', min.pct = 0.1, only.pos = T, logfc.threshold = 0.25, subset.ident = c('Young','Old sensitive'))# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
  deg_j_1 %>%
    slice_max(n = 1e5, order_by = avg_log2FC) -> deg_j_1
  deg_j_1$comparison <- 'OS_vs_YS'} else {deg_j_1 <- NULL}

  if(cellnumber[2] >0 & cellnumber[3] >0){
  deg_j_2 <- FindMarkers(  object = hm.integrated_j, ident.1 ='Old resistant', ident.2  = 'Old sensitive', min.pct = 0.1, only.pos = F, logfc.threshold = 0.25)# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
  # deg_j <- FindMarkers(  object = hm.integrated_j,  ident.1='Young', min.pct = 0.1, only.pos = T, logfc.threshold = 0.25, subset.ident = c('Young','Old sensitive'))# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
  deg_j_2 %>%
    slice_max(n = 1e5, order_by = avg_log2FC) -> deg_j_2
  deg_j_2$comparison <- 'OR_vs_OS' } else {deg_j_2 <- NULL}

  deg_j <- rbind(deg_j_1,deg_j_2)

  write.csv(deg_j,file = paste('DEG_',unique_cluster[j],'_1vs1.csv',sep = ''))
}


## DE peak across conditions, and then motif enrichment; correlation with RNA-seq
# hm.integrated <- readRDS('integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- 'ATAC'
unique_cluster <- 0:7
for(j in 1:length(unique_cluster)){
  hm.integrated_j <- subset(hm.integrated, subset = seurat_clusters == unique_cluster[j])
  Idents(hm.integrated_j)  <- ident_cluster[Cells(hm.integrated_j),'condition']
  hm.integrated_j@active.ident <- factor(x = hm.integrated_j@active.ident, levels = c('Young','Old sensitive','Old resistant'))

  cellnumber <- table(hm.integrated_j@active.ident)
  enriched.motifs.j.df <- NULL
  pos_neg_name <- c('positive','negative')

  if(cellnumber[1] >0 & cellnumber[2] >0){
    da_peaks_j_1 <- FindMarkers(  object = hm.integrated_j, ident.1 ='Old sensitive', ident.2  = 'Young', min.pct = 0.1, only.pos = F, logfc.threshold = 0.25)# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
    # deg_j <- FindMarkers(  object = hm.integrated_j,  ident.1='Young', min.pct = 0.1, only.pos = T, logfc.threshold = 0.25, subset.ident = c('Young','Old sensitive'))# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
    da_peaks_j_1 %>%
      slice_max(n = 1e5, order_by = avg_log2FC) -> da_peaks_j_1
    ClosestFeature_da_peaks_j_1 <- ClosestFeature(hm.integrated_j, regions = rownames(da_peaks_j_1))
    da_peaks_j_1$closest_gene <- ClosestFeature_da_peaks_j_1$gene_name[match(rownames(da_peaks_j_1),ClosestFeature_da_peaks_j_1$query_region)]
    da_peaks_j_1$comparison <- 'OS_vs_YS'

    idx_pos_neg <- list(which(da_peaks_j_1$avg_log2FC>0),which(da_peaks_j_1$avg_log2FC <0))
    for(jj in 1:2){
      if(length(rownames(da_peaks_j_1)[idx_pos_neg[[jj]]])==0){
        enriched.motifs.jj <- NULL
      } else{
        enriched.motifs.jj <- FindMotifs(
          object = hm.integrated_j,
          features = rownames(da_peaks_j_1)[idx_pos_neg[[jj]]]
        )
        enriched.motifs.jj %>%
          slice_max(n = 1e5, order_by = -pvalue) -> enriched.motifs.jj
        enriched.motifs.jj$category <- pos_neg_name[jj]
        enriched.motifs.jj$comparison <- 'OS_vs_YS'
        enriched.motifs.j.df <- rbind( enriched.motifs.j.df, enriched.motifs.jj)
      }
    }
    } else {da_peaks_j_1 <- NULL}

  if(cellnumber[2] >0 & cellnumber[3] >0){
    da_peaks_j_2 <- FindMarkers(  object = hm.integrated_j, ident.1 ='Old resistant', ident.2  = 'Old sensitive', min.pct = 0.1, only.pos = F, logfc.threshold = 0.25)# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
    # deg_j <- FindMarkers(  object = hm.integrated_j,  ident.1='Young', min.pct = 0.1, only.pos = T, logfc.threshold = 0.25, subset.ident = c('Young','Old sensitive'))# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
    da_peaks_j_2 %>%
      slice_max(n = 1e5, order_by = avg_log2FC) -> da_peaks_j_2
    ClosestFeature_da_peaks_j_2 <- ClosestFeature(hm.integrated_j, regions = rownames(da_peaks_j_2))
    da_peaks_j_2$closest_gene <- ClosestFeature_da_peaks_j_2$gene_name[match(rownames(da_peaks_j_2),ClosestFeature_da_peaks_j_2$query_region)]
    da_peaks_j_2$comparison <- 'OR_vs_OS'

    idx_pos_neg <- list(which(da_peaks_j_2$avg_log2FC>0),which(da_peaks_j_2$avg_log2FC <0))
    for(jj in 1:2){
      if(length(rownames(da_peaks_j_2)[idx_pos_neg[[jj]]])==0){
        enriched.motifs.jj <- NULL
      } else{
        enriched.motifs.jj <- FindMotifs(
          object = hm.integrated_j,
          features = rownames(da_peaks_j_2)[idx_pos_neg[[jj]]]
        )
        enriched.motifs.jj %>%
          slice_max(n = 1e5, order_by = -pvalue) -> enriched.motifs.jj
        enriched.motifs.jj$category <- pos_neg_name[jj]
        enriched.motifs.jj$comparison <- 'OR_vs_OS'
        enriched.motifs.j.df <- rbind( enriched.motifs.j.df, enriched.motifs.jj)
      }
    }
  } else {da_peaks_j_2 <- NULL}

  da_peaks_j <- rbind(da_peaks_j_1,da_peaks_j_2)

  write.csv(da_peaks_j,file=paste('ATAC_peak_closest_gene_',unique_cluster[j],'_DA_1vs1.csv',sep=''))
  write.csv(enriched.motifs.j.df,file = paste('enriched.motifs_ATAC_',unique_cluster[j],'_DA_1vs1.csv',sep = ''))
}

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
hm.integrated <- readRDS('integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- 'ATAC'

set.seed(1234)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
DefaultAssay(hm.integrated) <- 'ATAC'
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(hm.integrated)) %in% main.chroms)
hm.integrated <- hm.integrated[keep.peaks, ]

hm.integrated <- AddMotifs(
  object = hm.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
hm.integrated <- RunChromVAR(
  object = hm.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10
)
DefaultAssay(hm.integrated) <- 'chromvar'

differential.activity <- FindAllMarkers(
  object = hm.integrated,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",logfc.threshold = 0.1,
)
DefaultAssay(hm.integrated) <- 'ATAC'
differential.activity$TF <- ConvertMotifID(hm.integrated,id=differential.activity$gene)
differential.activity %>%
  group_by(cluster) %>%
  slice_max(n = 1e5, order_by = avg_diff) -> top10_TF
# print(top10_TF, n=160)
write.csv(top10_TF,file = 'ChromVAR_TF_all_by_cluster.csv')

# differential analysis across condition
setwd('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/')
dir.create('Diff_chromvar2')
hm.integrated <- readRDS('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/ATAC_integrated_chromvar.rds')
cluster_celltype_atac <- structure(c('zone 3','zone 2','zone 1','zone 1','EC','HSC','KC','Immune cells'),names=0:7)
hm.integrated$annotation <- cluster_celltype_atac[as.character(hm.integrated$seurat_clusters)]
tmp_atac <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(hm.integrated$dataset))
hm.integrated$condition <- factor(tmp_atac[hm.integrated$dataset],levels = c('Young','Old sensitive','Old resistant'))
hm.integrated$annotation <- factor(hm.integrated$annotation, levels=c('zone 1','zone 2','zone 3','EC','HSC','KC','Immune cells','Undefined'))

annotation_oi <- c('zone 1','zone 2','zone 3','EC','HSC','KC','Immune cells')
for(j in 1:length(annotation_oi)){
  hm.integrated.cv <- subset(hm.integrated, subset = annotation %in% annotation_oi[j])
  hm.integrated.cv@active.ident <- hm.integrated.cv$condition
  DefaultAssay(hm.integrated.cv) <- 'chromvar'
  differential.activity.cv <- FindAllMarkers(
  object = hm.integrated.cv,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0.1,
  )
  DefaultAssay(hm.integrated.cv) <- 'ATAC'
  differential.activity.cv$TF <- ConvertMotifID(hm.integrated.cv,id=differential.activity.cv$gene)
  differential.activity.cv %>%
    group_by(cluster) %>%
    slice_max(n = 1e5, order_by = avg_diff) -> top10_TF
  write.csv(top10_TF,file=paste('Diff_TF_fromATAC_',annotation_oi[j],'.csv',sep=''))

  # print(top10_TF, n=160)
  top10_TF %>% dplyr::filter(avg_diff>0.1 & p_val_adj < 0.001) -> top10_TF
  if(dim(top10_TF)[1]>0){
  avg_TF <- AverageExpression(subset(hm.integrated,subset = annotation ==annotation_oi[j] ), features = unique(top10_TF$gene), group.by = 'condition',assays = 'chromvar')
  avg_TF <- as.data.frame(avg_TF$chromvar)
  avg_TF$trend <- ifelse(avg_TF[,1] >= avg_TF[,2] & avg_TF[,2] >= avg_TF[,3],1,ifelse(avg_TF[,1] <= avg_TF[,2] & avg_TF[,2]  <= avg_TF[,3],-1,0))
  top10_TF$trend <- avg_TF$trend[match(top10_TF$gene,rownames(avg_TF))]
  top10_TF %>% dplyr::filter(trend!=0) -> top10_TF

  top10_TF %>%
    group_by(TF) %>%
    slice_max(n = 1, order_by = avg_diff) -> top10_TF
  top10_TF  %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_diff) -> top10_TF
  print(top10_TF, n=10*3)

  if(dim(top10_TF)[1]>0){
    DefaultAssay(hm.integrated.cv) <- 'chromvar'
    p_dot_chromvar <- DotPlot(hm.integrated.cv,features = unique(top10_TF$gene)) + coord_flip() +
      scale_x_discrete(breaks=top10_TF$gene,labels=top10_TF$TF, limits=rev) +
      ylab('condition') +  xlab('TF') + ggtitle(annotation_oi[j]) +
      theme(axis.text.y = element_text(size = 14), axis.text.x=element_text(size = 14,angle = 30, vjust = 0.65))
    ggsave(paste('p_dot_chromvar_',annotation_oi[j],'.pdf',sep=''),plot =  p_dot_chromvar, width=25,height=22,units='cm')
  }
}
}

seurat_sc@active.ident <- seurat_sc$annotation
DotPlot(seurat_sc,features=stringr::str_to_title(top10_TF$TF),group.by = 'condition',idents = 'zone 3') + coord_flip() +
  ylab('condition') +  xlab('TF') + ggtitle('Differential gene expression (scRNA-seq)') +
  theme(axis.text.y = element_text(size = 14), axis.text.x=element_text(size = 14,angle = 30, vjust = 0.65))


DefaultAssay(hm.integrated.cv) <- 'ATAC'
hm.integrated.cv <- subset(hm.integrated.cv, subset= annotation %in% 'zone 3')
cov_plot <- CoveragePlot(
  object = hm.integrated.cv, group.by = 'condition',
  region = c("Onecut1",'Onecut2','Onecut3'),
  annotation = T,
  peaks = FALSE,extend.upstream = 20000, extend.downstream = 20000, show.bulk	= F, window =1000,
) + theme(axis.text = element_text(size = 20)) + patchwork::plot_layout(ncol=1)
cov_plot

top10_TF  %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_diff) -> top10_TF
print(top10_TF, n=160)
MotifPlot(
  object = hm.integrated.cv,
  motifs = unique(top10_TF$gene)
) + patchwork::plot_layout(ncol=4)
