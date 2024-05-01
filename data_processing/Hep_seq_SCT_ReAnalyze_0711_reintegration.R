library(dplyr)
library(Seurat)
library(sctransform)
library(future)
library(ggplot2)
library(SoupX)
# check the current active plan
plan("multicore", workers = 12)
options(future.globals.maxSize= 5000*1024^2)
options(future.seed=TRUE)
set.seed(1234)
setwd('/Users/weizhao/Documents/Hep_Seq/')
pbmc_list <- list()
ind <- c(1,4,7,10,13,16)
filefolder_list <- paste('/Users/weizhao/Documents/Hep_Seq/Hep', ind, '_outs',sep='')
project_name <- c('Young_1','Old_resistant_1','Old_sensitive_1','Old_resistant_2','Old_sensitive_2','Young_2')
setwd(filefolder_list[[6]])

immune.combined <- readRDS('Hep_integration_SoupX_0712.rds')
DefaultAssay(immune.combined) <- "integrated"
FeatureScatter(immune.combined,'nCount_RNA','nFeature_RNA') + geom_function(fun = function(x) 2500*(x-300)/((x-300)+2000), colour = "red")
hill_fun <- function(x) 2500*(x-300)/((x-300)+2000)
immune.combined$nCount_over_nFeature <- immune.combined$nCount_RNA/immune.combined$nFeature_RNA
immune.combined$nCount_over_nFeature1 <- hill_fun(immune.combined$nCount_RNA)/immune.combined$nFeature_RNA
immune.combined <- subset(immune.combined,subset =  nFeature_RNA > 750 & nCount_over_nFeature1 < 1)

# immune.combined <- subset(immune.combined,subset = doublet == 'Singlet')
# immune.combined <- subset(immune.combined,subset=(Alb0 >  0 &  nFeature_RNA > 750) |  (Alb0 <=0 &  nFeature_RNA > 200))
# immune.combined$Alb0 <- immune.combined@assays$RNA@counts['Alb',]
# immune.combined <- subset(immune.combined,subset =  (Alb0 >  0 &  nFeature_RNA > 800) |  (Alb0 <=0 &  nFeature_RNA > 500) )
# Run the standard workflow for visualization and clustering
npcs = 15;
immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
ElbowPlot(immune.combined,ndims = npcs)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 3)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 2,
                           min.dist = 0.3, spread = 0.8, n.neighbors = 20)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs)
immune.combined <- FindClusters(immune.combined, resolution = 0.5) # common results 15 clusters
# immune.combined <- FindClusters(immune.combined, resolution = 0.15) # common results 15 clusters
DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 5,dims = c(1,2),group.by = 'seurat_clusters')

cluster_cellclass <- structure(c(rep('zone 2',4),rep('zone 3',2),'zone 3 Glul+',rep('zone 1',2),'Hep doublet','undefined 1','undefined 2','EC','EC/Hep doublet','KC','HSC/FB'),
                                                              names=as.character(c(0,1,2,6,3,7,10,5,8,4,11,13,9,14,12,15)))
immune.combined$annotation <- cluster_cellclass[as.character(immune.combined$seurat_clusters)]
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(immune.combined$orig.ident))
immune.combined$condition <- tmp[immune.combined$orig.ident]

table(immune.combined$orig.ident,immune.combined$annotation)/rowSums(table(immune.combined$orig.ident,immune.combined$annotation))*100
table(immune.combined$condition,immune.combined$annotation)/rowSums(table(immune.combined$condition,immune.combined$annotation))*100

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- FindSubCluster(object=immune.combined, cluster= 'KC', graph.name = 'integrated_snn', resolution = 0.1,subcluster.name = "sub.cluster.KC")
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE,group.by = 'sub.cluster.KC')
immune.combined$annotation[which(immune.combined@meta.data$sub.cluster.KC %in% c('KC_1'))] <- 'immune cells'
# immune.combined$annotation[which(immune.combined@meta.data$annotation %in% c('zone 3') & immune.combined@assays$RNA@data['Glul',] > 2.5)] <- 'zone 3 Glul+'
# immune.combined@active.ident <- factor(immune.combined$annotation,levels=c('zone 1','zone 3','zone 3 Glul+','zone 2', 'EC','HSC/FB','KC','immune cells','undefined 1','undefined 2','doublet'))
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE,pt.size = 0.8)

table(immune.combined$orig.ident,immune.combined$annotation)/rowSums(table(immune.combined$orig.ident,immune.combined$annotation))*100
table(immune.combined$condition,immune.combined$annotation)/rowSums(table(immune.combined$condition,immune.combined$annotation))*100


### re-do integration after subseting via Hill function
ifnb.list <- SplitObject(immune.combined, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- SCTransform(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",anchor.features = features, reference = c(1,6))
immune.combined.new <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.old <- immune.combined
immune.combined <- immune.combined.new
# Run the standard workflow for visualization and clustering
DefaultAssay(immune.combined) <- "integrated"
npcs = 10;
immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
ElbowPlot(immune.combined,ndims = npcs)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 2,
                           min.dist = 0.3, spread = 0.8, n.neighbors = 30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs)
immune.combined <- FindClusters(immune.combined, resolution = 0.5) # common results 15 clusters
DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 5,dims = c(1,2),group.by = 'seurat_clusters')

table(immune.combined$orig.ident,immune.combined$seurat_clusters)/rowSums(table(immune.combined$orig.ident,immune.combined$seurat_clusters))*100
table(immune.combined$condition,immune.combined$seurat_clusters)/rowSums(table(immune.combined$condition,immune.combined$seurat_clusters))*100

p1 <- DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 5,dims = c(1,2))
DefaultAssay(immune.combined) <- "RNA"
immune.combined.markers <- presto:::wilcoxauc.Seurat(immune.combined, group_by = "annotation", assay = 'data', seurat_assay = "RNA")
immune.combined.markers %>% filter(padj < 0.01, logFC > 0.25, pct_in > 25) %>%
  group_by(group) %>%
  slice_max(n = 1e5, order_by = abs(logFC)) -> immune.combined.markers

immune.combined.markers %>%
  group_by(group) %>%
  slice_max(n = 50, order_by = logFC) -> top10
print(top10,n=50*length(unique(immune.combined@active.ident)))

genes_to_plot <- c('Cyp2f2','Pck1','Uroc1','Pigr','Apoa4','Arg1','C3','G6pc','F2',   # zone 1
                   'Car3','Cyp1a2','Cyp2c29', 'Cyp2e1','Rgn', 'Gulo', 'Cldn2', 'Oat','Glul', # zone 3 'Lect2',  'Axin2',
                   'Hamp','Ttr','Fabp1','Apoa2','Apoc1','Apoe','Apoc3', 'Wfdc21',  # zone 2
                   'Dnase1l3','Ptprb','Clec4g','Igfbp7', # EC
                   'Dcn','Rgs5','Reln',                  # HSC/FB
                   'Ptgs1', 'Clec4f','Wfdc17','Vsig4','Cd74',   # KC
                   'Ptprc', 'Lsp1','Ccr2'
                   #,'Scd1','Eef1a1','Eef2','Klf2'
                   )
immune.combined@active.ident <-immune.combined$seurat_clusters
p_dot_cluster <- DotPlot(immune.combined,features = unique(genes_to_plot)) + coord_flip() +  scale_x_discrete(limits=rev) + ylab('cluster') +
  theme(axis.text.y = element_text(size = 14), axis.text.x=element_text(size = 14,angle = 30, vjust = 0.65))

cluster_cellclass <- structure(c(rep('zone 2',5),rep('zone 3',4),rep('zone 1',1),'EC','KC','HSC/FB'),
                               names=as.character(c(0,1,6,7,8,3,4,5,10,2,9,11,12)))
immune.combined$annotation <- cluster_cellclass[as.character(immune.combined$seurat_clusters)]
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(immune.combined$orig.ident))
immune.combined$condition <- tmp[immune.combined$orig.ident]
immune.combined@active.ident <- factor(immune.combined$annotation,levels=c('zone 1','zone 3','zone 3 Glul+','zone 2', 'EC','EC/Hep doublet','HSC/FB','KC','immune cells','undefined 1','undefined 2','Hep doublet'))

table(immune.combined$orig.ident,immune.combined$annotation)/rowSums(table(immune.combined$orig.ident,immune.combined$annotation))*100
table(immune.combined$condition,immune.combined$annotation)/rowSums(table(immune.combined$condition,immune.combined$annotation))*100

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- FindSubCluster(object=immune.combined, cluster= 'KC', graph.name = 'integrated_snn', resolution = 0.1,subcluster.name = "sub.cluster.KC")
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE,group.by = 'sub.cluster.KC')
immune.combined$annotation[which(immune.combined@meta.data$sub.cluster.KC %in% c('KC_1'))] <- 'immune cells'

immune.combined <- FindSubCluster(object=immune.combined, cluster= 'EC', graph.name = 'integrated_snn', resolution = 0.1,subcluster.name = "sub.cluster.EC")
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE,group.by = 'sub.cluster.EC')
immune.combined$annotation[which(immune.combined@meta.data$sub.cluster.EC %in% c('EC_1'))] <- 'EC/Hep doublet'
immune.combined$annotation[which(immune.combined@meta.data$annotation %in% c('zone 3') & immune.combined@assays$RNA@data['Glul',] > 2.5)] <- 'zone 3 Glul+'
immune.combined@active.ident <- factor(immune.combined$annotation,levels=c('zone 1','zone 3','zone 3 Glul+','zone 2', 'EC','EC/Hep doublet','HSC/FB','KC','immune cells','undefined 1','undefined 2','Hep doublet'))
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE,pt.size = 0.2,label.size = 10)

immune.combined$seurat_clusters1 <- as.character(immune.combined$seurat_clusters)
immune.combined$seurat_clusters1[which(immune.combined@meta.data$sub.cluster.EC %in% c('EC_1'))] <- '9-1'
immune.combined$seurat_clusters1[which(immune.combined@meta.data$seurat_clusters %in% '3' & immune.combined@assays$RNA@data['Glul',] > 2.5)] <- '3-1'
immune.combined$seurat_clusters1[which(immune.combined@meta.data$sub.cluster.KC %in% c('KC_1'))] <- '11-1'




immune.combined <- readRDS('/Users/weizhao/Documents/Hep_Seq/Hep16_outs/Hep_integration_SoupX_reintegration_0805.rds')
cluster_celltype <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                              names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
immune.combined$annotation <- cluster_celltype[as.character(immune.combined$seurat_clusters1)]
immune.combined$annotation <- factor(immune.combined$annotation,levels=c('zone 1','zone 3','zone 3 Glul+','zone 2', 'EC','EC/HSC chimera','HSC/FB','KC','immune cells'))
immune.combined$annotation[which(immune.combined@meta.data$sub.cluster.EC %in% c('EC_1'))] <- 'EC/HSC chimera'
immune.combined$seurat_clusters1 <- factor(immune.combined$seurat_clusters1,levels=c('2','10','3','4','5','3-1','6','1','0','7','8','9','9-1','12','11','11-1'))

dir.create('/Users/weizhao/Documents/Hep_Seq/scRNA_results')
setwd('/Users/weizhao/Documents/Hep_Seq/scRNA_results')
DefaultAssay(immune.combined) <- "RNA"
p_dot <- DotPlot(immune.combined,features = unique(genes_to_plot),group.by = 'annotation') + coord_flip() +  scale_x_discrete(limits=rev) + ylab('annotation') +
  theme(axis.text.y = element_text(size = 14), axis.text.x=element_text(size = 14,angle = 45, vjust = 0.9, hjust = 0.9))
p_dot_cluster <- DotPlot(immune.combined,features = unique(genes_to_plot),group.by = 'seurat_clusters1') + coord_flip() +  scale_x_discrete(limits=rev) + ylab('cluster') +
  theme(axis.text.y = element_text(size = 14), axis.text.x=element_text(size = 14,angle = 30, vjust = 0.65))
ggsave(p_dot, filename = 'dotplot_annotation.pdf',width = 9,height = 10.5)
ggsave(p_dot_cluster, filename = 'dotplot_cluster.pdf',width = 10,height = 12)

table(immune.combined$orig.ident,immune.combined$annotation)/rowSums(table(immune.combined$orig.ident,immune.combined$annotation))*100
table(immune.combined$condition,immune.combined$annotation)/rowSums(table(immune.combined$condition,immune.combined$annotation))*100


cluster_unique <- c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','HSC/FB','KC','immune cells')
color.use <- scales::hue_pal(l=75,c=150,h=c(0,360),h.start=10,direction = 1)(length(cluster_unique))
color.use <- structure(color.use, names=cluster_unique)
immune.combined$annotation <- factor(immune.combined$annotation, levels = cluster_unique)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident",shuffle=TRUE,pt.size = 0.5)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE, pt.size = 0.5,group.by = 'seurat_clusters1',label.size = 10)
p2_1 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE, pt.size = 0.5,group.by = 'annotation',label.size = 5) +
  scale_color_manual(name='Annotation',values = color.use) +
  theme(legend.text=element_text(size=14),legend.title = element_text(size=14)) +
  theme(text = element_text(size = 14), axis.text = element_text(size = 14))

ggsave('umap.pdf',plot =  p1 + p2 + p2_1, width=75,height=25,units='cm')
ggsave('umap_single.pdf',plot = p2_1, width=17,height=13,units='cm')

p3 <- DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident")
ggsave('umap_by_sample.pdf',plot =  p3, width=50,height=15,units='cm')

key_genes <- c('Glul','Gulo','Oat','Cyp2e1','Cyp2f2','Cyp1a2','Cyp2c29','Pigr','Pon1','Rgn','Arg1','Igfbp2','Sds','Fbp1','Pck1','Cps1','Alb','Serpina3k','Hamp','Hamp2')
pfeature <- FeaturePlot(immune.combined,features =key_genes )
ggsave('feature_plot.pdf',plot =  pfeature, width=50,height=50,units='cm')

p1 <- VlnPlot(immune.combined,features = 'percent.mt',group.by =  'seurat_clusters',pt.size = 0.1)
p2 <- VlnPlot(immune.combined,features = 'nCount_RNA',group.by =  'seurat_clusters',pt.size = 0.1)
p3 <- VlnPlot(immune.combined,features = 'nFeature_RNA',group.by =  'seurat_clusters',pt.size = 0.1)
p123 <- p1 + p2 +p3 + patchwork::plot_layout(ncol = 1)
ggsave('mt_percent_nCount_nFeature.pdf',plot =  p123, width=50,height=50,units='cm')

## statistics & table
ident_cluster <- data.frame(ident = immune.combined@meta.data$orig.ident,cluster=immune.combined@meta.data$annotation,value=1,row.names = colnames(immune.combined))
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(ident_cluster$ident))
ident_cluster$condition <- tmp[ident_cluster$ident]
ident_cluster$condition <- structure(c('YS','OS','OR'),names=c('Young','Old sensitive','Old resistant'))[ident_cluster$condition]
ident_cluster$condition <- factor(ident_cluster$condition, levels = c('YS','OS','OR'))
ident_cluster$cluster <- factor(ident_cluster$cluster, levels = c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','HSC/FB','KC','immune cells'))

ident_cluster_new <- ident_cluster %>%
  dplyr::group_by(condition, cluster) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(perc = count/sum(count))
ident_cluster_new2 <- ident_cluster %>%
  dplyr::group_by(ident, cluster) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(perc = count/sum(count))

library(ggplot2)
library("ggsci")
library("ggplot2")
library("gridExtra")
# pdistribution <- ggplot(ident_cluster_new[!(ident_cluster_new$cluster %in% c('zone 2','doublet')),], aes(fill=cluster, y=perc*100, x=condition)) +
pdistribution <- ggplot(ident_cluster_new, aes(fill=cluster, y=perc*100, x=condition)) +
  geom_bar(stat="identity")+labs(x = "Condition", y = "Percent")+
  # geom_text(aes(label = paste0(round(perc,4)*100,"%")),
  #           position = position_stack(vjust = 0.5), size = 4) +
  theme(text = element_text(size = 14), axis.text = element_text(size = 14), legend.text=element_text(size=14)) +
  scale_fill_discrete(name = "annotation") +
  scale_fill_manual(values = color.use) #+ scale_fill_simpsons()
ggsave('population_size_no_numbers.pdf',plot =  pdistribution, width=13,height=13,units='cm')
# ggsave('population_size_with_numbers.pdf',plot =  pdistribution, width=8,height=13,units='cm')

pdistribution_list <- list()
for(j in 1:length(cluster_unique)){
  pdistribution_list[[j]] <- ggplot(ident_cluster_new[ident_cluster_new$cluster==cluster_unique[j],], aes(fill=cluster, y=perc*100, x=condition)) +
    geom_bar(stat="identity")+labs(x = "Condition", y = "Percent")+ ggtitle(cluster_unique[j]) +
    geom_text(aes(label = paste0(round(perc,4)*100,"%")),
              position = position_stack(vjust = 0.5), size = 4) + scale_fill_manual(name = "annotation",values=color.use[j]) +
    theme(text = element_text(size = 14), axis.text = element_text(size = 14), legend.text=element_text(size=14)) +
    theme(legend.position="none")
}
ggsave('population_size_split.pdf',plot = wrap_plots(pdistribution_list,ncol = 3), width=27,height=25,units='cm')


#### mapping to spatial data not combined
sample_sc <- sort(unique(immune.combined$orig.ident))
sample_sp <- c('B1','C2','A2','C1','A1','B2')
sample_name <- c('A1','B2','A2','C1','B1','C2')
exp_condition <- structure(rep(c('Young','Old_sensitive','Old_resistant'),each=2),names=sample_name)

seurat_sc <- immune.combined#[,!immune.combined$annotation %in% c('Hep doublet','EC/Hep doublet','undefined 2')]
seurat_sp <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
DefaultAssay(seurat_sc) <- 'SCT'
seurat_sc <- PrepSCTFindMarkers(seurat_sc)
DefaultAssay(seurat_sp) <- 'SCT'
seurat_sp <- PrepSCTFindMarkers(seurat_sp)
seurat_sc_new <- CreateSeuratObject(seurat_sc@assays$SCT@counts)
seurat_sc_new@assays$RNA@data <- seurat_sc@assays$SCT@data
seurat_sc_new@meta.data <- seurat_sc@meta.data
seurat_sp_new <- CreateSeuratObject(seurat_sp@assays$SCT@counts)
seurat_sp_new@assays$RNA@data <- seurat_sp@assays$SCT@data
seurat_sp_new@meta.data <- seurat_sp@meta.data

anchors <- FindTransferAnchors(reference = seurat_sc_new , query = seurat_sp_new, reference.assay = 'RNA',
                               query.assay = 'RNA',normalization.method = 'LogNormalize', reduction = 'rpca', n.trees = 50, dims=1:30,features = rownames(seurat_sp_new),
                               k.anchor = 10,k.filter=200, k.score = 20) # 10, 200, 20
# anchors <- FindTransferAnchors(reference = seurat_sc , query = seurat_sp, reference.assay = 'integrated',
#                                query.assay = 'integrated',normalization.method = 'SCT', reduction = 'rpca', n.trees = 50, dims=1:30, #features = c(seurat_sc@assays$integrated@var.features,genes_to_plot),
#                                k.anchor = 5,k.filter=5, k.score = 20) # 5, 5, 20
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat_sc$seurat_clusters1, prediction.assay = TRUE,
                                   dims = 1:30, weight.reduction = seurat_sp[["pca"]],  #'cca'
                                  sd.weight = 1,
                                  n.trees=50,k.weight = 10) # 1：30, 1, 50, 10
seurat_sp[["predictions"]] <- predictions.assay

## Alluvial plot
sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
sp_to_sc$sp_cluster <- paste('sp-',seurat_sp$seurat_clusters[rownames(sp_to_sc)],sep='')
sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('sp-',c(0,1,2,5,3,4,6,7),sep=''))

df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
# df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  c('zone 1','zone 2', 'zone 3','zone 3 Glul+','EC','EC/Hep doublet','HSC/FB','KC','immune cells','undefined 1','undefined 2','Hep doublet'))
df_long$sc_cluster<- factor(df_long$sc_cluster,levels = c('2','0','1','6','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
sc_cluster_annotation <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                                   names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
df_long$sc_cluster2 <- factor(paste(sc_cluster_annotation[df_long$sc_cluster], ' (',df_long$sc_cluster,')',sep=''),
                              levels = paste(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'), ' (',c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),')',sep=''))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('sp-',c(0,1,2,5,3,4,7,6),sep=''))
df_long$sc_annotation <- factor(sc_cluster_annotation[df_long$sc_cluster],levels= c('zone 1','zone 2', 'zone 3','zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'))

df_long$Total <- table(seurat_sp$seurat_clusters[rownames(sp_to_sc)])[sub('sp-','',df_long$sp_cluster)]
df_long$ratio <- df_long$Freq/df_long$Total

##  heatmap
library(circlize)
col_fun = colorRamp2(c(0,  1), c("white","blue"))
sp_to_sc_hmp <- reshape2::dcast(df_long,  sc_annotation ~ sp_cluster,value.var = 'Freq',fun.aggregate=sum)
rownames(sp_to_sc_hmp) <- sp_to_sc_hmp[,1]
sp_to_sc_hmp <- sp_to_sc_hmp[,-1]
ComplexHeatmap::Heatmap((t(sp_to_sc_hmp)/colSums(sp_to_sc_hmp)),col = col_fun(seq(0, 1,by=0.1)),
                        cluster_rows = F, cluster_columns =F, name = 'ratio \n (row-normalized)',
                        row_names_side = 'left', column_names_rot = 45,
                        column_names_side = 'top')

df_long <- df_long[df_long$ratio >0.1 ,]
library(ggalluvial)
pmapping <- ggplot(data = df_long,
                   # aes(axis1 = sp_cluster, axis2 = sc_annotation , axis3=sc_cluster2, y= Freq)) +
                   aes(axis1 = sp_cluster, axis2 = sc_annotation, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/4) +
  geom_stratum(width = 1/4) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void()
pmapping
# saveRDS(immune.combined,'Hep_integration_SoupX_reintegration_0805.rds')
ggsave('sp_sc_mapping_new.pdf',plot =  pmapping, width=25,height=30,units='cm')

### bubble plot
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
## need to reformat the df_long (making each point unique)
df_long$ratio <- as.numeric(df_long$ratio)
xx = ggplot(df_long, aes(y = sp_cluster, x = sc_annotation)) +
  geom_point(aes(size = (Freq), fill = ratio), alpha = 0.5, shape = 21,color='black') +
  scale_size(limits = c(min(df_long$Freq),max(df_long$Freq)),range=c(1,17)) +
  labs( x= "sc-cluster", y = "sp-cluster", size = "# of spatial cells", fill = "proportion")  +
  scale_y_discrete(limits=rev) +
  # scale_x_discrete(limits=rev) +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 16, face = "bold", angle = 60, vjust = 0.8, hjust = 1),
        axis.text.y = element_text(colour = "black", face = "bold", size = 16),
        axis.title = element_text(colour = "black", face = "bold", size = 16),
        legend.text = element_text(size = 16, face ="bold", colour ="black"),
        legend.title = element_text(size = 16, face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "right") + scale_fill_gradient(low="gray", high="blue")

ggsave('sp_sc_mapping_bubble.pdf',plot =  xx, width=25,height=20,units='cm')

####### integration ATAC-seq with scRNA-seq as a whole
hm.integrated <- readRDS('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/integrated_ATAC.rds');


####### integration ATAC-seq with scRNA-seq,  dataset-wsie
project_name <- c('Young_1','Old_resistant_1','Old_sensitive_1','Old_resistant_2','Old_sensitive_2','Young_2')
hm.integrated <- readRDS('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/integrated_ATAC.rds');
DefaultAssay(immune.combined) <- "integrated"

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
  # sp_to_sc$sc_cluster <- paste('sc-',sp_to_sc$sc_cluster,sep = '')
  sp_to_sc$sp_cluster <- paste('ATAC-',atac_j$seurat_clusters[rownames(sp_to_sc)],sep='')
  sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))

  mapping_table_list[[j]] <- table(sp_to_sc[,c(1,2)])

  df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
  # df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  paste('sc-',c('0','1','4','6','9','2','2_1','3','5','8','7','10'),sep=''))
  df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('ATAC-',c(0,1,2,3,4,5,6,7),sep=''))
  df_long$sc_cluster<- factor(df_long$sc_cluster,levels = c('2','0','1','6','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
  sc_cluster_annotation <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/Hep doublet','HSC/FB','KC','immune cells'),
                                     names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
  df_long$sc_cluster2 <- factor(paste(sc_cluster_annotation[df_long$sc_cluster], ' (',df_long$sc_cluster,')',sep=''),
                                levels = paste(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/Hep doublet','HSC/FB','KC','immune cells'), ' (',c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),')',sep=''))
  df_long$sc_annotation <- factor(sc_cluster_annotation[df_long$sc_cluster],levels= c('zone 1','zone 2', 'zone 3','zone 3 Glul+','EC','EC/Hep doublet','HSC/FB','KC','immune cells'))

  df_long$Total <- table(atac_j$seurat_clusters[rownames(sp_to_sc)])[sub('ATAC-','',df_long$sp_cluster)] + 1e-6
  df_long$ratio <- df_long$Freq/df_long$Total

  df_long_list[[j]] <- df_long

  df_long <- df_long[df_long$ratio>0.05,]
  library(ggalluvial)
  pmapping <- ggplot(data = df_long,
                     aes(axis1 = sp_cluster, axis2 = sc_annotation , axis3=sc_cluster2, y= Freq)) +
    geom_alluvium(aes(fill = sp_cluster),width = 1/4) +
    geom_stratum(width = 1/4) +
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
mapping_table_sum <- matrix(0,nrow = length(unique(immune.combined$seurat_clusters1)),
                            ncol=length(unique(hm.integrated$seurat_clusters)),
                            dimnames = list(c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),paste('ATAC-',c(0,1,2,3,4,5,6,7), sep=''))
)
mapping_table_0 <- mapping_table_sum
for(j in 1:6){
  mapping_table_j <- mapping_table_0
  mapping_table_j[rownames(mapping_table_list[[j]]),colnames(mapping_table_list[[j]])] <- mapping_table_list[[j]]
  mapping_table_sum <- mapping_table_sum + mapping_table_j
}

## Alluvial plot
df_long <- reshape2::melt(mapping_table_sum); names(df_long) <- c('sc_cluster','sp_cluster','Freq')
# df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  c('zone 1','zone 2', 'zone 3','zone 3 Glul+','EC','EC/Hep doublet','HSC/FB','KC','immune cells','undefined 1','undefined 2','Hep doublet'))
df_long$sc_cluster<- factor(df_long$sc_cluster,levels = c('2','0','1','6','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
sc_cluster_annotation <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                                   names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
df_long$sc_cluster2 <- factor(paste(sc_cluster_annotation[df_long$sc_cluster], ' (',df_long$sc_cluster,')',sep=''),
                              levels = paste(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'), ' (',c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),')',sep=''))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels= paste('ATAC-',c(2,3,1,0,4,5,6,7),sep=''))
df_long$sc_annotation <- factor(sc_cluster_annotation[df_long$sc_cluster],levels= c('zone 1','zone 2', 'zone 3','zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'))

df_long$Total <- table(seurat_sp$seurat_clusters[rownames(sp_to_sc)])[sub('ATAC-','',df_long$sp_cluster)]
df_long$ratio <- df_long$Freq/df_long$Total

##  heatmap
library(circlize)
col_fun = colorRamp2(c(0,  1), c("white","blue"))
sp_to_sc_hmp <- reshape2::dcast(df_long,  sc_annotation ~ sp_cluster,value.var = 'Freq',fun.aggregate=sum)
rownames(sp_to_sc_hmp) <- sp_to_sc_hmp[,1]
sp_to_sc_hmp <- sp_to_sc_hmp[,-1]
ComplexHeatmap::Heatmap((t(sp_to_sc_hmp)/colSums(sp_to_sc_hmp)),col = col_fun(seq(0, 1,by=0.1)),
                        cluster_rows = F, cluster_columns =F, name = 'ratio \n (row-normalized)',
                        row_names_side = 'left', column_names_rot = 45,
                        column_names_side = 'top')

df_long <- df_long[df_long$ratio >0.1 ,]
library(ggalluvial)
pmapping_atac <- ggplot(data = df_long,
                   # aes(axis1 = sp_cluster, axis2 = sc_annotation , axis3=sc_cluster2, y= Freq)) +
                   aes(axis1 = sp_cluster, axis2 = sc_annotation, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/4) +
  geom_stratum(width = 1/4) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void()
pmapping_atac
ggsave('atac_sc_mapping.pdf',plot =  pmapping_atac, width=25,height=30,units='cm')

########
#### spatial plot
plot_sp_cluster <- function(sample_single,cluster_idx,sc_cluster_idx){
  # sample_single <- 'A1';
  # cluster_idx <- 0:2
  # sc_cluster_idx <- '3'
  sample_idx <- which(sample_name %in% sample_single)
  cell_names <- Cells(immune.combined.sct)[grepl(sample_single,Cells(immune.combined.sct))
                                           & (immune.combined.sct@meta.data$seurat_clusters) %in% cluster_idx
                                           & (sp_to_sc[Cells(immune.combined.sct),'sc_cluster'] %in% sc_cluster_idx)]
  cluster_id <- (immune.combined.sct@meta.data$seurat_clusters)[which( Cells(immune.combined.sct) %in% cell_names)]
  sc_cluster_id <- as.numeric(sp_to_sc[Cells(immune.combined.sct)[which( Cells(immune.combined.sct) %in% cell_names)],'sc_cluster'])
  cell_names <- sub(".*[A,B,C][1,2]-","",cell_names)
  roi_SpP <- roi_list[[sample_idx]]
  roi_SpP_subset <- roi_SpP[cell_names]
  # col_code <- c('red','yellow','blue','green','orange','purple') # 'white',
  col_code <- scales::hue_pal()(length(unique(immune.combined.sct@meta.data$seurat_clusters)))
  par(mar = c(2, 2, 2, 4))
  plot(roi_SpP_subset,col=col_code[as.numeric(sc_cluster_id)+1],lwd=0.2);
  legend('right',legend=as.numeric(sc_cluster_idx),fill=col_code[as.numeric(sc_cluster_idx)+1], cex=1.5,border = F,box.lwd=0,bty = "n",horiz = F,
         xpd=TRUE, inset=c(-.05,0))
}

plot_sp_cluster('A1',c(0,1,2,5))
plot_sp_cluster('A1',c(3,4,6,7))
plot_sp_cluster('C2',1,c('3','5','7'))
plot_sp_cluster('A1',3,c('10'))
plot_sp_cluster('A1',3,unique(seurat_sc$customer.cluster))

## heatmap for marker genes (RNA assay, Spatial)
immune.combined.sct <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
DefaultAssay(immune.combined.sct) <- "RNA"
immune.combined.sct <- NormalizeData(immune.combined.sct)
immune.combined.sct <- ScaleData(immune.combined.sct,  features = rownames(immune.combined.sct))
immune.combined.markers <- FindAllMarkers(immune.combined.sct, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
immune.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = abs(avg_log2FC)) -> top10
print(top10,n=80)
DoHeatmap(immune.combined.sct[,sample(1:dim(immune.combined.sct)[2],10000,replace=F)], features = top10$gene) + NoLegend()
## heatmap for marker genes (Integration assay, Spatial)
DefaultAssay(immune.combined.sct) <- "integrated"
DoHeatmap(immune.combined.sct[,sample(1:dim(immune.combined.sct)[2],10000,replace=F)], features = top10$gene) + NoLegend()
write.csv(immune.combined.markers,file='sp_markers.csv')

# GO analysis
library(org.Mm.eg.db)
mm <- org.Mm.eg.db
go <- lapply(unique(immune.combined.markers$cluster),FUN=function(x){
  cg <- immune.combined.markers$gene[immune.combined.markers$cluster==x]
  cg_entrezid <- select(mm,
                        keys = cg,
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL")
  library(limma)
  g <- limma::goana(cg_entrezid$ENTREZID,species='Mm')
  g <- g[order(g$P.DE,decreasing = F),]
})
names(go) <- unique(immune.combined.markers$cluster)
go_overlap <- matrix(0,ncol=length(go),nrow=length(go))
Ngo <- 50
for(i in 1:dim(go_overlap)[1]){
  for(j in 1:dim(go_overlap)[2]){
    go_overlap[i,j] <- length(which(rownames(go[[i]])[1:Ngo] %in% rownames(go[[j]])[1:Ngo]))
  }
}
rownames(go_overlap) <- paste('c',unique(immune.combined.markers$cluster),sep = '')
colnames(go_overlap) <- rownames(go_overlap)
library(circlize);
col_fun = circlize::colorRamp2(c(0, Ngo/2, Ngo), c("blue", "white", "red"))
ComplexHeatmap::Heatmap(go_overlap,cluster_rows = F,cluster_columns = F, name = '# of Overlap',
                        column_names_rot = 0,column_names_side = 'top',row_names_side = 'left', col=col_fun,
                        heatmap_legend_param=list(at = c(0,Ngo/2,Ngo), labels = c(0, Ngo/2, Ngo)))
## save go
df_go <- data.frame(go$'0'[1:50,],cluster=0)
for(j in 2:length(go)){
  df_go_tmp <- data.frame(go[[j]][1:50,],cluster=j-1)
  df_go <- rbind(df_go,df_go_tmp)
}
write.csv(df_go,file = 'GO_top50.csv')



# directly use gene from mouse database
go_senescence <- read.csv2('/Users/weizhao/Documents/Hep_Seq/go_senecense.txt',header = F, sep = '\t')
pathway_genes = go_senescence$V2 %>% str_to_title() %>% unique()  %>% .[. %in% rownames(immune.combined)]
immune.combined = AddModuleScore(immune.combined, list(pathway_genes), assay = "RNA", name = "go_senescence")
# plot
library(ggpubr)
p = ggboxplot(immune.combined@meta.data, x = "condition", y = "go_senescence1",
              color = "condition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              facet.by = "annotation", ncol = 5, ylim=c(-0.15,0.45))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

my_comparisons <- list(c("Young", "Old sensitive"), c("Old sensitive", "Old resistant"), c("Young", "Old resistant") )
p + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif", label.y = c(0.3, 0.35, 0.4))+ggtitle("go_senescence")+NoLegend()

go_apoptosis <- read.csv2('/Users/weizhao/Documents/Hep_Seq/go_apoptosis.txt',header = F, sep = '\t')
pathway_genes = go_apoptosis$V2 %>% str_to_title() %>% unique()  %>% .[. %in% rownames(immune.combined)]
immune.combined = AddModuleScore(immune.combined, list(pathway_genes), assay = "RNA", name = "go_apoptosis")
# plot
library(ggpubr)
p = ggboxplot(immune.combined@meta.data, x = "condition", y = "go_apoptosis1",
              color = "condition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              facet.by = "annotation", ncol = 5, ylim=c(-0.15,0.45))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

my_comparisons <- list(c("Young", "Old sensitive"), c("Old sensitive", "Old resistant"), c("Young", "Old resistant") )
p + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif", label.y = c(0.3, 0.35, 0.4))+ggtitle("go_apoptosis1")+NoLegend()


## DEG analysis
immune.combined <- readRDS('/Users/weizhao/Documents/Hep_Seq/Hep16_outs/Hep_integration_SoupX_reintegration_0805.rds')
cluster_celltype <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                              names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
immune.combined$annotation <- cluster_celltype[as.character(immune.combined$seurat_clusters1)]
immune.combined$annotation <- factor(immune.combined$annotation,levels=c('zone 1','zone 3','zone 3 Glul+','zone 2', 'EC','EC/HSC chimera','HSC/FB','KC','immune cells'))
immune.combined$annotation[which(immune.combined@meta.data$sub.cluster.EC %in% c('EC_1'))] <- 'EC/HSC chimera'
immune.combined$seurat_clusters1 <- factor(immune.combined$seurat_clusters1,levels=c('2','10','3','4','5','3-1','6','1','0','7','8','9','9-1','12','11','11-1'))

cluster_unique <- c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','EC periportal','HSC/FB','KC','immune cells')
cluster_unique_sc <- cluster_unique[cluster_unique %in% immune.combined$annotation]
library(ggplot2)
DE_list <- list()
p_list <- list()
pbmc_list <- list()
pbmc <- immune.combined;DefaultAssay(pbmc) <- 'RNA'

for(j in 1:length(cluster_unique_sc)){
  pbmc.j <- subset(pbmc,subset= (annotation==cluster_unique_sc[j]))
  Idents(pbmc.j) <- 'condition'
  pbmc.j@active.ident <- factor(x = pbmc.j@active.ident, levels = c('Young','Old sensitive','Old resistant'))
  pbmc_list[[j]] <- pbmc.j
  DE.j <- FindAllMarkers(pbmc.j, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5,
                         min.cells.group = 1)
  DE.j <- DE.j[grep('^MT',DE.j$gene,invert=TRUE),]
  DE.j <- DE.j[grep('^RPL',DE.j$gene,invert=TRUE),]
  DE.j <- DE.j[grep('^RPS',DE.j$gene,invert=TRUE),]
  DE.j <- DE.j[grep('^AY036118',DE.j$gene,invert = T),]
  DE.j <- DE.j[grep('^Gm42418',DE.j$gene,invert = T),]
  DE_list[[j]] <- DE.j
  names(DE_list)[j] <- cluster_unique_sc[j]
  if(dim(DE.j)[1]>0){
    DE.j %>%
      group_by(cluster) %>%
      slice_max(n = 30, order_by = avg_log2FC) -> top10
    pbmc.j <- ScaleData(pbmc.j, top10$gene)
    p_list[[j]] <- DoHeatmap(pbmc.j, features = top10$gene) + ggtitle(paste('cluster',cluster_unique_sc[j]))#+NoLegend()
  } else {
    p_list[[j]] <- ggplot() + theme_void()
  }
  ggsave(paste(paste('Heatmap_cluster_',sub('/','_',cluster_unique_sc[j]),sep=''),'.pdf',sep=''),
         plot =   p_list[[j]],
         width=25,height=40*length(top10$gene)/90,units='cm')
  write.csv(DE.j,paste(paste('sc_DEG_cluster_',sub('/','_',cluster_unique_sc[j]),sep=''),'.csv',sep=''))
}
p_list