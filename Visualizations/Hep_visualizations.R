genes_to_plot <- c('Nmnat1','Nap1l1','Zc3h13','Slc1a4','Rest')
genes_to_plot <- c('Wnt2','Wnt9b','Zbtb14','Zbtb33')
genes_to_plot <- c('Lrp5','Lrp6','Zbtb14')
genes_to_plot <- c(        'Tgfb1')

library(Seurat);library(Signac)
library(ggplot2)
library(dplyr)
setwd('/Users/weizhao/Documents/Hep_Seq/')

hm.integrated <- readRDS('ATAC-seq/integrated_ATAC.rds')
DefaultAssay(hm.integrated) <- 'RNA'
sample_name <- unique(hm.integrated$dataset)
exp_condition <- structure(c('Young','Old_resistant','Old_sensitive','Old_resistant','Old_sensitive','Young'),names=sample_name)
hm.integrated$condition <- exp_condition[hm.integrated$dataset]
VlnPlot(hm.integrated,genes_to_plot[1],group.by = 'seurat_clusters',pt.size = 0.1,split.by = 'condition')

sp_seurat <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
DefaultAssay(sp_seurat) <- 'RNA'
sample_name <- c('A1','B2','A2','C1','B1','C2')
exp_condition <- structure(rep(c('Young','Old_sensitive','Old_resistant'),each=2),names=sample_name)
sp_seurat$condition <- exp_condition[sp_seurat$orig.ident]
sp_seurat <- NormalizeData(sp_seurat)
VlnPlot(sp_seurat,genes_to_plot,group.by = 'seurat_clusters',pt.size = 0.01,split.by = 'condition')

vln_list <- list()
for(j in 1:length(genes_to_plot[genes_to_plot %in% rownames(sp_seurat)])){
  vln_list[[j]] <- VlnPlot(sp_seurat,genes_to_plot[genes_to_plot %in% rownames(sp_seurat)][j],group.by = 'seurat_clusters',pt.size = 0.01,split.by = 'condition')
  ggsave(paste(genes_to_plot[genes_to_plot %in% rownames(sp_seurat)][j],'_spatial.pdf',sep=''),plot =  vln_list[[j]], width=30,height=15,units='cm')
}


immune.combined <- readRDS('Hep_integration_NatMeta_SCT_low2000.rds')
# immune.combined <- readRDS('Hep_integration_NatMeta_SCT_low2000_filter_prior_integration.rds')
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- subset(immune.combined,subset =  percent.mt >=0.1)
# Run the standard workflow for visualization and clustering
npcs =30;
immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
ElbowPlot(immune.combined,ndims = npcs)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 3)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:npcs,n.components = 2)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs,k.param=10)
immune.combined <- FindClusters(immune.combined, resolution = 0.2) # common results 10 clusters
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- FindClusters(immune.combined, resolution = 0.5,method='matrix',n.start = 10, n.iter=10, algorithm=1) # common results 10 clusters
DimPlot(immune.combined, reduction = "umap",pt.size = 0.5,label=T, label.size = 5,dims = c(1,2))

# sub clustering
immune.combined <- FindSubCluster(object=immune.combined, cluster= 3, graph.name = 'integrated_snn', resolution = 0.1,subcluster.name = "sub.cluster3")
immune.combined@meta.data$customer.cluster <- as.character(immune.combined@meta.data$sub.cluster3)
immune.combined@meta.data$customer.cluster[which(immune.combined@meta.data$sub.cluster3 %in% c('3_0','3_1','3_2'))] <- '3'
immune.combined@meta.data$customer.cluster[which(immune.combined@meta.data$sub.cluster3 %in% c('3_3'))] <- '3_1'
immune.combined@meta.data$customer.cluster <- factor(immune.combined@meta.data$customer.cluster,levels = c('0','1','2','5','8','4','3','3_1','6','11','9','7','10','12','13'))
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE,group.by = 'sub.cluster3')
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE,group.by = 'customer.cluster',pt.size = 1)

DefaultAssay(immune.combined) <- 'RNA'
immune.combined <- NormalizeData(immune.combined)
sample_name <- unique(immune.combined$orig.ident)
exp_condition <- structure(c('Young','Old_resistant','Old_sensitive','Old_resistant','Old_sensitive','Young'),names=sample_name)
immune.combined$condition <- exp_condition[immune.combined$orig.ident]
vln_list <- list()
genes_to_plot <- c('Plg','Pard3','Nampt','Insr','F2','Agt','Agtr1a')
genes_to_plot <- stringr::str_to_title(top2_TF$TF)
genes_to_plot <- c('Hgf','Met')
for(j in 1:length(genes_to_plot)){
  vln_list[[j]] <- VlnPlot(immune.combined,genes_to_plot[j],group.by = 'seurat_clusters',pt.size = 0,split.by = 'condition')
  ggsave(paste(genes_to_plot[j],'.pdf'),plot =  vln_list[[j]], width=30,height=15,units='cm')
}


a <- AverageExpression(
  subset(immune.combined,subset = customer.cluster %in% c('7')),
  features = genes_to_plot[1],
  return.seurat = FALSE,
  assays = 'RNA',
  group.by = "condition"
)$RNA

b <- AverageExpression(
  subset(immune.combined,subset = customer.cluster %in% c('4','3')),
  features = genes_to_plot[2],
  return.seurat = FALSE,
  assays = 'RNA',
  group.by = "condition",
  slot = 'data'
)$RNA




library(Seurat)
library(ggplot2)
setwd('/Users/weizhao/Documents/Hep_Seq/')
seurat_sc <- readRDS('/Users/weizhao/Documents/Hep_Seq/Hep16_outs/Hep_integration_SoupX_reintegration_0805.rds')
cluster_celltype_sc <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                                 names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
seurat_sc$annotation <- cluster_celltype_sc[as.character(seurat_sc$seurat_clusters1)]
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(seurat_sc$orig.ident))
seurat_sc$condition <- tmp[seurat_sc$orig.ident]

seurat_sp <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
cluster_celltype_sp <- structure(c('zone 1','zone 2','zone 3','EC','HSC/FB','zone 3 Glul+','KC','EC periportal'),names=0:7)
seurat_sp$annotation <- cluster_celltype_sp[as.character(seurat_sp$seurat_clusters)]
tmp_sp <- structure(rep(c('Young','Old sensitive','Old resistant'),each=2),names=c('A1','B2','A2','C1','B1','C2'))
seurat_sp$condition <- tmp_sp[seurat_sp$orig.ident]

seurat_atac <- readRDS('ATAC-seq/integrated_ATAC.rds')
cluster_celltype_atac <- structure(c('zone 3','zone 2','zone 1','zone 1','EC','HSC/FB','KC','immune cells'),names=0:7)
seurat_atac$annotation <- cluster_celltype_atac[as.character(seurat_atac$seurat_clusters)]
tmp_atac <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(seurat_atac$dataset))
seurat_atac$condition <- tmp_atac[seurat_atac$dataset]

cluster_unique <- c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','EC periportal','HSC/FB','KC','immune cells')
color.use <- scales::hue_pal(l=75,c=150,h=c(0,360),h.start=10,direction = 1)(length(cluster_unique))
color.use <- structure(color.use, names=cluster_unique)

seurat_sc$annotation <- factor(seurat_sc$annotation, levels=cluster_unique)
seurat_sp$annotation <- factor(seurat_sp$annotation, levels=cluster_unique)
seurat_atac$annotation <- factor(seurat_atac$annotation, levels=cluster_unique)

seurat_sc$condition <- factor(seurat_sc$condition, levels=c('Young','Old sensitive','Old resistant'))

## statistics & table
seurat_list <- list(seurat_sc, seurat_sp,seurat_atac)
omics_names <- c('scRNA-seq','Spatial transcriptomics','scATAC-seq')
plist<- list()
ulist <- list()
for(j in 1:length(seurat_list)){
  ident_cluster <- data.frame(ident = colnames(seurat_list[[j]]),condition =seurat_list[[j]]$condition , cluster=seurat_list[[j]]$annotation,value=1)
  ident_cluster$condition <- structure(c('YS','OS','OR'),names=c('Young','Old sensitive','Old resistant'))[ident_cluster$condition]
  ident_cluster$condition <- factor(ident_cluster$condition, levels = c('YS','OS','OR'))
  ident_cluster$cluster <- factor(ident_cluster$cluster,levels=cluster_unique)
  ident_cluster_new <- ident_cluster %>%
    dplyr::group_by(condition, cluster) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(perc = count/sum(count))
  tmp <- data.frame(v1 = ifelse(ident_cluster_new$perc > 0.015,paste0(round(ident_cluster_new$perc,4)*100,"%"),''))
  ident_cluster_new$perc_label <- tmp$v1
  plist[[j]] <- ggplot(ident_cluster_new, aes(fill=cluster, y=perc*100, x=condition)) +
    geom_bar(stat="identity")+labs(x = "Condition", y = "Percent")+
    geom_text(aes(label = perc_label),
              position = position_stack(vjust = 0.5), size = 4) +
    theme(text = element_text(size = 14), axis.text = element_text(size = 14), legend.text=element_text(size=14)) +
    scale_fill_manual(name = "Annotation", values = color.use)
  # seurat_list[[j]]$annotation <- factor(seurat_list[[j]]$annotation,levels=c('zone 3','zone 2','zone 1','EC','HSC','KC','Immune cells','Undefined'))
  ulist[[j]]  <- DimPlot(seurat_list[[j]], reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE, pt.size = 0.5,group.by = 'annotation',label.size = 5, raster=FALSE) +
    scale_color_manual(name='Annotation',values = color.use) +
    theme(legend.text=element_text(size=14),legend.title = element_text(size=14)) +
    theme(text = element_text(size = 14), axis.text = element_text(size = 14))

  pdistribution_list <- list()
  idx <-  which(cluster_unique %in% unique(ident_cluster$cluster))
  condition_cluster <- table(ident_cluster$condition,ident_cluster$cluster)
  condition_cluster_prop <- data.frame((condition_cluster/rowSums(condition_cluster)))
  names(condition_cluster_prop) <- c('condition','cluster','perc')
  for(jj in 1:length(idx)){
    # pdistribution_list[[jj]] <- ggplot(ident_cluster_new[ident_cluster_new$cluster==cluster_unique[idx[jj]],], aes(fill=cluster, y=perc*100, x=condition)) +
    #   geom_bar(stat="identity")+labs(x = "Condition", y = "Percent")+ ggtitle(cluster_unique[idx[jj]]) +
    #   geom_text(aes(label = paste0(round(perc,4)*100,"%")),
    #             position = position_stack(vjust = 0.5), size = 4) + scale_fill_manual(name = "annotation",values=color.use[idx[jj]]) +
    #   theme(text = element_text(size = 14), axis.text = element_text(size = 14), legend.text=element_text(size=14)) +
    #   theme(legend.position="none")
    # https://stackoverflow.com/questions/55936572/add-significance-bars-to-proportions-plot-using-ggplot2
    library(ggsignif)
    library(broom)
    library(tidyverse)
    condition_cluster_prop_jj <- condition_cluster_prop %>% dplyr::filter(cluster ==cluster_unique[idx[jj]])
    res <- pairwise.prop.test(condition_cluster[,cluster_unique[idx[jj]]],table(ident_cluster$condition)) %>%
      tidy() %>% mutate_at(vars(contains("group")), ~factor(.x, rownames(condition_cluster)))
    res$annotation <- ifelse(res$p.value<0.001, "***",ifelse(res$p.value < 0.01,'**',ifelse(res$p.value < 0.05,"*","NS")))

    pdistribution_list[[jj]] <- ggplot(condition_cluster_prop_jj, aes(fill=cluster, y=perc*100, x=condition)) +
      geom_bar(stat="identity")+labs(x = "Condition", y = "Percent")+ ggtitle(cluster_unique[idx[jj]]) +
      geom_text(aes(label = paste0(round(perc,4)*100,"%")),
                position = position_stack(vjust = 0.5), size = 4) + scale_fill_manual(name = "annotation",values=color.use[idx[jj]]) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      theme(text = element_text(size = 14), axis.text = element_text(size = 14), legend.text=element_text(size=14)) +
      theme(legend.position="none") + ylim(c(0, max(condition_cluster_prop_jj$perc*100)*1.4)) +
      ggsignif::geom_signif( #comparisons = list(c("YS", "OS"), c("OS", "OR"), c("YS", "OR") ),
                             #map_signif_level = TRUE,
                            y_position = max(condition_cluster[,cluster_unique[idx[jj]]]/table(ident_cluster$condition)*100)*c(1.3,1.07,1.2),
                             xmax = as.integer(res$group1),
                             xmin = as.integer(res$group2),
                             annotation = res$annotation,
                             tip_length = 0.03, #min(max(condition_cluster_prop_jj$perc*100)*0.01,0.1)
                            size=0.2
                                )

  }
  ggsave(paste('umap',j, '.pdf',sep=''),plot = ulist[[j]], width=17,height=13,units='cm')
  ggsave(paste('prop',j, '.pdf',sep=''),plot = plist[[j]], width=13,height=13,units='cm')
  ggsave(paste('prop_split',j, '.pdf',sep=''),
         plot = wrap_plots(pdistribution_list,ncol = 3),
         width=27,height=25,units='cm')
}
# if(j %in% 1:2){
#   plist[[j]]  <- plist[[j]] + theme(legend.position = "none")
#   ulist[[j]] <-  ulist[[j]] + theme(legend.position = "none")
# } else{
# plist[[j]] <- plist[[j]] + scale_fill_manual(name='annotation',
#                                              limits=c('zone 3','zone 2','zone 1','EC','HSC','KC','Immune cells','Undefined'),
#                    values=structure(color.use,names=c('zone 3','zone 2','zone 1','EC','HSC','KC','Immune cells','Undefined')))
# ulist[[j]] <- ulist[[j]] + scale_color_manual(name='annotation',
#                                              limits=c('zone 3','zone 2','zone 1','EC','HSC','KC','Immune cells','Undefined'),
#                                              values=structure(color.use,names=c('zone 3','zone 2','zone 1','EC','HSC','KC','Immune cells','Undefined')))
# }
# ggsave('cell_props_multi.pdf',plot =  plist[[1]] + plist[[2]] + plist[[3]], width=32,height=10,units='cm')
# ggsave('umap_multi.pdf',plot =  ulist[[1]] + ulist[[2]] + ulist[[3]], width=32,height=10,units='cm')
seurat_sp$annotation2 <- factor(seurat_sp$annotation, levels = cluster_unique[cluster_unique %in% seurat_sp$annotation])

#### ST marker genes heatmap
DefaultAssay(seurat_sp) <- "RNA"
seurat_sp <- NormalizeData(seurat_sp)
seurat_sp <- ScaleData(seurat_sp,  features = rownames(seurat_sp))
sp.markers <- presto:::wilcoxauc.Seurat(seurat_sp, only.pos = TRUE,
                                        group_by = 'annotation2', seurat_assay = 'RNA')
DefaultAssay(seurat_sp) <- "SCT"
seurat_sp <- PrepSCTFindMarkers(seurat_sp)
sp.markers <- presto:::wilcoxauc.Seurat(seurat_sp, only.pos = TRUE,
                                        group_by = 'annotation2', seurat_assay = 'SCT')

sp.markers <- sp.markers %>% dplyr::filter(logFC > 0.25 & pct_in > 0.5 & padj < 0.001)
sp.markers %>%
  group_by(feature) %>%
  filter(max(pct_in)-min(pct_in) ==0 | (length(pct_in) >1 & max(pct_in)-sort(pct_in, TRUE)[2] > 30)| feature %in% c('Cyp2e1','Car3','Igfbp2','Txnip','G6pc')) %>%
  slice_max(n = 1, order_by = logFC) -> sp.markers
sp.markers %>%
  group_by(group) %>%
  slice_max(n = 10, order_by = logFC) %>%
  arrange(factor(group, levels = cluster_unique)) -> top10
print(top10,n=8*8)

# top10 <- top10[!(top10$feature %in% c('Ehd3',''))]

p_heat <-DoHeatmap(seurat_sp[,sample(1:dim(seurat_sp)[2],10000,replace=F)],
                   features = top10$feature,group.by = 'annotation2') +  guides(color="none") +
  theme(text = element_text(size = 17))

ggsave(p_heat, filename = 'heatmap_sp_markers.pdf',width = 9,height = 10.5)

#### ATAC marker genes heatmap
DefaultAssay(seurat_atac) <- "RNA"
seurat_atac <- NormalizeData(seurat_atac)
seurat_atac$annotation2 <- factor(seurat_atac$annotation, levels = cluster_unique[cluster_unique %in% seurat_atac$annotation])
atac.markers <- presto:::wilcoxauc.Seurat(seurat_atac, only.pos = TRUE,
                                        group_by = 'annotation2', seurat_assay = 'RNA')
atac.markers <- atac.markers %>% dplyr::filter(logFC > 0.2 & pct_in > 0.2 & padj < 0.001 & pct_in - pct_out > 0 & pct_out < 60)
atac.markers %>%
  group_by(feature) %>%
  slice_max(n = 1, order_by = logFC)   %>%
  group_by(group) %>%
  slice_max(n = 1e5, order_by = logFC) -> atac.markers
atac.markers$pct_diff <- atac.markers$pct_in - atac.markers$pct_out
atac.markers %>%
  group_by(group) %>%
  slice_max(n = 6, order_by = logFC+pct_diff) %>%
  arrange(factor(group, levels = cluster_unique),-pct_diff) -> top10
print(top10,n=6*10)
genes_to_plot <- c('Cyp2f2','Pck1','Uroc1','Pigr','Apoa4','Arg1','C3','G6pc','F2',   # zone 1
                   'Car3','Cyp1a2','Cyp2c29', 'Cyp2e1','Rgn', 'Gulo', 'Cldn2', 'Oat','Glul', # zone 3 'Lect2',  'Axin2',
                   'Hamp','Ttr','Fabp1','Apoa2','Apoc1','Apoe','Apoc3', 'Wfdc21', 'Aldob', # zone 2
                   'Dnase1l3','Ptprb','Clec4g','Igfbp7', # EC
                   'Dcn','Rgs5','Reln',                  # HSC/FB
                   'Ptgs1', 'Clec4f','Wfdc17','Cd74', "Dock10" ,     # KC
                   'Ptprc','Coro1a', 'Lsp1'
                   #,'Scd1','Eef1a1','Eef2','Klf2'
)
seurat_atac <- ScaleData(seurat_atac,features = c(top10$feature,genes_to_plot))
# seurat_atac_renamed <- RenameCells(seurat_atac, new.names= paste('cell',sample(1:dim(seurat_atac)[2],dim(seurat_atac)[2],replace=F),sep='_'))

genes_to_plot_atac <- c('Cyp2c29','Cyp2e1','Mug1',
                        'G6pc','Hal','Cyp2f2',
                        'Slc2a2','Lin7a','Kcnh7',
                        'Ptprb','Kdr','Plxnc1',
                        'Col4a1','Col5a1','Reln',
                        'Clec4f','Rgs1','Cd74',
                        'Ptprc','Fli1','Ebf1'
                        )
p_heat_atac <-myDoHeatmap(seurat_atac,
                          features =c(top10$feature),group.by = 'annotation2',disp.max = 1.5) + 
  theme(text = element_text(size = 17))
ggsave(p_heat_atac, filename = 'heatmap_atac_markers_raw.pdf',width = 9,height = 10.5)
p_heat_atac <-myDoHeatmap(seurat_atac,
                          features =genes_to_plot_atac,group.by = 'annotation2',disp.max = 1.5) +
  theme(text = element_text(size = 17))
ggsave(p_heat_atac, filename = 'heatmap_atac_markers.pdf',width = 9,height = 8)
p_dot_atac <- DotPlot(seurat_atac, features =genes_to_plot_atac,group.by = 'annotation2') + coord_flip() +  scale_x_discrete(limits=rev) + ylab('annotation') +
  theme(axis.text.y = element_text(size = 14), axis.text.x=element_text(size = 14,angle = 45, vjust = 0.9, hjust = 0.9))
ggsave(p_dot_atac, filename = 'dotplot_atac_markers.pdf',width = 9,height = 8)

DefaultAssay(seurat_atac) <- 'ATAC'
da_peaks <- presto:::wilcoxauc.Seurat(  seurat_atac, only.pos = TRUE,
                                        group_by = 'annotation2', seurat_assay = 'ATAC')# , test.use = 'LR', latent.vars = 'peak_region_fragments',)
da_peaks <- da_peaks %>% dplyr::filter(logFC > 0.1 & pct_in > 0.1 & padj < 0.001)
ClosestFeature_da_peaks <- ClosestFeature(seurat_atac, regions = (da_peaks$feature))
head(ClosestFeature_da_peaks)
da_peaks$closest_gene <- ClosestFeature_da_peaks$gene_name[match(da_peaks$feature,ClosestFeature_da_peaks$query_region)]
da_peaks$distance <- ClosestFeature_da_peaks$distance[match(da_peaks$feature,ClosestFeature_da_peaks$query_region)]

da_peaks <- da_peaks %>% dplyr::filter(!is.na(closest_gene))
da_peaks %>%
  group_by(feature) %>%
  slice_max(n = 1, order_by = logFC) -> da_peaks
da_peaks %>%
  group_by(group) %>%
  slice_max(n = 20, order_by = logFC) %>%
  arrange(factor(group, levels = cluster_unique)) -> top10_peak
print(top10_peak,n=8*8)
seurat_atac <- ScaleData(seurat_atac,  features = unique(top10_peak$feature))
chr2gene <- structure(paste(top10_peak$feature, ' (', top10_peak$closest_gene, ')',sep=''),names=top10_peak$feature)
p_heat_atac2 <-myDoHeatmap(seurat_atac,
                        features = top10_peak$feature,group.by = 'annotation2',assay = 'ATAC') + NoLegend() +
  scale_y_discrete(labels = chr2gene) +
  theme(text = element_text(size = 17))
ggsave(p_heat_atac2, filename = 'heatmap_atac_markers2.pdf',width = 12,height = 10.5)

hm.integrated <- readRDS('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/ATAC_integrated_chromvar.rds')
DefaultAssay(hm.integrated) <- 'chromvar'
cluster_celltype_atac <- structure(c('zone 3','zone 2','zone 1','zone 1','EC','HSC/FB','KC','Immune cells'),names=0:7)
hm.integrated$annotation <- cluster_celltype_atac[as.character(hm.integrated$seurat_clusters)]
tmp_atac <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(hm.integrated$dataset))
hm.integrated$condition <- factor(tmp_atac[hm.integrated$dataset],levels = c('Young','Old sensitive','Old resistant'))
hm.integrated$annotation <- factor(hm.integrated$annotation, levels=c('zone 3','zone 2','zone 1','EC','HSC/FB','KC','Immune cells'))
Idents(hm.integrated) <- 'annotation'
differential.activity <- FindAllMarkers(
  object = hm.integrated,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",logfc.threshold = 0.1,
)

DefaultAssay(hm.integrated) <- 'ATAC'
differential.activity$TF <- ConvertMotifID(hm.integrated,id=differential.activity$gene)
differential.activity %>%
  group_by(TF) %>%
  slice_max(n = 1, order_by = avg_diff) -> differential.activity
differential.activity  %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_diff) -> top10_TF
print(top10_TF, n=5*8)
motif2TF <- structure(top10_TF$TF,names=top10_TF$gene)
DefaultAssay(hm.integrated)  <-'chromvar'
hm.integrated <- ScaleData(hm.integrated, features = top10_TF$gene)
p_heat_atac3 <-myDoHeatmap(hm.integrated, disp.max = 2,
                           features = top10_TF$gene,group.by = 'annotation',assay = 'chromvar',slot = 'scale.data') + 
  scale_y_discrete(labels = motif2TF) +
  theme(text = element_text(size = 17))
ggsave(p_heat_atac3, filename = 'heatmap_atac_markers3.pdf',width = 12,height = 10.5)


DefaultAssay(seurat_atac) <- 'ATAC'
cov_plot <- CoveragePlot(
  object = seurat_atac, group.by = 'annotation2',
  region = c('Cyp2c29','Cyp2f2','Hamp','Ptprb','Reln','Clec4f','Ptprc'),
  annotation = T,
  peaks = FALSE,extend.upstream = 5000, extend.downstream = 0, show.bulk	= T, window = 500,
  ncol =3
) + theme(axis.title.x = element_text(size = 14))
ggsave(cov_plot, filename = 'coverage_plot.pdf',width = 16,height = 14)


p_cov1 <- CoveragePlot(
  object = seurat_atac, group.by = 'annotation2',
  region = 'chr1-72772703-72773616 ',
  annotation = T,
  peaks = FALSE,extend.upstream =0, extend.downstream = 60000, show.bulk	= T, window = 500,
  ncol =3
) + theme(axis.title.x = element_text(size = 14))
p_cov2 <- CoveragePlot(
  object = seurat_atac, group.by = 'annotation2',
  region = 'chr1-72772703-72773616 ',
  annotation = T,
  peaks = FALSE,extend.upstream =0, extend.downstream = 0, show.bulk	= T, window = 500,
  ncol =3
) + theme(axis.title.x = element_text(size = 14))

ggsave(p_cov1, filename = 'coverage_plot_Igfbp2_1.pdf',width = 15,height = 8)
ggsave(p_cov2, filename = 'coverage_plot_Igfbp2_2.pdf',width = 15,height = 8)


markers_CV <-  read.csv('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/Diff_chromvar/Diff_TF_fromATAC_zone 3.csv')
markers_CV %>% dplyr::filter(avg_diff>0.1 & p_val_adj < 0.001) -> markers_CV
hm.integrated <- readRDS('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/ATAC_integrated_chromvar.rds')
avg_TF <- AverageExpression(subset(hm.integrated,subset = annotation =='zone 3' ), features = unique(markers_CV$gene), group.by = 'condition',assays = 'chromvar')
avg_TF <- as.data.frame(avg_TF$chromvar)
avg_TF$trend <- ifelse(avg_TF[,1] >= avg_TF[,2] & avg_TF[,2] >= avg_TF[,3],1,ifelse(avg_TF[,1] <= avg_TF[,2] & avg_TF[,2]  <= avg_TF[,3],-1,0))
markers_CV$trend <- avg_TF$trend[match(markers_CV$gene,rownames(avg_TF))]
markers_CV %>% dplyr::filter(trend!=0) -> top10_TF
# markers_CV -> top10_TF
top10_TF %>%
  group_by(TF) %>%
  slice_max(n = 1, order_by = avg_diff) -> top10_TF
top10_TF  %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_diff) %>% arrange(desc(cluster)) -> top10_TF
print(top10_TF, n=10*3)

hm.integrated.cv <- subset(hm.integrated, subset = annotation %in% 'zone 3')
hm.integrated.cv@active.ident <- hm.integrated.cv$condition
DefaultAssay(hm.integrated.cv) <- 'chromvar'
p_dot_chromvar <- DotPlot(hm.integrated.cv,features = unique(top10_TF$gene)) + coord_flip() +
  scale_x_discrete(breaks=top10_TF$gene,labels=top10_TF$TF, limits=rev) +
  ylab('condition') +  xlab('TF') + ggtitle('Differential motif activity in zone 3 (ChromVAR)') +
  theme(axis.text.y = element_text(size = 14), axis.text.x=element_text(size = 14,angle = 30, vjust = 0.65))
ggsave(paste('p_dot_chromvar_','zone 3','.pdf',sep=''),plot =  p_dot_chromvar, width=25,height=22,units='cm')

# GO analysis on differential TF
library(org.Mm.eg.db)
mm <- org.Mm.eg.db
top10_TF$TF_mouse <- sub('\\(.*\\)','',top10_TF$TF)
motif_symbol <- c()
for(j in 1:dim(top10_TF)[1]){
  TF_gene_symbol <- unlist(strsplit(top10_TF$TF_mouse[j] ,split='\\:\\:')) %>% stringr::str_to_title()
  motif_symbol <- c(motif_symbol, structure(rep(top10_TF$trend[j],length(TF_gene_symbol)),names=TF_gene_symbol))
}
go <- lapply(c(1,-1),FUN=function(x){
  cg <- names(which(motif_symbol==x))
  cg_entrezid <- select(mm, keys = cg,
                        columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  library(limma)
  g <- limma::goana(cg_entrezid$ENTREZID,species='Mm')
  g <- g[order(g$P.DE,decreasing = F),]
})

## save go
df_go <- data.frame(go$'0'[1:50,],cluster=0)
for(j in 2:length(go)){
  df_go_tmp <- data.frame(go[[j]][1:50,],cluster=j-1)
  df_go <- rbind(df_go,df_go_tmp)
}
write.csv(df_go,file = 'GO_top50.csv')


#### mapping 
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
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat_sc$seurat_clusters1, prediction.assay = TRUE,
                                  dims = 1:30, weight.reduction = seurat_sp[["pca"]],  #1:30 pca 'cca'
                                  sd.weight = 1,
                                  n.trees=50,k.weight = 10) # 1：30, 1, 50, 10
seurat_sp[["predictions"]] <- predictions.assay

## Alluvial plot
sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
sp_to_sc$sp_cluster <-seurat_sp$annotation[rownames(sp_to_sc)]
sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels = cluster_unique)

cluster_unique2 <- c('zone 1','zone 2', 'zone 3','zone 3 Glul+','EC','EC/HSC chimera',"EC periportal", 'HSC/FB','KC','immune cells')
df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
df_long <- df_long[df_long$Freq >0,]
df_long$sc_cluster<- factor(df_long$sc_cluster,levels = c('2','0','1','6','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
sc_cluster_annotation <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                                   names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
df_long$sc_cluster2 <- factor(paste(sc_cluster_annotation[df_long$sc_cluster], ' (',df_long$sc_cluster,')',sep=''),
                              levels = paste(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'), ' (',c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),')',sep=''))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels=cluster_unique2)
df_long$sc_annotation <- factor(sc_cluster_annotation[df_long$sc_cluster],levels= cluster_unique2)

df_long$Total <- table(seurat_sp$annotation[rownames(sp_to_sc)])[as.character(df_long$sp_cluster)]
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
GeomAlluvium$default_aes$alpha <- 0.9
pmapping <- ggplot(data = df_long,
                   # aes(axis1 = sp_cluster, axis2 = sc_annotation , axis3=sc_cluster2, y= Freq)) +
                   aes(axis1 = sp_cluster, axis2 = sc_annotation, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/4) + scale_fill_manual(name='Annotation (ST)',values = color.use) +
  geom_stratum(width = 1/3) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),size=5) +
  theme_void() +   theme(text = element_text(size = 14), legend.text=element_text(size=14))
pmapping
# saveRDS(immune.combined,'Hep_integration_SoupX_reintegration_0805.rds')
ggsave('sp_sc_mapping_new.pdf',plot =  pmapping,  width=17,height=12,units='cm')

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
####### integration ATAC-seq with scRNA-seq,  dataset-wsie
project_name <- c('Young_1','Old_resistant_1','Old_sensitive_1','Old_resistant_2','Old_sensitive_2','Young_2')
pmapping_list <- list()
df_long_list <- list()
mapping_table_list <- list()
for(j in 1:length(project_name)) {
  dataset_j <- project_name[j]
  sc_j <- subset(seurat_sc, subset = orig.ident == dataset_j)
  DefaultAssay(sc_j) <- "RNA"
  sc_j <- NormalizeData(sc_j)
  sc_j <- ScaleData(sc_j,  features = rownames(sc_j))
  sc_j <- FindVariableFeatures(sc_j,nfeatures = 3000)

  atac_j <- subset(seurat_atac, subset = dataset == dataset_j)
  DefaultAssay(atac_j) <- "RNA"
  atac_j <- NormalizeData(atac_j)
  atac_j <- ScaleData(atac_j,  features = rownames(atac_j))
  atac_j <- FindVariableFeatures(atac_j,nfeatures = 3000)

  anchors <- FindTransferAnchors(reference = sc_j , query = atac_j, reference.assay = 'RNA', query.assay = 'RNA',normalization.method = 'LogNormalize',reduction = 'cca')
  predictions.assay <- TransferData(anchorset = anchors, refdata = sc_j$seurat_clusters1, prediction.assay = TRUE,
                                    weight.reduction = 'cca',dims = 1:20)
  # weight.reduction = atac_j[['lsi']],dims = 2:50)

  sp_to_sc <- apply(predictions.assay@data,2,FUN = function(cc){ind <- which(cc[1:(length(cc)-1)]==max(cc))[1]; names(cc)[ind]})
  sp_to_sc <- data.frame(sc_cluster=sp_to_sc)
  # sp_to_sc$sc_cluster <- paste('sc-',sp_to_sc$sc_cluster,sep = '')
  sp_to_sc$sp_cluster <-atac_j$annotation
  sp_to_sc$sp_cluster<- factor(sp_to_sc$sp_cluster,levels =cluster_unique2)

  mapping_table_list[[j]] <- table(sp_to_sc[,c(1,2)])

  df_long <- as.data.frame(table(sp_to_sc[,c(1,2)]))
  # df_long$sc_cluster<- factor(df_long$sc_cluster,levels =  paste('sc-',c('0','1','4','6','9','2','2_1','3','5','8','7','10'),sep=''))
  df_long$sp_cluster<- factor(df_long$sp_cluster,levels=cluster_unique2)
  df_long$sc_cluster<- factor(df_long$sc_cluster,levels = c('2','0','1','6','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
  sc_cluster_annotation <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                                     names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
  df_long$sc_cluster2 <- factor(paste(sc_cluster_annotation[df_long$sc_cluster], ' (',df_long$sc_cluster,')',sep=''),
                                levels = paste(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'), ' (',c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),')',sep=''))
  df_long$sc_annotation <- factor(sc_cluster_annotation[df_long$sc_cluster],levels= cluster_unique2)
  df_long <- df_long[df_long$Freq >0,]
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
p_sc <- DimPlot(seurat_sc, reduction = "umap",pt.size = 0.5,label=T, label.size = 10,dims = c(1,2))
pmapping_list_plus_sc <- list()
for(j in 1:6){
  pmapping_list_plus_sc[[j]] <- pmapping_list[[j]] + p_sc
}
mapping_table_sum <- matrix(0,nrow = length(unique(immune.combined$seurat_clusters1)),
                            ncol=length(cluster_unique2),
                            dimnames = list(c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),cluster_unique2))
mapping_table_0 <- mapping_table_sum
for(j in 1:6){
  mapping_table_j <- mapping_table_0
  mapping_table_j[rownames(mapping_table_list[[j]]),colnames(mapping_table_list[[j]])] <- mapping_table_list[[j]]
  mapping_table_sum <- mapping_table_sum + mapping_table_j
}

## Alluvial plot
df_long <- reshape2::melt(mapping_table_sum); names(df_long) <- c('sc_cluster','sp_cluster','Freq')
df_long <- df_long[df_long$Freq>0,]
df_long$sc_cluster<- factor(df_long$sc_cluster,levels = c('2','0','1','6','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
sc_cluster_annotation <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                                   names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
df_long$sc_cluster2 <- factor(paste(sc_cluster_annotation[df_long$sc_cluster], ' (',df_long$sc_cluster,')',sep=''),
                              levels = paste(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'), ' (',c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'),')',sep=''))
df_long$sp_cluster<- factor(df_long$sp_cluster,levels= cluster_unique2)
df_long$sc_annotation <- factor(sc_cluster_annotation[df_long$sc_cluster],levels= cluster_unique2)

df_long$Total <- table(seurat_sp$annotation[rownames(sp_to_sc)])[df_long$sp_cluster]
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
GeomAlluvium$default_aes$alpha <- 0.9
pmapping_atac  <- ggplot(data = df_long,
                   # aes(axis1 = sp_cluster, axis2 = sc_annotation , axis3=sc_cluster2, y= Freq)) +
                   aes(axis1 = sp_cluster, axis2 = sc_annotation, y= Freq)) +
  geom_alluvium(aes(fill = sp_cluster),width = 1/4) + scale_fill_manual(name='Annotation (ST)',values = color.use) +
  geom_stratum(width = 1/3) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),size=5) +
  theme_void() +   theme(text = element_text(size = 14), legend.text=element_text(size=14))
pmapping_atac
# saveRDS(immune.combined,'Hep_integration_SoupX_reintegration_0805.rds')
ggsave('atac_sc_mapping_new.pdf',plot =  pmapping_atac,  width=17,height=12,units='cm')


library(Seurat);library(Signac)
library(ggplot2)
setwd('/Users/weizhao/Documents/Hep_Seq/')
seurat_atac <- readRDS('ATAC-seq/integrated_ATAC.rds')
cluster_celltype_atac <- structure(c('zone 3','zone 2','zone 1','zone 1','EC','HSC','KC','Immune cells'),names=0:7)
seurat_atac$annotation <- cluster_celltype_atac[as.character(seurat_atac$seurat_clusters)]
tmp_atac <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(seurat_atac$dataset))
seurat_atac$condition <- tmp_atac[seurat_atac$dataset]
seurat_atac$annotation <- factor(seurat_atac$annotation, levels=c('zone 1','zone 2','zone 3','EC','HSC','KC','Immune cells','Undefined'))

DefaultAssay(seurat_atac) <- 'ATAC'
seurat_atac_cv <- subset(seurat_atac, subset= annotation %in% 'zone 3')
cov_plot <- CoveragePlot(
  object = seurat_atac_cv, group.by = 'condition',
  region = "Hnf4a",
  annotation = T,
  peaks = FALSE,extend.upstream = 00000, extend.downstream = 00000, show.bulk	= T,

)
cov_plot


library(sf);library(sp)
load(file='/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six.rda')
#### spatial plot
immune.combined.sct <- seurat_sp
plot_sp_cluster <- function(sample_single,cluster_idx,crop=c(0,0,1,1)){
  sample_single <- 'B2';
  cluster_idx <- 0:7
  crop=c(0,0,1,1)
  cluster_celltype_sp <- structure(c('zone 1','zone 2','zone 3','EC','HSC/FB','zone 3 Glul+','KC','EC periportal'),names=0:7)
  cluster_idx_anno <- cluster_celltype_sp[as.character(cluster_idx)]
  sample_idx <- which(sample_name %in% sample_single)
  cell_names <- Cells(immune.combined.sct)[grepl(sample_single,Cells(immune.combined.sct)) & (immune.combined.sct@meta.data$seurat_clusters) %in% cluster_idx]
  cluster_id <- (immune.combined.sct@meta.data$annotation)[which( Cells(immune.combined.sct) %in% cell_names)]
  cell_names <- sub(".*[A,B,C][1,2]-","",cell_names)
  roi_SpP <- roi_list[[sample_idx]]
  roi_SpP_subset <- roi_SpP[cell_names]
  cluster_unique <- c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','EC periportal','HSC/FB','KC','immune cells')
  color.use <- scales::hue_pal(l=75,c=150,h=c(0,360),h.start=10,direction = 1)(length(cluster_unique))
  color.use <- structure(color.use, names=cluster_unique)
  par(mar = c(2, 2, 2, 12))
  crop.margin.xmin <-  st_bbox(roi_SpP_subset)[1] + (st_bbox(roi_SpP_subset)[3]-st_bbox(roi_SpP_subset)[1])*crop[1]
  crop.margin.xmax <-  st_bbox(roi_SpP_subset)[1] + (st_bbox(roi_SpP_subset)[3]-st_bbox(roi_SpP_subset)[1])*crop[3]
  crop.margin.ymin <- st_bbox(roi_SpP_subset)[2] + (st_bbox(roi_SpP_subset)[4]-st_bbox(roi_SpP_subset)[2])*crop[2]
  crop.margin.ymax <- st_bbox(roi_SpP_subset)[2] + (st_bbox(roi_SpP_subset)[4]-st_bbox(roi_SpP_subset)[2])*crop[4]

  plot(roi_SpP_subset,col=color.use[cluster_id],main=sample_single,
       bg ='black',border=F,
       xlim = c(crop.margin.xmin, crop.margin.xmax),
       ylim = c(crop.margin.ymin, crop.margin.ymax));
  coord <- par("usr")
  legend(x = coord[2] * 1.02, y = coord[2]*0.4+coord[4]*0.6, legend=cluster_idx_anno,fill=color.use[cluster_idx_anno], cex=1.5,
         border = NA,box.lwd=0,bty = "n",horiz = F,
         xpd=TRUE, inset=c(-.05,0))
}
plot_sp_cluster('C2',0:7,crop=c(0.2,0.2,0.25,1))
plot_sp_cluster('B2',0:7,crop=c(0.3,0.57,0.55,0.58))

plot_sp_cluster('B2',0:7,crop=c(0.38,0.385,0.42,0.805)) # used to detemine the scale length 

plot_sp_cluster('C2',7,crop=c(0.,0.,1,1))

pdf(file="./CCI_analysis_mean/sp_B2.pdf",width = 8 ,height = 8)
plot_sp_cluster('B2',0:7,crop=c(0.28,0.35,0.57,0.85))
dev.off()
plot_sp_cluster('A2',0:7)
plot_sp_cluster('B1',0:7)
plot_sp_cluster('A1',c(0,2,6))
plot_sp_cluster('A1',c(5,6))
