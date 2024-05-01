#categorizing CCI both scRNA and spatial; write .csv
#####
load(file='/Users/weizhao/Documents/Hep_Seq/scRNA_results/CCI_analysis_mean/mean_lognorm/Hep_cellchat.rda')
target_names <- c('zone 3','zone 1','zone 2'); target_idx <- c(8,6,7)
# for(j in 1:3){
#   dir.create(paste('CCI_analysis_mean/',target_names[j],sep=''))
#   source('/Users/weizhao/Documents/Hep_Seq/subsetCommunication.R');
#   # excluding immune cells & EC/HSC chimera
#   my_netVisual_bubble(cellchat, sources.use = c(1,3,5,6:9), targets.use = target_idx[j],  comparison = c(1, 2, 3), angle.x = 45, font.size = 16,
#                       return.data = T,dir.to.save = paste('CCI_analysis_mean/',target_names[j],sep=''))
# }
j=1
gg12_j <- netVisual_bubble(cellchat, sources.use = c(1,3,5,6:9), targets.use = target_idx[j],  comparison = c(1, 2, 3), angle.x = 45, font.size = 16,
                           return.data = T)
df <- gg12_j$communication
pathway_name_j <- unique(df$pathway_name)
pathway_name_j <- pathway_name_j[!is.na(pathway_name_j)]
df_new <- df[!is.na(df$prob.original),]
df_tmp <- df_new[!is.na(df_new$pathway_name),]
df_tmp <- df_tmp[!duplicated(df_tmp$interaction_name_2),]
df_new$pathway_name <- df_tmp$pathway_name[match(df_new$interaction_name_2,df_tmp$interaction_name_2)]
comparison <- c(1,2,3)
df1_list <- list()
for(jj in 1:length(pathway_name_j)){
  df1 <- df_new[df_new$pathway_name == pathway_name_j[jj],]
  df1 <- df1[!is.na(df1$pathway_name),]
  interaction_name_2.jj <- unique(df1$interaction_name_2)
  for(i in 1:length(interaction_name_2.jj)){
    df <- df1[df1$interaction_name_2==interaction_name_2.jj[i],]
    group.names.jj <- unique(df$group.names)
    for(jjj in 1:length(group.names.jj)){
      dataset.jjj <- df$dataset[df$group.names==group.names.jj[jjj]]
      if(length(dataset.jjj) < length(comparison)){
        df_to_add <- df[which(df$group.names==group.names.jj[jjj])[rep(1,length(setdiff(names(cellchat@net),dataset.jjj)))],]
        df_to_add$dataset <-  setdiff(names(cellchat@net),dataset.jjj)
        df_to_add$source.target <- factor(paste(df_to_add$group.names,' (',df_to_add$dataset,')',sep=''),levels=levels(df$source.target))
        df_to_add$prob <- NA
        df_to_add$pval <- NA
        df_to_add$prob.original <- 0
        df1 <- rbind(df1,df_to_add)}
    }
  }
  df1_list[[jj]] <- df1
}
df2 <- do.call(rbind,df1_list)
df2 %>% arrange(interaction_name_2, group.names,dataset) -> df2
df2$interaction_cellpair <- paste(df2$interaction_name_2,df2$group.names)

cellpair_seq_all <- data.frame(interaction=df2$interaction_name_2,cellpair=df2$group.names, interaction_cellpair =df2$interaction_cellpair )
cellpair_seq_all <- cellpair_seq_all[!duplicated(cellpair_seq_all$interaction_cellpair),]
cellpair_seq_all$trend <- 0
cellpair_seq_all$category <- 'none'
tol.rate <- 0.05
len_interaction <- length(unique(cellpair_seq_all$interaction))
cellpair_seq_list <- list()
for(i in 1:len_interaction){
  interaction_i <- unique(cellpair_seq_all$interaction)[i]
  cellpair_tmp <- cellpair_seq_all[cellpair_seq_all$interaction==interaction_i,]
  len_cellpair <- length(unique(cellpair_tmp$cellpair))
  for(j in 1:len_cellpair){
    cellpair_j <- cellpair_tmp$cellpair[j]
    cellpair_seq <- cellpair_seq_all[cellpair_seq_all$interaction==interaction_i & cellpair_seq_all$cellpair==cellpair_j,]
    x_seq_j <- c(df2[df2$group.names==cellpair_j & df2$interaction_name_2 == interaction_i,'prob.original'])[c(3,2,1)]
    tmp1 <- diff(x_seq_j)
    for(jj in 1:2){
      tmp1[jj] <- tmp1[jj]/max(min(x_seq_j[jj:(jj+1)]),1e-16)
    }
    if(sum(tmp1 > tol.rate) >= 1 & sum(tmp1 < -tol.rate) == 0)
    {
      cellpair_seq$trend <- 1
      if(sum(tmp1 > tol.rate)>1) {cellpair_seq$category <- 'both' }
      else if ((tmp1 > tol.rate)[1] >0) {cellpair_seq$category <- 'aging'}
      else {cellpair_seq$category <- 'IR'}
    }
    if(sum(tmp1 < -tol.rate) >= 1 & sum(tmp1 > tol.rate) == 0)
    {
      cellpair_seq$trend <- -1
      if(sum(tmp1 < -tol.rate)  > 1) {cellpair_seq$category <- 'both' }
      else if ((tmp1 < -tol.rate)[1] >0) {cellpair_seq$category <- 'aging'}
      else {cellpair_seq$category <- 'IR'}
    }
    cellpair_seq_list[[length(cellpair_seq_list)+1]]  <- cellpair_seq
  }
}
# idx_plot <- which(cellpair_seq$trend_rank != 0 | cellpair_seq$trend != 0)
cellpair_seq_whole <- do.call(rbind, cellpair_seq_list)
cellpair_seq_whole$trend <- ifelse(cellpair_seq_whole$trend == 1, 'up', ifelse(cellpair_seq_whole$trend == -1, 'down','none'))
cellpair_seq_whole <- cellpair_seq_whole[ which( cellpair_seq_whole$trend != 'none'),]

cell_proximity_list <- readRDS("/Users/weizhao/Documents/Hep_Seq/CCI_analysis/Hep_cell_proximity.rds")
cellpair_seq <- data.frame(cellpair=cell_proximity_list[[1]]$CPScore_df$cellpair)
cellpair_seq$trend <- 0
cellpair_seq$trend_rank <- 0
cellpair_seq$category <- 'none'
cellpair_seq$category_rank <- 'none'
tol.rate <- 0.05
for(j in 1:dim(cellpair_seq)[1]){
  cellpair_j <- cellpair_seq$cellpair[j]
  x_seq_j <- c(cell_proximity_list[[1]]$CPScore_df$CPscore[cell_proximity_list[[1]]$CPScore_df$cellpair==cellpair_j],
               cell_proximity_list[[2]]$CPScore_df$CPscore[cell_proximity_list[[2]]$CPScore_df$cellpair==cellpair_j],
               cell_proximity_list[[3]]$CPScore_df$CPscore[cell_proximity_list[[3]]$CPScore_df$cellpair==cellpair_j])
  tmp1 <- diff(x_seq_j)
  for(jj in 1:2){
    tmp1[jj] <- tmp1[jj]/min(x_seq_j[jj:(jj+1)])
  }
  if(sum(tmp1 > tol.rate) >= 1 & sum(tmp1 < -tol.rate) ==0) {
    cellpair_seq$trend[j] <- 1
    if(sum(tmp1 > tol.rate)>1) {cellpair_seq$category[j]  <- 'both' }
    else if ((tmp1 > tol.rate)[1] >0) {cellpair_seq$category[j]  <- 'aging'}
    else {cellpair_seq$category[j]  <- 'IR'}
  }
  if(sum(tmp1 < -tol.rate) >= 1 & sum(tmp1 > tol.rate) == 0) {
    cellpair_seq$trend[j] <- -1
    if(sum(tmp1 < -tol.rate)  > 1) {cellpair_seq$category[j] <- 'both' }
    else if ((tmp1 < -tol.rate)[1] >0) {cellpair_seq$category[j] <- 'aging'}
    else {cellpair_seq$category[j] <- 'IR'}
  }
  x_seq_rank_j <- c(which(cell_proximity_list[[1]]$CPScore_df$cellpair==cellpair_j),
                    which(cell_proximity_list[[2]]$CPScore_df$cellpair==cellpair_j),
                    which(cell_proximity_list[[3]]$CPScore_df$cellpair==cellpair_j))
  tmp2 <- diff(x_seq_rank_j)
  for(jj in 1:2){
    tmp2[jj] <- tmp2[jj]/min(x_seq_rank_j[jj:(jj+1)])
  }
  if(sum(tmp2 > tol.rate) >= 1 & sum(tmp2 < -tol.rate) ==0) {
    cellpair_seq$trend_rank[j] <- 1
    if(sum(tmp2 > tol.rate)>1) {cellpair_seq$category_rank[j]  <- 'both' }
    else if ((tmp2 > tol.rate)[1] >0) {cellpair_seq$category_rank[j]  <- 'aging'}
    else {cellpair_seq$category_rank[j]  <- 'IR'}
  }
  if(sum(tmp2 < -tol.rate) >= 1 & sum(tmp2 > tol.rate) == 0) {
    cellpair_seq$trend_rank[j] <- -1
    if(sum(tmp2 < -tol.rate)  > 1) {cellpair_seq$category_rank[j] <- 'both' }
    else if ((tmp2 < -tol.rate)[1] >0) {cellpair_seq$category_rank[j] <- 'aging'}
    else {cellpair_seq$category_rank[j] <- 'IR'}
  }
}
# cellpair_seq_whole <- cellpair_seq_whole[!(cellpair_seq_whole$cellpair== "Undefined -> CV"),]
cluster_celltype_sp <- structure(c(0:6),names=c('zone 1','zone 2','zone 3','EC','HSC/FB','zone 3 Glul+','KC'))
cellpair_unique <- unique(cellpair_seq_whole$cellpair)
cellpair_idx <- sapply(stringr::str_split(cellpair_unique,' -> '), FUN = function(x){
  y <- paste(cluster_celltype_sp[x][1],cluster_celltype_sp[x][2], sep='--')
  z <- paste(cluster_celltype_sp[x][2],cluster_celltype_sp[x][1], sep='--')
  r <- unique(match(c(y,z),cellpair_seq$cellpair))
  r <- r[!is.na(r)]
})
tmp <- structure(cellpair_idx,names=cellpair_unique)
names(cellpair_seq) <- c( "sp_cellpair","sp_trend","sp_trend_rank","sp_category","sp_category_rank")
cellpair_seq_whole <- cbind(cellpair_seq_whole, cellpair_seq[tmp[cellpair_seq_whole$cellpair],])
cellpair_seq_whole$sp_trend <- ifelse(cellpair_seq_whole$sp_trend == 1, 'up', ifelse(cellpair_seq_whole$sp_trend == -1, 'down','none'))
cellpair_seq_whole$sp_trend_rank <- ifelse(cellpair_seq_whole$sp_trend_rank == 1, 'up', ifelse(cellpair_seq_whole$sp_trend_rank == -1, 'down','none'))
cellpair_seq_whole <- cellpair_seq_whole[,c(1:7,9,8,10)]
# write.csv(cellpair_seq_whole,file='CCI_analysis_mean/CCI_comparison_category.csv')
# CCI category visualization
#####
library(Seurat)
library(ggplot2)
seurat_sc <- readRDS(file='/Users/weizhao/Documents/Hep_Seq/Hep16_outs/Hep_integration_SoupX_reintegration_0805.rds');library(Seurat)
seurat_sc$annotation <- factor(seurat_sc$annotation,levels=c('zone 1','zone 3','zone 3 Glul+','zone 2', 'EC','EC/HSC chimera','HSC/FB','KC','immune cells'))
seurat_sc$annotation[which(seurat_sc@meta.data$sub.cluster.EC %in% c('EC_1'))] <- 'EC/HSC chimera'
seurat_sc$seurat_clusters1 <- factor(seurat_sc$seurat_clusters1,levels=c('2','10','3','4','5','3-1','6','1','0','7','8','9','9-1','12','11','11-1'))
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(seurat_sc$orig.ident))
seurat_sc$condition <- factor(tmp[seurat_sc$orig.ident],levels=c('Young','Old sensitive','Old resistant'))
tmp2 <- structure(c('YS','OR','OS','OR','OS','YS'),names=unique(seurat_sc$orig.ident))
seurat_sc$condition_abbr <- factor(tmp2[seurat_sc$orig.ident],levels=c('YS','OS','OR'))
Idents(seurat_sc) <- seurat_sc$annotation
DefaultAssay(seurat_sc) <- 'RNA'
# p1 <- Seurat::DotPlot(seurat_sc,features = 'Hgf',idents = 'HSC/FB', group.by  = 'condition') + ggtitle( c('Hgf in HSC'))
# p2 <- Seurat::DotPlot(seurat_sc,features = 'Met',idents = 'zone 3', group.by  = 'condition') + ggtitle( c('Met in zone 3'))
# p1 + p2
# 
# VlnPlot(seurat_sc,features = c('Foxe1', 'Foxd2', 'Hnf1a', 'Hnf1b' ),idents = 'zone 3', group.by  = 'condition',ncol = 2)

# Seurat::DotPlot(seurat_sc,features = 'Hgf', split.by  = 'condition',cols = c('blue','yellow','red'))

### visualization
# cellpair_seq_whole <- read.csv(file='CCI_analysis_mean/CCI_comparison_category.csv')
setwd('/Users/weizhao/Documents/Hep_Seq/scRNA_results/CCI_analysis_mean/')

# NicheNet analysis
# NicheNet analysis
#####
library(nichenetr)
library(tidyverse)
options(timeout = 600)
organism = "mouse"

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
}

lr_network = lr_network %>% distinct(from, to)

library(Seurat);library(ggplot2)
seurat_sc <- readRDS(file='/Users/weizhao/Documents/Hep_Seq/Hep16_outs/Hep_integration_SoupX_reintegration_0805.rds');library(Seurat)
cluster_celltype <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                              names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
seurat_sc$annotation <- cluster_celltype[as.character(seurat_sc$seurat_clusters1)]
seurat_sc$annotation <- factor(seurat_sc$annotation,levels=c('zone 1','zone 3','zone 3 Glul+','zone 2', 'EC','EC/HSC chimera','HSC/FB','KC','immune cells'))
seurat_sc$annotation[which(seurat_sc@meta.data$sub.cluster.EC %in% c('EC_1'))] <- 'EC/HSC chimera'
seurat_sc$seurat_clusters1 <- factor(seurat_sc$seurat_clusters1,levels=c('2','10','3','4','5','3-1','6','1','0','7','8','9','9-1','12','11','11-1'))
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(seurat_sc$orig.ident))
seurat_sc$condition <- factor(tmp[seurat_sc$orig.ident],levels=c('Young','Old sensitive','Old resistant'))
tmp2 <- structure(c('YS','OR','OS','OR','OS','YS'),names=unique(seurat_sc$orig.ident))
seurat_sc$condition_abbr <- factor(tmp2[seurat_sc$orig.ident],levels=c('YS','OS','OR'))
Idents(seurat_sc) <- seurat_sc$annotation
DefaultAssay(seurat_sc) <- 'RNA'

markers_CV <-  read.csv('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/Diff_chromvar/Diff_TF_fromATAC_zone 3.csv')
markers_CV %>% dplyr::filter(avg_diff>0.1 & cluster=='Young' & p_val_adj < 0.01) -> markers_CV
# hm.integrated <- readRDS('/Users/weizhao/Documents/Hep_Seq/ATAC-seq/ATAC_integrated_chromvar.rds')
avg_TF <- AverageExpression(subset(hm.integrated,subset = annotation =='zone 3' ), features = unique(markers_CV$gene), group.by = 'condition',assays = 'chromvar')
avg_TF <- as.data.frame(avg_TF$chromvar)
avg_TF$trend <- ifelse(avg_TF[,1] >= avg_TF[,2] & avg_TF[,2] >= avg_TF[,3],1,0)
markers_CV$trend <- avg_TF$trend[match(markers_CV$gene,rownames(avg_TF))]
markers_CV %>% dplyr::filter(trend==1) -> markers_CV
geneset_oi <- markers_CV$TF
geneset_oi <- gsub("\\(.*\\)",'',geneset_oi)
geneset_oi <- unlist(strsplit(geneset_oi,split='\\:\\:')) %>% str_to_title() %>% .[. %in% rownames(seurat_sc)]
# geneset_oi <- nichenetr::convert_human_to_mouse_symbols(geneset_oi)

ligands = lr_network %>% pull(from) %>% unique()
potential_ligands = intersect(ligands, unique(df2$ligand))

condition <- c('Young','Old sensitive','Old resistant')
ligand_activities_list <- list()
for(j in 1:3){
  seurat_CV_subset <- subset(seurat_sc, subset = annotation == 'zone 3' & seurat_sc$condition == condition[j])
  expressed_genes_receiver <- names(which(rowSums(seurat_CV_subset@assays$RNA@counts > 0)/dim(seurat_CV_subset)[2] > 0.3))
  print(length(expressed_genes_receiver))
  print(dim(seurat_CV_subset))
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] %>% .[. %in% rownames(seurat_sc)]
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities$condition <- condition[j]
  ligand_activities_list[[j]] <- ligand_activities
}
ligand_activities_all <- do.call(rbind, ligand_activities_list)
ligand_activities_metrics <- list()
for(jj in 1:3){
  metric_jj <- names(ligand_activities_all)[2:4][jj]
  ligand_activities_auroc <- reshape2::dcast(ligand_activities_all[,c('test_ligand',metric_jj,'condition')],formula = test_ligand ~ condition, value.var =metric_jj)
  ligand_activities_auroc$trend <-ifelse(ligand_activities_auroc[,2] <= ligand_activities_auroc[,3] & ligand_activities_auroc[,2] <= ligand_activities_auroc[,3] & ligand_activities_auroc[,3] < ligand_activities_auroc[,4],1,0)
  print(ligand_activities_auroc[ligand_activities_auroc$trend==1,])
  ligand_activities_metrics[[jj]] <- ligand_activities_auroc
}
# saveRDS(ligand_activities_metrics, file='./CCI_analysis_mean/mean_SCT/ligand_activities_nichenet_zone2.rds')
#ligand_activities_metrics <- readRDS('./CCI_analysis_mean/mean_SCT/ligand_activities_nichenet.rds')
library(ComplexHeatmap);library(circlize)
tol.rate <- 0.01
heat_mtx <- ligand_activities_metrics[[3]][,c(4,3,2)]
names(heat_mtx) <- c('YS','OS','OR')
rownames(heat_mtx) <- ligand_activities_metrics[[3]][,c(1)]
df_weight_cv_wide <- heat_mtx
df_weight_cv_wide$sign1 <- ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/abs(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))> tol.rate,1,
                                  ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/abs(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))< -tol.rate,-1,0))

df_weight_cv_wide$sign2 <- ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/abs(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))> tol.rate,1,
                                  ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/abs(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))< -tol.rate,-1,0))
df_weight_cv_wide$sign1_x_sign2 <- df_weight_cv_wide$sign1*df_weight_cv_wide$sign2
df_weight_cv_wide$sign1_plus_sign2 <- df_weight_cv_wide$sign1+df_weight_cv_wide$sign2
sign_category <- structure(c('up_both','up_one','down_both','down_one'),names=c(2,1,-2,-1))
df_weight_cv_wide$category <- 'others'
df_weight_cv_wide$category[df_weight_cv_wide$sign1_plus_sign2==2] <- 'up_both'
df_weight_cv_wide$category[df_weight_cv_wide$sign1_plus_sign2==-2] <- 'down_both'
df_weight_cv_wide$category[df_weight_cv_wide$sign1==1 & df_weight_cv_wide$sign2==0] <- 'up_aging'
df_weight_cv_wide$category[df_weight_cv_wide$sign1==0 & df_weight_cv_wide$sign2==1] <- 'up_IR'
df_weight_cv_wide$category[df_weight_cv_wide$sign1==-1 & df_weight_cv_wide$sign2==0] <- 'down_aging'
df_weight_cv_wide$category[df_weight_cv_wide$sign1==0 & df_weight_cv_wide$sign2==-1] <- 'down_IR'
df_weight_cv_wide$total_weight <- df_weight_cv_wide$YS+df_weight_cv_wide$OS+df_weight_cv_wide$OR
df_weight_cv_wide$category <- factor(df_weight_cv_wide$category,levels = c('up_both','up_aging','up_IR','down_IR','down_aging','down_both','others'))
df_weight_cv_wide %>% arrange(category,-total_weight) -> df_weight_cv_wide
labels_tmp <-  c('up_both','up_aging','up_IR','down_IR','down_aging','down_both','others')
labels_idx <- which(labels_tmp %in% df_weight_cv_wide$category )
col_tmp <- colorRampPalette(c("green",'grey', "red"))(7)[c(1:3,5:7,4)]
ht_opt$ROW_ANNO_PADDING= unit(1, "mm")
row_ha = rowAnnotation(trend= anno_block(gp = gpar(fill =col_tmp[labels_idx] ,alpha=0.8), labels = labels_tmp[labels_idx],
                                         labels_gp = gpar(col = "black", fontsize = 20, fontface=2),labels_rot=0,width=unit(45,'mm'),show_name=T),
                       total_activity = anno_barplot(df_weight_cv_wide$total_weight,  axis_param=list(gp=gpar(fontsize = 16)), 
                                                     width=unit(30,'mm')),
                       annotation_name_gp=gpar(fontsize = 16),gap=unit(1,'mm'))
col_fun = colorRamp2(c(-1.5, 1.5), c("white","orange"))
df_weight_cv_matrix <- as.matrix(df_weight_cv_wide[,c(1,2,3)])
rownames(df_weight_cv_matrix) <- rownames(df_weight_cv_wide)

h1 <- ComplexHeatmap::Heatmap(as.data.frame(t(scale(t(df_weight_cv_matrix)))),
                              #as.data.frame(t(scale(t(heat_mtx)))),
                              cluster_rows = F,cluster_columns = F,col = col_fun,
                              row_names_side = 'left', right_annotation = row_ha, #right_annotation = row_ha,
                              column_names_rot = 0,name = 'scaled ligand activity',
                              row_split=(df_weight_cv_wide$category),
                              row_title=NULL, row_names_gp = gpar(fontsize=20), column_names_gp = gpar(fontsize=20)
)
pdf(file="CCI_ligandactivity_zone3.pdf",width = 9 ,height = 25)
draw(h1,column_title='ligand signaling activity in zone 3',column_title_gp= gpar(fontsize=20))
dev.off()

# 
# library(nichenetr)
# library(tidyverse)
# 
# weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
# ligand_tf_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final_mouse.rds"))
# 
# lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
# sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_mouse_21122021.rds"))
# gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_mouse_21122021.rds"))
# 
# ligands_all = "Hgf" # this can be a list of multiple ligands if required
# targets_all = c('Onecut1','Onecut2','Onecut3')
# 
# active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)
# active_signaling_network$sig <- active_signaling_network$sig#[active_signaling_network$sig$weight>1.5,]
# active_signaling_network$gr <- active_signaling_network$gr#[active_signaling_network$gr$weight>1,]
# active_signaling_network$sig <- active_signaling_network$sig[-which(active_signaling_network$sig$from == 'Hgf' & active_signaling_network$sig$to != 'Met'),]
# active_signaling_network$gr <- active_signaling_network$gr[-which(active_signaling_network$gr$from == 'Hgf'),]
# 
# # For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
# active_signaling_network_min_max = active_signaling_network
# active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
# active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
# 
# graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
# graph_min_max$nodes_df$width <- 0.4
# graph_min_max$global_attrs$value[graph_min_max$global_attrs$attr == 'fontsize'][1] <- 6
# # To render the graph: uncomment following line of code
# DiagrammeR::render_graph(graph_min_max, layout = "tree")
# tmp <- DiagrammeR::render_graph(graph_min_max, layout = "kk",width = 1000, height=1000)#,output = "visNetwork")
# library(DiagrammeR)
# library(DiagrammeRsvg)
# library(magrittr)
# library(rsvg)
# tmp = DiagrammeRsvg::export_svg(tmp)
# tmp = charToRaw(tmp) # flatten
# rsvg::rsvg_pdf(tmp, "g.pdf") # saved graph as png in current working directory
# 
# seurat_sc = AddModuleScore(seurat_sc, list(unique(c(active_signaling_network$sig$to, active_signaling_network$sig$from[active_signaling_network$sig$from != 'Hgf'], active_signaling_network$gr$from, active_signaling_network$gr$to))), assay = "RNA", name = "Hgf_signaling_activity")
# VlnPlot(seurat_sc, 'Hgf_signaling_activity1', idents = 'zone 3',group.by = 'condition')
# VlnPlot(seurat_sc, c('Met','Stat3'), idents = 'zone 3',group.by = 'condition')
# 
# 
# VlnPlot(seurat_sc, c('Met','Stat3'), idents = 'zone 3',group.by = 'condition') + geom_boxplot(position=position_dodge(1))
# VlnPlot(subset(seurat_sc, subset = annotation == 'zone 3' & Met > 0),features = 'Stat3', group.by = 'condition',assay = 'RNA')
# 
# quantile(subset(seurat_sc,subset= annotation == 'zone 3' & condition=='Young')@assays$RNA@data['Stat3',])
# quantile(subset(seurat_sc,subset= annotation == 'zone 3' & condition=='Old sensitive')@assays$RNA@data['Stat3',])
# quantile(subset(seurat_sc,subset= annotation == 'zone 3' & condition=='Old resistant')@assays$RNA@data['Stat3',])
# DotPlot(subset(seurat_sc, subset = annotation == 'zone 3' & Met > 0.1),features = 'Stat3', group.by = 'condition',assay = 'RNA')
# 
# VlnPlot(seurat_sc, 'Met', idents = 'zone 3',group.by = 'seurat_clusters')
# 
# # seurat clusters for zone 3  -- '10','3','4','5'
# AverageExpression(subset(seurat_sc, subset = annotation == 'zone 3' & Met > 0),
#                   assays = 'RNA',group.by = 'condition', features = 'Stat3',slot = 'data')
# FeatureScatter(subset(seurat_sc, subset = annotation == 'zone 3'),'Met','Stat3')
# 
# 
# 
