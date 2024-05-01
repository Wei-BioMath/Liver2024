library(dplyr)
library(Seurat)
library(sctransform)
library(future)
library(ggplot2)
library(patchwork)
# check the current active plan
plan("multicore", workers = 12)
options(future.globals.maxSize= 5000*1024^2)
options(future.seed=TRUE)
set.seed(1234)
setwd('/Users/weizhao/Documents/Hep_Seq/scRNA_results/')
load(file='/Users/weizhao/Documents/Hep_Seq/scRNA_results/CCI_analysis_mean/mean_lognorm/Hep_cellchat.rda')
dir.create('CCI_analysis_mean/category_split_multiple')
dir.create('CCI_analysis_mean/category_multiple')
dir.create('CCI_analysis_mean/category_filter_combine')

seurat_sc <- readRDS('/Users/weizhao/Documents/Hep_Seq/Hep16_outs/Hep_integration_SoupX_reintegration_0805.rds')
cluster_celltype_sc <- structure(c('zone 1',rep('zone 2',5), rep('zone 3',4),'zone 3 Glul+','EC','EC/HSC chimera','HSC/FB','KC','immune cells'),
                                 names=c('2','6','1','0','7','8','10','3','4','5','3-1','9','9-1','12','11','11-1'))
seurat_sc$annotation <- cluster_celltype_sc[as.character(seurat_sc$seurat_clusters1)]
tmp <- structure(c('Young','Old resistant','Old sensitive','Old resistant','Old sensitive','Young'),names=unique(seurat_sc$orig.ident))
seurat_sc$condition <- tmp[seurat_sc$orig.ident]
cluster_unique <- c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','EC periportal','HSC/FB','KC','immune cells')
seurat_sc$annotation <- factor(seurat_sc$annotation, levels=cluster_unique)
seurat_sc$condition <- factor(seurat_sc$condition, levels=c('Young','Old sensitive','Old resistant'))
tmp2 <- structure(c('YS','OR','OS','OR','OS','YS'),names=unique(seurat_sc$orig.ident))
seurat_sc$condition_abbr <- factor(tmp2[seurat_sc$orig.ident],levels=c('YS','OS','OR'))
Idents(seurat_sc) <- seurat_sc$annotation

data_tmp <- seurat_sc@assays$RNA@data
max_condition <- sapply(1:3, function(x){
  return(max(cellchat_list[[x]]@data.signaling))
},simplify = T)
for(j in 1:3){
  condition_j <- names(cellchat_list)[j]
  idx_j <- which(seurat_sc$condition_abbr == condition_j)
  data_tmp[,idx_j] <- data_tmp[,idx_j]/max_condition[j]
}
seurat_sc[['RNA2']] <- CreateAssayObject(data = data_tmp)
DefaultAssay(seurat_sc) <- 'RNA2'

target_names <- rownames(cellchat_list$YS@netP$prob[,,1]); target_idx <- 1:9
for(j in 1:9){
gg12_j <- netVisual_bubble(cellchat, sources.use = c(1,3,5,6:9), targets.use = target_idx[j],  comparison = c(1, 2, 3), angle.x = 45, font.size = 16,
                           return.data = T)
  # gg12_j <- netVisual_bubble(cellchat, sources.use = 1:9, targets.use = target_idx[j],  comparison = c(1, 2, 3), angle.x = 45, font.size = 16,
  #                            return.data = T)
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
df2 <- df2[df2$interaction_name_2 != 'Nampt  - Insr',]

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
# cellpair_seq_whole <- cellpair_seq_whole[ which( cellpair_seq_whole$trend != 'none'),]

CCI_heatmap_interaction_pair_comparison <- function(cellpair_oi, top_n=5){
  # interaction_cellpair_oi <- cellpair_seq_whole[cellpair_seq_whole$cellpair %in% cellpair_oi,] #
  interaction_cellpair_oi <- cellpair_seq_whole[ cellpair_seq_whole$category!= 'none' & cellpair_seq_whole$cellpair %in% cellpair_oi,] #
  if(dim(interaction_cellpair_oi)[1]>0){
    df_weight_cv_wide <- reshape2::dcast(df2[df2$interaction_cellpair %in% interaction_cellpair_oi$interaction_cellpair,],formula = interaction_name_2 ~ dataset, value.var = 'prob.original')
    df_weight_cv_wide[is.na(df_weight_cv_wide)] <- 0
    df_weight_cv_wide <- df_weight_cv_wide[,c('interaction_name_2', 'YS', 'OS', 'OR')]
    colnames(df_weight_cv_wide) <- c('name', 'YS', 'OS', 'OR')
    tol.rate <- 0.05
    df_weight_cv_wide$sign1 <- ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))> tol.rate,1,
                                      ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))< -tol.rate,-1,0))

    df_weight_cv_wide$sign2 <- ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))> tol.rate,1,
                                      ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))< -tol.rate,-1,0))
    # df_weight_cv_wide <- df_weight_cv_wide[df_weight_cv_wide$sign1*df_weight_cv_wide$sign2 >=0,]
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
    df_weight_cv_wide %>% dplyr::filter(total_weight > 0) %>% arrange(category,-total_weight) -> df_weight_cv_wide
    if(!is.null(top_n) & top_n >0.5){
      df_weight_cv_wide %>% group_by(category)  %>%  slice_max(n = round(top_n), order_by = total_weight) -> df_weight_cv_wide
    }
    df_weight_cv_matrix <- as.matrix(df_weight_cv_wide[,c(2,3,4)])
    rownames(df_weight_cv_matrix) <- df_weight_cv_wide$name
    library(ComplexHeatmap);library(circlize)
    col_tmp <- colorRampPalette(c("green",'grey', "red"))(7)[c(1:3,5:7,4)]
    labels_tmp <-  c('up_both','up_aging','up_IR','down_IR','down_aging','down_both','others')
    labels_idx <- which(labels_tmp %in% df_weight_cv_wide$category )
    # row_ha = rowAnnotation(foo = anno_block(gp = gpar(fill =col_tmp[labels_idx] ,alpha=0.8), labels = labels_tmp[labels_idx],
    #                                         labels_gp = gpar(col = "black", fontsize = 20, fontface=2),labels_rot=0,width=unit(45,'mm')))
    ht_opt$ROW_ANNO_PADDING= unit(1, "mm")
    row_ha = rowAnnotation(trend= anno_block(gp = gpar(fill =col_tmp[labels_idx] ,alpha=0.8), labels = labels_tmp[labels_idx],
                                             labels_gp = gpar(col = "black", fontsize = 20, fontface=2),labels_rot=0,width=unit(45,'mm'),show_name=T),
                           total_strength = anno_barplot(df_weight_cv_wide$total_weight,  axis_param=list(gp=gpar(fontsize = 16)), 
                                                         width=unit(30,'mm')),
                           annotation_name_gp=gpar(fontsize = 16),gap=unit(1,'mm'))
    col_fun = colorRamp2(c(-1.5, 1.5), c("white","blue"))
    h1 <- ComplexHeatmap::Heatmap(as.data.frame(t(scale(t(df_weight_cv_matrix)))),
                                  cluster_rows = F,cluster_columns = F, col = col_fun,
                                  row_names_side = 'left', right_annotation = row_ha, #right_annotation = row_ha,
                                  column_names_rot = 0,name = 'scaled \n strength',
                                  row_split=(df_weight_cv_wide$category), 
                                  row_title=NULL,     row_names_gp = gpar(fontsize = 20),     column_names_gp = gpar(fontsize = 20),
                                  heatmap_legend_param = list(labels_gp = gpar(fontsize = 16),title_gp=gpar(fontsize = 20)),
    )
    filename_save <- sub(" -> ","_", cellpair_oi)
    filename_save <- gsub(" ","_", filename_save)
    filename_save <- gsub("/","-", filename_save)
    pdf(file= paste('CCI_analysis_mean/category_split_multiple/',filename_save,'.pdf',sep=''),width = 10 ,height =1.2 + 25*dim(df_weight_cv_matrix)[1]/55)
    draw(h1,column_title=cellpair_oi,column_title_gp = gpar(fontsize = 20))
    dev.off()
    ## plot violins
    idx_down <- (which(df_weight_cv_wide$sign1_plus_sign2 < 0))
    if(length(idx_down) >0 ){
      lr_list <- as.character(df_weight_cv_wide$name[idx_down])
      DB_subset <- CellChatDB.mouse$interaction[match(lr_list,CellChatDB.mouse$interaction$interaction_name_2) ,]
      gene_list <- lapply(lr_list,function(x){
        y <- strsplit(x,split=' *- *')
        ligand <- y[[1]][1]
        receptor <- unlist(strsplit( sub("\\)","", sub("\\(","",y[[1]][2])), split = '\\+'))
        receptor <- CellChatDB.mouse$interaction[match(x,CellChatDB.mouse$interaction$interaction_name_2) ,]$receptor
        if(sum(grep("_",receptor)!=0)>0){
          receptor22 <- paste(CellChatDB.mouse$complex[receptor,],collapse = '_')
          receptor22 <- sub("_*$",'',receptor22)
          receptor <- unlist(strsplit(receptor22, split = '_'))
        }
        return(list(ligand,receptor))
      })
      sender <- unlist(strsplit(cellpair_oi,split=' -> '))[1]
      receiver <- unlist(strsplit(cellpair_oi,split=' -> '))[2]
      for(ii in 1:length(idx_down)){
        p_ii_1 <- VlnPlot(seurat_sc,features = gene_list[[ii]][[1]],idents = sender, group.by = 'condition_abbr') + plot_annotation(title=paste('ligand in',sender)) + theme(legend.position="none")
        p_ii_2 <-  VlnPlot(seurat_sc,gene_list[[ii]][[2]],idents = receiver, group.by = 'condition_abbr') +  plot_annotation(title = paste('receptor in',receiver)) + theme(legend.position="none")
        if(length(gene_list[[ii]][[2]])>1) {
          p_ii_2 <- p_ii_2+ patchwork::plot_layout(nrow = 1)
        }
        p_ii <- p_ii_1 + p_ii_2  + plot_layout(ncol=2) + plot_annotation(title=paste('ligand in',sender, '&',paste('receptor in',receiver)) ) +
          theme(text=element_text(size=14), #change font size of all text
                axis.text=element_text(size=14,color = 'black'), #change font size of axis text
                axis.text.x = element_text(size=14,color = 'black'),
                axis.title=element_text(size=14), #change font size of axis titles
                legend.text=element_text(size=14), #change font size of legend text
                legend.title=element_text(size=14))#change font size of legend title
        ggsave(filename = paste(sub("/","_",sender),sub("/","_",receiver),gene_list[[ii]][[1]],paste( gene_list[[ii]][[2]],collapse = '_'),'.pdf',sep='_'),
               plot = p_ii, width = 3+3*length(gene_list[[ii]][[2]]) ,height = 5)
      }
    }
  }
}
for(cellpair_oi in unique(cellpair_seq_whole$cellpair)){
  CCI_heatmap_interaction_pair_comparison(cellpair_oi, top_n = 5)
}
CCI_heatmap_interaction_pair_comparison_filter_combine <- function(cellpair_oi, thresh=0.1, top_n=5){
  # cellpair_oi <- unique(cellpair_seq_whole$cellpair)
  # top_n= 5
  # thresh =0.1
  interaction_cellpair_oi <- cellpair_seq_whole[ cellpair_seq_whole$category!= 'none' & cellpair_seq_whole$cellpair %in% cellpair_oi,] #
  if(dim(interaction_cellpair_oi)[1]>0){
    df_tmp <- df2[df2$interaction_cellpair %in% interaction_cellpair_oi$interaction_cellpair,]
    df_tmp$interaction_name_3 <- paste("(",df_tmp$source,") ",df_tmp$interaction_name_2,sep="")
    df_weight_cv_wide <- reshape2::dcast(df_tmp,formula = interaction_name_3 ~ dataset, value.var = 'prob.original')
    df_weight_cv_wide[is.na(df_weight_cv_wide)] <- 0
    df_weight_cv_wide <- df_weight_cv_wide[,c('interaction_name_3', 'YS', 'OS', 'OR')]
    colnames(df_weight_cv_wide) <- c('name', 'YS', 'OS', 'OR')
    tol.rate <- 0.05
    df_weight_cv_wide$sign1 <- ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))> tol.rate,1,
                                      ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))< -tol.rate,-1,0))
    
    df_weight_cv_wide$sign2 <- ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))> tol.rate,1,
                                      ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))< -tol.rate,-1,0))
    # df_weight_cv_wide <- df_weight_cv_wide[df_weight_cv_wide$sign1*df_weight_cv_wide$sign2 >=0,]
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
    filename_save <- sub(" -> ","_", unique(df_tmp$target)[1])
    filename_save <- gsub(" ","_", filename_save)
    filename_save <- gsub("/","-", filename_save)
    write.csv(df_weight_cv_wide, file=paste('CCI_category_',filename_save,'.csv',sep=''))
    df_weight_cv_wide %>% dplyr::filter(total_weight > thresh) %>% arrange(category,-total_weight) -> df_weight_cv_wide
    if(dim(df_weight_cv_wide)[1]>0){
    if(!is.null(top_n) & top_n >0.5){
      df_weight_cv_wide %>% group_by(category)  %>%  slice_max(n = round(top_n), order_by = total_weight) -> df_weight_cv_wide
    }
    df_weight_cv_matrix <- as.matrix(df_weight_cv_wide[,c(2,3,4)])
    rownames(df_weight_cv_matrix) <- df_weight_cv_wide$name
    library(ComplexHeatmap);library(circlize)
    col_tmp <- colorRampPalette(c("green",'grey', "red"))(7)[c(1:3,5:7,4)]
    labels_tmp <-  c('up_both','up_aging','up_IR','down_IR','down_aging','down_both','others')
    labels_idx <- which(labels_tmp %in% df_weight_cv_wide$category )
    # row_ha = rowAnnotation(foo = anno_block(gp = gpar(fill =col_tmp[labels_idx] ,alpha=0.8), labels = labels_tmp[labels_idx],
    #                                         labels_gp = gpar(col = "black", fontsize = 20, fontface=2),labels_rot=0,width=unit(45,'mm')))
    ht_opt$ROW_ANNO_PADDING= unit(1, "mm")
    row_ha = rowAnnotation(trend= anno_block(gp = gpar(fill =col_tmp[labels_idx] ,alpha=0.8), labels = labels_tmp[labels_idx],
                                             labels_gp = gpar(col = "black", fontsize = 20, fontface=2),labels_rot=0,width=unit(45,'mm'),show_name=T),
                           total_strength = anno_barplot(df_weight_cv_wide$total_weight,  axis_param=list(gp=gpar(fontsize = 16)), 
                                                         width=unit(30,'mm')),
                           annotation_name_gp=gpar(fontsize = 16),gap=unit(1,'mm'))
    col_fun = colorRamp2(c(-1.5, 1.5), c("white","blue"))
    h1 <- ComplexHeatmap::Heatmap(as.data.frame(t(scale(t(df_weight_cv_matrix)))),
                                  cluster_rows = F,cluster_columns = F, col = col_fun,
                                  row_names_side = 'left', right_annotation = row_ha, #right_annotation = row_ha,
                                  column_names_rot = 0,name = 'scaled \n strength',
                                  row_split=(df_weight_cv_wide$category), 
                                  row_title=NULL,     row_names_gp = gpar(fontsize = 20),     column_names_gp = gpar(fontsize = 20),
                                  heatmap_legend_param = list(labels_gp = gpar(fontsize = 16),title_gp=gpar(fontsize = 20)),
    )
    filename_save <- sub(" -> ","_", unique(df_tmp$target)[1])
    filename_save <- gsub(" ","_", filename_save)
    filename_save <- gsub("/","-", filename_save)
    pdf(file= paste('CCI_analysis_mean/category_filter_combine/',filename_save,'.pdf',sep=''),width = 10 ,height =1.2 + 25*dim(df_weight_cv_matrix)[1]/55)
    draw(h1,column_title=unique(df_tmp$target)[1],column_title_gp = gpar(fontsize = 20))
    dev.off()
    }
  }
}
CCI_heatmap_interaction_pair_comparison_filter_combine(unique(cellpair_seq_whole$cellpair),thresh = 0.1,top_n = 1e5)
}

  for(jj in 1:9){
  cellchat_list_CV <- cellchat_list
  for(j in 1:length(cellchat_list_CV)){
    cellchat_list_CV_j <- cellchat_list_CV[[j]]
    cv_idx <- which(dimnames(cellchat_list_CV_j@netP$prob)[[2]]!=target_names[jj])
    cellchat_list_CV_j@netP$prob[,cv_idx,] <- 0
    cellchat_list_CV_j@net$prob[,cv_idx,] <- 0
    cellchat_list_CV_j@net$count[,cv_idx] <- 0
    cellchat_list_CV_j@net$weight[,cv_idx] <- 0
    cellchat_list_CV[[j]] <- cellchat_list_CV_j
  }
  cellchat_CV <- mergeCellChat(cellchat_list_CV, add.names = names(cellchat_list_CV))

  gg_list <- list()
  for(j in 1:2){
    gg_list[[j]]  <- rankNet(cellchat_CV, mode = "comparison",comparison = c(j,j+1), stacked = F, do.stat = F, font.size = 16,return.data = T)
  }
  df_weight_cv <- rbind(gg_list[[1]]$signaling.contribution,gg_list[[2]]$signaling.contribution[gg_list[[2]]$signaling.contribution$group=='OR',])
  df_weight_cv <- df_weight_cv[,c(1,2,4)]
  df_weight_cv_wide <- reshape2::dcast(df_weight_cv,formula = name ~ group, value.var = 'contribution')
  df_weight_cv_wide[is.na(df_weight_cv_wide)] <- 0
  df_weight_cv_wide <- df_weight_cv_wide[,c(1,3,2,4)]
  tol.rate <- 0.05
  top_n <- 5
  df_weight_cv_wide$sign1 <- ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))> tol.rate,1,
                                    ifelse((df_weight_cv_wide$OS-df_weight_cv_wide$YS)/(1e-16+apply(df_weight_cv_wide[,c('YS','OS')],MARGIN = 1,min))< -tol.rate,-1,0))

  df_weight_cv_wide$sign2 <- ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))> tol.rate,1,
                                    ifelse((df_weight_cv_wide$OR-df_weight_cv_wide$OS)/(1e-16+apply(df_weight_cv_wide[,c('OR','OS')],MARGIN = 1,min))< -tol.rate,-1,0))
  # df_weight_cv_wide <- df_weight_cv_wide[df_weight_cv_wide$sign1*df_weight_cv_wide$sign2 >=0,]
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
  df_weight_cv_wide %>% dplyr::filter(total_weight > 0) %>% arrange(category,-total_weight) -> df_weight_cv_wide
  if(!is.null(top_n) & top_n >0.5){
    df_weight_cv_wide %>% group_by(category)  %>%  slice_max(n = round(top_n), order_by = total_weight) -> df_weight_cv_wide
  }
  df_weight_cv_matrix <- as.matrix(df_weight_cv_wide[,c(2,3,4)])
  rownames(df_weight_cv_matrix) <- df_weight_cv_wide$name
  library(ComplexHeatmap);library(circlize)
  # row_ha = rowAnnotation(foo = anno_block(gp = gpar(fill = colorRampPalette(c("green", 'grey',"red"))(5)[c(1,2,4,5,3)],alpha=0.8), labels = (sort(unique(df_weight_cv_wide$category))),
  #                                         labels_gp = gpar(col = "white", fontsize = 14),labels_rot=0,width=unit(30,'mm')))
  col_tmp <- colorRampPalette(c("green",'grey', "red"))(7)[c(1:3,5:7,4)]
  labels_tmp <-  c('up_both','up_aging','up_IR','down_IR','down_aging','down_both','others')
  labels_idx <- which(labels_tmp %in% df_weight_cv_wide$category )
  # row_ha = rowAnnotation(foo = anno_block(gp = gpar(fill =col_tmp[labels_idx] ,alpha=0.8), labels = labels_tmp[labels_idx],
  #                                         labels_gp = gpar(col = "black", fontsize = 20, fontface=2),labels_rot=0,width=unit(45,'mm')))
  ht_opt$ROW_ANNO_PADDING= unit(1, "mm")
  #breaks = c(0.01, 0.1, 1, 10)
  row_ha = rowAnnotation(trend= anno_block(gp = gpar(fill =col_tmp[labels_idx] ,alpha=0.8), labels = labels_tmp[labels_idx],
                                labels_gp = gpar(col = "black", fontsize = 20, fontface=2),labels_rot=0,width=unit(45,'mm'),show_name=T),
                         total_strength = anno_barplot(log1p(df_weight_cv_wide$total_weight),  axis_param=list(gp=gpar(fontsize = 16)),#at = log10(breaks), labels =breaks,direction='normal'), 
                                                         width=unit(30,'mm')),
                        annotation_name_gp=gpar(fontsize = 16),gap=unit(1,'mm'))
  col_fun = colorRamp2(c(-1.5, 1.5), c("white","blue"))
  h1 <- ComplexHeatmap::Heatmap(as.data.frame(t(scale(t(df_weight_cv_matrix)))),
                                cluster_rows = F,cluster_columns = F,col = col_fun,
                                row_names_side = 'left', right_annotation = row_ha, #right_annotation = row_ha,
                                column_names_rot = 0,name = 'scaled \n strength',
                                row_split=(df_weight_cv_wide$category),
                                row_title=NULL, row_names_gp = gpar(fontsize=20), column_names_gp = gpar(fontsize=20),
                                heatmap_legend_param = list(labels_gp = gpar(fontsize = 16),title_gp=gpar(fontsize = 20)),
  )
  filename_save <- sub(" -> ","_", target_names[jj])
  filename_save <- gsub(" ","_", filename_save)
  filename_save <- sub("/","-", filename_save)
  pdf(file= paste('CCI_analysis_mean/category_multiple/',filename_save,'.pdf',sep=''),width = 10 ,height =1.2 + 25*dim(df_weight_cv_matrix)[1]/55)
  draw(h1,column_title=target_names[jj],column_title_gp = gpar(fontsize = 20))
  dev.off()
}
