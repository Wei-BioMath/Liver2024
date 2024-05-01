library(dplyr)
library(Seurat)
library(sctransform)
library(future)
library(ggplot2)
library(sp)
# check the current active plan
plan("multicore", workers = 12)
options(future.globals.maxSize= 5000*1024^2)
options(future.seed=TRUE)
set.seed(1234)
setwd('/Users/weizhao/Documents/Hep_Seq/')

sample_name <- c('A1','B2','A2','C1','B1','C2')
exp_condition <- structure(rep(c('Young','Old_sensitive','Old_resistance'),each=2),names=sample_name)

# meta_df <- read.csv(file='/Users/weizhao/Documents/Hep_Seq/metadata_sp.csv')
seurat_sp <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
load(file='/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six.rda')
centroid_list <- lapply(1:6,FUN = function(x){
  roi_j <- roi_list[[x]]@polygons
  centroid_j <- sapply(1:length(roi_j),FUN=function(y){
    roi_j[[names(roi_j)[y]]]@labpt
  })
  dimnames(centroid_j) <- list(c('centroid_x','centroid_y'),names(roi_j))
  centroid_j <- as.data.frame(t(centroid_j))
  rownames(centroid_j) <- paste(exp_condition[x],paste(sample_name[x], rownames(centroid_j),sep='-'), sep='_')
  centroid_j$batch <- x
  return(centroid_j)
})
meta_df <- do.call(rbind, centroid_list)
meta_df <- meta_df[rownames(meta_df) %in% colnames(seurat_sp),]
meta_df$cluster <- seurat_sp$seurat_clusters[match(rownames(meta_df),colnames(seurat_sp))]
meta_df$cell_id <- rownames(meta_df)
meta_df_list <- list()
for(j in 1:6){
  loc_centroid_subset_10K <- meta_df[meta_df$batch==j,]
  rownames(loc_centroid_subset_10K) <- NULL

  celltype <- sort(unique(as.character(loc_centroid_subset_10K[,'cluster'])),method='radix')
  n_celltype <- length(celltype)

  generate_new_cells <- function(idx_CV,pseudo_name){
    loc_centroid_subset_10K_CV <- loc_centroid_subset_10K[idx_CV,]
    loc_centroid_subset_10K_CV$idx_CV <- idx_CV
    coord_mtx <- as.matrix(loc_centroid_subset_10K_CV[,c('centroid_x','centroid_y')])
    Euclidean_dist <- dist(coord_mtx, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
    Euclidean_mtx <- as.matrix(Euclidean_dist)
    Euclidean_mtx[lower.tri(Euclidean_mtx,diag = T)] <- 0
    Euclidean_df <-  reshape2::melt(Euclidean_mtx)
    Euclidean_df <- Euclidean_df[Euclidean_df$value>0,]
    Euclidean_df_filter <- Euclidean_df[Euclidean_df$value <=400,]
    colnames(Euclidean_df_filter) <- c('from','to','dist')
    Euclidean_df_filter$link <- Euclidean_df_filter$dist >0
    
    library(igraph)
    g <- igraph::graph_from_data_frame(Euclidean_df_filter, directed = FALSE, vertices = unique(c(Euclidean_df_filter$from, Euclidean_df_filter$to)))
    coords = layout_with_fr(g)
    c1 = cluster_fast_greedy(g)
    # length(c1)
    # plot(g, layout=coords, vertex.label=NA, vertex.size=10)
    # plot(c1, g, layout=coords)
    
    idx_CV_use <- unique(c(Euclidean_df_filter$from, Euclidean_df_filter$to))
    cell_names <- cell_names <- sub(".*[A,B,C][1,2]-","",loc_centroid_subset_10K$cell_id[idx_CV_use])
    roi_SpP <- roi_list[[j]]
    roi_SpP_subset <- roi_SpP[cell_names]
    color.use <- scales::hue_pal(l=75,c=150,h=c(0,360),h.start=10,direction = 1)(length(c1))
    par(mar = c(2, 2, 2, 12))
    tmp <- structure(as.vector(membership(c1)),names= names(membership(c1)))
    plot(roi_SpP_subset,col=color.use[tmp[as.character(idx_CV_use)]],main='CV community detection',
         bg ='black',border=F)
    
    new_cells <- loc_centroid_subset_10K[idx_CV_use,]
    new_cells$membership <- tmp[as.character(idx_CV_use)]
    new_cells1 <- aggregate(x = list(centroid_x = new_cells$centroid_x,
                                     centroid_y = new_cells$centroid_y),
                            by = list(new_cells$membership),
                            FUN = mean)
    new_cells1 <- new_cells1[,2:3]
    new_cells1$batch <- j
    new_cells1$cluster <- pseudo_name
    new_cells1$cell_id <- paste('pseduo',pseudo_name, 1:dim(new_cells1)[1],sep='_')
    return(new_cells1)
}
  idx_CV <- which(loc_centroid_subset_10K$cluster %in% c(5))
  new_cells1 <- generate_new_cells(idx_CV, 999)
  idx_PV <- which(loc_centroid_subset_10K$cluster %in% c(7))
  new_cells2 <- generate_new_cells(idx_PV, 888)
  
  # library(sf)
  # crop <- c(0.1, 0.1, 0.9, 0.9)
  # crop.margin.xmin <-  st_bbox(roi_SpP)[1] + (st_bbox(roi_SpP)[3]-st_bbox(roi_SpP)[1])*crop[1]
  # crop.margin.xmax <-  st_bbox(roi_SpP)[1] + (st_bbox(roi_SpP)[3]-st_bbox(roi_SpP)[1])*crop[3]
  # crop.margin.ymin <- st_bbox(roi_SpP)[2] + (st_bbox(roi_SpP)[4]-st_bbox(roi_SpP)[2])*crop[2]
  # crop.margin.ymax <- st_bbox(roi_SpP)[2] + (st_bbox(roi_SpP)[4]-st_bbox(roi_SpP)[2])*crop[4]
  
  loc_centroid_subset_10K$cluster <- factor(loc_centroid_subset_10K$cluster, levels=c(0:7,888,999))
  loc_centroid_subset_10K <- rbind(loc_centroid_subset_10K, new_cells1, new_cells2)
  # loc_centroid_subset_10K %>% filter(centroid_x <= crop.margin.xmax & centroid_x >=crop.margin.xmin & 
  #                                      centroid_y <= crop.margin.ymax & centroid_y >=crop.margin.ymin) -> loc_centroid_subset_10K
  # 
  coord_mtx <- as.matrix(loc_centroid_subset_10K[,c('centroid_x','centroid_y')])
  Euclidean_dist <- dist(coord_mtx, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  Euclidean_mtx <- as.matrix(Euclidean_dist)
  idx_pseduo <- which(loc_centroid_subset_10K$cluster==999)
  Euclidean_mtx_CV <- Euclidean_mtx[,idx_pseduo]
  zonation_coord_CV <- apply(Euclidean_mtx_CV, MARGIN = 1, FUN = function(x){
    c(min(x), as.integer(names(which(x==min(x)))[1]))
  })
  zonation_coord_CV <- t(zonation_coord_CV)
  
  idx_pseduo_PV <- which(loc_centroid_subset_10K$cluster==888)
  Euclidean_mtx_PV <- Euclidean_mtx[,idx_pseduo_PV]
  zonation_coord_PV <- apply(Euclidean_mtx_PV, MARGIN = 1, FUN = function(x){
    c(min(x), as.integer(names(which(x==min(x)))[1]))
  })
  zonation_coord_PV <- t(zonation_coord_PV)
  
  loc_centroid_subset_10K <- cbind(loc_centroid_subset_10K, zonation_coord_CV, 
                                   loc_centroid_subset_10K[zonation_coord_CV[,2],c('centroid_x','centroid_y')],
                                   zonation_coord_PV,
                                   loc_centroid_subset_10K[zonation_coord_PV[,2],c('centroid_x','centroid_y')])
  names(loc_centroid_subset_10K)[6:13] <- c('zonation_coord','pesudo_CV_id','pesudo_CV_x','pesudo_CV_y', 
                                            'zonation_coord_PV','pesudo_PV_id','pesudo_PV_x','pesudo_PV_y') 
  loc_centroid_subset_10K$CV_PV_dist <- ((loc_centroid_subset_10K$pesudo_PV_x-loc_centroid_subset_10K$pesudo_CV_x)^2+
    (loc_centroid_subset_10K$pesudo_PV_y-loc_centroid_subset_10K$pesudo_CV_y)^2)^0.5
  
  meta_df_list[[j]] <- loc_centroid_subset_10K[loc_centroid_subset_10K$zonation_coord <=loc_centroid_subset_10K$CV_PV_dist,]
}
meta_df_new <- do.call(rbind, meta_df_list)
meta_df_new %>% filter(!(cluster %in% c(888,999))) -> meta_df_new
meta_df_new$batch <- as.factor(meta_df_new$batch )
meta_df_new$zonation_coord2 <- meta_df_new$zonation_coord/(meta_df_new$CV_PV_dist + 1e-12)
# meta_df_new <- meta_df_new[meta_df_new$zonation_coord <= 3000,]
# meta_df_new$zonation_coord2 <- meta_df_new$zonation_coord*138/1000 # unit micrometer
cluster_celltype_sp <- structure(c('zone 1','zone 2','zone 3','EC','HSC/FB','zone 3 Glul+','KC','EC periportal'),names=0:7)
meta_df_new$annotation <- cluster_celltype_sp[as.character(meta_df_new$cluster)]
cluster_unique <- c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','EC periportal','HSC/FB','KC','immune cells')
# meta_df_new$annotation <- factor(meta_df_new$annotation, levels=rev(cluster_unique))
tmp_sp_abbr <- structure(rep(c('YS','OS','OR'),each=2),names=c(1:6))
meta_df_new$condition <- tmp_sp_abbr[meta_df_new$batch]
meta_df_new$condition  <- factor(meta_df_new$condition, levels = c('YS','OS','OR'))
color.use <- scales::hue_pal(l=75,c=150,h=c(0,360),h.start=10,direction = 1)(length(cluster_unique))
color.use <- structure(color.use, names=cluster_unique)
cluster_unique_sp <- cluster_unique[cluster_unique %in% meta_df_new$annotation]

spatialDistribution_line <- function(meta_df_new, celltype_label='annotation', cluster_unique_sp= cluster_unique_sp, color.use=NULL, font.size=16,legend.symbol.size=5){
  coord_tmp <- seq(min(meta_df_new$zonation_coord2), max(meta_df_new$zonation_coord2),length.out = 21)
  coord_tmp[length(coord_tmp)] <- coord_tmp[length(coord_tmp)]+1e-6
  coord_bins <- data.frame(from=coord_tmp[-length(coord_tmp)],to=coord_tmp[-1])
  bin_prop <- matrix(0, nrow = length(cluster_unique_sp), ncol = dim(coord_bins)[1])
  for(j in 1:dim(bin_prop)[2]){
    for(i in 1:dim(bin_prop)[1]){
      bin_prop[i,j] <- length(which(meta_df_new[,celltype_label] == cluster_unique_sp[i] & 
                                      meta_df_new$zonation_coord2>= coord_bins[j,1] & 
                                      meta_df_new$zonation_coord2 < coord_bins[j,2]
      ))
      bin_prop[i,j]  <- 100*bin_prop[i,j]/ length(which( meta_df_new$zonation_coord2>= coord_bins[j,1] & 
                                                           meta_df_new$zonation_coord2 < coord_bins[j,2]))
    }
  }
  dimnames(bin_prop) <- list(cluster_unique_sp, paste('bin',1:dim(bin_prop)[2],sep=''))
  # bin_prop
  
  bin_prop_df <- reshape2::melt(bin_prop)
  names(bin_prop_df) <- c('annotation','bin','proportion')
  bin_prop_df$bin_coord <- structure((coord_bins[,1]+coord_bins[,2])/2,names=  colnames(bin_prop))[bin_prop_df$bin]
  p3 <- ggplot(data=bin_prop_df,aes(x=bin_coord, y=proportion, color=annotation)) +
    #geom_line(linewidth=1.5)+ #+ geom_point(size=0.1)
    stat_smooth(aes(x = bin_coord, y = proportion), method = "lm",
                formula = y ~ poly(x, 10), se = FALSE, linewidth=1.5)  + 
    labs(y ="Proportion of cells (%)", x = 'Zonation coordinate (along CV-PV axis)') +  guides(color = guide_legend(override.aes = list(size=legend.symbol.size))) + #coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(text=element_text(size=font.size), #change font size of all text
          axis.text=element_text(size=font.size,color = 'black'), #change font size of axis text
          axis.title=element_text(size=font.size), #change font size of axis titles
          legend.text=element_text(size=font.size), #change font size of legend text
          legend.title=element_text(size=font.size)) #change font size of legend title}
  if(!is.null(color.use)){
    p3 <- p3 + scale_color_manual(values = color.use)
  }
  return(p3)
}
# meta_df_new <- meta_df_new[meta_df_new$condition=='YS',]

# together
p_spaDistr_all <- spatialDistribution_line(meta_df_new, celltype_label='annotation', cluster_unique_sp= cluster_unique_sp, color.use = color.use)
ggsave('spaDistr_all.pdf',plot =  p_spaDistr_all, width=15,height=10,units='cm')

# separate by condition
p_spaDistr <- list()
cluster_unique_sp <- cluster_unique[cluster_unique %in% meta_df_new$annotation]
for(j in 1:3){
  idx_j <- which(meta_df_new$condition == c('YS','OS','OR')[j])
  p_spaDistr[[j]] <- spatialDistribution_line(meta_df_new[idx_j,], celltype_label='annotation', cluster_unique_sp= cluster_unique_sp) + xlim(c(0,1)) + coord_cartesian(ylim=c(0, 85))
    ggtitle(c('YS','OS','OR')[j]) 
 if(j %in% 1:2){
   p_spaDistr[[j]] <- p_spaDistr[[j]] + theme(legend.position="none")
 }
}
ggsave('spaDistr_seperate.pdf',plot =  patchwork::wrap_plots(p_spaDistr), width=40,height=10,units='cm')


# separate by annotation
bin_by_condition <- function(condition,celltype_label='annotation', cluster_unique_sp= cluster_unique_sp){
  meta_df_new_j <- meta_df_new[meta_df_new$condition==condition,]
  coord_tmp <- seq(min(meta_df_new$zonation_coord2), max(meta_df_new$zonation_coord2),length.out = 21)
  coord_tmp[length(coord_tmp)] <- coord_tmp[length(coord_tmp)]+1e-6
  coord_bins <- data.frame(from=coord_tmp[-length(coord_tmp)],to=coord_tmp[-1])
  bin_prop <- matrix(0, nrow = length(cluster_unique_sp), ncol = dim(coord_bins)[1])
  for(j in 1:dim(bin_prop)[2]){
    for(i in 1:dim(bin_prop)[1]){
      bin_prop[i,j] <- length(which(meta_df_new_j[,celltype_label] == cluster_unique_sp[i] & 
                                      meta_df_new_j$zonation_coord2>= coord_bins[j,1] & 
                                      meta_df_new_j$zonation_coord2 < coord_bins[j,2]
      ))
      bin_prop[i,j]  <- 100*bin_prop[i,j]/ length(which( meta_df_new_j$zonation_coord2>= coord_bins[j,1] & 
                                                           meta_df_new_j$zonation_coord2 < coord_bins[j,2]))
    }
  }
  dimnames(bin_prop) <- list(cluster_unique_sp, paste('bin',1:dim(bin_prop)[2],sep=''))
  # bin_prop
  
  bin_prop_df <- reshape2::melt(bin_prop)
  names(bin_prop_df) <- c('annotation','bin','proportion')
  bin_prop_df$bin_coord <- structure((coord_bins[,1]+coord_bins[,2])/2,names=  colnames(bin_prop))[bin_prop_df$bin]
  bin_prop_df$condition <- condition 
  return(bin_prop_df)
  
}
meta_df_by_condition <- do.call(rbind, lapply(c('YS','OS','OR'), function(x){
  bin_by_condition(x, celltype_label='annotation', cluster_unique_sp= cluster_unique_sp)
}))
meta_df_by_condition$condition <- factor(meta_df_by_condition$condition, levels=c('YS','OS','OR'))
plot_only <- function(bin_prop_df,celltype_label='annotation', color.use = NULL, font.size=16,legend.symbol.size=5){
p3 <- ggplot(data=bin_prop_df,aes_string(x='bin_coord', y='proportion', color=celltype_label)) +
  #geom_line(linewidth=1.5)+ #+ geom_point(size=0.1)
  stat_smooth(aes_string(x = 'bin_coord', y = 'proportion'), method = "lm",
              formula = y ~ poly(x, 10), se = FALSE, linewidth=1.5)  + 
  labs(y ="Proportion of cells (%)", x = 'Zonation coordinate') +  guides(color = guide_legend(override.aes = list(size=legend.symbol.size))) + #coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text=element_text(size=font.size), #change font size of all text
        axis.text=element_text(size=font.size,color = 'black'), #change font size of axis text
        axis.title=element_text(size=font.size), #change font size of axis titles
        legend.text=element_text(size=font.size), #change font size of legend text
        legend.title=element_text(size=font.size)) #change font size of legend title}
if(!is.null(color.use)){
  p3 <- p3 + scale_color_manual(values = color.use)
}
return(p3)
}
p_spaDistr_by_condition <- list()
cluster_unique_sp <- cluster_unique[cluster_unique %in% meta_df_by_condition$annotation]
for(j in 1:length(cluster_unique_sp)){
  idx_j <- which(meta_df_by_condition$annotation == cluster_unique_sp[j])
  p_spaDistr_by_condition[[j]] <- plot_only(meta_df_by_condition[idx_j,], celltype_label='condition') + xlim(c(0,1)) + 
  ggtitle(cluster_unique_sp[j]) 
  if(j %in% 1:(length(cluster_unique_sp)-1)){
    p_spaDistr_by_condition[[j]] <- p_spaDistr_by_condition[[j]] + theme(legend.position="none")
  }
}
ggsave('spaDistr_seperate_by_condition.pdf',plot =  patchwork::wrap_plots(p_spaDistr_by_condition), width=30,height=25,units='cm')

## distance between cell types
## using formula from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02783-y
t0 <- Sys.time()
celltype_cluster_sp <- structure(c('zone 1','zone 2','zone 3','EC','HSC/FB','zone 3 Glul+','KC','EC periportal'),names=as.character(0:7))
meta_df$annotation <- celltype_cluster_sp[as.character(meta_df$cluster)]
df_idx <- expand.grid(celltype_cluster_sp,celltype_cluster_sp)
df_idx <- cbind(df_idx, expand.grid(1:8,1:8))
colnames(df_idx) <- c('from','to','from_idx','to_idx')
df_idx$id <- 1:dim(df_idx)[1]
dist_list <- rep(list(1),dim(df_idx)[1])
for(j in 1:dim(df_idx)[1]){
  if(length(dist_list[[j]])==1){
  cell_from <- celltype_cluster_sp[df_idx$from_idx[j]]
  cell_to <- celltype_cluster_sp[df_idx$to_idx[j]]
  list_j1 <- list()
  list_j2 <- list()
  for(jj in 1:6){
    meta_df_tmp <- meta_df %>% filter(batch==jj & annotation %in% c(cell_from, cell_to))
    idx_from <- which(meta_df_tmp$annotation %in% cell_from)
    idx_to <- which(meta_df_tmp$annotation %in% cell_to)
    coord_mtx <- as.matrix(meta_df_tmp[,c('centroid_x','centroid_y')])
    Euclidean_dist <- dist(coord_mtx, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
    Euclidean_mtx <- as.matrix(Euclidean_dist)[idx_from, idx_to]
    diag(Euclidean_mtx) <- NA
    list_j1[[jj]] <- apply(Euclidean_mtx, 1,function(x){min(x, na.rm = T)})
    list_j2[[jj]] <- apply(Euclidean_mtx, 2,function(x){min(x, na.rm = T)})
    list_j1[[jj]] <- c(list_j1[[jj]], list_j2[[jj]])
    list_j2[[jj]] <- list_j1[[jj]]
  }
  dist_list[[j]] <- list_j1
  if(cell_from!=cell_to){
   idx_tmp <- which(df_idx$from==cell_to & df_idx$to==cell_from) 
   dist_list[[idx_tmp]] <- list_j2
  }
  }
}
save(dist_list,celltype_cluster_sp, df_idx, file='dist_list.rda')
Sys.time() - t0

# dist_list <- readRDS('/Users/weizhao/Documents/Hep_Seq/dist_list.rds')
load('/Users/weizhao/Documents/Hep_Seq/dist_list.rda')
dist_mtx <- matrix(0,8,8)
for(i in 1:8){
  for(j in 1:8){
    dist_mtx[i,j] <- mean(do.call(c, dist_list[[which(df_idx$from_idx==i & df_idx$to_idx==j) ]]))
  }
}
dimnames(dist_mtx) <- list(celltype_cluster_sp, celltype_cluster_sp)
dist_mtx_tmp0 <- (dist_mtx + t(dist_mtx))/2
dist_mtx_tmp0[lower.tri(dist_mtx_tmp0)] <- -1
dist_mtx_tmp0 <- reshape2::melt(dist_mtx_tmp0)
dist_mtx_tmp0$condition <- 'combined'
dist_mtx_df0 <- dist_mtx_tmp0[dist_mtx_tmp0$value>-1,]

library(ComplexHeatmap)
dist_mtx_list <- list()
dist_mtx_list1 <- list()
dist_mtx_df <- NULL
dist_mtx_df1 <- NULL
for(jj in 1:3){
  dist_mtx_tmp <- matrix(0,8,8)
  for(i in 1:8){
    for(j in 1:8){
      dist_mtx_tmp[i,j] <- mean(do.call(c, dist_list[[which(df_idx$from_idx==i & df_idx$to_idx==j) ]][(2*jj-1):2*jj]))
    }
  }
  dimnames(dist_mtx_tmp) <- list(celltype_cluster_sp, celltype_cluster_sp)
  dist_mtx_list[[jj]] <- dist_mtx_tmp
  dist_mtx_df_tmp <- reshape2::melt(dist_mtx_tmp)
  dist_mtx_df_tmp$condition <- c('YS','OS','OR')[jj]
  dist_mtx_df <- rbind(dist_mtx_df,dist_mtx_df_tmp )
  
  dist_mtx_tmp1 <- (dist_mtx_tmp + t(dist_mtx_tmp))/2
  dist_mtx_tmp1[lower.tri(dist_mtx_tmp1)] <- -1
  dist_mtx_list1[[jj]] <- dist_mtx_tmp1
  dist_mtx_df_tmp1 <- reshape2::melt(dist_mtx_tmp1)
  dist_mtx_df_tmp1$condition <- c('YS','OS','OR')[jj]
  dist_mtx_df1 <- rbind(dist_mtx_df1,dist_mtx_df_tmp1[dist_mtx_df_tmp1$value>-1,])
}

names(dist_mtx_df) <- c('from','to','distance','condition')
dist_mtx_df$cellpair <- paste(dist_mtx_df$from, dist_mtx_df$to, sep='->')
dist_mtx_mtx <- reshape2::dcast(dist_mtx_df, cellpair ~ condition, value.var='distance')
rownames(dist_mtx_mtx)  <-  dist_mtx_mtx$cellpair
dist_mtx_mtx <- dist_mtx_mtx[,-1]; dist_mtx_mtx <- dist_mtx_mtx[,c(3,2,1)]
# Heatmap((dist_mtx_mtx-rowSums(dist_mtx_mtx)/3)/(apply(dist_mtx_mtx, 1,function(x){max(x)-min(x)})),cluster_rows = F, cluster_columns = F)

names(dist_mtx_df1) <- c('from','to','distance','condition')
dist_mtx_df1$cellpair <- paste(dist_mtx_df1$from, dist_mtx_df1$to, sep='->')
dist_mtx_mtx1 <- reshape2::dcast(dist_mtx_df1, cellpair ~ condition, value.var='distance')
rownames(dist_mtx_mtx1)  <-  dist_mtx_mtx1$cellpair
dist_mtx_mtx1 <- dist_mtx_mtx1[,-1]; dist_mtx_mtx1 <- dist_mtx_mtx1[,c(3,2,1)]
# heat_mtx <- as.data.frame(t(scale(t(heat_mtx))))   # (dist_mtx_mtx1-rowSums(dist_mtx_mtx1)/3)/(apply(dist_mtx_mtx1, 1,function(x){max(x)-min(x)}))
# heat_mtx <- heat_mtx[!is.nan(heat_mtx$YS),]
# Heatmap(heat_mtx,name = 'scaled distance', cluster_rows = T, cluster_columns = F)
# Heatmap(dist_mtx_mtx1*138/1000,name = 'raw distance', cluster_rows = T, cluster_columns = F)

names(dist_mtx_df0) <- c('from','to','distance','condition')
rownames(dist_mtx_df0) <- paste(dist_mtx_df0$from, dist_mtx_df0$to, sep='->')
order_idx <- order(dist_mtx_df0$distance,decreasing = F)
rownames(dist_mtx_df0) <- sub('->','--',rownames(dist_mtx_df0))
rownames(dist_mtx_mtx1) <- sub('->','--',rownames(dist_mtx_mtx1))
row_ha = rowAnnotation(average_distance = anno_barplot(dist_mtx_df0[order_idx,'distance']*138/1000,  axis_param=list(gp=gpar(fontsize = 16))), 
                       width=unit(30,'mm'),annotation_name_gp=gpar(fontsize = 16))
ht <- Heatmap(log10(dist_mtx_mtx1[rownames(dist_mtx_df0)[order_idx],]*138/1000),name = 'scaled distance', cluster_rows = F, cluster_columns = F,
              row_names_side = 'left',column_names_rot = 0,right_annotation = row_ha,
              column_names_gp = grid::gpar(fontsize = 16), 
              row_names_gp = grid::gpar(fontsize = 16), 
              heatmap_legend_param = list(labels_gp = gpar(fontsize = 16), at = c(1,log10(50), 2,log10(250)), title = expression( paste('distance (',mu,'m)')),
                                          title_gp = gpar(fontsize = 16),  labels = c(10,50, 100,250),
                                          legend_direction = "horizontal"))
pdf(file="Spatial_distance.pdf",width =9.5,height = 12)
draw(ht, padding = unit(c(2, 20, 2, 20), "mm"),heatmap_legend_side="bottom") ## see right heatmap in following
dev.off()


### statistal test
dist_pairs_list <- list()
for(j in 1:dim(dist_mtx_mtx1)[1]){
  cell_from <- unlist(strsplit(rownames(dist_mtx_mtx1)[j],split = '--'))[1]
  cell_to <- unlist(strsplit(rownames(dist_mtx_mtx1)[j],split = '--'))[2]
  if(cell_from==cell_to){
  idx_tmp <- df_idx[which(df_idx$from==cell_from & df_idx$to==cell_to),'id']
  df_tmp <- NULL
  for(jj in 1:3){
    df_tmp_tmp <- do.call(c, dist_list[[idx_tmp]][(2*jj-1):2*jj])
    df_tmp_tmp <- data.frame(distance=df_tmp_tmp, condition = c('YS','OS','OR')[jj])
    df_tmp <- rbind(df_tmp, df_tmp_tmp)
  }
  dist_pairs_list[[j]] <- df_tmp
  } else {
    idx_tmp <- df_idx[which((df_idx$from==cell_from & df_idx$to==cell_to)|(df_idx$from==cell_to & df_idx$to==cell_from)),'id']
    df_tmp <- NULL
    for(jj in 1:3){
      for(jjj in 1:1){
      df_tmp_tmp <- do.call(c, dist_list[[idx_tmp[jjj]]][(2*jj-1):2*jj])
      df_tmp_tmp <- data.frame(distance=df_tmp_tmp, condition = c('YS','OS','OR')[jj])
      df_tmp <- rbind(df_tmp, df_tmp_tmp)}
    }
    dist_pairs_list[[j]] <- df_tmp
  }
}
names(dist_pairs_list) <- rownames(dist_mtx_mtx1)
dist_pairs_list <- dist_pairs_list[rownames(dist_mtx_df0)[order_idx]]
# plot
library(ggpubr)
dist_pairs_combined <- do.call(rbind, dist_pairs_list )
dist_pairs_combined$cellpair <- rep(names(dist_pairs_list),lapply(dist_pairs_list,function(x){ dim(x)[1]}))
dist_pairs_combined$distance <- dist_pairs_combined$distance*138/1000
dist_pairs_combined$cellpair  <- factor(dist_pairs_combined$cellpair , levels= rownames(dist_mtx_df0)[order_idx])
# p = ggboxplot(dist_pairs_combined, x = "condition", y = "distance",
#               color = "condition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#               facet.by = "cellpair", ncol = 6, ylim=c(10^0.5,10^3.2)
#               )+ yscale("log10", .format = TRUE) + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# my_comparisons <- list(c("YS", "OS"), c("OS", "OR"), c("YS", "OR") )
# p1 <- p + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif",label.y = c(2.4, 2.65, 2.9))+
#   NoLegend()
# ggsave('distance_significance.pdf',plot = p1 , width=40,height=50,units='cm')
# 
# 
# p = ggboxplot(dist_pairs_combined[grep('zone 3--|--zone 3$', dist_pairs_combined$cellpair),], x = "condition", y = "distance",
#               color = "condition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#               facet.by = "cellpair", ncol = 4, ylim=c(10^0.5,10^3.2)) + 
#   yscale("log10", .format = F ) +   
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   theme(text = element_text(size = 16))
# my_comparisons <- list(c("YS", "OS"), c("OS", "OR"), c("YS", "OR") )
# p2 <- p + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif",label.y = c(2.4, 2.65, 2.9))+
#   NoLegend()
# ggsave('distance_significance_zone3.pdf',plot = p2 , width=25,height=15,units='cm')


library(rstatix)
pairwise.test = dist_pairs_combined[grep('zone 3--|--zone 3$', dist_pairs_combined$cellpair),] %>% 
  group_by(cellpair) %>%
  t_test(distance ~ condition) %>% 
  adjust_pvalue(method = 'holm')

p2 <- ggboxplot(dist_pairs_combined[grep('zone 3--|--zone 3$', dist_pairs_combined$cellpair),], x = "condition", y = "distance",
          color = "condition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          facet.by = "cellpair", ncol = 4, ylim=c(10^0.5,10^3.2)) + 
  #yscale("log10", .format = F ) +   
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16)) +
  stat_pvalue_manual(
    pairwise.test, label = "p.adj.signif", 
    y.position =10^c(2.5, 2.75, 3)
  ) +   NoLegend() + ylab(expression( paste('distance (',mu,'m)'))) + coord_trans(y = "log10") + scale_y_continuous(breaks = 10^c(0.5,1,1.5,2,2.5,3),labels=scales::trans_format(log10, scales::math_format(10^.x))) #
ggsave('distance_significance_zone3_tmp.pdf',plot = p2 , width=25,height=15,units='cm')


pairwise.test = dist_pairs_combined[ dist_pairs_combined$cellpair %in% c('zone 3--HSC/FB','zone 3--EC','zone 2--HSC/FB','zone 2--EC'),] %>% 
  group_by(cellpair) %>%
  t_test(distance ~ condition) %>% 
  adjust_pvalue(method = 'holm')
ttt <- dist_pairs_combined %>% dplyr::group_by(condition, cellpair) %>% dplyr::summarize(Mean = mean(distance, na.rm=TRUE))
p2 <- ggboxplot(dist_pairs_combined[ dist_pairs_combined$cellpair %in% c('zone 3--HSC/FB','zone 3--EC','zone 2--HSC/FB','zone 2--EC'),], x = "condition", y = "distance",
                color = "condition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                facet.by = "cellpair", ncol = 1, ylim=c(10^0.5,10^2.85)) +  
  stat_summary(fun=mean, geom="point", shape=20, size=6,aes(color=condition)) +
  # stat_summary(aes(label=round(..y..,2)), fun=mean, geom="text", size=6,
  #              vjust = -0.5) + 
  # yscale("log10", .format = F ) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16)) +
  stat_pvalue_manual(
    pairwise.test, label = "p.adj.signif", 
    y.position =c(10^2.5, 10^2.65, 10^2.8)
  ) +   NoLegend() + ylab(expression( paste('distance (',mu,'m)'))) + coord_trans(y = "log10") + scale_y_continuous(breaks = 10^c(0.5,1,1.5,2,2.5,3),labels=scales::trans_format(log10, scales::math_format(10^.x))) #
# ggsave('distance_significance_examples.pdf',plot = p2 , width=25,height=12,units='cm') # ncol=4
# ggsave('distance_significance_examples.pdf',plot = p2 , width=8,height=30,units='cm') # ncol=1


pairwise.test = dist_pairs_combined %>% 
  group_by(cellpair) %>%
  t_test(distance ~ condition) %>% 
  adjust_pvalue(method = 'holm')

p1 <- ggboxplot(dist_pairs_combined, x = "condition", y = "distance",
          color = "condition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          facet.by = "cellpair", ncol = 6, ylim=c(10^0.5,10^3.2)) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6,aes(color=condition)) +
  # yscale("log10", .format = F ) +   
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16)) +
  stat_pvalue_manual(
    pairwise.test, label = "p.adj.signif", 
    y.position = 10^c(2.5, 2.75, 3)
  ) +   NoLegend() + ylab(expression( paste('distance (',mu,'m)'))) + coord_trans(y = "log10") + scale_y_continuous(breaks = 10^c(0.5,1,1.5,2,2.5,3),labels=scales::trans_format(log10, scales::math_format(10^.x))) #

ggsave('distance_significance.pdf',plot = p1 , width=40,height=50,units='cm')

unique_cellpair <- unique(pairwise.test$cellpair)
cellpair_trend <- data.frame(cellpair=unique_cellpair, trend = 0)
trend_category <- structure(c('down_both','down_aging','down_IR','none','up_IR','up_aging','up_both'),names=c(-3,-2,-1,0,1,2,3))
for(j in 1:length(unique_cellpair)){
  pairwise.test.tmp <- pairwise.test[pairwise.test$cellpair==unique_cellpair[j],]
  pairwise.test.tmp %>% arrange(group1, group2) -> pairwise.test.tmp
  x_seq <- (pairwise.test.tmp$p.adj< 0.05)*pairwise.test.tmp$statistic
  if((x_seq[1]*x_seq[3])>0){
      cellpair_trend[j,2] <-ifelse(x_seq[1]>0, '3','-3')
  } else if((x_seq[1]*x_seq[2])>0){
    cellpair_trend[j,2] <-ifelse(x_seq[1]>0 & x_seq[3]==0, '1','-1')
  } else if((x_seq[3]*x_seq[2])>0){
    cellpair_trend[j,2] <-ifelse(x_seq[2]>0 & x_seq[1]==0, '2','-2')
  } else{
    cellpair_trend[j,2] <- '0'
  }
}
cellpair_trend$category <- trend_category[cellpair_trend$trend]
rownames(cellpair_trend) <- cellpair_trend$cellpair
write.csv(cellpair_trend,file='spatial_trend_category.csv')

library(circlize)
col_fun = colorRamp2(-3:3,colorRampPalette(c("red",'grey', "green"))(7))

ht_opt$ROW_ANNO_PADDING= unit(5, "mm")
row_ha = rowAnnotation(trend=anno_simple(as.numeric(cellpair_trend[rownames(dist_mtx_df0)[order_idx],'trend']),col = col_fun, simple_anno_size=unit(10,'mm') ),
                       average_distance = anno_barplot(dist_mtx_df0[order_idx,'distance']*138/1000,  axis_param=list(gp=gpar(fontsize = 16)), 
                                                       width=unit(30,'mm') ), annotation_name_gp=gpar(fontsize = 16),gap=unit(5,'mm'))
ht <- Heatmap(log10(dist_mtx_mtx1[rownames(dist_mtx_df0)[order_idx],]*138/1000),name = 'scaled distance', cluster_rows = F, cluster_columns = F,
              row_names_side = 'left',column_names_rot = 0,right_annotation = row_ha,
              column_names_gp = grid::gpar(fontsize = 16), 
              row_names_gp = grid::gpar(fontsize = 16), 
              heatmap_legend_param = list(labels_gp = gpar(fontsize = 16), at = c(1,log10(50), 2,log10(250)), title = expression( paste('distance (',mu,'m)')),
                                          title_gp = gpar(fontsize = 16),  labels = c(10,50, 100,250),
                                          legend_direction = "horizontal"))
lgd = Legend(col_fun = col_fun, title = "trend category", at = -3:3, 
             labels = rev(c('up_both','up_aging','up_IR','others','down_IR','down_aging','down_both')))
lgd = Legend(at = -3:3, title = "trend category", legend_gp = gpar(fill = colorRampPalette(c("red",'grey', "green"))(7)),labels = rev(c('up_both','up_aging','up_IR','others','down_IR','down_aging','down_both')))
pdf(file="Spatial_distance_category.pdf",width =10,height = 12)
draw(ht, padding = unit(c(2, 20, 2, 20), "mm"),heatmap_legend_side="bottom", annotation_legend_list=list(lgd)) ## see right heatmap in following
dev.off()

# 
# p_spaDistr <- list()
# cluster_unique_sp <- cluster_unique[cluster_unique %in% meta_df_new$annotation]
# for(j in 1:8){
#   idx_j <- which(!(meta_df_new$cluster %in% c(999)) & meta_df_new$annotation %in% cluster_unique_sp[j])
#   p_spaDistr[[j]] <- spatialDistribution(meta_df_new[idx_j,],
#                                          celltype_label='condition',centroid='zonation_coord2') + xlim(c(0,3000)) +
#     ggtitle(cluster_unique_sp[j]) #+ theme(legend.position="none") + scale_fill_manual(values = color.use)
# }
# patchwork::wrap_plots(p_spaDistr)









# meta_df_new <- do.call(rbind, meta_df_list)
# meta_df_new$batch <- as.factor(meta_df_new$batch )
# meta_df_new <- meta_df_new[meta_df_new$zonation_coord <= 2500,]
# meta_df_new$zonation_coord2 <- meta_df_new$zonation_coord*138/1000 # unit micrometer
# cluster_celltype_sp <- structure(c('zone 1','zone 2','zone 3','EC','HSC/FB','zone 3 Glul+','KC','EC periportal'),names=0:7)
# meta_df_new$annotation <- cluster_celltype_sp[as.character(meta_df_new$cluster)]
# cluster_unique <- c('zone 3','zone 3 Glul+','zone 1','zone 2', 'EC','EC/HSC chimera','EC periportal','HSC/FB','KC','immune cells')
# # meta_df_new$annotation <- factor(meta_df_new$annotation, levels=rev(cluster_unique))
# tmp_sp_abbr <- structure(rep(c('YS','OS','OR'),each=2),names=c(1:6))
# meta_df_new$condition <- tmp_sp_abbr[meta_df_new$batch]
# meta_df_new$condition  <- factor(meta_df_new$condition, levels = c('YS','OS','OR'))
# color.use <- scales::hue_pal(l=75,c=150,h=c(0,360),h.start=10,direction = 1)(length(cluster_unique))
# color.use <- structure(color.use, names=cluster_unique)
# 
# library(scales)
# spatialDistribution <- function(meta,celltype_label='subclass',centroid='centroid_y', color.use= NULL, curve.alpha=0.6,font.size=20,legend.symbol.size=5){
#   attr_x <- celltype_label
#   p3 <- ggplot(data=meta[order(meta[,attr_x],decreasing = T),], aes_string(x=centroid, group=attr_x, fill=attr_x)) +
#     geom_density(adjust=1.5, alpha=curve.alpha) +
#     labs(y ="probability density", x = ~paste("zonation coordinate (",mu,'m)',sep=''),fill=attr_x) +  guides(color = guide_legend(override.aes = list(size=legend.symbol.size))) + #coord_flip() +
#     theme(text=element_text(size=font.size), #change font size of all text
#           axis.text=element_text(size=font.size,color = 'black'), #change font size of axis text
#           axis.title=element_text(size=font.size), #change font size of axis titles
#           legend.text=element_text(size=font.size), #change font size of legend text
#           legend.title=element_text(size=font.size)) #change font size of legend title
#   if(!is.null(color.use)){
#     p3 <- p3+ scale_fill_manual(values=color.use)
#   }
#   return(p3)
# }
# p_spaDistr_all <- spatialDistribution(meta_df_new[!(meta_df_new$cluster %in% c(999)) & meta_df_new$batch %in% c(1:6),],
#                     celltype_label='annotation',centroid='zonation_coord2',
#                     color.use = color.use, font.size = 16)
# ggsave('spaDistr_all.pdf',plot =  p_spaDistr_all, width=15,height=10,units='cm')
# p_spaDistr <- list()
# cluster_unique_sp <- cluster_unique[cluster_unique %in% meta_df_new$annotation]
# for(j in 1:8){
#   idx_j <- which(!(meta_df_new$cluster %in% c(999)) & meta_df_new$annotation %in% cluster_unique_sp[j])
#   p_spaDistr[[j]] <- spatialDistribution(meta_df_new[idx_j,],
#                                          celltype_label='annotation',centroid='zonation_coord2',font.size = 16) + xlim(c(0,3000*138/1000)) +
#     ggtitle(cluster_unique_sp[j]) + theme(legend.position="none") + scale_fill_manual(values = color.use)
# }
# patchwork::wrap_plots(p_spaDistr)
# ggsave('spaDistr_seperate.pdf',plot =  patchwork::wrap_plots(p_spaDistr), width=27,height=25,units='cm')
# 
# 
# p_spaDistr <- list()
# cluster_unique_sp <- cluster_unique[cluster_unique %in% meta_df_new$annotation]
# for(j in 1:8){
#   idx_j <- which(!(meta_df_new$cluster %in% c(999)) & meta_df_new$annotation %in% cluster_unique_sp[j])
#   p_spaDistr[[j]] <- spatialDistribution(meta_df_new[idx_j,],
#                                          celltype_label='condition',centroid='zonation_coord2') + xlim(c(0,3000)) +
#     ggtitle(cluster_unique_sp[j]) #+ theme(legend.position="none") + scale_fill_manual(values = color.use)
# }
# patchwork::wrap_plots(p_spaDistr)
# 
# 
# # spatialDistribution(meta_df_new[!(meta_df_new$cluster %in% c(999)) & meta_df_new$batch %in% c(1:6),],celltype_label='annotation',centroid='zonation_coord')
# # spatialDistribution(meta_df_new[!(meta_df_new$cluster %in%  c(5,999)) & meta_df_new$batch %in% c(1:6),],celltype_label='batch',centroid='zonation_coord')
# 
