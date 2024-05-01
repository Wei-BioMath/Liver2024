library(sp)
library(sf)
library(dplyr)
library(CellChat)
library(Seurat)
setwd('/Users/weizhao/Documents/Hep_Seq/R-spatial/')
sample_name <- c('A1','B2','A2','C1','B1','C2')
exp_condition <- structure(rep(c('Young','Old_sensitive','Old_resistance'),each=2),names=sample_name)
roi_path <- paste('~/Documents/Hep_Seq/32837-slide3_submission/RoiSet_',sample_name,'.zip',sep='')
sp_mtx_path <- paste('~/Documents/Hep_Seq/32837-slide3_submission/32837-slide3_',sample_name,'-1_results.txt',sep='')
rda_path <- paste('/Users/weizhao/Documents/Hep_Seq/R-spatial/',sample_name,'.rda',sep='')
roi_list <- list();sp_mtx_cell_by_gene_list <- list();seurat_list <- list()
library(foreach);library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
for(j in 1:6){
### load ROI
library(RImageJROI);c2_roi <- read.ijzip(roi_path[j])
roi_sp <- lapply(c2_roi, function(x){
  x_sr <- sp::Polygon(x$coords)
  x_sr_list <- sp::Polygons(list(x_sr),x$name)
  return(x_sr_list)
})
roi_SpP <- SpatialPolygons(roi_sp,1:length(roi_sp));
roi_list[[j]] <- roi_SpP

### import the expression data (count by location)
library(data.table);sp_mtx <- fread(sp_mtx_path[j])
### generate gene by cell matrix
sp_mtx_cell_by_gene <- matrix(0,nrow=length(c2_roi),ncol=length(unique(sp_mtx$V4)),dimnames=list(names(c2_roi),unique(sp_mtx$V4)))
count_single <- function(i){
  x <- c2_roi[[i]]
  x_polygon <- x$coords
  xmin <- min(x_polygon[,1])
  xmax <- max(x_polygon[,1])
  ymin <- min(x_polygon[,2])
  ymax <- max(x_polygon[,2])
  x_rectangle <- cbind(c(xmin,xmin,xmax,xmax,xmin),c(ymin,ymax,ymax,ymin,ymin))
  pos_idx <- which(sp_mtx$V1>=xmin & sp_mtx$V1<=xmax & sp_mtx$V2>=ymin & sp_mtx$V2<=ymax)
  pos_coords <- sp_mtx[pos_idx,1:2]
  pip <- sp::point.in.polygon(as.numeric(unlist(pos_coords[,1])),as.numeric(unlist(pos_coords[,2])),x_polygon[,1],x_polygon[,2])
  x_idx <- pos_idx[which(pip>0)]
  if(length(x_idx)>0){
    sp_mtx_x_sum <- aggregate(sp_mtx[x_idx,]$V3,list(sp_mtx[x_idx,]$V4), FUN=length)
    rownames(sp_mtx_x_sum) <- sp_mtx_x_sum[,1]
    sp_mtx_cell_by_gene[x$name,rownames(sp_mtx_x_sum)] <- sp_mtx_x_sum[,2]
  }
  return(sp_mtx_cell_by_gene[x$name,])
}
start.time <- Sys.time()
sp_mtx_cell_by_gene_tmp <- sp_mtx_cell_by_gene
sp_mtx_cell_by_gene_tmp <- foreach(i= 1:length(c2_roi),.combine = 'rbind') %dopar% { #length(c2_roi)
  .GlobalEnv$c2_roi <- c2_roi
  .GlobalEnv$sp_mtx_cell_by_gene <- sp_mtx_cell_by_gene
  count_single(i)
}
end.time <- Sys.time();time.taken <- end.time - start.time;
time.taken
rownames(sp_mtx_cell_by_gene_tmp) <- names(c2_roi)
sp_mtx_cell_by_gene <- sp_mtx_cell_by_gene_tmp
sp_mtx_cell_by_gene_list[[j]] <- sp_mtx_cell_by_gene
# load(file=rda_path[j])

### create seurat object
c2_sp_seurat <- CreateSeuratObject(count=t(sp_mtx_cell_by_gene),project = sample_name[j], min.cells = 20, min.features = 1)
c2_sp_seurat <- subset(c2_sp_seurat, subset = nCount_RNA >=20 & nFeature_RNA >=1)
seurat_list[[j]] <- c2_sp_seurat
# save(roi_SpP,sp_mtx_cell_by_gene,c2_sp_seurat,file=rda_path[j])
}
stopCluster(cl)
# save(roi_list,sp_mtx_cell_by_gene_list, seurat_list, file='sp_six.rda')
## method 1: SCT
for(j in 1:6){
  seurat_list[[j]] <- RenameCells(seurat_list[[j]], new.names=paste(sample_name[j],Cells(seurat_list[[j]]),sep='-'))
}
ifnb <-  merge(seurat_list[[1]], y = c(seurat_list[[2]],seurat_list[[3]],seurat_list[[4]],seurat_list[[5]],seurat_list[[6]] ), add.cell.ids = exp_condition, project = "Hep_all")
ifnb.list <- SplitObject(ifnb, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- SCTransform(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
# saveRDS(immune.combined.sct,file='sp_six_integrationSCT.rds')

load(file='sp_six.rda')
## method 2: Normalization
for(j in 1:6){
  seurat_list[[j]] <- RenameCells(seurat_list[[j]], new.names=paste(sample_name[j],Cells(seurat_list[[j]]),sep='-'))
}
ifnb <-  merge(seurat_list[[1]], y = c(seurat_list[[2]],seurat_list[[3]],seurat_list[[4]],seurat_list[[5]],seurat_list[[6]] ), add.cell.ids = exp_condition, project = "Hep_all")
ifnb.list <- SplitObject(ifnb, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = dim(x)[1])
})
features <- SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = 2000) # nfeatures = 2000 is useless
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
# saveRDS(immune.combined,file='sp_six_integrationNorm.rds')

immune.combined.sct <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationSCT.rds')
DefaultAssay(immune.combined.sct) <- 'integrated'
immune.combined.sct <- RunPCA(immune.combined.sct, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined.sct,ndims = 30)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:15) # 1:30 gives the similar results
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:15) # 1:30 gives the similar results
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.2)
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "orig.ident",shuffle=TRUE)
p2 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE)
p1 + p2
DimPlot(immune.combined.sct, reduction = "umap", split.by = "orig.ident",pt.size = 1)
DimPlot(immune.combined.sct,raster = F)

# immune.combined.sct <- readRDS('/Users/weizhao/Documents/Hep_Seq/R-spatial/sp_six_integrationNorm.rds')
DefaultAssay(immune.combined.sct) <- "RNA"
immune.combined.sct <- NormalizeData(immune.combined.sct)
immune.combined.sct <- ScaleData(immune.combined.sct,  features = rownames(immune.combined.sct))#,vars.to.regress='nCount_RNA') # vars.to.regress = 'percent.mt')
immune.combined.sct <- RunPCA(immune.combined.sct, npcs = 30, verbose = FALSE,features = rownames(immune.combined.sct))
ElbowPlot(immune.combined.sct,ndims = 30)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:15) # 1:30 gives the similar results
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:15) # 1:30 gives the similar results
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.2)
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "orig.ident",shuffle=TRUE)
p2 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = TRUE,shuffle=TRUE)
p1 + p2
DimPlot(immune.combined.sct, reduction = "umap", split.by = "orig.ident",pt.size = 1)
DimPlot(immune.combined.sct)

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
# VlnPlot(immune.combined.sct,features = 'Cyp2e1')
# FeaturePlot(immune.combined.sct,features = 'Cyp2e1')
# FeaturePlot(immune.combined.sct,features = 'Pck1')

#### cluster distribution
ident_cluster <- data.frame(ident = immune.combined.sct@meta.data$orig.ident,cluster=immune.combined.sct@meta.data$seurat_clusters,value=1,row.names = colnames(immune.combined.sct))
ident_cluster$condition <- sub("-.*","",rownames(ident_cluster))
ident_cluster$condition <- sub("_[A,B,C].*","",ident_cluster$condition)
ident_cluster$condition <- factor(ident_cluster$condition, levels = (c('Young','Old_sensitive','Old_resistance')))
ident_cluster_new <- ident_cluster %>%
  group_by(condition, cluster) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))
ggplot(ident_cluster_new, aes(fill=cluster, y=perc*100, x=condition)) + geom_bar(stat="identity")+labs(x = "Condition", y = "Percent")+
  geom_text(aes(label = paste0(round(perc,4)*100,"%")), position = position_stack(vjust = 0.5), size = 4)

#### spatial plot
plot_sp_cluster <- function(sample_single,cluster_idx){
  # sample_single <- 'A1';
  # cluster_idx <- 0:2
sample_idx <- which(sample_name %in% sample_single)
cell_names <- Cells(immune.combined.sct)[grepl(sample_single,Cells(immune.combined.sct)) & (immune.combined.sct@meta.data$seurat_clusters) %in% cluster_idx]
cluster_id <- (immune.combined.sct@meta.data$seurat_clusters)[which( Cells(immune.combined.sct) %in% cell_names)]
cell_names <- sub(".*[A,B,C][1,2]-","",cell_names)
roi_SpP <- roi_list[[sample_idx]]
roi_SpP_subset <- roi_SpP[cell_names]
col_code <- c('red','yellow','blue','green','orange','purple') # 'white',
col_code <- scales::hue_pal()(length(unique(immune.combined.sct@meta.data$seurat_clusters)))
par(mar = c(2, 2, 2, 4))
plot(roi_SpP_subset,col=col_code[cluster_id],border=F,main=sample_single);
legend('right',legend=cluster_idx,fill=col_code[cluster_idx+1], cex=1.5,border = F,box.lwd=0,bty = "n",horiz = F,
       xpd=TRUE, inset=c(-.05,0))
}
plot_sp_cluster('A1',0:7)
plot_sp_cluster('A2',0:7)
plot_sp_cluster('B1',0:7)
plot_sp_cluster('A1',c(0,2,6))
plot_sp_cluster('A1',c(5,6))

sample_name <- c('A1','B2','A2','C1','B1','C2')
for(j in 1:6){
  plot_sp_cluster(sample_name[j],c(2,6))
}

# saveRDS(immune.combined.sct,file='sp_six_integrationSCT.rds')

# plot in spatial -- continuous variable
plot_var <- immune.combined.sct@assays$RNA['Pck1',]
# plot_var <- c2_sp_seurat@assays$SCT@scale.data['Pck1',]
colmap <- circlize::colorRamp2(c(-2,0,2), c("white","slategray1","blue"), transparency = 0, space = "LAB")
roi_SpP_subset <- roi_SpP[which(names(roi_SpP) %in% colnames(c2_sp_seurat))];roi_SpP_subset <- roi_SpP_subset[colnames(c2_sp_seurat)]
# plot(roi_SpP_subset,col=colmap(as.matrix((plot_var-mean(plot_var))/sd(plot_var))),lwd=0.1);
color=(t(as.matrix((plot_var-mean(plot_var))/sd(plot_var))));color[color>2] <- 2;color[color<-2] <- -2
attr = data.frame(ID=names(roi_SpP_subset), color=color,row.names=names(roi_SpP_subset))
# roi_SpP_subset_df = SpatialPolygonsDataFrame(roi_SpP_subset, attr)
# spplot(roi_SpP_subset_df,zcol=2,col.regions = colmap(seq(-2,2,length.out=100)),key.space = "right",lwd=0.1,par.settings = list(axis.line = list(col = 'transparent')), colorkey=list(title='Expression level'))
par(mar = c(2, 2, 2, 4))
plot(roi_SpP_subset,col=colmap(color),lwd=0.2);
bbox <- roi_SpP_subset@bbox
library(plotrix);color.legend(xl=bbox[1,2]*1.02,yb=bbox[2,1],xr=bbox[1,2]*1.05,yt=bbox[2,2],legend=c('Low','','High'),rect.col=colmap(seq(-2,2,length.out=100)),cex=1,gradient="y",align='rb')




