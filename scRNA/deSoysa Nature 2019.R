
#############
# Begin
#############
#load Seurat
library(Seurat)
library(patchwork)

#load generic libraries
library(stringi)
library(tidyverse)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(readxl)

#load topGO
#library(BiocManager)
#BiocManager::install("topGO")
library(topGO)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(GOplot)


#############
# Preferences
#############
seurat_umap_alpha_hex_string = "d0"
original_beach_colors <- c("#1a748e","#f0df99","#ecdfcf","#55c4d7","#5f6c24","#992915","#d38e31","#81c7f8","#bbc4af","#393430","#92dccd","#63fba9","#1155d4","#abb2ba","#eccd16","#7b5c52","#063581")
beach_colors <- c("#1a748e","#f0df99","#bbc4af","#ecdfcf","#5f6c24","#992915","#81c7f8","#63fba9","#55c4d7","#d38e31","#abb2ba","#eccd16","#1155d4","#7b5c52")
paired_colors <- brewer.pal(name = "Set1", n=12)
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
paired_colors_alpha <- paste0( paired_colors, seurat_umap_alpha_hex_string )

library(scales)
show_col(original_beach_colors)
show_col(beach_colors)
show_col(paired_colors)

col2gray <- function(cols){
  rgb <- col2rgb(cols)
  gry <- rbind( c(0.3, 0.59, 0.11) ) %*% rgb
  rgb(gry,gry,gry, maxColorValue=255)
}

beach_colors_gray <- paste0( col2gray(beach_colors), "70" )
show_col(beach_colors_gray)

#############
# Data download, read, and create Seurat object
#############
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126128&format=file&file=GSE126128_E775E825_WTKO_10X.csv.gz",method="curl",destfile="GSE126128_E775E825_WTKO_10X.csv.gz")
untar("GSE126128_E775E825_WTKO_10X.csv.gz")
deSoysaMatrix <- read.csv(file = "GSE126128_E775E825_WTKO_10X.csv",header=TRUE,row.name=1,sep=",")
file.remove("GSE126128_E775E825_WTKO_10X.csv")
deSoysaSeurat <- CreateSeuratObject(counts = deSoysaMatrix, project = "deSoysaSeurat")


#############
# Add metadata to Seurat Object
#############
deSoysaSamples <- as.data.frame(matrix(c('1'='E7.75','2'='E7.75','3'='E7.75','4'='E7.75','5'='E7.75','6'='E7.75KO','7'='E7.75KO','8'='E7.75KO','9'='E8.25','10'='E8.25','11'='E8.25KO','12'='E8.25KO')))
colnames(deSoysaSamples)[1] <- "sample"
newtable <- colnames(deSoysaMatrix)
newtable1 <- stri_split_fixed(newtable,".",2)
metadata <- as.data.frame(cbind(newtable,matrix(unlist(newtable1),ncol=2,byrow=TRUE,dimnames = NULL)[,2],deparse.level = 0))
metadata2 <- matrix(metadata[,-1])
metadata3 <- matrix(metadata[,1])
rownames(metadata2) <- metadata3
colnames(metadata2)[1] <- colnames(deSoysaSamples)[1]
metadata2 <- as.data.frame(metadata2)
metadata2$stage <- deSoysaSamples$sample[match(metadata2$sample, rownames(deSoysaSamples))]
deSoysaSeurat <- AddMetaData(object = deSoysaSeurat, metadata = metadata2, col.name=c("sample","stage") )
remove(metadata3)
remove(metadata)
remove(newtable)
remove(newtable1)


#############
# QC and Sample subsetting
#############
deSoysaSeurat[["percent.mt"]] <- PercentageFeatureSet(deSoysaSeurat, pattern = "^mt-")
VlnPlot(deSoysaSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","Epcam","Cdh5","Hbb-b2","Sox2"), ncol = 3)
deSoysaSeurat[["old.ident"]] <- Idents(deSoysaSeurat)
Idents(deSoysaSeurat) <- "stage"
#we remove erythrocytes (Gypa), epiblast/neuroepithelium and derivatives (Sox2/Pax6), endoderm (Epcam), and endothelium (Cdh5)
deSoysaSeurat <- subset(deSoysaSeurat, idents=c("E7.75","E8.25"), subset = nCount_RNA > 2000 & nFeature_RNA < 7000 & percent.mt > 0.1 & Epcam < 0.5 & Pax6 < 1 & Sox2 < 0.5 & Gypa < 1 ) #Cdh5 < 0.25
VlnPlot(deSoysaSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","Hand2"), ncol = 4, group.by="sample")
#save(deSoysaSeurat,file="deSoysaSeurat_QCed_2.Rds")
#load(file="deSoysaSeurat_QCed_2.Rds")


#############
# Split object by sample, and integrate
#############
deSoysa.list <- SplitObject(deSoysaSeurat, split.by = "sample")
deSoysa.list <- lapply(X = deSoysa.list, FUN = function(x) {
  #x <- NormalizeData(x) # the data is provided pre-normalized
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = deSoysa.list)
deSoysa.anchors <- FindIntegrationAnchors(object.list = deSoysa.list, anchor.features = features)
deSoysa.combined <- IntegrateData(anchorset = deSoysa.anchors)
DefaultAssay(deSoysa.combined) <- "integrated"


#############
#Run the standard workflow for visualization and clustering
#############
deSoysa.combined <- ScaleData(deSoysa.combined, verbose = FALSE)
deSoysa.combined <- RunPCA(deSoysa.combined, npcs = 30, verbose = FALSE)
ElbowPlot( deSoysa.combined, ndims=30)
deSoysa.combined <- RunUMAP(deSoysa.combined, reduction = "pca", dims = 1:15)
deSoysa.combined <- FindNeighbors(deSoysa.combined, nn.method="rann", reduction = "pca", dims = 1:15)
deSoysa.combined <- FindClusters(deSoysa.combined, resolution = 1)

# Visualization
p1 <- DimPlot(deSoysa.combined, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.25) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(deSoysa.combined,  reduction = "umap", group.by = "seurat_clusters", pt.size=0.25, label=TRUE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
#ggsave( plot = p1 + p2, height=200, width=400, filename = "deSoysa_original_UMAP.png", units = "mm",dpi = 300 )

#save(deSoysa.combined,file="deSoysaSeurat_integrated.Rds")
#load(file="deSoysaSeurat_integrated.Rds")




#############
#Eliminate non-cardiac clusters
#############
#DefaultAssay(deSoysa.combined) <- "RNA"
#FeaturePlot(deSoysa.combined, features = c("Isl1", "Mef2c", "Nkx2-5","Tnnt2","Tbx5","Hand1","Hand2","Smarcd3","Mab21l2","Msx1","Twist1","Tdgf1","Tek","Foxc2","Cdh5"), cols = c("#FEFEFE", "black"))
#FeaturePlot(deSoysa.combined, features = c("Tcf15","Mab21l2","Msx1","Isl1","Hand2","Tdgf1","Hand1","Six2","Meox1","Lrrn4","Hoxb6","Postn"), cols = c("#FEFEFE", "black"))
#FeaturePlot(deSoysa.combined, features = c("Tcf21","Kdr","Hoxa13","Tdo2","Foxf1","Smarcd3","Mki67","Bmp5","Msx1","Postn","Irx1","Twist1"), cols = c("#FEFEFE", "black"))
#VlnPlot(deSoysa.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","Wt1","Cdh5"), ncol = 3)
Idents(deSoysa.combined) <- "seurat_clusters"
deSoysa.combined <- subset( deSoysa.combined, idents=c(0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,19,20) )



#############
#Re-cluster, now with bad cells gone, and fewer dims needed for biology-recapitulating UMAP
#############
DefaultAssay(deSoysa.combined) <- "integrated"
deSoysa.combined <- ScaleData(deSoysa.combined, verbose = FALSE)
deSoysa.combined <- RunPCA(deSoysa.combined, npcs = 30, verbose = FALSE)
ElbowPlot( deSoysa.combined, ndims=30)
deSoysa.combined <- RunUMAP(deSoysa.combined, reduction = "pca", dims = 1:10)
deSoysa.combined <- FindNeighbors(deSoysa.combined, reduction = "pca", nn.method="rann", dims = 1:10)
deSoysa.combined <- FindClusters(deSoysa.combined, resolution = 0.75)

# Visualization
p1 <- DimPlot(deSoysa.combined, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.25) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(deSoysa.combined,  reduction = "umap", cols=beach_colors_alpha, group.by = "seurat_clusters", pt.size=0.25, label=TRUE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
#ggsave( plot = p1 + p2, height=200, width=400, filename = "deSoysa_complete_UMAP.png", units = "mm",dpi = 300 )

save(deSoysa.combined,file="deSoysaSeurat_integrated_final.Rds")
#load("deSoysaSeurat_integrated_final.Rds")
#table(Idents(deSoysa.combined))

DefaultAssay(deSoysa.combined) <- "RNA"
deSoysa.combined <- ScaleData(deSoysa.combined, verbose = FALSE)

#Cluster identification
#FeaturePlot( deSoysa.combined, features=c("Osr1","Tbx5","Hoxb1","Hoxa1","Tbx1","Cdh5","Foxc2","Pax3","Nr2f2","Etv2","Foxc1","Isl1","Postn","Krt8","Lhx2","Msc"))
#VlnPlot( deSoysa.combined, features=c("Osr1","Hoxb1","Hoxa1","Tbx5"))
#VlnPlot( deSoysa.combined, features=c("Aurkb","Cdc20","Cdkn1b","Mab21l2","Foxc2"))

#Cluster naming
Idents(deSoysa.combined) <- "seurat_clusters"
deSoysa.combined[["old.ident"]] <- Idents(deSoysa.combined)
deSoysa.combined <- RenameIdents(object = deSoysa.combined, 
                                 '0' = "Foxc2+ mesoderm", #this may include cranial mesoderm as well
                                 '1' = "Cardiomyocytes (CMs)",  #these are quiescent non-proliferative cardyomyocytes
                                 '2' = "Somitic mesoderm",
                                 '3' = "Posterior second heart field (pSHF)",
                                 '4' = "Anterior second heart field (aSHF)", #this is the ground state of aSHF
                                 '5' = "Differentiating SHF CMs",  #this becomes the Tdgf1+ cluster
                                 '6' = "Cardiomyocytes (CMs)",  #these are proliferative cardiomyocytes
                                 '7' = "Endothelial cells (ECs)",
                                 '8' = "Cardiopharyngeal progenitors",
                                 '9' = "First heart field (FHF)",
                                 '10' = "Juxtacardiac field (JCF)",   #this is the pro-epicardial portion of the JCF
                                 '11' = "Anterior second heart field (aSHF)",  #these are intermediate progenitors
                                 '12' = "First heart field (FHF)",  #this is the myocardial portion of the JCF
                                 '13' = "Endothelial progenitors",
                                 '14' = "Msx1+ mesoderm",
                                 '15' = "Mesothelium",
                                 '16' = "Proepicardium / S. transversum"
)
deSoysa.combined[["cluster_names"]] <- Idents(deSoysa.combined)


#############
#Visualize final object
#############
# UMAP
p1 <- DimPlot(deSoysa.combined, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.3) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(deSoysa.combined,  reduction = "umap", cols=beach_colors_alpha, group.by = "cluster_names", pt.size=0.3, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
ggsave( plot = p1 + p2, height=200, width=500, filename = "deSoysa_complete_UMAP.png", units = "mm",dpi = 300 )

# CLuster identity
all.markers <- FindAllMarkers(object = deSoysa.combined)
all.markers$ranking = all.markers$avg_log2FC * all.markers$pct.1 * all.markers$pct.1 / all.markers$pct.2
save(all.markers,file="all.markers.Rds")
#load("all.markers.Rds")

top10_uniqueness <- all.markers %>% group_by(cluster) %>% top_n(8, avg_log2FC * pct.1 / pct.2 )
HM <- DoHeatmap(object = deSoysa.combined, features = top10_uniqueness$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15),legend.key.size = unit(1.5,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=550, width=500, filename = "deSoysa_complete_clusterHM1.png", units = "mm",dpi = 300 )

top10_uniqueness <- all.markers %>% group_by(cluster) %>% top_n(8, avg_log2FC )
HM <- DoHeatmap(object = deSoysa.combined, features = top10_uniqueness$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15),legend.key.size = unit(1.5,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=550, width=500, filename = "deSoysa_complete_clusterHM2.png", units = "mm",dpi = 300 )



#############
# Visualize final object feature plots
#############
genes_feature_plot = c("Nkx2-5","Isl1","Tbx5","Mab21l2","Msx1","Foxc2","Mef2c","Foxc1","Twist1","Snai1","Tdgf1","Tnnt2","Tbx1")
#FeaturePlot(new_object, features=genes_feature_plot )

vln <- FeaturePlot(deSoysa.combined, features=genes_feature_plot,cols=c("#000000E0","#FFFFFFE0"), pt.size=0.5,ncol=3,combine=FALSE)
for(i in 1:length(vln)) {
  vln[[i]] <- vln[[i]] + NoLegend() + NoAxes()
  #vln[[i]]$layers[[1]]$aes_params$size = 0
  ggsave( plot = vln[[i]], height=200, width=200, filename = paste("deSoysa_FtPlt_", genes_feature_plot[[i]],".png",sep=""), units = "mm",dpi = 300)
}
#ggsave( plot = wrap_plots(vln,nrow=1), height=100, width=250, filename = "deSoysa_FtPlts.png", units = "mm",dpi = 300)
#wrap_plots(vln,nrow=1)


#############
# Analyze characteristics of JCF cells
#############
# UMAP to highlight JCF
show_colors =  beach_colors_gray
show_colors[10] = "#d38e31"
p1 <- DimPlot(deSoysa.combined, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.3) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(deSoysa.combined,  reduction = "umap", cols=show_colors, group.by = "cluster_names", pt.size=0.3, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
ggsave( plot = p1 + p2, height=200, width=500, filename = "deSoysa_complete_UMAP_JCF.png", units = "mm",dpi = 300 )

Idents(deSoysa.combined) <- "cluster_names"
JCF.markers <- FindMarkers(object = deSoysa.combined, ident.1 = c("Juxtacardiac field (JCF)"), logfc.threshold = 0)
save(JCF.markers,file="JCF.markers.Rds")
#load("JCF.markers.Rds")

top15_log2FC <- JCF.markers %>% top_n(20, avg_log2FC)
top15_neglog2FC <- JCF.markers %>% top_n(20, -avg_log2FC)
top30_total <- rbind( top15_log2FC,top15_neglog2FC )



#############
# Analyze characteristics of Tdgf1 cells
#############
# UMAP to highlight Tdgf1 cluster
show_colors =  beach_colors_gray
show_colors[6] = "#992915"
p1 <- DimPlot(deSoysa.combined, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.3) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(deSoysa.combined,  reduction = "umap", cols=show_colors, group.by = "cluster_names", pt.size=0.3, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
ggsave( plot = p1 + p2, height=200, width=500, filename = "deSoysa_complete_UMAP_Tdgf1cluster.png", units = "mm",dpi = 300 )

Idents(deSoysa.combined) <- "cluster_names"
Tdgf1vsNeighbor.markers <- FindMarkers(object = deSoysa.combined, ident.1 = c("Differentiating SHF CMs"), ident.2 = c("Anterior second heart field (aSHF)","Cardiomyocytes (CMs)"), logfc.threshold = 0)
save(Tdgf1vsNeighbor.markers,file="Tdgf1vsNeighbor.markers.Rds")
#load("Tdgf1vsNeighbor.markers.Rds")

top40_log2FC <- Tdgf1vsNeighbor.markers %>% top_n(40, abs(avg_log2FC))
top40_log2FC <- top40_log2FC[order(-top40_log2FC$avg_log2FC),] 

#subset data so we can properly re-order the clusters for heatmap ease of reading
Tdgf1_analyze_subset <- subset(deSoysa.combined,idents=c("Differentiating SHF CMs","Anterior second heart field (aSHF)","Cardiomyocytes (CMs)"))
Tdgf1_analyze_subset <- RenameIdents(object = Tdgf1_analyze_subset, 
                            "Anterior second heart field (aSHF)" = "Anterior second heart field (aSHF)",
                            "Differentiating SHF CMs" = "Differentiating SHF CMs",
                            "Cardiomyocytes (CMs)" = "Cardiomyocytes (CMs)"
)
Tdgf1_analyze_subset.markers <- FindAllMarkers( object = Tdgf1_analyze_subset )

top10_uniqueness <- Tdgf1_analyze_subset.markers %>% group_by(cluster) %>% top_n(12, avg_log2FC * pct.1 / pct.2 )
HM <- DoHeatmap(object = Tdgf1_analyze_subset, features = top10_uniqueness$gene, label=FALSE, group.colors =c(beach_colors[5],beach_colors[6],beach_colors[2]))  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15),legend.key.size = unit(1.5,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=220, width=400, filename = "deSoysa_Tdgf1vsNeighbor_HM1.png", units = "mm",dpi = 300 )

top10_uniqueness <- Tdgf1_analyze_subset.markers %>% group_by(cluster) %>% top_n(12, avg_log2FC )
HM <- DoHeatmap(object = Tdgf1_analyze_subset, features = top10_uniqueness$gene, label=FALSE, group.colors = c(beach_colors[5],beach_colors[6],beach_colors[2]))  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15),legend.key.size = unit(1.5,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=220, width=400, filename = "deSoysa_Tdgf1vsNeighbor_HM2.png", units = "mm",dpi = 300 )

remove(Tdgf1_analyze_subset)

#take just positives for GO analysis, negatives seem to be mostly related to CM-related contraction
Tdgf1vsNeighbor.markers <- subset( Tdgf1vsNeighbor.markers, subset = avg_log2FC > 0 )


#############
# Analyze characteristics of Tdgf1 cells, E7.75 vs. E8.25
#############
Idents(deSoysa.combined) <- "cluster_names"
Tdgf1_analyze_subset <- subset(deSoysa.combined,idents=c("Differentiating SHF CMs"))
Idents(Tdgf1_analyze_subset) <- "stage"
Tdgf1_analyze_subset.stagemarkers <- FindMarkers(object = Tdgf1_analyze_subset, ident.1 = "E7.75", ident.2 = "E8.25", logfc.threshold = 0)
save(Tdgf1_analyze_subset.stagemarkers,file="Tdgf1clusterStage.markers.Rds")

top15_log2FC <- Tdgf1_analyze_subset.stagemarkers %>% top_n(15, avg_log2FC)
top15_neglog2FC <- Tdgf1_analyze_subset.stagemarkers %>% top_n(15, -avg_log2FC)
top30_total <- rbind( top15_log2FC,top15_neglog2FC )

#top40_log2FC <- Tdgf1_analyze_subset.stagemarkers %>% top_n(-40, p_val_adj)
#top40_log2FC <- top40_log2FC[order(-top40_log2FC$avg_log2FC),] 

HM <- DoHeatmap(object = subset(Tdgf1_analyze_subset,idents=c("E7.75","E8.25")), group.colors = c(paired_colors[1],paired_colors[2]), features = rownames(top30_total), label=FALSE)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 16),legend.key.size = unit(1.55,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=190, width=150, filename = "deSoysa_Tdgf1Stage_HM.png", units = "mm",dpi = 300 )



#############
# Download Gao Cell Res 2019 and Quaranta eLife 2018 data
#############
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126406/suppl/GSE126406%5FE8%5F75%5FIsl1%5FKO%5Fquanty%5Ftranscripts%5Fdiff%5Fexpr%2Etxt%2Egz",method="curl",destfile = "GSE126406_E8_75_Isl1_KO_quanty_transcripts_diff_expr.txt.gz")
Gao.matrix <- fread("GSE126406_E8_75_Isl1_KO_quanty_transcripts_diff_expr.txt.gz",sep="\t",header=TRUE,data.table = FALSE)
#Gao.matrix <- column_to_rownames(Gao.matrix, "SYMBOL")
Gao.matrix.normalize = Gao.matrix %>% group_by(SYMBOL) %>% summarise_all(mean)
Gao.matrix.normalize <- data.frame(Gao.matrix.normalize$CTR_1_raw,Gao.matrix.normalize$CTR_2_raw,Gao.matrix.normalize$CTR_3_raw,Gao.matrix.normalize$CTR_4_raw,Gao.matrix.normalize$Isl1_KO_1_raw,Gao.matrix.normalize$Isl1_KO_2_raw,Gao.matrix.normalize$Isl1_KO_3_raw,Gao.matrix.normalize$Isl1_KO_4_raw, Gao.matrix.normalize$log2FoldChange ,Gao.matrix.normalize$pvalue ,row.names= Gao.matrix.normalize$SYMBOL )
names(Gao.matrix.normalize) <- c("Control 1","Control 2","Control 3","Control 4","Isl1KO 1","Isl1KO 2","Isl1KO 3","Isl1KO 4","avg_log2FC","p_val")
#Gao.matrix.normalize.bak = Gao.matrix.normalize
Gao.matrix.normalize[1:8] <- lapply(X = Gao.matrix.normalize[1:8], FUN = function(x) {
  x <- 1000000 * x/sum(x)
})
#rename "3632451O06Rik" / C14orf37 to name "Armh4"
feature_names <- rownames(Gao.matrix.normalize)
feature_names <- replace(feature_names,feature_names=="3632451O06Rik","Armh4")
rownames(Gao.matrix.normalize) <- feature_names
which(feature_names %in% "Armh4")


download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5770158/bin/elife-31706-supp1.xls",method="curl",destfile = "elife-31706-supp1.xls")
Quaranta.matrix <- as.data.frame(read_xls("elife-31706-supp1.xls",sheet=1,skip=0))
#Quaranta.matrix <- column_to_rownames(Quaranta.matrix, "Symbol")
Quaranta.matrix.normalize = Quaranta.matrix %>% filter( Symbol != "" ) %>% group_by(Symbol) %>% summarise_all(mean)
Quaranta.matrix.normalize <- data.frame(Quaranta.matrix.normalize$`ILSL1 d0 Signal`,Quaranta.matrix.normalize$`ILSL1 d1 Signal`,Quaranta.matrix.normalize$`ILSL1 d2 Signal`,Quaranta.matrix.normalize$`ILSL1 d3 Signal`,Quaranta.matrix.normalize$`ILSL1 d4 Signal`,Quaranta.matrix.normalize$`ILSL1 d5 Signal`,Quaranta.matrix.normalize$`ILSL1 d6 Signal`,Quaranta.matrix.normalize$`ILSL1 d8 Signal`,Quaranta.matrix.normalize$`no ILSL1 d4 Signal`,Quaranta.matrix.normalize$`no ILSL1 d5 Signal`,Quaranta.matrix.normalize$`no ILSL1 d6 Signal`,Quaranta.matrix.normalize$`no ILSL1 d8 Signal`,
                                        Quaranta.matrix.normalize$`RA d4 Signal`,Quaranta.matrix.normalize$`RA d5 Signal`,Quaranta.matrix.normalize$`RA d6 Signal`,Quaranta.matrix.normalize$`RA d8 Signal`,Quaranta.matrix.normalize$`no RA d4 Signal`,Quaranta.matrix.normalize$`no RA d5 Signal`,Quaranta.matrix.normalize$`no RA d6 Signal`,Quaranta.matrix.normalize$`no RA d8 Signal`,
                                        row.names= Quaranta.matrix.normalize$Symbol )
names(Quaranta.matrix.normalize) <- c("Control d0","Control d1","Control d2","Control d3","Control d4","Control d5","Control d6","Control d7","Isl1KO d4","Isl1KO d5","Isl1KO d6","Isl1KO d7",
                                      "RA d4","RA d5","RA d6","RA d7","Control d4","Control d5","Control d6","Control d7"
                                      ) #,"avg_log2FC","p_val_adj")
#rename "3632451O06Rik" / C14orf37 to name "Armh4"
feature_names <- rownames(Quaranta.matrix.normalize)
feature_names <- replace(feature_names,feature_names=="C14orf37","ARMH4")
rownames(Quaranta.matrix.normalize) <- feature_names
which(feature_names %in% "ARMH4")



#############
# Analyze Isl1 KO condition: Compare SHF clusters with Gao Cell Res 2019 and Quaranta eLife 2018
#############
SHF_analyze_subset <- subset(deSoysa.combined,idents=c("Foxc2+ mesoderm","Cardiopharyngeal progenitors","Anterior second heart field (aSHF)","Differentiating SHF CMs","Cardiomyocytes (CMs)"))
SHF_analyze_subset <- RenameIdents(object = SHF_analyze_subset, 
                                   "Foxc2+ mesoderm" = "Foxc2+ mesoderm",
                                   "Cardiopharyngeal progenitors" = "Cardiopharyngeal progenitors",
                                   "Anterior second heart field (aSHF)" = "Anterior second heart field (aSHF)",
                                     "Differentiating SHF CMs" = "Differentiating SHF CMs",
                                     "Cardiomyocytes (CMs)" = "Cardiomyocytes (CMs)"
)
SHF_analyze_subset[["cluster_names2"]] <- Idents(SHF_analyze_subset)
#rename "3632451O06Rik" / C14orf37 to name "Armh4"
feature_names <- rownames(SHF_analyze_subset@assays$RNA@data)
feature_names <- replace(feature_names,feature_names=="3632451O06Rik","Armh4")
rownames(SHF_analyze_subset@assays$RNA@data) <- feature_names
which(feature_names %in% "Armh4")

#p1 <- DimPlot(SHF_analyze_subset, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.25) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
#p2 <- DimPlot(SHF_analyze_subset,  reduction = "umap", cols=beach_colors_alpha, group.by = "cluster_names", pt.size=0.25, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
#p1 + p2


BG <- ggplot(SHF_analyze_subset@meta.data, aes(cluster_names2, fill=stage)) + geom_bar(stat="count",position=position_dodge2(reverse=TRUE),width=0.7) + scale_fill_manual(values=paired_colors) + coord_flip() + scale_x_discrete(limits=rev) + xlab("") + ylab("Cells") + theme(panel.grid.major.x = element_line(colour="#C0C0C0",size=0.5, linetype=1),legend.title = element_text(size = 14),rect = element_blank(),axis.ticks = element_blank(),legend.position = c(0.65, 0.65),legend.text = element_text(size = 12),axis.text.y = element_text(size = 12),axis.text = element_text(size = 12),axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))
#p5a <- ggplot(FB.Seurat@meta.data, aes(cluster_names, fill=experiment_names))+geom_bar(stat="count",position="fill") + scale_fill_manual(values=c(paired_colors[1],paired_colors[3],paired_colors[4])) + RotatedAxis() + xlab("Cluster") + ylab("Cell %") + NoLegend() + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15)) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 14),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))
ggsave( plot = BG, height=50, width=220, filename = "deSoysa_SHF_cell_numbers.png", units = "mm",dpi = 300 )


SHF.cluster.markers <- FindAllMarkers(object = SHF_analyze_subset)
SHF.cluster.markers$ranking = SHF.cluster.markers$avg_log2FC * SHF.cluster.markers$pct.1 * SHF.cluster.markers$pct.1 / SHF.cluster.markers$pct.2
save(SHF.cluster.markers,file="SHF.cluster.markers.Rds")
#load("SHF.cluster.markers.Rds")
#top10_logFC <- all.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10_uniqueness <- SHF.cluster.markers %>% group_by(cluster) %>% top_n(15, ranking)
#DoHeatmap(object = SHF_analyze_subset, features = top10_logFC$gene, label=FALSE, group.colors = beach_colors_alpha)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) #slim.col.label = TRUE, remove.key = TRUE)

#DotPlot(object = SHF_analyze_subset, features =c("Sox4","Gata4","Dlk1")) +scale_y_discrete(limits=rev) + scale_color_viridis() +RotatedAxis() + xlab("") + ylab("")

HM <- DoHeatmap(object = SHF_analyze_subset, features = top10_uniqueness$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15),legend.key.size = unit(1.5,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=550, width=500, filename = "deSoysa_SHF_clusterHM.png", units = "mm",dpi = 300 )

DefaultAssay(SHF_analyze_subset) <- "RNA"
list_SHF_timeline <- c("Lum","Snai1","Foxc2","Twist1","Cdh11","Tbx1","Id3","Sox4","Kdr","Dkk1","Isl1","Bmp4","Hand2","Rgs5","Wnt5a","Abra","Tdgf1","Fhl1","Armh4","Mef2c","Smyd1","Myh7","Fermt2","Gata6","Nkx2-5","Nexn","Ankrd1","Tnni3")
#also include "3632451O06Rik" / Armh4 / C14orf37
#also include Gja1 ?
DP <- DotPlot(object = SHF_analyze_subset, features =list_SHF_timeline) +scale_y_discrete(limits=rev) + scale_color_viridis() + theme(axis.text.x = element_text(angle=90, vjust=0.4, hjust=1)) + xlab("") + ylab("")
ggsave( plot = DP, height=55, width=280, filename = "deSoysa_SHF_clusterDP.png", units = "mm",dpi = 300 )

Quaranta.matrix.subset <-Quaranta.matrix.normalize[  rownames(Quaranta.matrix.normalize) %in% toupper(list_SHF_timeline), ]
Quaranta.matrix.subset <-Quaranta.matrix.subset[match(toupper(list_SHF_timeline), rownames(Quaranta.matrix.subset)),seq(1:12)]
Deltad4 <- Quaranta.matrix.subset[9]/Quaranta.matrix.subset[5]
Deltad5 <- Quaranta.matrix.subset[10]/Quaranta.matrix.subset[6]
Deltad6 <- Quaranta.matrix.subset[11]/Quaranta.matrix.subset[7]
Deltad7 <- Quaranta.matrix.subset[12]/Quaranta.matrix.subset[8]
Quaranta.matrix.subset$expr_change <- log2(rowMeans(data.frame(Deltad4,Deltad5,Deltad6,Deltad7)))

Gao.matrix.subset <- Gao.matrix.normalize[ rownames(Gao.matrix.normalize) %in% list_SHF_timeline, ]    
Gao.matrix.subset <- Gao.matrix.subset[match(list_SHF_timeline, rownames(Gao.matrix.subset)),]
Gao.matrix.subset$log_p <- -log(Gao.matrix.subset$p_val)

#############
# Gao Isl1KO heatmap
#############
plot.new()
png(filename="Gao_SHF_factors_HM.png", width=1600, height=1280, bg="white")
heatmap(as.matrix(t(apply(log2(Gao.matrix.subset[1:8]), 1, FUN= rescale))), margins = c(20,60), col = plasma(256),cexRow=3,cexCol=3)
dev.off()



#############
# Gao Isl1KO SHF features bar graph
#############
bilat_sqrt_trans = function() {
  bilat_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  inv_bilat_sqrt <- function(x){x^2*sign(x)}
  trans_new("bilat_sqrt",bilat_sqrt,inv_bilat_sqrt)
}

BG <- Gao.matrix.subset %>%
ggplot(aes(fill=log_p, y=avg_log2FC, x=rownames(Gao.matrix.subset))) + 
  geom_bar(position="dodge", stat="identity") +
  scale_x_discrete(limits=rownames(Gao.matrix.subset)) +
  scale_y_continuous(trans="bilat_sqrt") +
  scale_fill_gradient2( low="Black",mid="Red",high="#FFFF20", midpoint=5, limits=c(0,77), trans='sqrt', breaks=c(0,2,25,50,75) ) + #1.3)
  #scale_fill_gradient( palette=scales::area_pal(), breaks=c(0,1.3,10), limits=c(0,67) )
  ylab("Average Log2FC") +
  xlab("") +
  #NoLegend() +
  guides(fill=guide_legend(title="-log(p):")) +
  theme(
    #legend.position="none",
    legend.position = c(0.1, 0.3),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_text(size=16,color="black",family="Liberation Sans", angle=90, vjust=0.55),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    #axis.ticks = element_line(color="#444444",size=3),
    axis.ticks = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    legend.title = element_text(size=16,color="black",family="Liberation Sans", angle=0, vjust=0.4),
    legend.text = element_text(size=13,color="black",family="Liberation Sans", angle=0, vjust=0.4),
    axis.text.x = element_text(size=16,color="black",family="Liberation Sans", angle=90, vjust=0.4, hjust=1),
    axis.text.y = element_text(size=16,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#C0C0C0",size=0.2, linetype=1),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=0.5, linetype=2),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(0, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill=c("Light Blue","White"))
  ) 
ggsave( plot = BG, height=100, width=200, filename = "Gao_SHF_cluster_Isl1KOlogFC_BG.png", units = "mm",dpi = 300 )



#############
# Quaranta Isl1KO SHF features bar graph
#############
#quickly calculate 4-observation paired t-tests for dataset
Quaranta.matrix.subset$p_val <- (apply(X = Quaranta.matrix.subset[5:12], MARGIN=1, FUN = function(x) {
  #print( x[1:4] )
  #print( '\n')
  t.test(x[1:4],x[5:8],paired=TRUE)$p.value
}))
Quaranta.matrix.subset$p_val[is.na(Quaranta.matrix.subset$p_val)] <- 1
Quaranta.matrix.subset$log_p <- -log(Quaranta.matrix.subset$p_val)

BG <- Quaranta.matrix.subset %>%
  ggplot(aes(fill=log_p, y=expr_change, x=rownames(Quaranta.matrix.subset))) + 
  geom_bar(position="dodge", stat="identity") +
  scale_x_discrete(limits=rownames(Quaranta.matrix.subset)) +
  scale_y_continuous(trans="bilat_sqrt") +
  scale_fill_gradient2( low="Black",mid="Red",high="#FFFF20", midpoint=5, limits=c(0,8), trans='sqrt', breaks=c(0,2,4,8) ) + #1.3)
  #scale_fill_gradient( palette=scales::area_pal(), breaks=c(0,1.3,10), limits=c(0,67) )
  ylab("Average Log2FC") +
  xlab("") +
  #NoLegend() +
  guides(fill=guide_legend(title="-log(p):")) +
  theme(
    #legend.position="none",
    legend.position = c(0.05, 0.27),
    legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
    plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
    #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
    axis.title.y = element_text(size=16,color="black",family="Liberation Sans", angle=90, vjust=0.55),
    #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
    #axis.ticks = element_line(color="#444444",size=3),
    axis.ticks = element_blank(),
    #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
    legend.title = element_text(size=16,color="black",family="Liberation Sans", angle=0, vjust=0.4),
    legend.text = element_text(size=13,color="black",family="Liberation Sans", angle=0, vjust=0.4),
    axis.text.x = element_text(size=16,color="black",family="Liberation Sans", angle=90, vjust=0.4, hjust=1),
    axis.text.y = element_text(size=16,color="black",family="Liberation Sans", angle=0, vjust=0.65),
    #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
    strip.text = element_text(size=72,color="black",family="Liberation Sans"),
    #panel.grid = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour="#C0C0C0",size=0.2, linetype=1),
    panel.grid.major.x = element_line(colour="#A0A0A0",size=0.5, linetype=2),
    #strip.text = element_blank(),
    #strip.background = element_blank(),
    #rect = element_rect(fill = "#FF0000"),
    #plot.margin = margin(0, 30, 0, 0),
    rect = element_blank(),
    panel.spacing = unit(4,"lines"),
    #rect = element_rect( color="white", size=1)   
    strip.background=element_rect(fill=c("Light Blue","White"))
  ) 
ggsave( plot = BG, height=110, width=200, filename = "Quaranta_SHF_cluster_Isl1KOexprchange_BG.png", units = "mm",dpi = 300 )


#############
# Quaranta differentation SHF features
#############
Quaranta.matrix.subset2 <- Quaranta.matrix.subset[!(row.names(Quaranta.matrix.subset) %in% c("TDGF1","ABRA")),8:1]
names(Quaranta.matrix.subset2) <- paste0("d",as.character(seq(7,0)))
# BG <- melt(as.matrix(t(apply(Quaranta.matrix.subset2, 1, FUN= rescale)))) %>%
#   #ggplot(aes(fill=rep(c("d0","d1","d2","d3","d4","d5","d6","d7"),2),y=names(Quaranta.matrix.subset[1:8]),x=rownames(Quaranta.matrix.subset))) +
#   ggplot(aes(y=value,x=Var1,fill=Var2)) + 
#   geom_bar(position="dodge", stat="identity") +
#   scale_x_discrete(limits=rownames(Quaranta.matrix.subset)) +
#   scale_y_continuous(trans="bilat_sqrt") +
#   #scale_fill_gradient2( low="Black",mid="Red",high="#FFFF20", midpoint=5, limits=c(0,8), trans='sqrt', breaks=c(0,2,4,8) ) + #1.3)
#   scale_fill_brewer(palette = "Spectral", guide="Colourbar") +
#   #scale_fill_gradient( palette=scales::area_pal(), breaks=c(0,1.3,10), limits=c(0,67) )
#   ylab("Average Log2FC") +
#   xlab("") +
#   #NoLegend() +
#   guides(fill=guide_legend(title="Day")) +
#   theme(
#     #legend.position="none",
#     legend.position = c(0.05, 0.27),
#     legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
#     plot.title = element_text(size=32,vjust=1,hjust=0.5,color="black",family="Liberation Sans", face="bold"),
#     #axis.text.y = element_text(size=32,color="black",family="Liberation Sans"),
#     axis.title.y = element_text(size=16,color="black",family="Liberation Sans", angle=90, vjust=0.55),
#     #axis.ticks.y = element_line(colour="#DDDDDD",size=0.5),
#     #axis.ticks = element_line(color="#444444",size=3),
#     axis.ticks = element_blank(),
#     #axis.title.y = element_text(size=32,color="black",family="Liberation Sans", face="bold"),
#     legend.title = element_text(size=16,color="black",family="Liberation Sans", angle=0, vjust=0.4),
#     legend.text = element_text(size=13,color="black",family="Liberation Sans", angle=0, vjust=0.4),
#     axis.text.x = element_text(size=16,color="black",family="Liberation Sans", angle=90, vjust=0.4, hjust=1),
#     axis.text.y = element_text(size=16,color="black",family="Liberation Sans", angle=0, vjust=0.65),
#     #axis.title.x = element_text(size=32,color="black",family="Liberation Sans", face="bold" ),
#     strip.text = element_text(size=72,color="black",family="Liberation Sans"),
#     #panel.grid = element_blank(),
#     #panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_line(colour="#C0C0C0",size=0.2, linetype=1),
#     panel.grid.major.x = element_line(colour="#A0A0A0",size=0.5, linetype=2),
#     #strip.text = element_blank(),
#     #strip.background = element_blank(),
#     #rect = element_rect(fill = "#FF0000"),
#     #plot.margin = margin(0, 30, 0, 0),
#     rect = element_blank(),
#     panel.spacing = unit(4,"lines"),
#     #rect = element_rect( color="white", size=1)   
#     strip.background=element_rect(fill=c("Light Blue","White"))
#   ) 
# ggsave( plot = BG, height=110, width=200, filename = "Quaranta_SHF_cluster_differentation_BG.png", units = "mm",dpi = 300 )

plot.new()
png(filename="Quaranta_SHF_factors_diff_HM.png", width=1600, height=960, bg="white")
heatmap(as.matrix(apply((Quaranta.matrix.subset2), 1, FUN= rescale)), Colv=NA, Rowv=NA, margins = c(60,00), col = plasma(256),cexRow=3,cexCol=3, scale="column")
dev.off()




#############
# Perform GO analysis for CC on differentially expressed genes in JCF
#############
GOdata <- new("topGOdata",
              ontology = "CC", # use biological process ontology
              allGenes = setNames(JCF.markers$p_val_adj,rownames(JCF.markers)),
              geneSel = function(p) p < 0.01,
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol", nodeSize = 10)

resultWeight01 <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allGO <- usedGO(GOdata)

allRes_orig <- GenTable(GOdata, classicFisher = resultFisher,weightFisher=resultWeight01,
                        classicKS = resultKS, elimKS = resultKS.elim,
                        orderBy = "elimKS", ranksOf = "weightFisher", topNodes=length(allGO), numChar=1000)
write.csv(allRes_orig,file="deSoysa_JCF_GOres_CC.csv")


#now, make GO p-value heatmap for top XXX terms
allRes <- head(allRes_orig,  10 )
allRes_ <- data.frame(-log10(as.numeric(allRes$classicFisher)),-log10(as.numeric(allRes$weightFisher)),-log10(as.numeric(allRes$classicKS)),-log10(as.numeric(allRes$elimKS)))
rownames(allRes_) <- allRes$Term
allRes_heatmap <- as.matrix(allRes_)
allRes_heatmap <- allRes_heatmap[nrow(allRes_heatmap):1, ] #avoid using reorderfun, since that is more complicated


#############
# GO CC Heatmap on differentially expressed genes in JCF
#############
plot.new()
png(filename="deSoysa_JCF_GOres_topCC.png", width=1600, height=600, bg="white")
heatmap(allRes_heatmap, Colv = NA, Rowv = NA, margins = c(20,40), col=heat.colors(108)[16:96],cexRow=3.5,cexCol=3,labCol=c("Fisher","Weight","KS","elimKS"), scale="none")

lgd_ = rep(NA, 16)
lgd_[c(1,16)] = c(as.character(as.integer(max(allRes_$X.log10.as.numeric.allRes.classicKS..)+0.5)),as.character(as.integer(min(allRes_$X.log10.as.numeric.allRes.weightFisher..)+0.5)))

legend(x = 0.02, y = 0.95,
       legend = lgd_,
       fill = rev(heat.colors(25)[1:17]), #, direction = -1),
       border = NA,
       bty = "n",
       y.intersp = 0.5,
       cex = 3, text.font = 1,
       title = "-log(p)    \n")
dev.off()



#############
# Perform GO analysis for BP on differentially expressed genes in JCF
#############
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = setNames(JCF.markers$p_val_adj,rownames(JCF.markers)),
              geneSel = function(p) p < 0.01,
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol", nodeSize = 10)

resultWeight01 <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allGO <- usedGO(GOdata)

allRes_orig <- GenTable(GOdata, classicFisher = resultFisher,weightFisher=resultWeight01,
                        classicKS = resultKS, elimKS = resultKS.elim,
                        orderBy = "elimKS", ranksOf = "weightFisher", topNodes=length(allGO), numChar=1000)
write.csv(allRes_orig,file="deSoysa_JCF_GOres_BP.csv")


#now, make GO p-value heatmap for top XXX terms
allRes <- head(allRes_orig,  12 )
allRes_ <- data.frame(-log10(as.numeric(allRes$classicFisher)),-log10(as.numeric(allRes$weightFisher)),-log10(as.numeric(allRes$classicKS)),-log10(as.numeric(allRes$elimKS)))
rownames(allRes_) <- allRes$Term
allRes_heatmap <- as.matrix(allRes_)
allRes_heatmap <- allRes_heatmap[nrow(allRes_heatmap):1, ] #avoid using reorderfun, since that is more complicated


#############
# GO BP Heatmap on differentially expressed genes in JCF
#############
plot.new()
png(filename="deSoysa_JCF_GOres_topBP.png", width=1680, height=680, bg="white")
heatmap(allRes_heatmap, Colv = NA, Rowv = NA, margins = c(20,46), col=heat.colors(108)[16:96],cexRow=3.5,cexCol=3,labCol=c("Fisher","Weight","KS","elimKS"), scale="none")

lgd_ = rep(NA, 16)
lgd_[c(1,16)] = c(as.character(as.integer(max(allRes_$X.log10.as.numeric.allRes.classicKS..)+0.5)),as.character(as.integer(min(allRes_$X.log10.as.numeric.allRes.weightFisher..)+0.5)))

legend(x = 0.02, y = 0.95,
       legend = lgd_,
       fill = rev(heat.colors(25)[1:17]), #, direction = -1),
       border = NA,
       bty = "n",
       y.intersp = 0.5,
       cex = 3, text.font = 1,
       title = "-log(p)    \n")
dev.off()



#############
# GO Chord plot on differentially expressed genes in JCF
#############
#courtesy of code at https://github.com/W-L/ADA_coursework
JCF.markers <- arrange(JCF.markers, desc((avg_log2FC)/p_val_adj))
#JCF.markers <- arrange(JCF.markers, desc(abs(avg_log2FC)))
gt_BP <- genesInTerm(GOdata)
candidates <- data.frame( rownames(JCF.markers), (JCF.markers$pct.1 + JCF.markers$pct.2) / 2, JCF.markers$avg_log2FC, JCF.markers$p_val_adj )
names(candidates) <- c("ID", "baseMean", "logFC", "adj_pval")

add_genes_to_terms <- function(genes2GO, results){
  for (i in seq(1, nrow(results))){
    GOterm <- results$GO.ID[i]
    results[i, "genes"] <- paste(genes2GO[GOterm][[1]], collapse=', ')
  }
  return(results)
}
results_BP <- add_genes_to_terms(genes2GO = gt_BP, results = allRes_orig)
results_BP$Category <- rep("BP", nrow(results_BP))
#full_res <- rbind(results_BP, results_MF, results_CC)
full_res <- subset( results_BP, select = -c(`classicKS`, `weightFisher`,  `classicFisher`, `Rank in weightFisher`) )
full_res$genes <- as.factor(full_res$genes)
names(full_res) <- c("ID", "term", "Annotated", "Significant", "Expected", "adj_pval", "genes", "category")
full_res$adj_pval <- as.numeric(full_res$adj_pval)

#this function is modified from GOplot so gene names don't get capitalized
circle_dat_mod <- function(terms, genes){
  colnames(terms) <- tolower(colnames(terms))
  #terms$genes <- toupper(terms$genes)
  #genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    value <- value[!is.na(value)]
    zsc <- c(zsc, sum(value) / sqrt(count[c]))
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

GO_plot <- circle_dat_mod(full_res, candidates)
#GOCircle(GO_plot)




sub_IDs <- c('GO:0007411','GO:0030198','GO:0007156','GO:0030335')
#sub_IDs <- c('axon guidance', 'neural crest migration', 'positive regulation of epithelial to mesenchymal transition', 'cell migration involved in heart development', 'heterotypic cell-cell adhesion' )
#GOCircle(GO_plot, nsub = sub_IDs)
chord_genes = data.frame( rownames(JCF.markers), JCF.markers$avg_log2FC )
names(chord_genes) <- c("ID", "logFC" )
chord <- chord_dat( data = GO_plot, genes = chord_genes[1:450,], process = sub_IDs ) #
#head(chord)
chord_plot <- GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5, border.size = 0, ribbon.col = (brewer.pal (n=12,name="Set2"))[3:(2+length(sub_IDs))] )
ggsave( plot = chord_plot, height=242, width=200, filename = "deSoysa_JCF_GOBPchord.png", units = "mm",dpi = 300 )
#GOBubble(data=GO_plot, labels = 3, display="single", ID = T, table.legend = T, table.col = T)






#############
# Perform GO analysis for CC on differentially expressed genes in Tdgf1 domain
#############
GOdata <- new("topGOdata",
              ontology = "CC", # use biological process ontology
              allGenes = setNames(Tdgf1vsNeighbor.markers$p_val_adj,rownames(Tdgf1vsNeighbor.markers)),
              geneSel = function(p) p <  0.01,
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol", nodeSize = 10)

resultWeight01 <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allGO <- usedGO(GOdata)

allRes_orig <- GenTable(GOdata, classicFisher = resultFisher,weightFisher=resultWeight01,
                        classicKS = resultKS, elimKS = resultKS.elim,
                        orderBy = "elimKS", ranksOf = "weightFisher", topNodes=length(allGO), numChar=1000)
write.csv(allRes_orig,file="deSoysa_Tdgf1vsNeighbor_GOres_CC.csv")


#now, make GO p-value heatmap for top XXX terms
allRes <- head(allRes_orig,  10 )
allRes_ <- data.frame(-log10(as.numeric(allRes$classicFisher)),-log10(as.numeric(allRes$weightFisher)),-log10(as.numeric(allRes$classicKS)),-log10(as.numeric(allRes$elimKS)))
rownames(allRes_) <- allRes$Term
allRes_heatmap <- as.matrix(allRes_)
allRes_heatmap <- allRes_heatmap[nrow(allRes_heatmap):1, ] #avoid using reorderfun, since that is more complicated


#############
# GO CC Heatmap on differentially expressed genes in Tdgf1 domain
#############
plot.new()
png(filename="deSoysa_Tdgf1vsNeighbor_GOres_topCC.png", width=1600, height=600, bg="white")
heatmap(allRes_heatmap, Colv = NA, Rowv = NA, margins = c(20,40), col=heat.colors(108)[16:96],cexRow=3.5,cexCol=3,labCol=c("Fisher","Weight","KS","elimKS"), scale="none")

lgd_ = rep(NA, 16)
lgd_[c(1,16)] = c(as.character(as.integer(max(allRes_$X.log10.as.numeric.allRes.classicKS..)+0.5)),as.character(as.integer(min(allRes_$X.log10.as.numeric.allRes.weightFisher..)+0.5)))

legend(x = 0.02, y = 0.95,
       legend = lgd_,
       fill = rev(heat.colors(25)[1:17]), #, direction = -1),
       border = NA,
       bty = "n",
       y.intersp = 0.5,
       cex = 3, text.font = 1,
       title = "-log(p)    \n")
dev.off()



#############
# Perform GO analysis for BP on differentially expressed genes in Tdgf1 domain
#############
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = setNames(Tdgf1vsNeighbor.markers$p_val_adj,rownames(Tdgf1vsNeighbor.markers)),
              geneSel = function(p) p < 0.01,
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol", nodeSize = 10)

resultWeight01 <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allGO <- usedGO(GOdata)

allRes_orig <- GenTable(GOdata, classicFisher = resultFisher,weightFisher=resultWeight01,
                        classicKS = resultKS, elimKS = resultKS.elim,
                        orderBy = "elimKS", ranksOf = "weightFisher", topNodes=length(allGO), numChar=1000)
write.csv(allRes_orig,file="deSoysa_Tdgf1vsNeighbor_GOres_BP.csv")


#now, make GO p-value heatmap for top XXX terms
allRes <- head(allRes_orig,  12 )
allRes_ <- data.frame(-log10(as.numeric(allRes$classicFisher)),-log10(as.numeric(allRes$weightFisher)),-log10(as.numeric(allRes$classicKS)),-log10(as.numeric(allRes$elimKS)))
rownames(allRes_) <- allRes$Term
allRes_heatmap <- as.matrix(allRes_)
allRes_heatmap <- allRes_heatmap[nrow(allRes_heatmap):1, ] #avoid using reorderfun, since that is more complicated


#############
# GO BP Heatmap on differentially expressed genes in Tdgf1 domain
#############
plot.new()
png(filename="deSoysa_Tdgf1vsNeighbor_GOres_topBP.png", width=1600, height=680, bg="white")
heatmap(allRes_heatmap, Colv = NA, Rowv = NA, margins = c(20,46), col=heat.colors(108)[16:96],cexRow=3.5,cexCol=3,labCol=c("Fisher","Weight","KS","elimKS"), scale="none")

lgd_ = rep(NA, 16)
lgd_[c(1,16)] = c(as.character(as.integer(max(allRes_$X.log10.as.numeric.allRes.classicKS..)+0.5)),as.character(as.integer(min(allRes_$X.log10.as.numeric.allRes.weightFisher..)+0.5)))

legend(x = 0.02, y = 0.95,
       legend = lgd_,
       fill = rev(heat.colors(25)[1:17]), #, direction = -1),
       border = NA,
       bty = "n",
       y.intersp = 0.5,
       cex = 3, text.font = 1,
       title = "-log(p)    \n")
dev.off()



#############
# GO Chord plot on differentially expressed genes in Tdgf1 domain
#############
Tdgf1vsNeighbor.markers <- arrange(Tdgf1vsNeighbor.markers, desc((avg_log2FC)/p_val_adj))
gt_BP <- genesInTerm(GOdata)
candidates <- data.frame( rownames(Tdgf1vsNeighbor.markers), (Tdgf1vsNeighbor.markers$pct.1 + Tdgf1vsNeighbor.markers$pct.2) / 2, Tdgf1vsNeighbor.markers$avg_log2FC, Tdgf1vsNeighbor.markers$p_val_adj )
names(candidates) <- c("ID", "baseMean", "logFC", "adj_pval")
results_BP <- add_genes_to_terms(genes2GO = gt_BP, results = allRes_orig)
results_BP$Category <- rep("BP", nrow(results_BP))
#full_res <- rbind(results_BP, results_MF, results_CC)
full_res <- subset( results_BP, select = -c(`classicKS`, `weightFisher`,  `classicFisher`, `Rank in weightFisher`) )
full_res$genes <- as.factor(full_res$genes)
names(full_res) <- c("ID", "term", "Annotated", "Significant", "Expected", "adj_pval", "genes", "category")
full_res$adj_pval <- as.numeric(full_res$adj_pval)

GO_plot <- circle_dat_mod(full_res, candidates)
#GOCircle(GO_plot)

sub_IDs <- c('GO:0051017', 'GO:0030335', 'GO:0007411', 'GO:0008038', 'GO:0007229', 'GO:0048675' )
#sub_IDs <- c('actin filament bundle assembly', 'positive regulation of cell migration', 'axon guidance', 'neuron recognition', 'integrin-mediated signaling pathway', 'axon extension' )
#GOCircle(GO_plot, nsub = sub_IDs)
chord_genes = data.frame( rownames(Tdgf1vsNeighbor.markers), Tdgf1vsNeighbor.markers$avg_log2FC )
names(chord_genes) <- c("ID", "logFC" )

chord <- chord_dat( data = GO_plot, genes = chord_genes[1:350,], process = sub_IDs ) #

#head(chord)
chord_plot <- GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5, border.size = 0, ribbon.col = (brewer.pal (n=12,name="Set2"))[1:(0+length(sub_IDs))] )
ggsave( plot = chord_plot, height=260, width=220, filename = "deSoysa_Tdgf1vsNeighbor_GOBPchord.png", units = "mm",dpi = 300 )
#GOBubble(data=GO_plot, labels = 3, display="single", ID = T, table.legend = T, table.col = T)


