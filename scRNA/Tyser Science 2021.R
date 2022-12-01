

#############
#Begin
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
library(viridis)
library(reshape2)
library(igraph)

#load topGO
#library(BiocManager)
#BiocManager::install("topGO")
library(topGO)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)


#############
#Preferences
#############
seurat_umap_alpha_hex_string = "c0"
original_beach_colors <- c("#1a748e","#f0df99","#ecdfcf","#55c4d7","#5f6c24","#992915","#d38e31","#81c7f8","#bbc4af","#393430","#92dccd","#63fba9","#1155d4","#abb2ba","#eccd16","#7b5c52","#063581")
beach_colors <- c("#1a748e","#f0df99","#bbc4af","#ecdfcf","#5f6c24","#992915","#81c7f8","#63fba9","#55c4d7","#7b5c52","#abb2ba","#eccd16","#1155d4","#393430")
paired_colors <- brewer.pal(name = "Dark2", n=6)
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
paired_colors_alpha <- paste0( paired_colors, seurat_umap_alpha_hex_string )


#############
#Data download and read
#############
download.file("https://github.com/MarioniLab/mouseHeart_2020/blob/master/data/Mus_musculus.GRCm38.87.tsv")
download.file("https://www.science.org/doi/suppl/10.1126/science.abb2986/suppl_file/abb2986_datas3.csv")
#file.rename("abb2986_datas3.csv", "DataS3.csv")
download.file("https://content.cruk.cam.ac.uk/jmlab/mouseEmbryonicHeartAtlas/heartData_unbiased.RAW.Rds")

#read Tyser data
TyserData <- readRDS("heartData_unbiased.RAW.Rds")


#############
#Process raw Tyser data, and use SCT to integrate different batches
#############
#read sample metadata
meta <- read.csv("abb2986_datas3.csv", header = TRUE, stringsAsFactors = FALSE)
stopifnot(identical(colnames(TyserData), meta$cell))
meta$batch <- as.factor(paste0("batch_", meta$batch)) ## make 'batch' categorical
meta$stage <- as.factor(meta$stage) ## make 'stage' categorical
colnames(meta)[colnames(meta) == "cellID"] = "cells"
rownames(meta) <- meta$cells
meta <- meta[ , !(names(meta) %in% c("cells"))]

#replace Ensembl IDs with gene name
ann <- read.table("Mus_musculus.GRCm38.87.tsv", sep="\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
colnames(ann) <- c("gene", "chr", "start", "end", "strand")
TyserData$geneID <- ann$gene[match(rownames(TyserData), rownames(ann))]
TyserData = TyserData %>% group_by(geneID) %>% summarise_all(sum)
TyserData = head(TyserData,-1)
geneIDlist <- TyserData$geneID
TyserData <- TyserData[ , !(names(TyserData) %in% c("geneID"))]
rownames(TyserData) <- geneIDlist

#create Seurat Object and QC cleanup low-quality cells
TyserSeurat <- CreateSeuratObject(counts = TyserData, project = "TyserSeurat", meta.data = meta )
TyserSeurat[["percent.mt"]] <- PercentageFeatureSet(TyserSeurat, pattern = "^mt-")
TyserSeurat <- subset(TyserSeurat, subset = nFeature_RNA > 500 & nCount_RNA < 4e6 & percent.mt < 10)
VlnPlot(TyserSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#use SCT for batch correction and re-integrate batches
TyserSeurat.list <- SplitObject(TyserSeurat, split.by = "batch")
TyserSeurat.list <- lapply(X = TyserSeurat.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = TyserSeurat.list, nfeatures = 2500)
TyserSeurat.list <- PrepSCTIntegration( object.list = TyserSeurat.list[1:7], anchor.features = features, verbose = T )
Tyser.anchors <- FindIntegrationAnchors(object.list = TyserSeurat.list, normalization.method = "SCT", anchor.features = features)
TyserSeurat.corrected <- IntegrateData(anchorset = Tyser.anchors, normalization.method = "SCT")

#visualize batch and stage effects
TyserSeurat.corrected <- RunPCA(TyserSeurat.corrected, npcs = 30, verbose = FALSE)
TyserSeurat.corrected <- RunUMAP(TyserSeurat.corrected, reduction = "pca", dims = 1:20)
TyserSeurat.corrected <- FindNeighbors(TyserSeurat.corrected, nn.method="rann", reduction = "pca", dims = 1:20, nn.eps = 0.5)
TyserSeurat.corrected <- FindClusters(TyserSeurat.corrected, resolution = 0.2)


#############
#Cluster identification -- this is for reference only, not for making figure panels
#############
p1 <- DimPlot(TyserSeurat.corrected, reduction = "umap", group.by = "batch")
p2 <- DimPlot(TyserSeurat.corrected, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.25) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(TyserSeurat.corrected,  reduction = "umap", cols=beach_colors_alpha, group.by = "seurat_clusters", pt.size=0.25, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p2 + p3
#ggsave( plot = p2 + p3, height=200, width=400, filename = "Tyser_complete_UMAP.png", units = "mm",dpi = 300 )

DefaultAssay(TyserSeurat.corrected) <- "RNA"
table(Idents(TyserSeurat.corrected))


#FeaturePlot(TyserSeurat.corrected, features=c("Mef2c","Isl1","Tbx5","Pou5f1"),cols=c("#F2F2F2C0","#91391CC0"), pt.size=1.5,ncol=3,max.cutoff = 500)
#FeaturePlot(TyserSeurat.corrected, features=c("Pou3f1","Sox2","Epcam","Tek","Foxa2","Tubb3"),cols=c("#F2F2F2C0","#91391CC0"), pt.size=1.5,ncol=3,max.cutoff = 500)
#FeaturePlot(TyserSeurat.corrected, features = c("Isl1", "Mef2c", "Nkx2-5","Tnnt2","Tbx5","Hand1","Hand2","Smarcd3","Mab21l2"), cols = c("#FEFEFE", "black"))
#VlnPlot(TyserSeurat.corrected, features = c("Smarcd3","Mef2c","Isl1","Nkx2-5","Twist1","Snai1"))#, group.by = "stage")
#table(TyserSeurat.corrected@meta.data$seurat_clusters, TyserSeurat.corrected@meta.data$stage)
#VlnPlot(TyserSeurat.corrected,features=c("percent.mt","nFeature_RNA","nCount_RNA"))
#FindMarkers(TyserSeurat.corrected,ident.1 = 4)

#############
#Store integrated object
#############
save(TyserSeurat.corrected,file="TyserSeurat.corrected.Rds")
#load("TyserSeurat.corrected.Rds")
table(TyserSeurat.corrected@meta.data$seurat_clusters, TyserSeurat.corrected@meta.data$stage)


#############
#Subset clusters 0,1,3 (i.e. cardiac progenitors and early CMs) within timepoints -1 to 2, 
#############
Idents(TyserSeurat.corrected) <- "seurat_clusters"
keep_clusters <- c(0,1,3)
new_object = subset(TyserSeurat.corrected, idents = keep_clusters ) #, subset=Mef2c>0 )
Idents(new_object) <- "stage"
new_object = subset(new_object, idents = c( "-1","0","1","2" )) #, subset = Smarcd3 > 5 )
Idents(new_object) <- factor(x = Idents(new_object), levels = sort(levels(new_object)))


#############
#Now, for visualization, need to re-cluster using SCT data
#############
DefaultAssay(new_object) <- "SCT"
new_object <- FindVariableFeatures(new_object, selection.method = "vst", nfeatures = 2000)

new_object <- RunPCA(new_object, features = VariableFeatures(object = new_object))
new_object <- FindNeighbors(new_object, nn.method="rann", reduction="pca", dims = 1:10, nn.eps = 0.5)
new_object <- FindClusters(new_object, resolution = 0.5)
new_object <- RunUMAP(new_object, dims = 1:10)


#############
#Rename identities
#############
#Rename stages
Idents(new_object) <- "stage"
Idents(new_object) <- factor(x = Idents(new_object), levels = sort(levels(new_object)))
new_object[["old.ident"]] <- Idents(new_object)
new_object <- RenameIdents(object = new_object, 
                           '-1' = "-1 (~ EHF)", 
                           '0' = "0 (~ LHF)", 
                           '1' = "1 (~ 1-2 somites)", 
                           '2' = "2 (~ 3-4 somites)"
)
Idents(new_object) <- factor(x = Idents(new_object), levels = sort(levels(new_object)))
new_object[["stage_names"]] <- Idents(new_object)


#Rename clusters and colorize to match deSoysa
#FeaturePlot(new_object,features=c("Mef2c","Isl1","Osr1","Hoxb1","Tdgf1","Foxc1","Foxc2","Dkk1","Isl1","Meox1","Six1","Prdm1"))

Idents(new_object) <- "seurat_clusters"
new_object[["old.ident"]] <- Idents(new_object)
new_object <- RenameIdents(object = new_object, 
                           '2' = "Msx1+ mesoderm",
                           '6' = "Foxc2+ mesoderm",
                           '3' = "First heart field (FHF)", 
                           '1' = "Anterior second heart field (aSHF)",
                           '5' = "Posterior second heart field (pSHF)",
                           '4' = "Differentiating SHF CMs",
                           '0' = "Cardiomyocytes (CMs)"
                           )
new_object[["cluster_names"]] <- Idents(new_object)

library(scales)
show_col(beach_colors)

#re-order the cluster colors to jibe with deSoysa
beach_colors <- c(beach_colors[12],beach_colors[1],beach_colors[9],beach_colors[5],beach_colors[3],beach_colors[6],beach_colors[2])
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
show_col(beach_colors)



#############
#Visualize final object
#############
p2 <- DimPlot(new_object, reduction = "umap", cols=paired_colors_alpha, group.by = "stage_names", pt.size=0.6) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=16))
p3 <- DimPlot(new_object, reduction = "umap", cols=beach_colors_alpha, group.by="cluster_names", label = FALSE, repel = FALSE, pt.size=0.6) + NoAxes() + labs(title = "Cluster")  + theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=16))
p2 + p3
ggsave( plot = p2 + p3, height=100, width=400, filename = "Tyser_final_UMAP.png", units = "mm",dpi = 300)

# run sctransform on the data subset
new_object <- SCTransform(new_object, vars.to.regress = c("percent.mt"), verbose = FALSE)




#############
#Pseudo-bulk differential expression between earliest timepoints -1 and 0, with SCT data (scaled RNA assay yields similar results)
#############
Idents(new_object) <- "stage"

diffex_early_vs_mid <- FindMarkers(subset(new_object,subset=Mef2c>0), ident.1 = c("-1"), ident.2 = c("0"), base=2, logfc.threshold = 0)
save( diffex_early_vs_mid, file="Mef2c_-1vs0.Rds" )
#load("Mef2c_-1vs0.Rds")
diffex_early_vs_mid <- diffex_early_vs_mid[grep("^Gm", row.names(diffex_early_vs_mid),invert=TRUE),]
diffex_early_vs_mid <- diffex_early_vs_mid[grep("^RP", row.names(diffex_early_vs_mid),invert=TRUE),]
diffex_early_vs_mid <- diffex_early_vs_mid[grep("^Rp", row.names(diffex_early_vs_mid),invert=TRUE),]
diffex_early_vs_mid <- diffex_early_vs_mid[grep("-ps\\d?$", row.names(diffex_early_vs_mid),invert=TRUE),]
write.csv(diffex_early_vs_mid,file="Tyser_-10vs12_HM.csv" )

#make list for heatmap -- top 20 up and down between stages -1 and 0
top15_log2FC <- diffex_early_vs_mid %>% top_n(20, avg_log2FC ) #* (pct.1-pct.2) ) #* pct.1/(pct.2) )
top15_neglog2FC <- diffex_early_vs_mid %>% top_n(20, -avg_log2FC  ) #* (pct.1-pct.2)  ) #*pct.2/(pct.1) )
#top30_total <- rbind( top15_log2FC,top15_neglog2FC )
top30_total <- diffex_early_vs_mid %>% top_n(40, avg_log2FC * (pct.1-pct.2) ) %>% arrange(desc(avg_log2FC))

#Idents(new_object) <- "stage_names"
HM <- DoHeatmap(object = subset(new_object,idents=c("-1","0")), group.by="stage", slot="data", group.colors = paired_colors[1:2], features = rownames(top30_total), label=FALSE)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 16),legend.key.size = unit(1.55,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=300, width=200, filename = "Tyser_-1vs0_HM.png", units = "mm",dpi = 100 )


#############
#Perform GO analysis for BP on differentially expressed genes
#############
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = setNames(diffex_early_vs_mid$p_val_adj,rownames(diffex_early_vs_mid)),
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
write.csv(allRes_orig,file="GOres_clusters_early_mid.csv")

#now, make GO p-value heatmap for top XXX terms
allRes <- head(allRes_orig,  12 )
allRes_ <- data.frame(-log10(as.numeric(allRes$classicFisher)),-log10(as.numeric(allRes$weightFisher)),-log10(as.numeric(allRes$classicKS)),-log10(as.numeric(allRes$elimKS)))
rownames(allRes_) <- allRes$Term
allRes_heatmap <- as.matrix(allRes_)
allRes_heatmap <- allRes_heatmap[nrow(allRes_heatmap):1, ] #avoid using reorderfun, since that is more complicated


#############
#GO Heatmap
#############
plot.new()
png(filename="GOres013_top.png", width=1400, height=680, bg="white")
heatmap(allRes_heatmap, Colv = NA, Rowv = NA, margins = c(20, 45), col=heat.colors(108)[16:96],cexRow=3.5,cexCol=3,labCol=c("Fisher","Weight","KS","elimKS"), scale="none")

lgd_ = rep(NA, 16)
lgd_[c(1,16)] = c(as.character(as.integer(max(allRes_$X.log10.as.numeric.allRes.classicKS..)+0.5)),as.character(as.integer(min(allRes_$X.log10.as.numeric.allRes.weightFisher..)+0.5)))

legend(x = 0.05, y = 0.95,
       legend = lgd_,
       fill = rev(heat.colors(25)[1:17]), #, direction = -1),
       border = NA,
       bty = "n",
       y.intersp = 0.5,
       cex = 3, text.font = 1,
       title = "-log(p)    \n")
dev.off()



#############
# Dotplots GO_0010718 and 0007043
#############
Idents(new_object) <- "stage"
DP <- DotPlot(object = subset(new_object,idents=c("-1","0","1","2")), cols = c("lightgrey", "blue"), group.by="stage_names", features = genesInTerm(GOdata)["GO:0010718"] ) +coord_flip() +RotatedAxis()  + scale_color_viridis() + theme(axis.text.x = element_text(angle=90, vjust=0.4, hjust=1)) + xlab("") + ylab("")
ggsave( plot = DP, height=250, width=105, filename = "Tyser_013_Dotplt_GO0010718.png", units = "mm",dpi = 300 )

DP <- DotPlot(object = subset(new_object,idents=c("-1","0","1","2")), cols = c("lightgrey", "blue"), group.by="stage_names", features = genesInTerm(GOdata)["GO:0007043"] ) +coord_flip() +RotatedAxis()  + scale_color_viridis() + theme(axis.text.x = element_text(angle=90, vjust=0.4, hjust=1)) + xlab("") + ylab("")
ggsave( plot = DP, height=460, width=105, filename = "Tyser_013_Dotplt_GO0007043.png", units = "mm",dpi = 300 )



#############
# Vlns Down GO_0010718 and 0007043
#############
#Idents(new_object) <- factor(x = Idents(new_object), levels = sort(levels(new_object)))
Idents(new_object) <- "stage_names"
DefaultAssay(new_object) <- "SCT"
#vln <- VlnPlot(new_object, cols=paired_colors, features = c("Foxc1","Snai1","Twist1","Notch1","Alx1","Zfp703","Smad3"), ncol = 2, log=TRUE, combine = FALSE, pt.size=0)
vln <- VlnPlot(subset(new_object,subset=Mef2c>0), cols=paired_colors, features = c("Smad3","Foxc1","Snai1","Twist1","Notch1","Epha2"), ncol = 3, log=FALSE, combine = FALSE, pt.size=0)
for(i in 1:length(vln)) {
  if ( i == length(vln) ) {
    vln[[i]] <- vln[[i]] + NoAxes() + NoLegend()  + theme(plot.margin = margin(0, 0, 0, 2, "mm"), plot.title = element_text(size=12,face="plain" ), legend.key.size = unit(3,'mm'), legend.text = element_text(size=9) )
  } else {
    vln[[i]] <- vln[[i]] + NoAxes() + NoLegend() + theme(plot.margin = margin(0, 8, 0, 0, "mm"), plot.title = element_text(size=12,face="plain" ) )
  }
  vln[[i]]$layers[[1]]$aes_params$size = 0
}
#wrap_plots(vln,nrow=vln_rows)
ggsave( plot = wrap_plots(vln,nrow=1) , height=50, width=140, filename = "Tyser_early_mid_Mef2c_Vlns_Down.pdf", units = "mm",dpi = 300,useDingbats = FALSE)

#############
# Vlns Up GO_0010718 and 0007043
#############
#Idents(new_object) <- factor(x = Idents(new_object), levels = sort(levels(new_object)))
Idents(new_object) <- "stage_names"
DefaultAssay(new_object) <- "SCT"
#vln <- VlnPlot(new_object, cols=paired_colors, features = c("Foxc1","Snai1","Twist1","Notch1","Alx1","Zfp703","Smad3"), ncol = 2, log=TRUE, combine = FALSE, pt.size=0)
vln <- VlnPlot(subset(new_object,subset=Mef2c>0), cols=paired_colors, features = c("Cdh2","Vcl","Jup","Afdn","Crb2"), ncol = 3, log=FALSE, combine = FALSE, pt.size=0) #"Ephb2"
for(i in 1:length(vln)) {
  if ( i == length(vln) ) {
    vln[[i]] <- vln[[i]] + NoAxes() + theme(plot.margin = margin(0, 0, 0, 2, "mm"), plot.title = element_text(size=12,face="plain" ), legend.key.size = unit(3,'mm'), legend.text = element_text(size=9) )
  } else {
    vln[[i]] <- vln[[i]] + NoAxes() + NoLegend() + theme(plot.margin = margin(0, 8, 0, 0, "mm"), plot.title = element_text(size=12,face="plain" ) )
  }
  vln[[i]]$layers[[1]]$aes_params$size = 0
}
#wrap_plots(vln,nrow=vln_rows)
ggsave( plot = wrap_plots(vln,nrow=1) , height=50, width=148, filename = "Tyser_early_mid_Mef2c_Vlns_Up.pdf", units = "mm",dpi = 300,useDingbats = FALSE)


#############
#Visualize final object feature plots
#############
genes_feature_plot = c("Nkx2-5","Isl1","Tbx5","Mab21l2","Msx1","Foxc2","Mef2c","Foxc1","Twist1","Snai1","Tdgf1","Tnnt2","Tbx1")
#FeaturePlot(new_object, features=genes_feature_plot )

vln <- FeaturePlot(new_object, features=genes_feature_plot,cols=c("#000000E0","#FFFFFFE0"), pt.size=0.6,ncol=3,combine=FALSE)
for(i in 1:length(vln)) {
  vln[[i]] <- vln[[i]] + NoLegend() + NoAxes()
  #vln[[i]]$layers[[1]]$aes_params$size = 0
  ggsave( plot = vln[[i]], height=100, width=125, filename = paste("Tyser_013_", genes_feature_plot[[i]],".png",sep=""), units = "mm",dpi = 300)
}
#ggsave( plot = wrap_plots(vln,nrow=1), height=100, width=250, filename = "Tyser_013_FeatPlts.png", units = "mm",dpi = 300)
#wrap_plots(vln,nrow=1)

#############
#Visualize final object cluster heatmap
#############
#Idents(new_object) <- "seurat_clusters"
Idents(new_object) <- "cluster_names"
all.markers <- FindAllMarkers(object = new_object)
all.markers$ranking = all.markers$avg_log2FC * all.markers$pct.1 * all.markers$pct.1 / all.markers$pct.2
top15_logFC <- all.markers %>% group_by(cluster) %>% top_n(12, avg_log2FC)
top15_uniqueness <- all.markers %>% group_by(cluster) %>% top_n(12, ranking)

HM <- DoHeatmap(object = new_object, features = top15_logFC$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 16),legend.key.size = unit(1.55,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=450, width=300, filename = "Tyser_013_clusterHM.png", units = "mm",dpi = 100 )




#############
#Analyze reverse EMT process by investigating potential Snai1/Twist1 GRNs -- Correlation analysis
#############
Idents(new_object) <- "stage"
DefaultAssay(new_object) <- "RNA"
new_object2 <- subset( new_object, idents = c("-1","0") )
new_object2 <- NormalizeData(new_object2)
new_object2 <- ScaleData(new_object2)

#correlation analysis
new_object2.matrix  <- new_object2@assays$RNA@scale.data#begins as sparse matrix)

new_object2.matrix <- new_object2.matrix[grep("^Gm", row.names(new_object2.matrix),invert=TRUE),]
new_object2.matrix <- new_object2.matrix[grep("^RP", row.names(new_object2.matrix),invert=TRUE),]
new_object2.matrix <- new_object2.matrix[grep("^Rp", row.names(new_object2.matrix),invert=TRUE),]
new_object2.matrix <- new_object2.matrix[grep("-ps\\d?$", row.names(new_object2.matrix),invert=TRUE),]


new_object2.matrix.abs <- abs(new_object2.matrix)
keep <- which( rowSums(new_object2.matrix.abs) > ncol(new_object2.matrix)/3 ) #convert to dense matrix

for ( ii in 1:length(cor_genes)) {
  keep <- append(keep,which(row.names(new_object2.matrix)==cor_genes[ii]))  
}
keep_mtx <- t(new_object2.matrix[keep,])
cor <- cor(keep_mtx, method="pearson")

save(cor,file="Tyser_013_-1and0_cor_matrix.Rds")
#load("Tyser_013_-1and0_cor_matrix.Rds")


#############
#Analyze reverse EMT process by investigating potential Snai1/Twist1 GRNs -- Heatmap drawing
#############
cor_genes = c( "Twist1", "Snai1" )

#functions courtesy of http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

combined_gene_names = c()
orientation = "horizontal"
#now, make heatmaps to unveil GRN
for ( ii in 1:length(cor_genes)) {
  Thiscor = cor[c(cor_genes[ii]),]
  Thiscor <- sort(Thiscor)
  combined_gene_names <- append( combined_gene_names, names(tail( Thiscor, n=70)) )
  Thiscor <- sort(abs(Thiscor))
  #print(cor_genes[ii])
  #gene_names <- append(names(tail( Thiscor, n=25)),names(head( Thiscor, n=25)))
  #gene_names <- append(gene_names,names(tail( Thiscor, n=30)))
  gene_names <- names(tail( Thiscor, n=50))
  
  if ( orientation == "horizontal" ) {
    final_cor <- cor[c(cor_genes[ii]),gene_names,drop=FALSE]
    tri <- get_upper_tri(final_cor)
  } else if ( orientation == "vertical") {
    final_cor <- cor[gene_names,c(cor_genes[ii]),drop=FALSE]
    tri <- get_lower_tri(final_cor)
  } else {
    final_cor <- cor[gene_names,gene_names]
    if ( (ii %% 2) == 0 ) { #even, do bottom triangle
      tri <- get_lower_tri(final_cor)
    } else { #odd, do top triangle
      tri <- get_upper_tri(final_cor)
    }
  }
  melted_cormat <- melt(tri, na.rm = TRUE)

  
  HM <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_brewer(palette = "Accent",type="div") +
    scale_fill_gradient2(low = "#466baf", high = "#e34e33", mid = "#ecf8de",
                         midpoint = 0, limit = c(-1,1),
                         name="Pearson\nCorrelation") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = if((ii %% 2) == 0) 270 else 90, vjust = if((ii %% 2) == 0) -10 else 0.5,
                                     hjust = 1 ), panel.grid.major = element_blank()) +
    scale_y_discrete(position = if((ii %% 2) == 0) "left" else "right") +
    scale_x_discrete(position = if((ii %% 2) == 0) "top" else "bottom") +
    ylab("") + 
    xlab("") +
    coord_fixed()
  
  ggsave( plot = HM, height=180, width=180, filename = paste("Tyser_013_-1and0_cor_HM_", cor_genes[ii],".png",sep=""), units = "mm",dpi = 300 )
  
}
combined_gene_names <- distinct(as.data.frame(combined_gene_names))$combined_gene_names

#############
#Draw network connectivity for Twist1/Snai1 GRN
#############
g <- graph.adjacency(
  cor[combined_gene_names,combined_gene_names],
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

# Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "blue"

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "red"

# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

# Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.25)])


# Remove any vertices remaining that have no edges
g <- delete.vertices(g, igraph::degree(g)==0)

# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

## Change shape of graph vertices
V(g)$shape <- "circle"

# Change colour of graph vertices
V(g)$color <- "skyblue"
V(g)$font

# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(cor[combined_gene_names,combined_gene_names], 1, mean)) + 1.0) * 5

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 15

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="unweighted")

mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
#mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

# Plot the tree object
plot.new()
png(filename="Twist1-Snai1_GNR.png", width=1980, height=1980, bg="white")
plot.igraph(
  mst,
  #layout=layout.fruchterman.reingold,
  layout=layout_with_lgl,
  edge.curved=FALSE,
  vertex.size=vSizes,
  vertex.label.dist=0, #-0.5,
  vertex.label.degree = 1,,
  vertex.label.color="black",
  vertex.label.cex=3,   
  vertex.label.family="Arial",
  vertex.frame.color= "white",
  asp=FALSE,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="")
dev.off()

