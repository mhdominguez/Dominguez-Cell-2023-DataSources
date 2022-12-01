

#############
#Begin
#############
#load Seurat
library(Seurat)
library(patchwork)

#library(BiocManager)
#BiocManager::install("biomaRt")
library(biomaRt)

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

#load BioMart
#library(BiocManager)
#BiocManager::install("biomaRt")
library(biomaRt)

#load topGO
#library(BiocManager)
#BiocManager::install("topGO")
library(topGO)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(GOplot)

#############
#Preferences
#############
seurat_umap_alpha_hex_string = "c0"
beach_colors <- c("#1a748e","#f0df99","#ecdfcf","#55c4d7","#5f6c24","#992915","#d38e31","#81c7f8","#bbc4af","#393430","#92dccd","#63fba9","#1155d4","#abb2ba","#eccd16","#7b5c52","#063581")
paired_colors <- brewer.pal(name = "Paired", n=12)
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
paired_colors_alpha <- paste0( paired_colors, seurat_umap_alpha_hex_string )


#############
#Data download and read
#############
download.file("http://singlecell.stemcells.cam.ac.uk/data/HTSeq_counts_genesonly_combined.txt")
download.file("http://singlecell.stemcells.cam.ac.uk/data/metadata_combined.txt")

#read Lescroart data
Lescroart.matrix <- read.csv("HTSeq_counts_genesonly_combined.txt", header = TRUE, sep="\t", check.names = FALSE)
meta <- read.csv("metadata_combined.txt", header = TRUE, sep="\t", stringsAsFactors=TRUE)
rownames(meta) <- str_replace(meta$ID, "\\+", "") #get rid of plus signs in cell name
#rownames(meta) <- meta$ID
meta = meta[,-1]



#############
#Convert useless ENSEMBL ids
#############
m_ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset="mmusculus_gene_ensembl")
G_list_Lescroart <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description"),rownames(Lescroart.matrix),mart= m_ensembl,useCache = FALSE)

#Create data frames ready for Seurat
Lescroart.df <- as.data.frame(Lescroart.matrix)
Lescroart.df$geneID <- G_list_Lescroart$external_gene_name[match(rownames(Lescroart.df), G_list_Lescroart$ensembl_gene_id)]
Lescroart.df = Lescroart.df %>% group_by(geneID) %>% summarise_all(sum)
geneIDlist <- Lescroart.df$geneID
#Lescroart.df1 = head(Lescroart.df,-1)
Lescroart.df <- Lescroart.df[ , !(names(Lescroart.df) %in% c("geneID"))]
geneIDlist[1] = "EmptyGene"
geneIDlist[length(geneIDlist)] = "UnknownGene"
rownames(Lescroart.df) <-geneIDlist
names(Lescroart.df) <- str_replace(names(Lescroart.df), "\\+", "") #get rid of plus signs in cell name



#############
#Process raw Lescroart data, and use SCT to integrate different batches
#############
#create Seurat Object and QC cleanup low-quality cells
LescroartSeurat <- CreateSeuratObject(counts = Lescroart.df, project = "LescroartSeurat", meta.data = meta )
LescroartSeurat[["percent.mt"]] <- PercentageFeatureSet(LescroartSeurat, pattern = "^mt-")
Idents(LescroartSeurat) <- "Study"
VlnPlot(LescroartSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LescroartSeurat <- subset(LescroartSeurat, subset = nFeature_RNA > 2500 & nCount_RNA < 5e6 & percent.mt < 10)

#use SCT for batch correction and re-integrate batches
LescroartSeurat.list <- SplitObject(LescroartSeurat, split.by = "Study")
LescroartSeurat.list <- lapply(X = LescroartSeurat.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = LescroartSeurat.list, nfeatures = 2500)
LescroartSeurat.list <- PrepSCTIntegration( object.list = LescroartSeurat.list, anchor.features = features, verbose = T )
Lescroart.anchors <- FindIntegrationAnchors(object.list = LescroartSeurat.list, normalization.method = "SCT", anchor.features = features)
LescroartSeurat.corrected <- IntegrateData(anchorset = Lescroart.anchors, normalization.method = "SCT")

#visualize batch and stage effects
LescroartSeurat.corrected <- RunPCA(LescroartSeurat.corrected, npcs = 30, verbose = FALSE)
ElbowPlot(LescroartSeurat.corrected)
LescroartSeurat.corrected <- RunUMAP(LescroartSeurat.corrected, reduction = "pca", dims = 1:15)
LescroartSeurat.corrected <- FindNeighbors(LescroartSeurat.corrected, nn.method="rann", reduction = "pca", dims = 1:15, nn.eps = 0.5)
LescroartSeurat.corrected <- FindClusters(LescroartSeurat.corrected, resolution = 0.5)


#############
#Cluster identification -- this is for reference only, not for making figure panels
#############
p1 <- DimPlot(LescroartSeurat.corrected, reduction = "umap", cols=brewer.pal(3,"Set2"), group.by = "Study", pt.size=1.25) + NoAxes() + labs(title = "Batch") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(LescroartSeurat.corrected, reduction = "umap", cols=paired_colors_alpha, group.by = "Condition", pt.size=1.25) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(LescroartSeurat.corrected,  reduction = "umap", cols=beach_colors_alpha, group.by = "seurat_clusters", pt.size=1.25, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2 + p3
#ggsave( plot = p2 + p3, height=200, width=400, filename = "Lescroart_complete_UMAP.png", units = "mm",dpi = 300 )
table(LescroartSeurat.corrected@meta.data$seurat_clusters, LescroartSeurat.corrected@meta.data$Condition)


#############
#Store integrated object
#############
save(LescroartSeurat.corrected,file="LescroartSeurat.corrected.Rds")
#load("LescroartSeurat.corrected.Rds")


#############
#Identify clusters
#############
DefaultAssay(LescroartSeurat.corrected) <- "SCT"
#LescroartSeurat.corrected <- NormalizeData(LescroartSeurat.corrected)
#LescroartSeurat.corrected <- ScaleData(LescroartSeurat.corrected)
LescroartSeurat.corrected <- PrepSCTFindMarkers(LescroartSeurat.corrected )
#FeaturePlot(LescroartSeurat.corrected, features=c("Msx1","Foxc2","Mesp1","Sox7","Wnt2b","Eomes","Pou5f1","Nodal","Foxa2","Mef2c","Smarcd3"))

Idents(LescroartSeurat.corrected) <- "seurat_clusters"
all.markers <- FindAllMarkers(object = LescroartSeurat.corrected)

all.markers$ranking = all.markers$avg_log2FC * all.markers$pct.1 * all.markers$pct.1 / all.markers$pct.2
top15_logFC <- all.markers %>% group_by(cluster) %>% top_n(12, avg_log2FC/p_val)
top15_uniqueness <- all.markers %>% group_by(cluster) %>% top_n(12, ranking)

HM <- DoHeatmap(object = LescroartSeurat.corrected, features = top15_uniqueness$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 16),legend.key.size = unit(1.55,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=450, width=300, filename = "Lescroart_clusterHM.png", units = "mm",dpi = 100 )



#############
#Rename identities
#############
#Rename clusters and colorize to match deSoysa
Idents(LescroartSeurat.corrected) <- "seurat_clusters"
LescroartSeurat.corrected[["old.ident"]] <- Idents(LescroartSeurat.corrected)
LescroartSeurat.corrected <- RenameIdents(object = LescroartSeurat.corrected, 
                           '0' = "Msx1+ mesoderm",
                           '1' = "Foxc2+ mesoderm",
                           '2' = "Primitive Streak", 
                           '3' = "Extraembryonic Mesoderm",
                           '4' = "Lateral mesoderm",
                           '5' = "Epiblast"

                           )
LescroartSeurat.corrected[["cluster_names"]] <- Idents(LescroartSeurat.corrected)
#Idents(LescroartSeurat.corrected) <- "cluster_names"

library(scales)
show_col(beach_colors)

#re-order the cluster colors to jibe with deSoysa/Tyser
beach_colors <- c(beach_colors[15],beach_colors[1],beach_colors[9],beach_colors[4],beach_colors[10],beach_colors[5])
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
show_col(beach_colors)


#############
#Complete object visualization -- FYI only
#############
p3 <- DimPlot(LescroartSeurat.corrected,  reduction = "umap", cols=beach_colors_alpha, group.by = "cluster_names", pt.size=1.25, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))

#p4 <- ggplot(LescroartSeurat.corrected@meta.data, aes(orig.ident, fill=cluster_names))+geom_bar(stat="count") + scale_fill_manual(values=beach_colors) + RotatedAxis() + xlab("Sample") + ylab("Cell count") + NoLegend()
p4a <- ggplot(LescroartSeurat.corrected@meta.data, aes(Condition, fill=cluster_names))+geom_bar(stat="count",position="fill") + scale_fill_manual(values=beach_colors) + RotatedAxis() + xlab("Sample") + ylab("Cell %") + NoLegend()
#p5 <- ggplot(LescroartSeurat.corrected@meta.data, aes(seurat_clusters, fill=orig.ident))+geom_bar(stat="count") + scale_fill_manual(values=paired_colors)+ xlab("Cluster") + RotatedAxis() + ylab("Cell count") + NoLegend()
p5a <- ggplot(LescroartSeurat.corrected@meta.data, aes(cluster_names, fill=Condition))+geom_bar(stat="count",position="fill") + scale_fill_manual(values=paired_colors) + RotatedAxis() + xlab("Cluster") + ylab("Cell %") + NoLegend()
p6a <- ggplot(LescroartSeurat.corrected@meta.data, aes(Study, fill=cluster_names))+geom_bar(stat="count",position="fill") + scale_fill_manual(values=beach_colors) + RotatedAxis() + xlab("Sample") + ylab("Cell %") + NoLegend()
p6b <- ggplot(LescroartSeurat.corrected@meta.data, aes(Study, fill=Condition))+geom_bar(stat="count",position="fill") + scale_fill_manual(values=paired_colors) + RotatedAxis() + xlab("Sample") + ylab("Cell %") + NoLegend()

#all.markers <- FindAllMarkers(object = LescroartSeurat.corrected)
#top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)

HM <- DoHeatmap(object = LescroartSeurat.corrected, features = top15_uniqueness$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15),legend.key.size = unit(1.5,"line")) + NoLegend()#slim.col.label = TRUE, remove.key = TRUE) + coord

#pPlot <- (p1/p2/p3) | ((p6a/p6b)/(p5a|p4a))
pPlot <- (p1|p2) / (p3|(p5a|p4a))
ggsave( plot = (pPlot/HM) + plot_layout(heights=c(1,1,2.5)), height=750, width=400, filename = "Lescroart 2018 clustering.pdf", units = "mm",dpi = 300 )




#############
#Visualize final object
#############
p2 <- DimPlot(LescroartSeurat.corrected, reduction = "umap", cols=paired_colors_alpha, group.by = "Condition", pt.size=1) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=16))
p3 <- DimPlot(LescroartSeurat.corrected, reduction = "umap", cols=beach_colors_alpha, group.by="cluster_names", label = FALSE, repel = FALSE, pt.size=1) + NoAxes() + labs(title = "Cluster")  + theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=16))
p2 + p3
ggsave( plot = p2 + p3, height=120, width=300, filename = "Lescroart_UMAP.png", units = "mm",dpi = 300)


#############
#Pseudo-bulk differential expression between KO and WT
#############
Idents(LescroartSeurat.corrected) <- "Condition"

Mesp1KOvsWT <- FindMarkers(LescroartSeurat.corrected, ident.2 = c("E6.75_KO"), ident.1 = c("E6.75_WT"), base=2, logfc.threshold = 0)
#Mesp1KOvsWT <- Mesp1KOvsWT[grep("^Gm", row.names(Mesp1KOvsWT),invert=TRUE),]
#Mesp1KOvsWT <- Mesp1KOvsWT[grep("^RP", row.names(Mesp1KOvsWT),invert=TRUE),]
#Mesp1KOvsWT <- Mesp1KOvsWT[grep("^Rp", row.names(Mesp1KOvsWT),invert=TRUE),]
#Mesp1KOvsWT <- Mesp1KOvsWT[grep("-ps\\d?$", row.names(Mesp1KOvsWT),invert=TRUE),]

#make list for heatmap -- top 20 up and down between stages -1 and 0
top15_log2FC <- Mesp1KOvsWT %>% top_n(35, avg_log2FC)
top15_neglog2FC <- Mesp1KOvsWT %>% top_n(5, -avg_log2FC)
top30_total <- rbind( top15_log2FC,top15_neglog2FC )
Mesp1KOvsWT$ranking = abs(Mesp1KOvsWT$avg_log2FC) * max(Mesp1KOvsWT$pct.1 , Mesp1KOvsWT$pct.2)
top20_uniqueness <- Mesp1KOvsWT %>% top_n(40, ranking)

VlnPlot(LescroartSeurat.corrected, group.by="Condition", features=c("Mesp1","Fgf3","Fgf5","Pdgfra","Fgfr1","Fgfr2","Wnt5a","Kdr") )
#Idents(new_object) <- "stage_names"
HM <- DoHeatmap(object = subset(LescroartSeurat.corrected,idents=c("E6.75_KO","E6.75_WT")), group.by="Condition", slot="data", group.colors = paired_colors[2:3], features = rev(rownames(top30_total)), label=FALSE)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 16),legend.key.size = unit(1.55,"line")) #slim.col.label = TRUE, remove.key = TRUE)
ggsave( plot = HM, height=300, width=200, filename = "Lescroart_KOvsWT_HM.png", units = "mm",dpi = 100 )






#############
#Perform GO analysis for CC on differentially expressed genes in KO/WT
#############
GOdata <- new("topGOdata",
              ontology = "CC", # use biological process ontology
              allGenes = setNames(Mesp1KOvsWT$p_val_adj,rownames(Mesp1KOvsWT)),
              geneSel = function(p) p < 0.05,
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol", nodeSize = 10)

resultWeight01 <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allGO <- usedGO(GOdata)

allRes_orig <- GenTable(GOdata, classicFisher = resultFisher,weightFisher=resultWeight01,
                        classicKS = resultKS, elimKS = resultKS.elim,
                        orderBy = "elimKS", ranksOf = "weightFisher", topNodes=length(allGO), numChar=1000)
write.csv(allRes_orig,file="Lescroart_KOvsWT_GOres_CC.csv")


#now, make GO p-value heatmap for top XXX terms
allRes <- head(allRes_orig,  15 )
allRes_ <- data.frame(-log10(as.numeric(allRes$classicFisher)),-log10(as.numeric(allRes$weightFisher)),-log10(as.numeric(allRes$classicKS)),-log10(as.numeric(allRes$elimKS)))
rownames(allRes_) <- allRes$Term
allRes_heatmap <- as.matrix(allRes_)
allRes_heatmap <- allRes_heatmap[nrow(allRes_heatmap):1, ] #avoid using reorderfun, since that is more complicated


#############
#GO CC Heatmap on differentially expressed genes in KO/WT
#############
plot.new()
png(filename="Lescroart_KOvsWT_GOres_topCC.png", width=1600, height=800, bg="white")
heatmap(allRes_heatmap, Colv = NA, Rowv = NA, margins = c(20,40), col=heat.colors(108)[16:96],cexRow=3.5,cexCol=2.8,labCol=c("Fisher","Weight","KS","elimKS"), scale="none")

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
#Perform GO analysis for BP on differentially expressed genes in KO/WT
#############
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = setNames(Mesp1KOvsWT$p_val_adj,rownames(Mesp1KOvsWT)),
              geneSel = function(p) p < 0.05,
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol", nodeSize = 10)

resultWeight01 <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allGO <- usedGO(GOdata)

allRes_orig <- GenTable(GOdata, classicFisher = resultFisher,weightFisher=resultWeight01,
                        classicKS = resultKS, elimKS = resultKS.elim,
                        orderBy = "elimKS", ranksOf = "weightFisher", topNodes=length(allGO), numChar=1000)
write.csv(allRes_orig,file="Lescroart_KOvsWT_GOres_BP.csv")


#now, make GO p-value heatmap for top XXX terms
allRes <- head(allRes_orig,  10 )
allRes_ <- data.frame(-log10(as.numeric(allRes$classicFisher)),-log10(as.numeric(allRes$weightFisher)),-log10(as.numeric(allRes$classicKS)),-log10(as.numeric(allRes$elimKS)))
rownames(allRes_) <- allRes$Term
allRes_heatmap <- as.matrix(allRes_)
allRes_heatmap <- allRes_heatmap[nrow(allRes_heatmap):1, ] #avoid using reorderfun, since that is more complicated


#############
#GO BP Heatmap on differentially expressed genes in KO/WT
#############
plot.new()
png(filename="Lescroart_KOvsWT_GOres_topBP.png", width=1600, height=600, bg="white")
heatmap(allRes_heatmap, Colv = NA, Rowv = NA, margins = c(20,40), col=heat.colors(108)[16:96],cexRow=3.5,cexCol=2.8,labCol=c("Fisher","Weight","KS","elimKS"), scale="none")

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
#GO:0010640 dotplot
#############
DP <- DotPlot(object = subset(LescroartSeurat.corrected,idents=c("E6.75_WT","E6.75_KO")), cols = c("lightgrey", "blue"), group.by="Condition", features = genesInTerm(GOdata)["GO:0008543"] ) +RotatedAxis() # +coord_flip()  #, cols = paired_colors[1:2],  )
ggsave( plot = DP, height=60, width=500, filename = "Lescroart_KOvsWT_Dotplt_GO0008543.png", units = "mm",dpi = 300 )


#############
#GO:0030178 dotplot
#############
DP <- DotPlot(object = subset(LescroartSeurat.corrected,idents=c("E6.75_WT","E6.75_KO")), cols = c("lightgrey", "blue"), group.by="Condition", features = genesInTerm(GOdata)["GO:0030178"] ) +RotatedAxis() # +coord_flip()  #, cols = paired_colors[1:2],  )
ggsave( plot = DP, height=60, width=800, filename = "Lescroart_KOvsWT_Dotplt_GO0030178.png", units = "mm",dpi = 300 )


#############
#GO:1900746 dotplot
#############
DP <- DotPlot(object = subset(LescroartSeurat.corrected,idents=c("E6.75_WT","E6.75_KO")), cols = c("lightgrey", "blue"), group.by="Condition", features = genesInTerm(GOdata)["GO:1900746"] ) +RotatedAxis() # +coord_flip()  #, cols = paired_colors[1:2],  )
ggsave( plot = DP, height=60, width=200, filename = "Lescroart_KOvsWT_Dotplt_GO1900746.png", units = "mm",dpi = 300 )

#############
#GO:1900746 dotplot
#############
DP <- DotPlot(object = subset(LescroartSeurat.corrected,idents=c("E6.75_WT","E6.75_KO")), cols = c("lightgrey", "blue"), group.by="Condition", features = genesInTerm(GOdata)["GO:0035025"] ) +RotatedAxis() # +coord_flip()  #, cols = paired_colors[1:2],  )
ggsave( plot = DP, height=60, width=200, filename = "Lescroart_KOvsWT_Dotplt_GO0035025.png", units = "mm",dpi = 300 )






