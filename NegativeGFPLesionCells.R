library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)
library(Matrix)
library(umap)


endodata.data <- Read10X_h5('RM4_filtered_feature_bc_matrix.h5')
endodata.data
endodata <- CreateSeuratObject(counts = endodata.data, project = "endodata3k", min.genes = 0, min.cells = 0, min.features = 200)

eGFP_expression = GetAssayData(object = endodata, assay = "RNA", slot = "data")["eGFP",]
neg_ids = names(which(eGFP_expression==0))
neg_cells = subset(endodata,cells=neg_ids)

neg_cells[["percent.mt"]] <- PercentageFeatureSet(neg_cells, pattern = "mt-")
VlnPlot(neg_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

neg_cells <- subset(neg_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

neg_cells <- NormalizeData(neg_cells, normalization.method = "LogNormalize", scale.factor = 10000)
neg_cells <- NormalizeData(neg_cells)

neg_cells <- FindVariableFeatures(neg_cells, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(neg_cells)
neg_cells <- ScaleData(neg_cells, features = all.genes)

neg_cells <- RunPCA(neg_cells, features = VariableFeatures(object = neg_cells))

DimPlot(neg_cells, reduction = "pca")
FeaturePlot(neg_cells, features = c("eGFP"), cols = c("grey","grey"))

#Angptl1_expression1 = GetAssayData(object = neg_cells, assay = "RNA", slot = "data")["Angptl1",]
#pos_ids1 = names(which(Angptl1_expression1>0))
#neg_ids1 = names(which(Angptl1_expression1==0))
#pos_exp1 = subset(neg_cells,cells=pos_ids1)
#neg_exp1 = subset(neg_cells,cells=neg_ids1)

#Cells(neg_cells)
#Cells(pos_exp1)

DimHeatmap(neg_cells, dims = 1, cells = 500, balanced = TRUE)

neg_cells <- FindNeighbors(neg_cells, dims = 1:10)
neg_cells <- FindClusters(neg_cells, resolution = c(0.75))

neg_cells <- RunUMAP(neg_cells, dims = 1:10)

DimPlot(neg_cells, reduction = "umap", group.by = 'seurat_clusters')
DimPlot(neg_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


FeaturePlot(neg_cells, features = c("eGFP"), cols = c("grey","grey"))

FeaturePlot(neg_cells, features = c("Angptl1"), cols = c("grey","blue"))

VlnPlot(neg_cells, features = "Hoxa10")


Gene_Expression1 <- FindMarkers(neg_cells, ident.1 = 11, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
  
FindMarkers(neg_cells, ident.1 = "0", logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
  
write.table(Gene_Expression1, file = "my_data.csv", append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)

new.cluster.ids <- c('Natural Killer Cells','B Cells', 'Macrophages', 'B Cells',
                     'T Cells', 'T Cells', 'Fibroblasts', 'Neutrophils', 
                     'Macrophages', 'Fibroblasts', 'B Cells', 'Macrophages',
                     'Neutrophils', 'Endothelial Cells', 'Neutrophils')

names(new.cluster.ids) <- levels(neg_cells)
new.cluster.ids
neg_cells <- RenameIdents(neg_cells, new.cluster.ids)

DimPlot(neg_cells, reduction = "umap", label = TRUE, pt.size = 0.25, label.size = 7.5) + theme_classic(base_size = 20)

F_cluster_neg <- subset(neg_cells, idents = "Fibroblasts")


cell.num2 <- table(Idents(neg_cells))
cell.num2

g1_untreat <- WhichCells(neg_cells, idents = c("B Cells",
              "Endothelial Cells","Natural Killer Cells",
              "T Cells", "B Cells","B Cells","Neutrophils","Macrophages"))
g1_treat <- WhichCells(neg_cells, idents = c( "Fibroblasts"))
DimPlot(neg_cells, reduction = "umap", label = TRUE, label.size = 4.0, 
        pt.size = 0.5, cells.highlight= list(g1_treat, g1_untreat), 
        cols.highlight = "turquoise 3", cols = "grey") + NoLegend()

FeaturePlot(F_cluster_neg, features = c("Angptl1"), cols = c("grey","blue"))

#negcells.markers <- FindAllMarkers(neg_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#negcells.markers %>% group_by(cluster) %>% top_n(n = 2)
#top10_neg <- negcells.markers %>% group_by(cluster) %>% top_n(n = 10)
#top10_neg

FeaturePlot(neg_cells, features = c("Ctla2a"), cols = c("grey","blue"))
VlnPlot(neg_cells, features = c("Igf1"))

AverageExpression(neg_cells, features = "Ctla2a")

#Angptl1_expression3 = GetAssayData(object = F_cluster, assay = "RNA", slot = "data")["Angptl1",]
#pos_ids3 = names(which(Angptl1_expression3>0))
#neg_ids3 = names(which(Angptl1_expression3==0))
#pos_exp3 = subset(F_cluster,cells=pos_ids3)
#neg_exp3 = subset(F_cluster,cells=neg_ids3)

#Cells(pos_exp3)

VlnPlot(neg_cells, features = c("Cd14", "Cd68","Cd86","Cd80"))


FindMarkers(neg_cells, ident.1 = "Endothelial Cells", logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)

VlnPlot(neg_cells,  ncol = 2, features = c("Egfl7","Mmrn1","Cldn5","Tie1","Ccl21a","Cdh5",
                                "Flt4","Sdpr","Ptprb","Reln"))

VlnPlot(neg_cells, features = c("Sf1"))

VlnPlot(neg_cells, ncol = 2, features = c("Il1b","Kras","Igf1","Mmp2"))

VlnPlot(neg_cells, ncol = 2, features = c("Ccl8", "Ccl7", "Ccl12", "Ccl2", "Ccl5",
                                "Ccl4", "Cxcl12", "Cxcl16","Ccl21a"))

micro_genes <- c("Ccl8", "Ccl7", "Ccl12")

VlnPlot(neg_cells, features = micro_genes[1:6]) & NoLegend()

ggplot(neg_cells, aes(x = "Ccl8", "Ccl7", "Ccl12")) + geom_violin()

Gene_Expression_B_Cells_neg <- FindMarkers(neg_cells, ident.1 = "B Cells", logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)

Endothelial_Cells <- subset(neg_cells, idents = "Endothelial Cells")


AverageExpression(neg_cells, features = "Hoxa10")
mean_vector <- c(0.007663928, 0.005740001, 0.07586678, 4.064982,
                 0.02164355, 0.01686745)
mean(mean_vector)

write.table(Gene_Expression_B_Cells_neg, file = "B_cells_Non-GFP_markers.csv", append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

#vector_B_Cells_Expression_NonGFP <- list()
#for(z in 1:423)
#{
#  val3 = as.character(vector_matching_genes[z])
#  expression_pos_cells <- AverageExpression(BC_cluster_neg, features = val3)
#  vector_B_Cells_Expression_NonGFP = append(vector_B_Cells_Expression_NonGFP,expression_pos_cells)
#}

#vector_B_Cells_Expression_NonGFP <- as.data.frame(vector_B_Cells_Expression_NonGFP)
#write.csv(vector_B_Cells_Expression_NonGFP, file = 'B_cells_NonGFP_expression.csv', append = FALSE, sep = "\t", dec = ".",
#          row.names = TRUE, col.names = TRUE)


M_cluster_neg <- FindNeighbors(M_cluster_neg, dims = 1:10)
M_cluster_neg <- FindClusters(M_cluster_neg, resolution = c(0.5))

M_cluster_neg <- RunUMAP(M_cluster_neg, dims = 1:10)

DimPlot(M_cluster_neg, reduction = "umap", group.by = 'seurat_clusters')
DimPlot(M_cluster_neg, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(M_cluster_neg, features = c("Ccrl2"), cols = c("grey","blue"))

FindMarkers(M_cluster_neg, ident.1 = 3, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)

new.cluster.ids.macrophages <- c('0', '1', 'M2 Macrophages', 
                                 'M1 Macrophages', '4')

names(new.cluster.ids.macrophages) <- levels(M_cluster_neg)
new.cluster.ids.macrophages
M_cluster_neg <- RenameIdents(M_cluster_neg, new.cluster.ids.macrophages)

M2_Macrophages_cellcount <- subset(M_cluster_neg, idents = "M2 Macrophages")

DimPlot(M_cluster_neg, reduction = "umap", label = TRUE, pt.size = 0.5)
