# Vizgen MERFISH liver data ======================

setwd("D:/ETH/Merscope/liver_cancer_merscope")

library(Seurat)
library(ggplot2)
library(reshape2)
library(gplots)
library(vroom)
library(scCustomize)

#read in data
countTable <- vroom("Liver1Slice1_cell_by_gene.csv")

input<-countTable
row.names(input)<-input$cell
input$cell<-NULL

analysis1 <- CreateSeuratObject(counts =input, 
                                project = "merscope", 
                                min.cells = 100, 
                                min.features = 8, 
                                meta.data = input$seurat_clusters)

#defining the graph.name parameter
Graphs(analysis)

#QC
analysis <- subset(analysis1, subset = nCount_RNA > 25 & nCount_RNA < 1000) 

#Normalizing
analysis <- NormalizeData(analysis1, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 1)

#discovering the most variable features (genes)
analysis <- FindVariableFeatures(analysis, 
                                 selection.method = "vst", 
                                 nfeatures = 100)
plot1 <- VariableFeaturePlot(analysis)
plot1

#identifying the 10 most highly variable genes
top10 <- head(VariableFeatures(analysis), 10)
top10

#converting normalized gene expression to Z-score 
#(values centered at 0 and with variance of 1)
all.genes <- rownames(analysis)
analysis <- ScaleData(analysis, 
                      features = all.genes)

#perfom PCA
#common way of linear dimensionality reduction
analysis <- RunPCA(analysis, 
                   features = VariableFeatures(object = analysis), 
                   approx=FALSE)

#Prinicpal component “loadings” 
#should match markers of distinct populations 
#for well-behaved datasets
VizDimLoadings (analysis, 
                dims = 1:9, 
                reduction = 'pca') & 
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=8,
                                face="bold"))

# visualizing all reduced representations
DimPlot_scCustom(analysis, 
                 reduction = 'pca', 
                 label = T, 
                 figure_plot = TRUE,
                 repel = TRUE)

#clustering the cells
#Computing nearest neighbor graph
analysis <- FindNeighbors(analysis, 
                          dims = dimensions, 
                          reduction = "pca",
                          graph.name = "RNA_snn") 

#Computing SNN (shared nearest neighbor)
analysis <- FindClusters(analysis, 
                         resolution = 2, 
                         random.seed = 2, 
                         algorithm = 1, 
                         graph.name = "RNA_snn",
                         verbose = TRUE)

#generating UMAP reduced dimensionality representation
analysis <- RunUMAP(analysis, 
                    dims = 1:10,
                    spread = 2, 
                    seed.use = 5)

table1 <- table(analysis@meta.data$seurat_clusters)
table1

#UMAP by default, with Seurat clusters as identity
DimPlot_scCustom(seurat_object = analysis, 
                 reduction = "umap", 
                 repel = T, 
                 label = 6)

#Looking for the marker genes in every cluster and grouping them
analysis.markers <- FindAllMarkers(analysis, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.01) %>%
  Add_Pct_Diff()


analysis.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, 
            order_by = avg_log2FC)%>% 
  print(n=60)

#defining clusters 
#0       Fgb              Hepatocytes PV
# 1       Cyp2e1       Hepatocytes CV 
# 2       Fgb              Hepatocytes PV
# 3       Cyp2e1       Hepatocytes CV 
# 4       App              Cholangiocytes (Met)
# 5       Tgfbi             Monocytes 
# 6       Cyp2e1        Hepatocytes CV (PV)
# 7       Pglyrp1        Metastasis (Fib)
# 8       Cyp2e1       Hepatocytes CV 
# 9        Pglyrp1        Metastasis (Chol)
# 10      Fgb              Hepatocytes PV
# 11     Tgfbi             Monocytes (Stellate)
# 12     Spp1             Cholangiocytes 
# 13     Cd5l               Kupffer
# 14     Fgb                Hepatocytes PV (B)
# 15     Fgb                Hepatocytes PV
# 16     Fgb                Hepatocytes PV
# 17     Cyp2e1       Hepatocytes CV
# 18     Cyp2e1       Hepatocytes CV
# 19     App              Cholangiocytes (Fib)
# 20     Cyp2e1       Hepatocytes CV
# 21     Fn1              Fibroblasts 
# 22     Spp1             Cholangiocytes 
# 23     Pglyrp1        Metastasis (Chol)
# 24     Il2rb              T   (Fibr)
# 25     Cald1            Fibroblasts 


#cell type annotation
library(readxl)
annotation <-  read_excel('cluster_Annotation1.xlsx')
annotation_df <- data.frame(annotation)

annotation_info <- Pull_Cluster_Annotation(annotation = annotation_df)

#assigning cluster labels from marker
analysis_renamed <- Rename_Clusters(seurat_object = analysis, 
                                    new_idents = annotation_info$new_cluster_idents)
#plot with clusters names
DimPlot_scCustom(seurat_object = analysis_renamed, 
                 reduction = "umap", 
                 label = TRUE,
                 figure_plot = TRUE,
                 pt.size = 1) + NoLegend()