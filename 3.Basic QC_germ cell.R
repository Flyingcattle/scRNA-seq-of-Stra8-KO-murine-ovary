### load required packages 
library(Seurat)
library(SeuratData)
library(dplyr)
library(cowplot)
library(DoubletFinder)
library(ggplot2)

### read inputdata 
seurat_obj<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/seurat_obj.rds")

### Extract the germ cells  
Germcell <- subset(seurat_obj, subset = celltype == 'Germ cell')

Germcell <- NormalizeData(Germcell, normalization.method = "LogNormalize", scale.factor = 10000)
Germcell <- FindVariableFeatures(Germcell, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Germcell)
Germcell <- ScaleData(Germcell, features = all.genes)
Germcell <- RunPCA(Germcell,features = VariableFeatures(object = Germcell))

ElbowPlot(Germcell)

Germcell <- FindNeighbors(Germcell, dims = 1:6)
Germcell <- FindClusters(Germcell, resolution = 1.1)
Germcell <- RunUMAP(Germcell,dims = 1:6)

my_color=c("#1F8A42", "#832688", "#F47D2C","#6E4C9F",
                     "#0C737B", "#282E6C", "#91D5E5","#89C75F",
                     "#E4989C", "#D44C26", "#D7A764","#DB2228",
                     "#87A0D1", "#C06DAD", "#3DBBA7","#2CA029","#FD7F0D")

table(Germcell$seurat_clusters)

DimPlot(Germcell,label = T)
DimPlot(Germcell,label = T,cols = my_color)
DimPlot(Germcell,label = T,group.by = "orig.ident",cols = my_color)
DimPlot(Germcell,label = T,group.by = "time",cols = my_color)

saveRDS(Germcell, file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/Germcell.rds")

Germcell<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/Germcell.rds")

FeaturePlot(Germcell,features = c("Nanog","Ccnb1","Aurka","Cenpa"))#FGC_mitotic
FeaturePlot(Germcell,features = c("Rec8","Sct","Ccnd1","Stra8"))#Oogonia_STRA8
FeaturePlot(Germcell,features = c("Prdm9","Inca1","Rad51ap2","Spdya"))#Oogonia_meiotic
FeaturePlot(Germcell,features = c("Syce3","Dmrtc2","Figla","Syce2"))#Pre_oocyte

FeaturePlot(Germcell,features = c("Nanog","Ccnb1","Aurka"))#FGC_mitotic_KO
FeaturePlot(Germcell,features = c("Rec8","Sct","Ccnd1"))#Oogonia_STRA8_KO
FeaturePlot(Germcell,features = c("Syce3","Lhx8","Figla"))#Pre_oocyte_KO

Germcell@meta.data$celltype <- as.numeric(as.character(Germcell@meta.data$seurat_clusters))

Germcell@meta.data <- within(Germcell@meta.data, celltype[celltype==1] <- c("FGC_mitotic"))
Germcell@meta.data <- within(Germcell@meta.data, celltype[celltype==6|celltype==3|celltype==5|celltype==4|celltype==11|celltype==12|celltype==14] <- c("Oogonia_STRA8"))
Germcell@meta.data <- within(Germcell@meta.data, celltype[celltype==2|celltype==9|celltype==15] <- c("Oogonia_meiotic"))
Germcell@meta.data <- within(Germcell@meta.data, celltype[celltype==8] <- c("Pre_oocyte"))
Germcell@meta.data <- within(Germcell@meta.data, celltype[celltype==7] <- c("FGC_mitotic_KO"))
Germcell@meta.data <- within(Germcell@meta.data, celltype[celltype==13|celltype==0] <- c("Oogonia_STRA8_KO"))
Germcell@meta.data <- within(Germcell@meta.data, celltype[celltype==10] <- c("Oocyte_like"))

saveRDS(Germcell, file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/Germcell.rds")





