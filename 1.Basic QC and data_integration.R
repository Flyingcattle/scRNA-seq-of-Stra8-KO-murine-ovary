### load required packages 
library(Seurat)
library(dplyr)
library(cowplot)
library(DoubletFinder)
library(ggplot2)

### read inputdata 
Stra8_KO_K14.data <- Read10X('/home/lab401_G/10xGenomics/scRNA_Stra8_KO/K14/outs/filtered_feature_bc_matrix')
Stra8_KO_E16.5.data <- Read10X('/home/lab401_G/10xGenomics/scRNA_Stra8_KO/K16/outs/filtered_feature_bc_matrix/')
E14.data <- Read10X('/home/lab401_G/10xGenomics/scRNA_Stra8_KO/GR14/filtered_feature_bc_matrix')
E16.data <- Read10X('/home/lab401_G/10xGenomics/scRNA_Stra8_KO/GR16/filtered_feature_bc_matrix')


## Construct Seurat object
Stra8_KO_E14.5 <- CreateSeuratObject(counts = Stra8_KO_K14.data, project = "Stra8_KO_E14.5", min.cells = 3, min.features = 200)
Stra8_KO_E14.5[["percent.mt"]] <- PercentageFeatureSet(Stra8_KO_E14.5, pattern = "^mt-")
VlnPlot(Stra8_KO_E14.5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
Stra8_KO_E14.5 <- subset(Stra8_KO_E14.5, subset = nFeature_RNA >2500 & nFeature_RNA < 5000 & percent.mt < 15)
Stra8_KO_E14.5 <- NormalizeData(Stra8_KO_E14.5, normalization.method = "LogNormalize", scale.factor = 10000)
Stra8_KO_E14.5 <- FindVariableFeatures(Stra8_KO_E14.5, selection.method = "vst", nfeatures = 2000)
Stra8_KO_E14.5$group<-"Stra8_KO_E14.5"
Stra8_KO_E14.5@meta.data$time <- "E14.5"
all.genes <- rownames(Stra8_KO_E14.5)
Stra8_KO_E14.5 <- ScaleData(Stra8_KO_E14.5, features = all.genes)  


Stra8_KO_E16.5 <- CreateSeuratObject(counts = Stra8_KO_E16.5.data, project = "Stra8_KO_E16.5", min.cells = 3, min.features = 200)
Stra8_KO_E16.5[["percent.mt"]]<-PercentageFeatureSet(Stra8_KO_E16.5,pattern = "^mt-")
VlnPlot(Stra8_KO_E16.5, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
Stra8_KO_E16.5 <- subset(Stra8_KO_E16.5, subset = nFeature_RNA > 2500 & nFeature_RNA < 5000 & percent.mt < 10)
Stra8_KO_E16.5 <- NormalizeData(Stra8_KO_E16.5, normalization.method = "LogNormalize", scale.factor = 10000)
Stra8_KO_E16.5 <- FindVariableFeatures(Stra8_KO_E16.5, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
Stra8_KO_E16.5$group<-"Stra8_KO_E16.5"
Stra8_KO_E16.5@meta.data$time <- "E16.5"
all.genes <- rownames(Stra8_KO_E16.5)
Stra8_KO_E16.5 <- ScaleData(Stra8_KO_E16.5, features = all.genes)

E14 <- CreateSeuratObject(counts = E14.data, project = "E14.5", min.cells = 3, min.features = 200)
E14[["percent.mt"]]<-PercentageFeatureSet(E14,pattern = "^mt-")
VlnPlot(E14, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
E14 <- subset(E14, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & percent.mt < 7.5)
E14 <- NormalizeData(E14, normalization.method = "LogNormalize", scale.factor = 10000)
E14 <- FindVariableFeatures(E14, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
E14$group<-"E14.5"
E14@meta.data$time <- "E14.5"
all.genes <- rownames(E14)
E14 <- ScaleData(E14, features = all.genes)

E16 <- CreateSeuratObject(counts = E16.data, project = "E16.5", min.cells = 3, min.features = 200)
E16[["percent.mt"]]<-PercentageFeatureSet(E16,pattern = "^mt-")
VlnPlot(E16, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
E16 <- subset(E16, subset = nFeature_RNA > 1800  & nFeature_RNA < 3500 & percent.mt < 10)# > 1000 < 5000   > 1800 < 3500
E16 <- NormalizeData(E16, normalization.method = "LogNormalize", scale.factor = 10000)
E16 <- FindVariableFeatures(E16, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
E16$group<-"E16.5"
E16@meta.data$time <- "E16.5"
all.genes <- rownames(E16)
E16 <- ScaleData(E16, features = all.genes)

#################################### remove douplets cells with Doublefinder #################################


############################################### Module 1 ##################################################
Stra8_KO_E14.5<-RunPCA(Stra8_KO_E14.5)
Stra8_KO_E14.5<-RunUMAP(Stra8_KO_E14.5,dims = 1:10)
## pK Identification (no ground-truth)
sweep.res.list_Stra8_KO_E14.5 <- paramSweep_v3(Stra8_KO_E14.5, PCs = 1:20, sct = FALSE)
head(sweep.res.list_Stra8_KO_E14.5)
sweep.stats_Stra8_KO_E14.5 <- summarizeSweep(sweep.res.list_Stra8_KO_E14.5, GT = FALSE)
bcmvn_CON11 <- find.pK(sweep.stats_Stra8_KO_E14.5)
## Homotypic Doublet Proportion Estimate
Stra8_KO_E14.5<-FindNeighbors(Stra8_KO_E14.5,reduction ="pca",dims = 1:20)
Stra8_KO_E14.5<-FindClusters(Stra8_KO_E14.5,resolution = 0.5)
head(Stra8_KO_E14.5@meta.data)
annotationscon_Stra8_KO_E14.5<-Stra8_KO_E14.5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_Stra8_KO_E14.5)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.075*(length(Stra8_KO_E14.5@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_Stra8_KO_E14.5 <- doubletFinder_v3(Stra8_KO_E14.5, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_Stra8_KO_E14.5@meta.data)

table(seu_Stra8_KO_E14.5$ DF.classifications_0.25_0.09_359)
# Doublet Singlet 
#  359    4428

seu_Stra8_KO_E14.5@meta.data$cellfilter <- seu_Stra8_KO_E14.5@meta.data$DF.classifications_0.25_0.09_359
seu_Stra8_KO_E14.5@meta.data <-seu_Stra8_KO_E14.5@meta.data[,-10]

seu_Stra8_KO_E14.5@meta.data$time<- "E14.5"

DimPlot(seu_Stra8_KO_E14.5, group.by = "cellfilter",
        cols = c("#92D713","#832688")) + ggtitle("Stra8_KO_E14.5")

############################################### Module 1 ##################################################

############################################### Module 2 ##################################################
Stra8_KO_E16.5<-RunPCA(Stra8_KO_E16.5)
Stra8_KO_E16.5<-RunUMAP(Stra8_KO_E16.5,dims = 1:10)
## pK Identification (no ground-truth)
sweep.res.list_Stra8_KO_E16.5 <- paramSweep_v3(Stra8_KO_E16.5, PCs = 1:20, sct = FALSE)
head(sweep.res.list_Stra8_KO_E16.5)
sweep.stats_Stra8_KO_E16.5 <- summarizeSweep(sweep.res.list_Stra8_KO_E16.5, GT = FALSE)
bcmvn_CON11  <- find.pK(sweep.stats_Stra8_KO_E16.5)
## Homotypic Doublet Proportion Estimate
Stra8_KO_E16.5<-FindNeighbors(Stra8_KO_E16.5,reduction ="pca",dims = 1:20)
Stra8_KO_E16.5<-FindClusters(Stra8_KO_E16.5,resolution = 0.5)
head(Stra8_KO_E16.5@meta.data)
annotationscon_Stra8_KO_E16.5<-Stra8_KO_E16.5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_Stra8_KO_E16.5)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.075*(length(Stra8_KO_E16.5@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_Stra8_KO_E16.5 <- doubletFinder_v3(Stra8_KO_E16.5, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_Stra8_KO_E16.5@meta.data)

table(seu_Stra8_KO_E16.5$DF.classifications_0.25_0.09_271)
#Doublet Singlet
# 271    3343

seu_Stra8_KO_E16.5@meta.data$cellfilter <- seu_Stra8_KO_E16.5@meta.data$DF.classifications_0.25_0.09_271
seu_Stra8_KO_E16.5@meta.data <-seu_Stra8_KO_E16.5@meta.data[,-10]

seu_Stra8_KO_E16.5@meta.data$time<- "E16.5"

DimPlot(seu_Stra8_KO_E16.5, group.by = "cellfilter",
        cols = c("#92D713","#832688"))  + ggtitle("Stra8_KO_E16.5")

############################################### Module 2 ##################################################

############################################### Module 3 ##################################################
E14<-RunPCA(E14)
E14<-RunUMAP(E14,dims = 1:10)
## pK Identification (no ground-truth)
sweep.res.list_E14 <- paramSweep_v3(E14, PCs = 1:20, sct = FALSE)
head(sweep.res.list_E14)
sweep.stats_E14 <- summarizeSweep(sweep.res.list_E14, GT = FALSE)
bcmvn_CON11 <- find.pK(sweep.stats_E14)
## Homotypic Doublet Proportion Estimate
E14<-FindNeighbors(E14,reduction ="pca",dims = 1:20)
E14<-FindClusters(E14,resolution = 0.5)
head(E14@meta.data)
annotationscon_E14<-E14@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E14)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.09*(length(E14@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_E14 <- doubletFinder_v3(E14, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_E14@meta.data)

table(seu_E14$DF.classifications_0.25_0.09_563)
#Doublet Singlet
#563   5688

seu_E14@meta.data$cellfilter <- seu_E14@meta.data$DF.classifications_0.25_0.09_563
seu_E14@meta.data <-seu_E14@meta.data[,-10]

seu_E14@meta.data$time<- "E14.5"

DimPlot(seu_E14, group.by = "cellfilter",
        cols = c("#92D713","#832688"))  + ggtitle("E14.5")

############################################### Module 3 ##################################################

############################################### Module 4 ##################################################
E16<-RunPCA(E16)
E16<-RunUMAP(E16,dims = 1:10)
## pK Identification (no ground-truth)
sweep.res.list_E16 <- paramSweep_v3(E16, PCs = 1:20, sct = FALSE)
head(sweep.res.list_E16)
sweep.stats_E16 <- summarizeSweep(sweep.res.list_E16, GT = FALSE)
bcmvn_CON11 <- find.pK(sweep.stats_E16)
## Homotypic Doublet Proportion Estimate
E16<-FindNeighbors(E16,reduction ="pca",dims = 1:20)
E16<-FindClusters(E16,resolution = 0.5)
head(E16@meta.data)
annotationscon_E16<-E16@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E16)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.09*(length(E16@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_E16 <- doubletFinder_v3(E16, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_E16@meta.data)

table(seu_E16$ DF.classifications_0.25_0.09_431)
#Doublet Singlet
# 431    4356

seu_E16@meta.data$cellfilter <- seu_E16@meta.data$ DF.classifications_0.25_0.09_431
seu_E16@meta.data <-seu_E16@meta.data[,-10]

seu_E16@meta.data$time<- "E16.5"

DimPlot(seu_E16, group.by = "cellfilter",
        cols = c("#92D713","#832688"))  + ggtitle("E16.5")

############################################### Module 4 ##################################################


#########################################  Extract the singlets  ############################################
Stra8_KO_E14.5.singlet <- subset(seu_Stra8_KO_E14.5, subset = cellfilter == 'Singlet')
Stra8_KO_E16.5.singlet <- subset(seu_Stra8_KO_E16.5, subset = cellfilter == 'Singlet')
E14.singlet <- subset(seu_E14, subset = cellfilter == 'Singlet')
E16.singlet <- subset(seu_E16, subset = cellfilter == 'Singlet')

## Here, add some UMAP plots to show the doublets dstribution is desired 
saveRDS(Stra8_KO_E14.5.singlet, file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/K14/Stra8_KO_E14.5.singlet.rds")
saveRDS(Stra8_KO_E16.5.singlet, file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/K16/Stra8_KO_E16.5.singlet.rds")
saveRDS(E14.singlet, file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/GR14/E14.singlet.rds")
saveRDS(E16.singlet, file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/GR16/E16.singlet.rds")

########################################### end doublefinder  ################################################

###########################################  harmony ##########################################
Stra8_KO_E14.5.singlet<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/K14/Stra8_KO_E14.5.singlet.rds")
Stra8_KO_E16.5.singlet<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/K16/Stra8_KO_E16.5.singlet.rds")
E14.5.singlet<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/GR14/E14.singlet.rds")
E16.5.singlet<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/GR16/E16.singlet.rds")

library(harmony)

seurat_list <- list(E14.5.singlet,new_E16.singlet,Stra8_KO_E14.5.singlet,Stra8_KO_E16.5.singlet)

# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
head(seurat_obj@meta.data)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# evaluating batch effect 
options(repr.plot.height = 5, repr.plot.width = 18)
p1 <- DimPlot(object = seurat_obj, reduction = "pca", pt.size = .1, group.by = "time")
p2 <- VlnPlot(object = seurat_obj, features = "PC_1", group.by = "time", pt.size = .1)
plot_grid(p1,p2)

seurat_obj=RunHarmony(seurat_obj,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(seurat_obj, 'harmony')
harmony_embeddings[1:5, 1:5]


options(repr.plot.height = 5, repr.plot.width = 18)
p1 <- DimPlot(object = seurat_obj, reduction = "harmony", pt.size = .1, group.by = "time")
p2 <- VlnPlot(object = seurat_obj, features = "harmony_1", group.by = "time", pt.size = .1)
plot_grid(p1,p2)


seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:20)    ## 23
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20,reduction = "harmony") #20
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

options(repr.plot.width=24,repr.plot.height=4)

my_color=("#1F8A42", "#832688", "#F47D2C","#6E4C9F",
          "#0C737B", "#282E6C", "#91D5E5","#89C75F",
          "#E4989C", "#D44C26", "#D7A764","#DB2228",
          "#87A0D1", "#C06DAD", "#3DBBA7","#2CA029","#FD7F0D")

DimPlot(seurat_obj,reduction = "umap",label = T,cols =my_color)

head(seurat_obj@meta.data)
DimPlot(seurat_obj,reduction = "umap",label = T,cols = my_color)
DimPlot(seurat_obj,reduction = "umap",label = F,cols = my_color,group.by = "seurat_clusters")
DimPlot(seurat_obj,reduction = "umap", group.by = "group")

seurat_obj@meta.data$celltype <- as.numeric(as.character(seurat_obj@meta.data$seurat_clusters))

seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==5|celltype==2|celltype==12|celltype==13|celltype==14|celltype==18] <- c("Germ cell"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==3|celltype==9|celltype==1|celltype==4] <- c("Pregranulosa cell"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==11|celltype==15] <- c("Mesothelial cell"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==7|celltype==0|celltype==6|celltype==17] <- c("Interstitial cell"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==10] <- c("Endothelial cell"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==8] <- c("Erythroid cell"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==16] <- c("Immune cell"))

DimPlot(seurat_obj,reduction = "umap",label = T,cols =my_color,group.by = "celltype")

FeaturePlot(seurat_obj,features = c('Sycp3','Ddx4','Dazl'),cols = c("lightgrey" ,"#DE1F1F"),slot = "data") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F))#Germ cell
FeaturePlot(seurat_obj,features = c('Wnt4','Wnt6','Fst','Kctd14'),cols = c("lightgrey" ,"#DE1F1F"),slot = "data") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F))#Pregranulosa cell
FeaturePlot(seurat_obj,features = c('Upk3b','Tm4sf5','Krt19','Lhx9'),cols = c("lightgrey" ,"#DE1F1F"),slot = "data") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F))#Mesothelial cell
FeaturePlot(seurat_obj,features = c('Bgn','Col1a2','Ptn','Pdgfra'),cols = c("lightgrey" ,"#DE1F1F"),slot = "data") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F))#Interstitial cell
FeaturePlot(seurat_obj,features = c('Pecam1','Kdr','Cldn5','Cd34'),cols = c("lightgrey" ,"#DE1F1F"),slot = "data") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F))#Endothelial cell
FeaturePlot(seurat_obj,features = c('Alas2','Slc25a37','Gypa','Alad'),cols = c("lightgrey" ,"#DE1F1F"),slot = "data") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F))#Erythroid cell
FeaturePlot(seurat_obj,features = c('Cd52','C1qb','Lyz2','Car2'),cols = c("lightgrey" ,"#DE1F1F"),slot = "data") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F))#Immune cell

saveRDS(seurat_obj, file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/seurat_obj.rds")


