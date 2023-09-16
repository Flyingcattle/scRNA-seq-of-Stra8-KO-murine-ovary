### load required packages 
library(Seurat)
library(SeuratData)
library(dplyr)
library(cowplot)
library(DoubletFinder)
library(ggplot2)

### read inputdata 
seurat_obj<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/seurat_obj.rds")

###Dot plot of the gonadal cells
marker <- c('Sycp3','Ddx4','Dazl',
                     'Wnt4','Wnt6','Fst',
                     'Upk3b','Tm4sf5','Krt19',
                     'Bgn','Col1a2','Ptn',
                     'Pecam1','Cd34','Cldn5',
                     'Alas2','Slc25a37','Gypa',
                     'Cd52','C1qb','Lyz2')

DotPlot(seurat_obj, features = marker,group.by='seurat_clusters')+coord_flip()+theme_bw()+
               theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
               labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
               scale_color_gradientn(values = seq(0,1,0.2),colours = c('#DADAEB','#9E9AC8','#6A51A3','#54278F'))