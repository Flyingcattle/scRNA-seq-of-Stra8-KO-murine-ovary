### load required packages 
library(Seurat)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

### read inputdata 
Germcell<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/Germcell.rds")
DimPlot(Germcell,label = T)

###############################GO enrichment analysis of DEGs from different gene sets.########################
######################################## FGC_mitotic vs FGC_mitotic_KO ########################################
DEGs_Germcell_FGC_mitotic <- FindMarkers(Germcell,
                                    ident.1 ="FGC_mitotic",
                                    ident.2 ='FGC_mitotic_KO',
                                    group.by = "celltype",
                                    min.pct = 0.25)
write.csv(DEGs_Germcell_FGC_mitotic,file = "DEGs_Germcell_FGC_mitotic.csv",row.names = TRUE, quote = F)
diff_FGC_mitotic <- read.table("DEGs_Germcell_FGC_mitotic.csv", sep = ",", header = T)

gene.df_diff_FGC_mitotic <- bitr(diff_FGC_mitotic$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb =org.Mm.eg.db)
gene_FGC_mitotic <- gene.df_diff_FGC_mitotic$ENTREZID
gene_ALL_FGC_mitotic <- enrichGO(gene = gene_FGC_mitotic,
                                 OrgDb = org.Mm.eg.db,
                                 keyType = 'ENTREZID',
                                 ont = 'ALL',
                                 pAdjustMethod = 'BH',
                                 minGSSize = 1,
                                 pvalueCutoff = 0.01,
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)

gene_ALL_1 <- as.data.frame(gene_ALL_FGC_mitotic)

ggplot(gene_ALL_FGC_mitotic, aes(y=Description,x=GeneRatio)) +
  geom_point(aes(size=Count,color=p.adjust))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
######################################## FGC_mitotic vs FGC_mitotic_KO ########################################

################## Oogonia_STRA8 vs Oogonia_STRA8_KO #####################################################
DEGs_Germcell_Oogonia_STRA8 <- FindMarkers(Germcell,
                                         ident.1 ="Oogonia_STRA8",
                                         ident.2 ='Oogonia_STRA8_KO',
                                         group.by = "celltype",
                                         min.pct = 0.25)
write.csv(DEGs_Germcell_Oogonia_STRA8,file = "DEGs_Germcell_Oogonia_STRA8.csv",row.names = TRUE, quote = F)
diff_Oogonia_STRA8 <- read.table("DEGs_Germcell_Oogonia_STRA8.csv", sep = ",", header = T)

gene.df_diff_Oogonia_STRA8 <- bitr(diff_Oogonia_STRA8$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb =org.Mm.eg.db)
gene_Oogonia_STRA8 <- gene.df_diff_Oogonia_STRA8$ENTREZID
gene_ALL_Oogonia_STRA8 <- enrichGO(gene = gene_Oogonia_STRA8,
                                   OrgDb = org.Mm.eg.db,
                                   keyType = 'ENTREZID',
                                   ont = 'ALL',
                                   pAdjustMethod = 'BH',
                                   minGSSize = 1,
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)

gene_ALL_1 <- as.data.frame(gene_ALL_Oogonia_STRA8)

ggplot(gene_ALL_Oogonia_STRA8, aes(y=Description,x=GeneRatio)) +
  geom_point(aes(size=Count,color=p.adjust))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

################## Oogonia_STRA8 vs Oogonia_STRA8_KO #####################################################

################## Pre_oocyte vs Pre_oocyte_KO ###########################################################
DEGs_Germcell_Pre_oocyte <- FindMarkers(Germcell,
                                        ident.1 ="Pre_oocyte",
                                        ident.2 ='Pre_oocyte_KO',
                                        group.by = "celltype",
                                        min.pct = 0.25)
write.csv(DEGs_Germcell_Pre_oocyte,file = "DEGs_Germcell_Pre_oocyte.csv",row.names = TRUE, quote = F)
diff_Pre_oocyte <- read.table("DEGs_Germcell_Pre_oocyte.csv", sep = ",", header = T)

gene.df_diff_Pre_oocyte <- bitr(diff_Pre_oocyte$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb =org.Mm.eg.db)
gene_Pre_oocyte <- gene.df_diff_Pre_oocyte$ENTREZID
gene_ALL_Pre_oocyte <- enrichGO(gene = gene_Pre_oocyte,
                                OrgDb = org.Mm.eg.db,
                                keyType = 'ENTREZID',
                                ont = 'ALL',
                                pAdjustMethod = 'BH',
                                minGSSize = 1,
                                pvalueCutoff = 0.01,
                                qvalueCutoff = 0.05,
                                readable = TRUE)

gene_ALL_1 <- as.data.frame(gene_ALL_Pre_oocyte)

ggplot(gene_ALL_Pre_oocyte, aes(y=Description,x=GeneRatio)) +
  geom_point(aes(size=Count,color=p.adjust))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
################## Pre_oocyte vs Pre_oocyte_KO ###########################################################