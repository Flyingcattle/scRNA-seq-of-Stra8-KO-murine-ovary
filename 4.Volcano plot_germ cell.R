### load required packages 
library(Seurat)
library(ggplot2)
library(ggrepel)

### read inputdata 
Germcell<-readRDS(file = "/home/lab401_G/10xGenomics/scRNA_Stra8_KO/outcomes/Germcell.rds")
DimPlot(Germcell,label = T)

### Volcano plot demonstrating DEGs in different germ cell stages.
################## FGC_mitotic vs FGC_mitotic_KO #####################################################
DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO <- FindMarkers(Germcell,
                                                        ident.1 ="FGC_mitotic",
                                                        ident.2 ='FGC_mitotic_KO',
                                                        group.by = "celltype",
                                                        logfc.threshold = 0,
                                                        min.pct = 0)

DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$difference <- DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$pct.1 - DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$pct.2
DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO_sig <- DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO[which(DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$p_val_adj<0.05 & abs(DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$avg_log2FC) >0.25),]
DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO_sig$label <- rownames(DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO_sig)

ggplot(DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=0.5, color="grey60") + 
  geom_text_repel(data=DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO_sig, aes(label=label),
                  color="black",fontface="italic")+
  geom_point(data=DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO[which(DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$p_val_adj<0.05 & DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$avg_log2FC>0.1),],
             aes(x=difference, y=avg_log2FC),
             size=0.5, color="red")+
  geom_point(data=DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO[which(DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$p_val_adj<0.05 & DEGs_Germcell_FGC_mitotic_FGC_mitotic_KO$avg_log2FC< -0.1),],
             aes(x=difference, y=avg_log2FC),
             size=0.5, color="blue")+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
################## FGC_mitotic vs FGC_mitotic_KO #########################################################

################## Oogonia_STRA8 vs Oogonia_STRA8_KO #####################################################
DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO <- FindMarkers(Germcell,
                                                            ident.1 ="Oogonia_STRA8",
                                                            ident.2 ='Oogonia_STRA8_KO',
                                                            group.by = "celltype",
                                                            logfc.threshold = 0,
                                                            min.pct = 0)

DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$difference <- DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$pct.1 - DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$pct.2
DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO_sig <- DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO[which(DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$p_val_adj<0.05 & abs(DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$avg_log2FC) >0.25),]
DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO_sig$label <- rownames(DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO_sig)

ggplot(DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=0.5, color="grey60") + 
  geom_text_repel(data=DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO_sig, aes(label=label),
                  color="black",fontface="italic")+
  geom_point(data=DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO[which(DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$p_val_adj<0.05 & DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$avg_log2FC>0.1),],
             aes(x=difference, y=avg_log2FC),
             size=0.5, color="red")+
  geom_point(data=DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO[which(DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$p_val_adj<0.05 & DEGs_Germcell_Oogonia_STRA8_Oogonia_STRA8_KO$avg_log2FC< -0.1),],
             aes(x=difference, y=avg_log2FC),
             size=0.5, color="blue")+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
################## Oogonia_STRA8 vs Oogonia_STRA8_KO #####################################################

################## Pre_oocyte vs Pre_oocyte_KO ###########################################################
DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO <- FindMarkers(Germcell,
                                                      ident.1 ="Pre_oocyte",
                                                      ident.2 ='Pre_oocyte_KO',
                                                      group.by = "celltype",
                                                      logfc.threshold = 0,
                                                      min.pct = 0)

DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$difference <- DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$pct.1 - DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$pct.2
DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO_sig <- DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO[which(DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$p_val_adj<0.05 & abs(DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$avg_log2FC) >0.25),]
DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO_sig$label <- rownames(DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO_sig)

ggplot(DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=0.5, color="grey60") + 
  geom_text_repel(data=DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO_sig, aes(label=label),
                  color="black",fontface="italic")+
  geom_point(data=DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO[which(DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$p_val_adj<0.05 & DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$avg_log2FC>0.1),],
             aes(x=difference, y=avg_log2FC),
             size=0.5, color="red")+
  geom_point(data=DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO[which(DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$p_val_adj<0.05 & DEGs_Germcell_Pre_oocyte_Pre_oocyte_KO$avg_log2FC< -0.1),],
             aes(x=difference, y=avg_log2FC),
             size=0.5, color="blue")+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
################## Pre_oocyte vs Pre_oocyte_KO ###########################################################