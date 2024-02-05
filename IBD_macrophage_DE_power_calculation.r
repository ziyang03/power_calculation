rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
library(Seurat)
library(Signac)
library(ggplot2)
library(magrittr)
library(sva)
library(RColorBrewer)
library(scales)
library(dplyr)
library(reshape2)
library(MAST)
library(edgeR)
library(RNASeqPower)
})
source('/data/miraldiNB/wayman/scripts/scRNA_utils.R')

####### INPUTS #######

# Output directory
dir_out <- '/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/gene_expression/DE_analysis/signature_gene/EdgeR'
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

# scRNA data
file_in_data <- '/data/GastroAI/genomics/analysis/multiome/IBD/outs/202401/annotation_final/All/multiome_Multi_All_harmony_annotated.rds'
data_assay <- 'RNA' # data assay
order_disease<-c('CD','UC')
order_tissue<-c('R','TI')
order_status<-c('I','NI')
order_disease_tissue<-c('CD_R','CD_TI','UC_R')
order_disease_tissue_status<-c('CD_R_I','CD_R_NI','CD_TI_I','CD_TI_NI','UC_R_I','UC_R_NI')
# Params
file_save <- 'macrophage' # file save base
meta_celltype <- 'CellTypeLevel1' # celltype metadata name
min_frac <- 0.05 # expressed in minimum fraction of cells per cluster
min_cells <- 50 # min cells to keep gene scrna
cutoff_log2fc <- 0.58 # log2FC cutoff for celltype signature genes
cutoff_fdr <- 0.1# FDR cutoff for celltype signature gene
type_de <- 'EdgeR' # type of DE test, options: 'DESeq2','wilcox', 'MAST', 'LR'
gene_y <- readLines('/data/GastroAI/genomics/analysis/multiome/liver/data/gene_filter/gene_set_refdata-cellranger-arc-GRCh38-2020-A-2.0.0_chrY.txt') # chrY genes (removing these from list)
######################
# Output directory
label_meta_celltype <- meta_celltype
if (min_frac > 0){
    label_meta_celltype <- paste(label_meta_celltype,paste0('pct',100*min_frac),sep='_')
}
dir_out <- file.path(dir_out, label_meta_celltype, type_de)
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

# Load scRNA
print('Load scRNA')
scrna <- readRDS(file_in_data)
DefaultAssay(scrna) <- data_assay
Idents(scrna) <- meta_celltype
# filter genes
print('Filter genes')
keep_gene <- rownames(CreateSeuratObject(scrna[["RNA"]]$counts, min.cells=min_cells))
keep_gene <- setdiff(keep_gene,gene_y)
scrna<-subset(scrna,features=keep_gene)
scrna[['ClusterSubtype']] <- paste(as.character(scrna@meta.data[,paste0(meta_celltype)]),as.character(scrna@meta.data[,'DiseaseTissueStatus']),
                                   as.character(scrna@meta.data[,'Sample']), sep='_')
Idents(scrna) <- meta_celltype

# DE genes - cell type signatures
print('DE: cell type signatures')
sig_gene <- NULL
Idents(scrna)<-'ClusterSubtype'
keep_idents_scrna <- names(which(table(Idents(scrna)) >= 5))
scrna <- subset(scrna, idents=keep_idents_scrna)
order_celltype<-unique(scrna$CellTypeLevel1)
order_sample<-unique(scrna$Sample)
file_counts_bulk<- '/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/gene_expression/DE_analysis/signature_gene/Status_within_celltype_disease_tissue/CellTypeLevel1_pct5/DESeq2/raw_counts_pesudobulk.tsv'
file_meta_bulk<-'/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/gene_expression/DE_analysis/signature_gene/Status_within_celltype_disease_tissue/CellTypeLevel1_pct5/DESeq2/meta_bulk.txt'
counts_bulk<-read.table(file_counts_bulk,sep='\t',header=T,row.names=1)
meta_bulk<-read.table(file_meta_bulk,sep='\t',header=T,row.names=1)
order_clustersubtype<-rownames(meta_bulk)
counts_bulk<-counts_bulk[,order_clustersubtype]
order_celltype<-c('CD4_T_Eff','CD4_T_cell','CD8_T_cell','gdT_cell','Naive_T_cell','B_cell','Plasma_cell',
'Macrophage','DC','NK_cell','Fibroblast','Enterocyte','Colonocyte','Goblet')
df_cv<-meta_bulk[,c('CellType','DiseaseTissueStatus')]
rownames(df_cv)<-NULL
df_cv<-unique(df_cv)
df_cv<-df_cv[df_cv$CellType%in%order_celltype,]
df_cv$disp<-NA
df_cv$BCV<-NA
for (ix in order_celltype){
    for (jx in order_disease_tissue_status){
    order_context_sample<-NULL
    for (kx in order_sample){
    order_context_sample<-c(order_context_sample,paste(ix,jx,kx,sep='_'))
    }    
    order_context_sample<-intersect(order_context_sample,rownames(meta_bulk))
    if(length(order_context_sample)>0){    
    counts_bulk_sub<-counts_bulk[,order_context_sample]
    counts_bulk_sub<-counts_bulk_sub[rowMeans(counts_bulk_sub)>=5,]
    meta_sub<-meta_bulk[order_context_sample,]                       
    obj<-DGEList(counts=counts_bulk_sub,group=meta_sub[,'DiseaseTissueStatus'])
    obj<-calcNormFactors(obj)
    obj <- estimateCommonDisp(obj,verbose=T)
    obj <- estimateTagwiseDisp(obj)
    file_out<-file.path(dir_out,paste0('CV_estimation_',ix,'_',jx,'.pdf'))
    pdf(file_out,height=6,width=6)
    plotBCV(obj)
    dev.off()
    df_cv[df_cv$CellType%in%ix&df_cv$DiseaseTissueStatus%in%jx,'disp']<-round(obj$common.dispersion,2)
    df_cv[df_cv$CellType%in%ix&df_cv$DiseaseTissueStatus%in%jx,'BCV']<-round(sqrt(obj$common.dispersion),2)
    }
  }
}
file_out<-file.path(dir_out,paste('Dispersion_BCV_estimation_filtered.txt'))
write.table(df_cv,file_out,sep='\t',col.names=T,row.names=F,quote=F)
#Power calculation for DE analysis
#CD_R
CD_R_FC_2<-rnapower(5, cv=0.36, cv2 = 0.35, effect=2, alpha=0.1, power=seq(0.2,0.95,by=0.05))
CD_R_FC_2<-ceiling(CD_R_FC_2)
CD_R_FC_2<-melt(CD_R_FC_2)
colnames(CD_R_FC_2)<-'N'
CD_R_FC_2$Power<-rownames(CD_R_FC_2)
CD_R_FC_2$Foldchange<-2
CD_R_FC_2<-rbind(c(1,0,2),CD_R_FC_2)
CD_R_FC_1.5<-rnapower(5, cv=0.36, cv2 = 0.35, effect=1.5, alpha=0.1, power=seq(0.2,0.95,by=0.05))
CD_R_FC_1.5<-ceiling(CD_R_FC_1.5)
CD_R_FC_1.5<-melt(CD_R_FC_1.5)
colnames(CD_R_FC_1.5)<-'N'
CD_R_FC_1.5$Power<-rownames(CD_R_FC_1.5)
CD_R_FC_1.5$Foldchange<-1.5
CD_R_FC_1.5<-rbind(c(1,0,1.5),CD_R_FC_1.5)
CD_R<-rbind(CD_R_FC_2,CD_R_FC_1.5)
CD_R<-as.data.frame(CD_R)
CD_R$N<-as.numeric(CD_R$N)
CD_R$Power<-as.numeric(CD_R$Power)
CD_R$Foldchange<-as.character(CD_R$Foldchange)
CD_R$Foldchange<-factor(CD_R$Foldchange,levels=c(2,1.5))
ggplot(CD_R, aes(x=N, y=Power, group=Foldchange)) +
  geom_line(aes(color=Foldchange))+
  scale_color_manual(values=c("purple", "green")) +
  labs(title = "Crohn's disease rectum samples macrophage RNA-seq",x="n",
       y = "Power") +
  theme_minimal() +
  geom_hline(yintercept = c(0.8), linetype = "dashed", color = "purple")+
  geom_vline(xintercept = c(9), linetype = "dashed", color = "purple")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold",hjust=0.5),
        panel.grid.major.x = element_line(linetype=2,color='grey',linewidth=0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype=2,color='grey',linewidth=0.3),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(2,2,2,2),"cm"),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,by=0.1)) +
  scale_x_continuous(expand = c(0, 0),breaks=c(1,5,9,15,20,25,30),limits=c(0,30))
file_out<-file.path(dir_out,paste('CD_R_power_calculation.pdf'))
ggsave(file_out,height=4.5,width=5.5)

#CD_TI
CD_TI_FC_2<-rnapower(5, cv=0.39, cv2 = 0.35, effect=2, alpha=0.1, power=seq(0.2,0.95,by=0.05))
CD_TI_FC_2<-ceiling(CD_TI_FC_2)
CD_TI_FC_2<-melt(CD_TI_FC_2)
colnames(CD_TI_FC_2)<-'N'
CD_TI_FC_2$Power<-rownames(CD_TI_FC_2)
CD_TI_FC_2$Foldchange<-2
CD_TI_FC_2<-rbind(c(1,0,2),CD_TI_FC_2)
CD_TI_FC_1.5<-rnapower(5, cv=0.39, cv2 = 0.35, effect=1.5, alpha=0.1, power=seq(0.2,0.95,by=0.05))
CD_TI_FC_1.5<-ceiling(CD_TI_FC_1.5)
CD_TI_FC_1.5<-melt(CD_TI_FC_1.5)
colnames(CD_TI_FC_1.5)<-'N'
CD_TI_FC_1.5$Power<-rownames(CD_TI_FC_1.5)
CD_TI_FC_1.5$Foldchange<-1.5
CD_TI_FC_1.5<-rbind(c(1,0,1.5),CD_TI_FC_1.5)
CD_TI<-rbind(CD_TI_FC_2,CD_TI_FC_1.5)
CD_TI<-as.data.frame(CD_TI)
CD_TI$N<-as.numeric(CD_TI$N)
CD_TI$Power<-as.numeric(CD_TI$Power)
CD_TI$Foldchange<-as.character(CD_TI$Foldchange)
CD_TI$Foldchange<-factor(CD_TI$Foldchange,levels=c(2,1.5))
ggplot(CD_TI, aes(x=N, y=Power, group=Foldchange)) +
  geom_line(aes(color=Foldchange))+
  scale_color_manual(values=c("purple", "green")) +
  labs(title = "Crohn's disease ileum samples macrophage RNA-seq",x="n",
       y = "Power") +
  theme_minimal() +
  geom_hline(yintercept = c(0.8), linetype = "dashed", color = "purple")+
  geom_vline(xintercept = c(9), linetype = "dashed", color = "purple")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold",hjust=0.5),
        panel.grid.major.x = element_line(linetype=2,color='grey',linewidth=0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype=2,color='grey',linewidth=0.3),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(2,2,2,2),"cm"),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,by=0.1)) +
  scale_x_continuous(expand = c(0, 0),breaks=c(1,5,9,15,20,25,30),limits=c(0,30))
file_out<-file.path(dir_out,paste('CD_TI_power_calculation.pdf'))
ggsave(file_out,height=4.5,width=5.5)

#UC_R
UC_R_FC_2<-rnapower(5, cv=0.4, cv2 = 0.28, effect=2, alpha=0.1, power=seq(0.2,0.95,by=0.05))
UC_R_FC_2<-ceiling(UC_R_FC_2)
UC_R_FC_2<-melt(UC_R_FC_2)
colnames(UC_R_FC_2)<-'N'
UC_R_FC_2$Power<-rownames(UC_R_FC_2)
UC_R_FC_2$Foldchange<-2
UC_R_FC_2<-rbind(c(1,0,2),UC_R_FC_2)
UC_R_FC_1.5<-rnapower(5, cv=0.4, cv2 = 0.28, effect=1.5, alpha=0.1, power=seq(0.2,0.95,by=0.05))
UC_R_FC_1.5<-ceiling(UC_R_FC_1.5)
UC_R_FC_1.5<-melt(UC_R_FC_1.5)
colnames(UC_R_FC_1.5)<-'N'
UC_R_FC_1.5$Power<-rownames(UC_R_FC_1.5)
UC_R_FC_1.5$Foldchange<-1.5
UC_R_FC_1.5<-rbind(c(1,0,1.5),UC_R_FC_1.5)
UC_R<-rbind(UC_R_FC_2,UC_R_FC_1.5)
UC_R<-as.data.frame(UC_R)
UC_R$N<-as.numeric(UC_R$N)
UC_R$Power<-as.numeric(UC_R$Power)
UC_R$Foldchange<-as.character(UC_R$Foldchange)
UC_R$Foldchange<-factor(UC_R$Foldchange,levels=c(2,1.5))
ggplot(UC_R, aes(x=N, y=Power, group=Foldchange)) +
  geom_line(aes(color=Foldchange))+
  scale_color_manual(values=c("purple", "green")) +
  labs(title = "Ulcerative colitis rectum samples macrophage RNA-seq",x="n",
       y = "Power") +
  theme_minimal() +
  geom_hline(yintercept = c(0.8), linetype = "dashed", color = "purple")+
  geom_vline(xintercept = c(9), linetype = "dashed", color = "purple")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold",hjust=0.5),
        panel.grid.major.x = element_line(linetype=2,color='grey',linewidth=0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype=2,color='grey',linewidth=0.3),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(2,2,2,2),"cm"),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
  scale_y_continuous(expand=c(0,0),breaks=seq(0,1,by=0.1)) +
  scale_x_continuous(expand = c(0, 0),breaks=c(1,5,9,15,20,25,30),limits=c(0,30))
file_out<-file.path(dir_out,paste('UC_R_power_calculation.pdf'))
ggsave(file_out,height=4.5,width=5.5)