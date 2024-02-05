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
library(DESeq2)
})
source('/data/miraldiNB/wayman/scripts/scRNA_utils.R')

####### INPUTS #######

# Output directory
dir_out <- '/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/gene_expression/DE_analysis/viz'
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

# raw count matrix and meta data
file_counts_bulk<- '/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/gene_expression/DE_analysis/signature_gene/Status_within_celltype_disease_tissue/CellTypeLevel1_pct5/DESeq2/raw_counts_pesudobulk.tsv'
file_meta_bulk<-'/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/gene_expression/DE_analysis/signature_gene/Status_within_celltype_disease_tissue/CellTypeLevel1_pct5/DESeq2/meta_bulk.txt'
# Params
data_assay <- 'RNA' # data assay
order_celltype<-c('Macrophage')
order_disease<-c('CD','UC')
order_disease_tissue<-c('CD_R','CD_TI','UC_R')
file_save <- 'macrophage' # file save base
#load data
counts_bulk<-read.table(file_counts_bulk,sep='\t',header=T,row.names=1)
meta_bulk<-read.table(file_meta_bulk,sep='\t',header=T,row.names=1)
order_sample<-unique(meta_bulk$Sample)
order_celltype_disease_tissue_status_sample<-NULL
for (ix in order_celltype){
    for (jx in order_disease_tissue_status){
        for (kx in order_sample){
            order_celltype_disease_tissue_status_sample<-c(order_celltype_disease_tissue_status_sample,paste(ix,jx,kx,sep='_'))
        }
    }
}
order_celltype_disease_tissue_status_sample<-intersect(order_celltype_disease_tissue_status_sample,rownames(meta_bulk))

counts_bulk<-counts_bulk[,order_celltype_disease_tissue_status_sample]
counts_bulk<-counts_bulk[rowMeans(counts_bulk)>=5,]
meta_bulk<-meta_bulk[order_celltype_disease_tissue_status_sample,]
#Calculate mean counts distribution for each disease tissue status
counts_mean_all<-NULL
for (ix in order_disease_tissue){
    counts_bulk_sub<-counts_bulk[,grepl(ix,order_celltype_disease_tissue_status_sample)]
    counts_mean<-as.data.frame(rowMeans(counts_bulk_sub))
    colnames(counts_mean)<-ix
    if (ix=='CD_R_I'){
        counts_mean_all<-counts_mean
    } else {
    counts_mean_all<-cbind(counts_mean_all,counts_mean)
    }
}
    df_counts <- counts_mean_all[,'UC_R',drop=FALSE]
    df_counts <- df_counts[df_counts >= 5,,drop=FALSE ]
    ggplot(df_counts, aes(x =UC_R)) +
    labs(title = paste0("Macrophage Counts distribution under UC Rectum"),
       y = "Frequency",x="Mean counts") +
    geom_histogram(position="identity",fill='grey',binwidth=100)+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
        scale_y_continuous(expand=c(0,0),limits=c(0,5000)) +
        scale_x_continuous(expand = c(0, 0))
        file_out<-file.path(dir_out,paste('Mean_raw_counts_within_macrophage_UC_R.pdf'))
        ggsave(file_out,height=5,width=6)
    breaks<-c(0,4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,84,89,94,99)
    ggplot(df_counts, aes(x = UC_R)) +
    labs(title = paste0("Macrophage Counts distribution under UC Rectum"),
       y = "Frequency",x="Mean counts") +
    geom_histogram(position="identity",fill='grey',breaks=breaks)+
    geom_hline(yintercept = c(250,500,750,1000,1250), linetype = "dashed", color = "gray51")+
    geom_vline(xintercept = seq(4,94,by=5), linetype = "dashed", color = "gray51")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
        scale_y_continuous(expand=c(0,0),limits=c(0,1500),breaks=seq(0,1500,by=250)) +
        scale_x_continuous(expand = c(0, 0),limits=c(0,100),breaks=seq(0,100,by=10))
        file_out<-file.path(dir_out,paste('Mean_raw_counts_within_macrophage_UC_R_zoomed_in_100.pdf'))
        ggsave(file_out,height=5,width=6)
        quantile(df_counts$UC_R,prob=seq(0,1,by=0.1))
#CD_R mean counts distribution
    df_counts <- counts_mean_all[,'CD_R',drop=FALSE]
    df_counts <- df_counts[df_counts >= 5,,drop=FALSE ]
    ggplot(df_counts, aes(x =CD_R)) +
    labs(title = paste0("Macrophage Counts distribution under CD Rectum"),
       y = "Frequency",x="Mean counts") +
    geom_histogram(position="identity",fill='grey',binwidth=100)+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(expand = c(0, 0))
        file_out<-file.path(dir_out,paste('Mean_raw_counts_within_macrophage_CD_R.pdf'))
        ggsave(file_out,height=5,width=6)
    breaks<-c(0,4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,84,89,94,99)
    ggplot(df_counts, aes(x = CD_R)) +
    labs(title = paste0("Macrophage Counts distribution under CD Rectum"),
       y = "Frequency",x="Mean counts") +
    geom_histogram(position="identity",fill='grey',breaks=breaks)+
    geom_hline(yintercept = c(250,500,750,1000,1250,1500,1750), linetype = "dashed", color = "gray51")+
    geom_vline(xintercept = seq(4,94,by=5), linetype = "dashed", color = "gray51")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
        scale_y_continuous(expand=c(0,0),limits=c(0,2000),breaks=seq(0,2000,by=250)) +
        scale_x_continuous(expand = c(0, 0),limits=c(0,100),breaks=seq(0,100,by=10))
        file_out<-file.path(dir_out,paste('Mean_raw_counts_within_macrophage_CD_R_zoomed_in_100.pdf'))
        ggsave(file_out,height=5,width=6)
        quantile(df_counts$CD_R,prob=seq(0,1,by=0.1))
#CD TI mean counts distribution
df_counts <- counts_mean_all[,'CD_TI',drop=FALSE]
df_counts <- df_counts[df_counts >= 5,,drop=FALSE ]
    ggplot(df_counts, aes(x =CD_TI)) +
    labs(title = paste0("Macrophage Counts distribution under CD Ileum"),
       y = "Frequency",x="Mean counts") +
    geom_histogram(position="identity",fill='grey',binwidth=100)+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(expand = c(0, 0))
        file_out<-file.path(dir_out,paste('Mean_raw_counts_within_macrophage_CD_TI.pdf'))
        ggsave(file_out,height=5,width=6)
    breaks<-c(0,4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,84,89,94,99)
    ggplot(df_counts, aes(x = CD_TI)) +
    labs(title = paste0("Macrophage Counts distribution under CD Ileum"),
       y = "Frequency",x="Mean counts") +
    geom_histogram(position="identity",fill='grey',breaks=breaks)+
    geom_hline(yintercept = c(250,500,750,1000,1250,1500), linetype = "dashed", color = "gray51")+
    geom_vline(xintercept = seq(4,94,by=5), linetype = "dashed", color = "gray51")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold",hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 10,face="bold"),  # Adjust legend text size
        legend.title = element_text(size = 12, face = "bold")) +  # Adjust legend title size
        scale_y_continuous(expand=c(0,0),limits=c(0,1750),breaks=seq(0,1750,by=250)) +
        scale_x_continuous(expand = c(0, 0),limits=c(0,100),breaks=seq(0,100,by=10))
        file_out<-file.path(dir_out,paste('Mean_raw_counts_within_macrophage_CD_TI_zoomed_in_100.pdf'))
        ggsave(file_out,height=5,width=6)
        round(quantile(df_counts$CD_TI,prob=seq(0,1,by=0.1)),2)