#!/usr/local/bin/Rscript
# title: "Papillary Thyroid Carcinoma Landscape and Its Immunological Link With Hashimoto Thyroiditis at Single-Cell Resolution"
# author: "Ruihan Luo"
# date: "Apr 24th,2023"
options(bitmapType = 'cairo')
library(Rcpp)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR, lib.loc='/home/rluo4/R/x86_64-pc-linux-gnu-library/4.1')
library(viridis)
library(DoubletFinder)
data_path <- '/data/rluo4/All/Output/Cluster/'
setwd(data_path)
# load(paste0(data_path,'../MP_new.rds'))
organ_all <- list.files(data_path)
cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
for (organ in organ_all[-(1:7)]) {
  print(organ)
  # save(Trans_G,DEG_list, Sig_G, Dat_sig, rdata_filter, file =paste0(cohort_directory, 'Epi/Results/','Malignant-Transformation-',organ, '.rds' ))
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds')
  load(file)
  print(table(rdata_filter$cell_class, rdata_filter$Trajectory))
  options(stringsAsFactors = F)
  dat = rdata_filter
  DefaultAssay(dat)="RNA"
  table(dat$cell_class)
  # dat$cell_class <- gsub('C20','C2',dat$cell_class) #
  Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
  dat.markers = FindAllMarkers(dat, only.pos = TRUE,test.use="MAST",min.pct = 0.25, logfc.threshold = 0.25)
  dat.markers$mean.1=NA
  dat.markers$mean.2=NA
  dat.markers$median.1=NA
  dat.markers$median.2=NA
  for(i in 1:nrow(dat.markers)){
    case.label=rownames(dat@meta.data)[dat$cell_class==as.character(dat.markers$cluster[i])]
    control.label=rownames(dat@meta.data)[dat$cell_class!=as.character(dat.markers$cluster[i])]
    dat.markers$mean.1[i]=mean(dat@assays$RNA@data[dat.markers$gene[i],case.label])
    dat.markers$mean.2[i]=mean(dat@assays$RNA@data[dat.markers$gene[i],control.label])
    dat.markers$median.1[i]=median(dat@assays$RNA@data[dat.markers$gene[i],case.label])
    dat.markers$median.2[i]=median(dat@assays$RNA@data[dat.markers$gene[i],control.label]) }
  f <-  paste0(cohort_directory, organ, '_deg_all.rds') 
  save(dat.markers, file = f)
}
