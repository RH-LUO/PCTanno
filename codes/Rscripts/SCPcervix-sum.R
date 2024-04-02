#!/usr/local/bin/Rscript
# title: "Integrated single-cell transcriptome analysis of the tumor ecosystems underlying cervical cancer metastasis"
# author: "Ruihan Luo"
# date: "Apr 25th,2023"
##############################################################################################################################
# SCP_CC data
##############################################################################################################################
# mv Pelvic.cavity/ Cervix
# SCP.obs <- read.csv("/public/home/lorihan/lrh/database/Cervix/SCP/SCP1950/documentation/AllGene.avg_exp.annot.xls",sep = '\t')
# SCP.obs <- read.csv("/public/home/lorihan/lrh/database/Cervix/SCP/SCP1950/cluster/tSNE_clustering_file.txt",sep = '\t')
options(bitmapType = 'cairo')
SCP.obs <- read.csv("/public/home/lorihan/lrh/database/Cervix/SCP/SCP1950/metadata/metadata.txt",sep = '\t')
SCP.obs <- SCP.obs[-1,]
table(SCP.obs$cell_type__ontology_label)
# Import libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR, lib.loc='/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(viridis)
library(DoubletFinder)
library(data.table)
library(stringr)
library(tidyverse) 
library(GEOquery)
cohort_directory <- '/public/home/lorihan/lrh/database/Cervix/SCP/'
path = '/public/home/lorihan/lrh/database/Cervix/SCP/SCP1950/expression'
setwd(cohort_directory)
# library(reshape2)
# library(Matrix)
# library(Seurat)
# solution <- list.dirs(path)
# solution
# dir <- solution[grep('62f', solution)]
# dir
# solution <- list.files(dir,pattern='*_barcodes_mod.tsv')
# batch <- str_split(solution,'_barcodes',simplify = T)[,1]
# batch
# library(DropletUtils)
# SCP.meta<- as.data.frame(batch)
# SCP.meta$sample <- dir
# # path = '/public/home/lorihan/lrh/database/Cervix/SCP/SCP1950'
# # for (i in unique(SCP.meta$batch)) {
# #   setwd(path)
# #   dir = SCP.meta$sample[match(i,SCP.meta$batch)]
# #   path_new = paste(path,i,sep = '/')
# #   dir.create(path_new)
# #   setwd(paste(dir,sep = ''))
# #   cmd = paste('mv ./',i,'*barcodes*.tsv ',path_new,'/barcodes.tsv',sep = '')
# #   system(cmd)
# #   cmd = paste('mv ./',i,'*genes*.tsv ',path_new,'/features.tsv',sep = '')
# #   system(cmd)
# #   cmd = paste('mv ./',i,'*matrix*.mtx ',path_new,'/matrix.mtx',sep = '')
# #   system(cmd)
# #   
# #   setwd(path_new)
# #   cmd = paste('gzip ','./barcodes.tsv',sep = '')
# #   system(cmd)
# #   cmd = paste('gzip ','./features.tsv',sep = '')
# #   system(cmd)
# #   cmd = paste('gzip ','./matrix.mtx',sep = '')
# #   system(cmd)
# #   }
# setwd(cohort_directory)
# path
# solution <- list.files(path)
# solution <- solution[grep('CCI', solution)]
# solution
# n <- solution
# df_empty <- NULL
# allDat <- lapply(setNames(n, n), function(nameindex) {
#   x <- Read10X(paste(path, nameindex, sep = '/'))
#   # y <- CreateSeuratObject(counts = x, project = nameindex)
#   #colnames(x) <- paste(nameindex,colnames(x),sep = '.')
#   df_empty <<- cbind(df_empty,x) 
#   return(df_empty)
# })
# dim(df_empty)#36784 80861
# remove(allDat)
# dt <- CreateSeuratObject(counts = df_empty)
# dt$batch <- str_split(colnames(dt), '[_]',simplify = T)[,2]
# table(dt$batch)
# dt$orig.ident <- paste(dt$orig.ident,dt$batch,sep = '_')
# table(dt$orig.ident)
# dt$batch <- dt$orig.ident
# erccs <- grep(pattern = "^ERCC-", x = rownames(x = dt), value = TRUE)
# percent.ercc <- Matrix::colSums(dt[erccs, ])/Matrix::colSums(dt)
# ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = dt), value = FALSE)
# table(ercc.index)
# # dt <- dt[-ercc.index,]
# dim(dt)#36784 80861
# library(stringr)
# SCP_CC <- read.csv(paste(cohort_directory,'SCP1950/metadata/metadata.txt',sep = "/"),fill = T,header = T,sep = '\t')
# SCP_CC <- SCP_CC[-1,]
# table(SCP_CC$cell_type__ontology_label)
# dt$cell_id <- str_split(rownames(dt@meta.data),'-',simplify = T)[,1]
# table(str_split(rownames(dt@meta.data),'-',simplify = T)[,2])
# SCP_CC$cell_id <- paste(SCP_CC$Sample, str_split(SCP_CC$Cell,'@',simplify = T)[,2],sep = ".")
# table(dt$cell_id %in% SCP_CC$NAME)
# table(SCP_CC$NAME %in% dt$cell_id)
# dt <- subset(dt,cell_id %in% SCP_CC$NAME)
# dt_SCP <- dt
# colnames(SCP_CC)[1] = 'cell_id'
# SCP.metadata <- dt_SCP@meta.data
# SCP.metadata <- left_join(SCP.metadata, SCP_CC,by='cell_id')#merge(meta.abn, meta.2, by='barcode')
# rownames(SCP.metadata) <- paste0(SCP.metadata$cell_id,'-1')
# # Create the Seurat object with all the data (unfiltered)
# dt_SCP <- AddMetaData(object = dt_SCP, metadata = SCP.metadata)
# table(dt_SCP$cell_type__ontology_label)
# dt_SCP$batch <- dt_SCP$orig.ident
# dt_SCP$sample_name <- rep('SCP', ncol(dt_SCP))
# dt_SCP$Tissue <- 'CC'
# epi.cells <- c('basal cell','epithelial cell','keratinocyte','malignant cell')
# imm.cells <- c('regulatory T cell','plasmacytoid dendritic cell','natural killer cell',
#                'naive B cell','macrophage','granulocyte monocyte progenitor cell','gamma-delta T cell',
#                'CD8-positive, alpha-beta T cell','CD141-positive myeloid dendritic cell','B cell')
# SCP_s <- DietSeurat(subset(dt_SCP, subset = cell_type__ontology_label %ni% c(imm.cells, epi.cells)))
# SCP_i <- DietSeurat(subset(dt_SCP, subset = cell_type__ontology_label %in% imm.cells))
# SCP_e <- DietSeurat(subset(dt_SCP, subset = cell_type__ontology_label %in% epi.cells))
# setwd('initial_clustering')
# saveRDS(SCP_i, file = "SCP_immune_initial.rds")
# saveRDS(SCP_s, file = "SCP_stromal_initial.rds")
# saveRDS(SCP_e, file = "SCP_epithelial_initial.rds")
setwd('/public/home/lorihan/lrh/database/All/')
library(SeuratDisk, lib.loc='/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(Seurat)
library(data.table)
library(stringr)
library(dplyr)
All_organ <- LoadH5Seurat("AdultCervix.h5seurat",meta.data=F)
sdata.obs <- openxlsx::read.xlsx("AdultCervix.obs.xlsx", rowNames = T)
path = '/public/home/lorihan/lrh/database/All/annotation_rmbatch_data_revised417/'
solution <- list.files('annotation_rmbatch_data_revised417/',pattern='Adult-Cervix')
filePath <- sapply(solution, function(x){ 
  paste(path,x,sep='/')})  
solutions <- lapply(filePath, function(x){
  data.table::fread(x)})
df_empty <-NULL
All.metadata <- lapply(solutions,function(y){
  df_empty <<- rbind(df_empty, y)
  return(df_empty)
})#read.table('GSE134355_RAW/HU_0246_Stomach_GSE134355_meta.txt',sep = '\t',fill = T,header = T)
All.metadata <- df_empty
table(colnames(All_organ) %in% All.metadata$Cell_id)
table(str_split(colnames(All_organ),'[-]',simplify = T)[,2])
length(unique(str_split(colnames(All_organ),'[-]',simplify = T)[,1]))
#All.metadata$Cell_id <- substr(All.metadata$Cell_id,1,33)
length(unique(All.metadata$Cell_id))#14521-->14553?
table(str_split(colnames(All_organ),'[-]',simplify = T)[,1] %in% All.metadata$Cell_id)
sdata.obs$Cell_id <- str_split(rownames(sdata.obs),'[-]',simplify = T)[,1]
sdata.obs <- left_join(sdata.obs, All.metadata[,-c(2:3)])
rownames(sdata.obs) <- colnames(All_organ)
# write.csv(sdata.obs,'/public/home/lorihan/lrh/database/All/Stomach_cell_type.csv')
All_organ@meta.data <- sdata.obs
table(sdata.obs$Celltype)
All_organ$Tissue <- All_organ$tissue
All_organ$sample_name <- 'Adult'
table(All_organ$batch)
All_organ$orig.ident <- All_organ$batch
Imm.cells <- sdata.obs$Celltype[c(grep('B cell',sdata.obs$Celltype),grep('Mast cell',sdata.obs$Celltype),
                                  grep('CD8_T cell',sdata.obs$Celltype), grep('M1 Macrophage',sdata.obs$Celltype))]#9848
table(Imm.cells)
Str.cells <- sdata.obs$Celltype[c(grep('Endothelial', sdata.obs$Celltype),
                                  grep('Smooth',sdata.obs$Celltype),
                                  # grep('muscle',sdata.obs$Celltype),
                                  grep('Stromal',sdata.obs$Celltype)
)]
table(Str.cells)
# ------
All_s <- subset(All_organ, Celltype %in% Str.cells)
All_i <- subset(All_organ, Celltype %in% Imm.cells)
All_e <- subset(All_organ, Celltype %ni% c(Str.cells,Imm.cells))
All_e
#---------snRNA_all_epithelial_analysis.R#------
# Import libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR, lib.loc='/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(viridis)
library(DoubletFinder)
set.seed(1)
packageVersion('Rcpp')#
# assuming that you have some Seurat object cWiliamsed colon_Wiliams_epi:
SCP_e <- readRDS('/public/home/lorihan/lrh/database/Cervix/SCP/SCP_epithelial_initial.rds')
table(SCP_e$Tissue)
Hua1_e <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua1/initial_clustering/Hua1_epithelial_initial.rds')
Hua2_e <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua2/initial_clustering/Hua2_epithelial_initial.rds')
Hua3_e <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua3/initial_clustering/Hua3_epithelial_initial.rds')
Hua4_e <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua4/initial_clustering/Hua4_epithelial_initial.rds')
Guo_e <- readRDS('/public/home/lorihan/lrh/database/Cervix/Guo/initial_clustering/Guo_epithelial_initial.rds')
table(Guo_e$batch,Guo_e$Tissue)
#Hua1_e$Tissue <- 'ESCC'#substr(Hua1_e$orig.ident,1,3)
table(Hua3_e$Tissue)
table(Guo_e$Tissue)
colon_full <- merge(subset(Guo_e, Tissue!='NO_HPV'),c(SCP_e,Hua1_e,Hua2_e,Hua3_e,Hua4_e))
colon_full
table(colon_full$batch)
table(colon_full$Tissue)
All_e$orig.ident <- All_e$batch
table(All_e$orig.ident)
table(rownames(All_e) %in% rownames(colon_full))
packageVersion('irlba')
###############################################################################################################################
###############################################################################################################################
# Import libraries
# lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ArchR, lib.loc='/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(future, lib.loc = "/public/home/lorihan/miniconda3/lib/R/library")
library(Matrix, lib.loc = "/public/apps/R-4.1.0/library")
set.seed(1)
cohort_directory <- '/public/home/lorihan/lrh/database/Cervix/SCP'
###############################################################################################################################
###############################################################################################################################
# Define variables
execute_steps <- 1:10
# Define variables
sample_name <- "all_samples"
# folder to save the results
analysis_parent_folder <- "./epithelial_results/"
# load seurat object containing all epithelial cells
# load metadata
# metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
# LSI Params
resolution <- c(0.1,0.2,0.4,0.8)
umapNeighbors <- 50
umapMinDist <- 0.45
umapDistMetric <- "cosine"
nTop <- 3200
initial_nPCs <- 1:32#1:8
final_nPCs <- 1:32#c(1:4,6:8)

# nTop <- 1600#for only Fillio Healthy
# initial_nPCs <- 1:8
# final_nPCs <- 1:8#c(1:4,6:8)

nPCs <- initial_nPCs
# Create directory
setwd(cohort_directory)
if (!dir.exists(paste0(analysis_parent_folder))){
  dir.create(paste0(analysis_parent_folder))
}
setwd(analysis_parent_folder)
getwd()

###############################################################################################################################
###############################################################################################################################
# Define Functions--as noted above, multiple functions copied from Granja et al 2019
sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  using namespace std;
  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)          

sparseRowVariances <- function (m){
  #Compute Fast Sparse Row Variances--From Granja et al 2019
  rM <- Matrix::rowMeans(m)
  rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
  return(rV)
}

getVarGenesFilterBlacklist <- function(mat, nTopGenes = 2000, blacklist = NULL){
  # Get the top nTopGenes variable genes present in mat (a gene x sample/cell matrix)
  # If blacklist is present, do not return any genes in the blacklist
  if(!is.null(blacklist)){
    mat <- mat[!rownames(mat) %in% blacklist,]
  }
  # compute row variances and return the most variable genes for eith matrix or sparse matrix
  if (class(mat) == "matrix"){
    rownames(mat)[head(order(rowVars(mat), decreasing = TRUE), nTopGenes)]
  } else {
    rownames(mat)[head(order(sparseRowVariances(mat), decreasing = TRUE), nTopGenes)]
  }
}

calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
  #From Granja et al 2019
  set.seed(1)
  
  #TF IDF LSI adapted from flyATAC
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1
  }
  
  #Calc RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm)
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  #Return Object
  out <- list(
    matSVD = matSVD, 
    rowSm = rowSm, 
    colSm = colSm, 
    svd = svd, 
    binarize = binarize, 
    nComponents = nComponents,
    seed = 1)
  
  out
  
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

projectLSI <- function(mat, lsi, binarize){   
  # print(length(lsi$varFeatures))
  # lsi$svd$u <- lsi$svd$u[names(lsi$rowSm) %in% rownames(mat)]  
  # lsi$varFeatures <- lsi$varFeatures[lsi$varFeatures %in% rownames(mat)] 
  # lsi$rowSm <- lsi$rowSm[names(lsi$rowSm) %in% rownames(mat)] 
  #Get Same Features
  mat <- mat[lsi$varFeatures,]
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1       
  }
  
  #Calc TF IDF
  rowsToZero <- which(lsi$rowSm == 0)
  setToZero <- which((mat@i + 1) %in% rowsToZero)
  if(length(setToZero) > 0){
    mat@x[setToZero] <- 0
  }
  
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
    tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
  }
  
  #Calc V
  V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
  diag(svdDiag) <- lsi$svd$d
  matSVD <- t(svdDiag %*% t(V))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  return(matSVD)
}
seurat_feature_plot <- function(sample_name, reduction, cell_type, markers){
  # function to make grid layout of seurat feature plots
  p1 <- FeaturePlot(colon, features = markers, reduction = reduction, sort.cell = TRUE, combine = FALSE, pt.size = 2)
  if (length(p1)==1){
    width <- 4
    height <- 4
  } else if (length(p1)==2){
    width <- 8
    height <- 4
  } else if (length(p1)<5){
    width <- 8
    height <- 8
  } else if (length(p1)<7){
    width <- 12
    height <- 8
  } else if (length(p1)<10){
    width <- 12
    height <- 12
  } else if (length(p1)<13){
    width <- 16
    height <- 12
  } else if (length(p1)<17){
    width <- 16
    height <- 16
  }
  else if (length(p1)>18){
    width <- 20
    height <- 18
  }
  pdf(paste0("./", reduction, "_feature_plot_", sample_name, "_", cell_type ,".pdf"), width = width, height = height)
  print(CombinePlots(lapply(p1, function (x) AugmentPlot(x))))
  dev.off()
}

vln_plot <- function(features, save_name){
  pdf(save_name, width = 20, onefile=F)
  print(VlnPlot(colon, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
          geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
          scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
  dev.off()
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) Subset to create normal project
All_e$Tissue <- All_e$tissue
table(All_e$Tissue)
All_e$sample_name <- 'AdultCervix'
table(All_e$Celltype)
if (1 %in% execute_steps){
  # define normal project containing all normal samples except B001-A-104, which was much lower quality and was not used in reference as a result
   colon <- DietSeurat(subset(Guo_e, subset = Tissue == "NO_HPV"))
   Cervix.Adult <- DietSeurat(subset(All_e, subset = Celltype != "Unknown"))
   # colon <- merge(All_e, y=colon)
   colon <- merge(Cervix.Adult, y=colon)
}
unique(colon$batch)#5
unique(colon$sample_name)#5
colon#6979 -->5579
table(colon$Tissue)
table(Idents(colon))
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) QC
if (2 %in% execute_steps){
  # create and set working directory to save qc plots
  setwd(cohort_directory)
  if (!dir.exists(paste0(analysis_parent_folder, "all_samples_qc_plots"))){
    dir.create(paste0(analysis_parent_folder, "all_samples_qc_plots"))
  }
  setwd(paste0(analysis_parent_folder, "all_samples_qc_plots"))
  
  # already done, but can add here if you got here a different way
  colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")
  # Now subset the project and make some nice qc plots for just the epithelial cells--this was already done in the first step
  colon <- subset(colon, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 25 & nCount_RNA < 10000)
  vln_plot(c("nFeature_RNA"), paste0("./n_genes_violin_", sample_name, ".pdf"))
  
table(colon_full$orig.ident)
table(rownames(colon) %in% rownames(colon_full))
colon <- colon[rownames(colon) %in% rownames(colon_full),]# this step was critical
colon #17637 features across 4404 samples 
vln_plot(c("nCount_RNA"), paste0("./n_counts_violin_", sample_name, ".pdf"))
vln_plot(c("percent.mt"), paste0("./pMT_violin_", sample_name, ".pdf"))
setwd(cohort_directory)
setwd(analysis_parent_folder)
}
colon
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Dim reduction with LSI  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(colon)
colon <- ScaleData(colon, features = all.genes)#, vars.to.regress = c("orig.ident"))#c("nFeature_RNA", "percent.mito","percent.ribo",

##############################################################################################################################
if (3 %in% execute_steps){
  #Initialize list for storing iterative LSI output
  lsiOut <- list()
  nPCs <- initial_nPCs
  # Get the raw counts
  rawCounts <- GetAssayData(object = colon, slot = "counts")
  
  # Identify genes (mitochondrial, ribosomal, and HLA) to blacklist for dimensionality reduction
  mt_genes <- grep("^MT-", rownames(colon), value = TRUE)
  rps_genes <- grep("^RPS", rownames(colon), value = TRUE)
  rpl_genes <- grep("^RPL", rownames(colon), value = TRUE)
  hla_genes <- grep("^HLA-", rownames(colon), value = TRUE)
  blacklist <- c(mt_genes, rps_genes, rpl_genes, hla_genes)
  
  # Here we will use the normalized data in the data slot from seurat, which is ln rather than log2 normalized
  log2CP10k <- GetAssayData(object = colon, slot = "data")
  # log2CP10k <- cca.combined@assays$integrated@data#GetAssayData(object = colon, slot = "data")
  # alternatively could do:
  # log2CP10k <- as(log2(t(t(rawCounts)/Matrix::colSums(rawCounts)) * 10000 + 1), "matrix")
  
  for(i in seq_along(resolution)){
    if(i == 1){
      #Initial LSI uses variances that are across all single cells and will have larger batch relationships
      message("Running initial LSI...")
      varGenes <- getVarGenesFilterBlacklist(log2CP10k, nTopGenes = nTop, blacklist = blacklist)
    }else{
      message(sprintf("Running LSI %s of %s...", i,  length(resolution)))
      # Calculate variable genes using clusters defined in previous LSI iteration
      clusterMat <- edgeR::cpm(groupSums(rawCounts, clusters, sparse = TRUE), log=TRUE, prior.count = 3)
      varGenes <- getVarGenesFilterBlacklist(clusterMat, nTopGenes = nTop, blacklist = blacklist)
    }
    
    # Now run LSI and find clusters
    LSIi <- calcLSI(rawCounts[varGenes,], nComponents = max(nPCs), binarize = FALSE)
    colon[[paste0("LSI_iter",i)]] <- CreateDimReducObject(embeddings = LSIi$matSVD, key = sprintf("LSI%s_", i), assay = "RNA")
    colon <- FindNeighbors(object = colon, reduction = paste0("LSI_iter",i), dims = nPCs, force.recalc = TRUE)
    colon <- FindClusters(object = colon, resolution = resolution[i])
    clusters <- Idents(colon)
    table(Idents(colon))
    # Save LSI iteration
    lsiOut[[paste0("LSI_iter",i)]] <- list(
      lsiMat = LSIi$matSVD, 
      varFeatures = varGenes, 
      clusters = clusters,
      colSm = LSIi$colSm,
      rowSm = LSIi$rowSm,
      svd = LSIi$svd, 
      binarize = LSIi$binarize, 
      nComponents = LSIi$nComponents
    )
  }
  
  # Run uwot to compute the UMAP
  nPCs <- final_nPCs
  uwotUmap <- uwot::umap(
    LSIi$matSVD[,nPCs], 
    n_neighbors = umapNeighbors, 
    min_dist = umapMinDist, 
    metric = umapDistMetric, 
    n_threads = 1, 
    verbose = TRUE, 
    ret_nn = TRUE,
    ret_model = TRUE
  )
  
  # Add to the seurat object and plot
  umap_vals <- uwotUmap[[1]][,1:2]
  rownames(umap_vals) <- rownames(LSIi$matSVD[,nPCs])
  colon[["uwot_UMAP"]] <- CreateDimReducObject(embeddings = umap_vals, key = "uwot_UMAP", assay = "RNA")
  
  # pdf(paste0("./uwot_UMAP_samples_", nTop, "_vargenes_all_", max(nPCs), "PCs.pdf"), width = 8)
  # print(DimPlot(colon, reduction = "uwot_UMAP", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
  # dev.off()
  
  # Recluster with final nPCs
  colon <- FindNeighbors(object = colon, reduction = paste0("LSI_iter",i), dims = nPCs, force.recalc = TRUE)
  colon <- FindClusters(object = colon, resolution = 1.0)
  table(Idents(colon))
}
saveRDS(colon, "clustered_Healthy_colon_proj_seurat.rds")
save(nPCs, nTop, LSIi,lsiOut,umap_vals,varGenes, file="lsiOut.RData")
uwot::save_uwot(uwotUmap, file='colon_umap', unload = T)
pdf(paste0("./UMAPclustering" , ".pdf"))
DimPlot(colon,reduction = "uwot_UMAP", label = T, pt.size = 1.5,label.size = 8,label.color = 'black',
        cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))+
  NoLegend()+labs(x = "UMAP 1", y = "UMAP 2",title =  "Resolution.1.uwot_UMAP") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
# print(DimPlot(colon, reduction = "uwot_UMAP", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
dev.off()

pdf(paste0("./UMAP_cell_type_orig.pdf"), height=5,width = 5,onefile = F)
DimPlot(colon, group.by = "Celltype", reduction = "uwot_UMAP", label = T, pt.size = 1.5,label.size = 2.5,label.color = 'black',
        cols = paletteDiscrete(values = unique(colon@meta.data$Celltype), set = "stallion", reverse = FALSE))+
  NoLegend()+labs(x = "UMAP 1", y = "UMAP 2",title =  "Epithelial cell type") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()
# options(bitmapType = 'cairo')
# DimPlot(colon, reduction = "uwot_UMAP", label=T,
#         cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))+NoLegend()+
#   labs(x = "UMAP1", y = "UMAP2") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#                                          axis.text.x = element_blank(),axis.ticks.x = element_blank())+
#   theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #加边框
table(Idents(colon))# first version 1:32 3200  1:8 1600
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
# 644 640 537 457 432 387 340 333 329 299 234 221 210 181 144 141 120 115 111 
# 19  20  21  22  23  24 
# 82  65  54  52  49  24 
# final version (not nohup but Rstudio)
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
# 467 439 387 383 367 340 326 317 254 223 192 188 181 176 166 
# 15  16  17  18  19  20  21  22  23  24  25 
# 152 150 142 120 114 111 107 100  98  55  24 
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) ID Cluster Markers to help with manual annotation
# if (4 %in% execute_steps){
#   colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
# }
table(colon$Celltype)
table(SCP_e$cell_type__ontology_label)
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) Plot Some Markers on the UMAPs to help with annotation
if (5 %in% execute_steps){
  reductions_to_plot <- c("uwot_UMAP")
  for (reduction in reductions_to_plot){
    ##  Squamous epithelial cells: These are the most abundant cells in the cervix and are flat, scale-like cells that form the outer layer of the cervix. Squamous epithelial cells are stratified, meaning that they are arranged in layers, and are responsible for providing a protective barrier to the underlying tissues.
    # In the cervix, the squamous epithelial layer consists mainly of keratinocytes, which are derived from the basal layer of the epithelium.
    #Keratin genes, such as KRT5 and KRT14
    seurat_feature_plot(sample_name, reduction, 'BAS', c("KRT5", "KRT14", 'KRT6A',"KRT15","TP63","ITGA6","ITGB4")) #basal layer/differentiating keratinocytes
    seurat_feature_plot(sample_name, reduction, 'MES', c("WT1", "CALB1", "VIM", "CDH1", "MSLN", "BAP1"))
    seurat_feature_plot(sample_name, reduction, 'COL', c('SOX17',"CLDN15","MUC17",'MUC16',"MUC5B", "MUC6", "PAX8")) #mucinous‐?
    # seurat_feature_plot(sample_name, reduction, 'SUP', c("KRT4", "KRT13", "KRT18", "FLG", "TGM1"))
    # seurat_feature_plot(sample_name, reduction, 'GOB', c("MUC2", 'MUC6', "AGR2","FOXA3", 'ATOH1','TFF3',"HES6"))
    # seurat_feature_plot(sample_name, reduction, 'QUIE', c("DST", 'DLK2', "KRT15"))
    # seurat_feature_plot(sample_name, reduction, 'PRO', c("HMGB2", "TOP2A"))  #keratinization and cornification (SPRR3, CRNN, CNFN) typically encountered in the outmost layers of the epitheSCPm
    seurat_feature_plot(sample_name, reduction, 'KER', c('KRT5', 'KRT13','KRT14', 'KRT4','KRT18', 'KRT17', 'IVL','LOR','SPRR3'))  #keratinocytes cell (differentiated)
    seurat_feature_plot(sample_name, reduction, 'TRANS', c("UPK1A",'DSC', 'DSP', 'SERPINB3')) #
    # seurat_feature_plot(sample_name, reduction, 'DIFF', c('PADI1', 'TGM3', 'MT1G', 'RNASE7', 'KRT16', 'RNASE7'))  #
    seurat_feature_plot(sample_name, reduction, 'Stem', c('KRT15','KRT5', 'KRT14','CD44','PROM1','ABCG2','KRT19','TP63',
                                                          'ITGB1', 'ITGA6', 'ALDH1A1', 'ABCB5','ABCG2',"ALCAM"))#stem cells
    seurat_feature_plot(sample_name, reduction, 'STM', c("ALDH3A1", "ALDH1A1", 'BMI1', 'SFRP1', 'ITGB1',
                                                         'LGR5','SOX2','POU5F1','NANOG','CD44','TP63','ABCG2', 'PROM1', 'ALCAM'))#from ChatGPT #'CD133'(PROM1),'CD166'(ALCAM)
  }#POU5F1 (also known as OCT4)  ITGB1 (also known as CD29)
}
table(colon$seurat_clusters,colon$Celltype)
# biomarkers from 'Single‐Cell Landscape Highlights Heterogenous Microenvironment, Novel Immune Reaction Patterns, Potential Biomarkers and Unique Therapeutic Strategies of Cervical Squamous Carcinoma, Human Papillomavirus‐Associated (HPVA) and Non‐HPVA Adenocarcinoma'
# 1) stratified squamous epithelial cells (SSECs; KRT5);
# 2) gastrointestinal glandular epithelial cells (GGECs; ANXA10, CLDN18);
# 3) mucinous glandular epithelial cells (MGECs; MUC17 and CLDN15);
# 4) usual‐type glandular epithelial cells (UGECs; MUC16); and 5) NMECs (SOX17). To elucidate the relationship between HPV infection and cell malignancy, we compared the distribution of HPV‐infected epithelial cell
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) Define clusters
if (6 %in% execute_steps){
  # cluster identification
  new.cluster.ids <- c(  "BAS", #0
                         "COL", #1
                         "STM", #2
                         "KER", #3 #
                         "TRANS", #4 #
                         "STM", #5 ?KER
                         "COL", #6
                         "MES", #7
                         "STM", #8 #?KER
                         "KER", #9
                         "KER", #10
                         "COL", #11
                         "KER",  #12
                         "MES",  #13
                         "STM", #14 ?STM
                         "KER", #15
                         "TRANS", #16
                         "STM", #17
                         "STM", #18
                         "STM", #19
                         "KER", #20
                         "STM", #21
                         "BAS", #22
                         "KER",  #23
                         "KER", #24
                         "BAS") #25
  # new.cluster.ids <- c("COL", #0
  #                     "STM", #1
  #                     "BAS", #2
  #                     "TRANS", #3 #?KER
  #                     "COL", #4 #?KER
  #                     "STM", #5
  #                     "STM", #6 ?
  #                     "MES", #7
  #                     "TRANS", #8 #?KER
  #                     "TRANS", #9
  #                     "BAS", #10
  #                     "COL", #11 ?STM
  #                     "STM",  #12
  #                     "KER",  #13
  #                     "KER", #14 ?STM
  #                     "MES", #15
  #                     "STM", #16
  #                     "STM", #17
  #                     "KER", #18
  #                     "KER", #19
  #                     "BAS", #20
  #                     "KER", #21
  #                     "BAS", #22 ?STM
  #                     "STM",  #23
  #                     "BAS") #24
  # Set the cell types in the project
  identities <- data.frame(colon[['seurat_clusters']])
  identities$seurat_clusters <- as.character(colon[['seurat_clusters']]$seurat_clusters)
  for (i in 0:(length(new.cluster.ids)-1)){
    identities[identities$seurat_clusters==as.character(i),] <- new.cluster.ids[i+1]
  }
  # plot
  colon <- AddMetaData(colon, identities$seurat_clusters, col.name = "Cell_Type")
  #   pdf(paste0("./UMAP_cell_type_id.pdf"), width = 6,onefile = F)
  #   DimPlot(colon, group.by = "Cell_Type", reduction = "uwot_UMAP", cols = paletteDiscrete(values = unique(colon@meta.data$CellType), set = "stallion", reverse = FALSE)) + theme_ArchR()
  #   dev.off()
}
pdf(paste0("./UMAP_cell_type_id.pdf"), height=5.5,width = 5.5,onefile = F)
DimPlot(colon, group.by = "Cell_Type", reduction = "uwot_UMAP", label = T, pt.size = 1.5,label.size = 6,label.color = 'black',
        cols = paletteDiscrete(values = unique(colon@meta.data$Cell_Type), set = "stallion", reverse = FALSE))+
  NoLegend()+labs(x = "UMAP 1", y = "UMAP 2",title =  "Epithelial Cell Type") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()
pdf(paste0("./UMAP_cell_type_id.nolabel.pdf"), height=5.5,width = 5.5,onefile = F)
DimPlot(colon, group.by = "Cell_Type", reduction = "uwot_UMAP", label = F, pt.size = 1.5,label.size = 6,label.color = 'black',
        cols = paletteDiscrete(values = unique(colon@meta.data$Cell_Type), set = "stallion", reverse = FALSE))+
  NoLegend()+labs(x = "UMAP 1", y = "UMAP 2",title =  "Epithelial Cell Type") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) Project all samples into normal manifold
load("lsiOut.RData")
uwotUmap <- uwot::load_uwot('colon_umap')
colon_full#
table(colon_full$orig.ident)
t <- data.frame(table(colon_full$orig.ident))
colon_full <- colon_full[,!colon_full$orig.ident %in% t$Var1[t$Freq==1]]# this sample with one cell was removed
#colon_full <- colon_full[rownames(colon_full) %in% rownames(colon),]
paletteDiscrete(values = unique(colon@meta.data$Cell_Type))
if (7 %in% execute_steps){
  # Subset data for a sample
  umap_df <- NULL
  cell_type_df <- DataFrame()#data.frame()
  for (sampleName in unique(colon_full@meta.data$orig.ident)){
    # Subset seurat project to get sample only and get the normalized data
    temp <- subset(colon_full, subset = orig.ident == sampleName)
    rawCounts <- GetAssayData(object = temp, slot = "counts")

    # project in the previously defined lsi dimensionality reduction and add to temp seurat project--assumes lsiOut is still around, save it and load it if its not
    lsiProjectionMat <- projectLSI(rawCounts,lsiOut[[paste0("LSI_iter",length(resolution))]], binarize = FALSE)
    temp[["LSI_project"]] <- CreateDimReducObject(embeddings = as.matrix(lsiProjectionMat), key = sprintf("LSI_project"), assay = "RNA")

    # project into umap and add to temp seurat project
    umapProjection <- uwot::umap_transform(as.matrix(Embeddings(temp, reduction = "LSI_project"))[,nPCs], uwotUmap, verbose = TRUE)
    proDF <- data.frame(X1 = umapProjection[,1], X2 = umapProjection[,2])
    rownames(proDF) <- rownames(Embeddings(temp, reduction = "LSI_project"))
    temp[["uwot_UMAP_projection"]] <- CreateDimReducObject(embeddings = as.matrix(proDF), key = "uwot_UMAP_projection", assay = "RNA")

    # now plot the projection and the reference
    refDF <- data.frame(X1 = uwotUmap$embedding[,1], X2 = uwotUmap$embedding[,2], Type = "Reference")
    proDF <- data.frame(X1 = umapProjection[,1], X2 = umapProjection[,2], Type = "Sample")
    projectionDF <- rbind(refDF, proDF)
    p <- ggplot(projectionDF, aes(X1,X2,color=Type)) +
      geom_point(size = 0.1, alpha = 1) +
      xlab("UMAP Dimension 1") +
      ylab("UMAP Dimension 2") +
      theme_ArchR(baseSize = 10) + ggtitle("projection") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
      scale_color_manual(values=c("Reference"="#D5D5D5","Sample"="#008F95"))
    pdf(paste0("./uwot_projection_", sampleName, "_", nTop, "_vargenes_all_", max(nPCs), "PCs.pdf"), width = 12)
    print(p)
    dev.off()

    # Now id cell based on nearest neghbors
    library(FNN)
    input_knn <- 25
    nPCs <- final_nPCs
    svdReference <- as.data.frame(lsiOut[[paste0("LSI_iter",length(resolution))]]$lsiMat)
    svdDisease <- as.data.frame(as.matrix(lsiProjectionMat))
    knnDisease <- get.knnx(
      data = svdReference[,nPCs],
      query = svdDisease[,nPCs],
      k = input_knn)

    cellTypes <- c()
    for (j in 1:length(knnDisease$nn.index[,1])){
      types <- colon@meta.data$Cell_Type[knnDisease$nn.index[j,]]
      cellTypes <- c(cellTypes,labels(sort(table(types),decreasing=TRUE)[1]))
    }
    refDF$CellType <- "Reference"
    proDF$CellType <- cellTypes
    projectionDF2 <- rbind(refDF, proDF)

    pal <- paletteDiscrete(values = unique(proDF$CellType), set = "stallion", reverse = FALSE)
    pal["Reference"] <- "#D5D5D5"
    p <- ggplot(projectionDF2, aes(X1,X2,color=CellType)) +
      geom_point(size = 2, alpha = 1) +
      xlab("UMAP Dimension 1") +
      ylab("UMAP Dimension 2") +
      theme_ArchR(baseSize = 10) + ggtitle(sampleName) +
      # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
      scale_color_manual(values=c("Reference"="#D5D5D5", 'BAS'='#D51F26','COL'='#272E6A','KER'='#208A42','MES'='#89288F', 'STM'='#F47D2B', 'TRANS'='#FEE500'))+ #,'TSK'='#D8A767'
      theme(
        plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
        axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
        legend.title  = element_text(color = 'black', size  = 16),
        legend.text   = element_text(color = 'black', size   = 16),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
        panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
      )
    pdf(paste0("./uwot_projection_", sampleName, "_", nTop, "_vargenes_all_", max(nPCs), "PCs.pdf"), width = 6.5)
    print(p)
    dev.off()

    current_cell_types <- DataFrame("Cell" = rownames(svdDisease), "CellType" = cellTypes)
    cell_type_df <- rbind(cell_type_df, current_cell_types)
    umap_df <- rbind(umap_df, proDF)
  }
  write.csv(cell_type_df, "allEpithelialCellTypes.csv")
  colon_full <- AddMetaData(colon_full, cell_type_df$CellType, col.name = "Cell_Type")
  write.csv(umap_df[,1:2], file = "cell_embeddings.csv")#,row.names = F)
  
  # save seurat objects
  saveRDS(colon, "clustered_normal_colon_proj_seurat.rds")#cluster 10 == TFC
  saveRDS(DietSeurat(colon), "diet_clustered_normal_colon_proj_seurat.rds")#cluster 10 = iATC

  saveRDS(colon_full, "clustered_all_samples_epithelial_colon_proj_seurat.rds")
  saveRDS(DietSeurat(colon_full), "diet_clustered_all_samples_epithelial_colon_proj_seurat.rds")
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 8) Compute differential tests of stem cells relative to normal and unaffected stem cells
# if (8 %in% execute_steps){
colon_full <- readRDS("diet_clustered_all_samples_epithelial_colon_proj_seurat.rds")
colon <- readRDS("diet_clustered_normal_colon_proj_seurat.rds")
table(colon$batch)
table(colon_full@meta.data$orig.ident )#<- colon_full@meta.data$HTAN.Parent.Data.File.ID
colon_full <- merge(colon_full,y=colon)
# backgrounds <- c("Unaffected", "NormalColon")
table(colon_full$Tissue)
# table(colon_full$tissue)
table(colon_full$Cell_Type)
background <- 'NormalColon'
colon_full@meta.data$SimplifiedSampleName <- colon_full@meta.data$orig.ident
STM <- data.frame(table(colon_full$SimplifiedSampleName,colon_full$Cell_Type))
STM <- STM[STM$Var2=='STM',]
STM$Tissue <- colon_full$Tissue[match(STM$Var1,colon_full$SimplifiedSampleName)]
STM$sample_name <- colon_full$sample_name[match(STM$Var1,colon_full$SimplifiedSampleName)]
summary(STM$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    61.5   134.0   311.2   473.5  1822.0
STM <- STM[STM$Tissue %ni% c("AdultCervix",'NO_HPV'),]
#ggplot(STM, aes(x = (Freq))) +   geom_histogram()
summary(STM$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.0    62.0   129.0   287.2   455.0  1822.0
table(colon_full$sample_name, colon_full$Tissue)
table(colon$Tissue)
table(colon_full$Tissue, colon_full$SimplifiedSampleName)
table(colon$Tissue)
colon_full@meta.data$SimplifiedSampleName[colon_full@meta.data$Tissue %in% c('AdultCervix','NO_HPV')]='NormalColon'
table(colon_full@meta.data$SimplifiedSampleName)
table(colon_full@meta.data$batch)
# first set up a colon test object where the idents include test and background groups
# for (background in backgrounds){
background = "NormalColon"
if (background == "NormalColon"){
  sample_name_list <- paste0(unique(colon_full@meta.data$SimplifiedSampleName), "STM")
  sample_name_list <- sample_name_list[sample_name_list != "NormalColonSTM"]
  colon_full <- AddMetaData(colon_full, paste0(colon_full@meta.data$SimplifiedSampleName, colon_full@meta.data$Cell_Type), col.name = "CellTypeSimplifiedSample")
  colon_test_object <- colon_full
  Idents(colon_test_object) <- "CellTypeSimplifiedSample"
  background_sample <- "NormalColonSTM"
}
table(colon_full$orig.ident,colon_full$Cell_Type)
table(colon_full$SimplifiedSampleName,colon_full$Cell_Type)
# identify samples with at least 100 cells
cutoff <- ifelse(summary(STM$Freq)[2]>100,100,
                 ifelse(summary(STM$Freq)[2]<10,10,summary(STM$Freq)[2]))
cutoff
new_sample_list <- c()
for (i in 1:length(sample_name_list)){
  current_sample <- sample_name_list[i]
  if (sum(Idents(colon_test_object) == current_sample)>= cutoff){
    new_sample_list <- c(new_sample_list, current_sample)
  }
}
setdiff(sample_name_list, new_sample_list)#7
sample_name_list <- new_sample_list#
# compute differential genes using seurats findMarkers, merge into single object
# for (i in 1:length(sample_name_list)){
#   current_sample <- sample_name_list[i]
#   message(paste0("Computing differential genes for ", current_sample))
#   # added use.test = "DESeq2" and got similar results
#   differential_test <- FindMarkers(colon_test_object, ident.1 = current_sample, ident.2 = background_sample, verbose = FALSE, min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, max.cells.per.ident = 300, test.use = "MAST")
#   colnames(differential_test) <- paste0(colnames(differential_test), current_sample)
#   if (i == 1){
#     colon_full_diff_test <- differential_test
#   } else {
#     message(paste0("Merging differential expression data for ", current_sample))
#     colon_full_diff_test <- merge(colon_full_diff_test, differential_test, by=0, all=TRUE)
#     rownames(colon_full_diff_test) <- colon_full_diff_test$Row.names
#     #    colon_full_diff_test <- colon_full_diff_test[,colnames(colon_full_diff_test) %ni% c("Row.names")]
#     colon_full_diff_test <- colon_full_diff_test[,! colnames(colon_full_diff_test) %in% c("Row.names")]
#   }
# }
# saveRDS(colon_full_diff_test, paste0("./colon_full_diff_test_", background_sample, "_background.rds"))
# }}
table(colon_full$sample_name, colon_full$Tissue)
table(colon_full$SimplifiedSampleName,colon_full$Tissue)
tab <- data.frame(table(colon_full$SimplifiedSampleName,colon_full$Tissue))
tab <- tab[tab$Freq!=0,]
normal <- paste0(unique(tab$Var1[tab$Var2 %in% c('normal')]),'STM')
# new_sample_list <- new_sample_list[!grepl('NORM',new_sample_list)]#
new_sample_list <- new_sample_list[!new_sample_list %in% normal]#
sample_name_list <- new_sample_list#70 in 107 for >20
sample_name_list #77
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 9) Compute malignancy continuum
if (9 %in% execute_steps){
  # load the differential test relative to normal stem cells
  colon_full_diff_test <- readRDS("./Cervix_full_diff_test_NormalColonSTM_background.rds")#"./colon_full_diff_test_NormalColonSTM_background.rds")
  background <- "Healthy"
  typeA <- 'STM'#"STM"

  # select just the logfc and p values
  colon_full_diff_test_avg_logFC <- colon_full_diff_test[,grepl("avg_log2FC", colnames(colon_full_diff_test))]
  # colon_full_diff_test_avg_logFC <- colon_full_diff_test_avg_logFC/0.69314718056
  colon_full_diff_test_p_val_adj <- colon_full_diff_test[,grepl("p_val_adj", colnames(colon_full_diff_test))]
  colon_full_diff_test_avg_logFC[is.na(colon_full_diff_test_avg_logFC)] <- 0
  colon_full_diff_test_p_val_adj[is.na(colon_full_diff_test_p_val_adj)] <- 1

  # select significant genes with the following cutoffs
  colon_full_diff_test_avg_logFC_significant <- colon_full_diff_test_avg_logFC[(rowSums(colon_full_diff_test_avg_logFC>0.5 & colon_full_diff_test_p_val_adj<0.05)>1 | rowSums(colon_full_diff_test_avg_logFC<(-0.5) & colon_full_diff_test_p_val_adj<0.05)>1), ]
  colon_full_diff_test_p_val_adj_significant <- colon_full_diff_test_p_val_adj[(rowSums(colon_full_diff_test_avg_logFC>0.5 & colon_full_diff_test_p_val_adj<0.05)>1 | rowSums(colon_full_diff_test_avg_logFC<(-0.5) & colon_full_diff_test_p_val_adj<0.05)>1), ]

  sig_sample <- substr(colnames(colon_full_diff_test_avg_logFC_significant),11,nchar(colon_full_diff_test_avg_logFC_significant))
  table(sig_sample %in% new_sample_list)
  setdiff(new_sample_list,sig_sample)
  setdiff(sig_sample, new_sample_list)
  table(!is.na(match(new_sample_list,sig_sample)))
  new_sample_list <- new_sample_list[new_sample_list %in% sig_sample]
  colon_full_diff_test_avg_logFC_significant <- colon_full_diff_test_avg_logFC_significant[,match(new_sample_list,sig_sample)]
  colon_full_diff_test_p_val_adj_significant <- colon_full_diff_test_p_val_adj_significant[,match(new_sample_list,sig_sample)]
  # compute the pcs on the logfold changes
  pcs <- prcomp(colon_full_diff_test_avg_logFC_significant)
  pc_df <- data.frame(pcs$rotation)
  pc_df$sample <- rownames(pc_df)
  pc_df$sample
  pc_df <- pc_df[ substr(rownames(pc_df),11,nchar(rownames(pc_df))) %in% new_sample_list,]
  # remove atac column from metadata
  # load metadata
  gsub('STM','',new_sample_list)
  sample_names <- substr(rownames(pc_df),11,nchar(rownames(pc_df))-(nchar(typeA)))
  metadata <- colon_full@meta.data[match(sample_names,colon_full$orig.ident),]#read.table("../hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
  # metadata <- metadata[,colnames(metadata)[2:ncol(metadata)]]#[2:28]]
  colnames(metadata) <- c("Sample", colnames(metadata)[2:ncol(metadata)])#[2:27])
  # metadata$Sample <- str_split(metadata$Sample, '_cd45',simplify = T)[,1]
  metadata <- metadata[metadata$Sample != "",]
  sample_names <- substr(rownames(pc_df),11,nchar(rownames(pc_df))-(nchar(typeA)))
  # metadata <- metadata[metadata$SimplifiedSampleName %in% sample_names,]
  metadata <- metadata[!duplicated(metadata),]
  unique(metadata$Sample)#27
  rownames(metadata) <- metadata$Sample
  metadata <- metadata[sample_names,]
  table(metadata$Tissue)
  # metadata$Tissue <- str_split(metadata$Tissue,' ',simplify = T)[,2]
  # combine the pcs and the metadata
  pc_df <- cbind(pc_df, metadata)
  scalef <- 1
  pc_df$PC1 <- pc_df$PC1*scalef
  #
  table(pc_df$Tissue,pc_df$CellTypeSimplifiedSample)
  table(pc_df$Tissue)
  # pc_df$Tissue <- str_split(pc_df$Tissue,' ', simplify = T)[,1]
  pc_df$Tissue <- gsub('cervical cancer tissue','CC',pc_df$Tissue)
  pc_df$Tissue <- gsub('CA_HPV','CC',pc_df$Tissue)
  pc_df$Tissue <- gsub('neoplasm','CC',pc_df$Tissue)
  pc_df$Tissue <- gsub('high','HSIL_HPV',pc_df$Tissue)
  pc_df$Tissue <- gsub('metastatic','CC',pc_df$Tissue)
  # pc_df$Tissue[grepl("carcinoma",(pc_df$Tissue))] <- 'CC'
  color <-  c(brewer.pal(9, "Set1"),
              "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
              "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
              "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")
  table(pc_df$Tissue)

  # plot the first two pcs
  pc_df <- pc_df[,!grepl('tissue',colnames(pc_df))]
  p <- ggplot(pc_df, aes(x=PC1, y=PC2, color=Tissue)) +
    geom_point(size=2) + scale_color_manual(values=color) +
    theme_ArchR()
  p
  # ggsave(plot=p,height=5,width=5, filename=paste0("All_Cancer_Polyp_", typeA, "_vs_", background, "_", typeA, "_pca_on_diff_genes.pdf"), useDingbats=FALSE)

  table(colon_full$Tissue,colon_full$sample_name)
  scalef <- (-1)
  pc_df$PC1 <- pc_df$PC1*scalef
  pc_df$PC2 <- pc_df$PC2*scalef
  color <-  c(brewer.pal(9, "Set1"),
              "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
              "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
              "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")
  # combine the pcs and the metadata
  # plot the first two pcs
  #######-version1 -- y=PC1-#######
  # # plot the first two pcs(y=PC1)
  # table(pc_df$Tissue)
  # #pc_df$DiseaseState <- gsub('ADJ_','',pc_df$DiseaseState)
  # # plot the first two pcs
  # p <- ggplot(pc_df, aes(x=PC2, y=PC1, color=Tissue)) +
  #   geom_point(size=2) + scale_color_manual(values=color) +
  #   theme_ArchR()
  # p
  # # fit a spline (0.11 knot selected based on previous plot)
  # require(splines)
  # r <- summary(pc_df$PC2)
  # r
  # x <- seq(round(as.numeric(r[1]),2)-0.01,round(as.numeric(r[6]),2), by = 0.01)
  # spline_fit<-vector("list",length(x))
  # names(spline_fit) = x
  # A<-vector("list",length(x))
  # names(A) = x
  # for (i in 1:length(x)) {
  #   # set.seed(321)
  #   spline_fit[[i]] <- lm(PC1 ~ bs(PC2,knots = x[i]),data = pc_df )
  #   A[[i]] <-  AIC(spline_fit[[i]])
  #   # B <-  BIC(spline_fit[[i]])
  #   print(A[[i]]);
  #   #print(B)
  # }
  # unlist(A)
  # a <- which.min(unlist(A))
  # k <- as.numeric(names(a))
  # k
  # spline_fit <- lm(PC1 ~ bs(PC2,knots = k),data = pc_df )
  # fit<-lm(PC1 ~ bs(PC2,knots = 0.22),data = pc_df )
  # AIC(fit, spline_fit)
  # BIC(fit,spline_fit)
  # age.grid<-seq(from=(round(as.numeric(r[1]),2)-0.01)*10000*scalef, to = round(as.numeric(r[6]),2)*10000*scalef)/10000
  # splinefit <- data.frame(age.grid,predict(spline_fit,newdata = list(PC2=age.grid)))
  # colnames(splinefit) <- c("PC2", "PC1")
  # # plot the pcs with the spline fit
  # p <- ggplot(pc_df, aes(x=PC2, y=PC1, color=Tissue)) +
  #   geom_point(size=2) + scale_color_manual(values=color) +
  #   xlab(paste0("PC2 (",100*summary(pcs)$importance["Proportion of Variance","PC2"], "%)")) + ylab(paste0("PC1 (",100*summary(pcs)$importance["Proportion of Variance","PC1"], "%)")) + theme_ArchR()
  # p <- p+ geom_line(data=splinefit, colour="#CC0000") + theme(
  #   plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
  #   plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
  #   plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
  #   axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
  #   axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
  #   axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
  #   axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
  #   legend.title  = element_text(color = 'black', size  = 16),
  #   legend.text   = element_text(color = 'black', size   = 16),
  #   axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
  #   axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
  #   #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  # )
  # p
  # ggsave(plot=p,height=6,width=7, filename=paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_pca_on_diff_genes.pdf"), useDingbats=FALSE)
  # #get the nearest spline point for each sample
  # # get x values of nearest points--this will order from "least to most cancerous"
  # #get the nearest spline point for each sample
  # x_vals <- c()
  # for (i in 1:length(pc_df[,1])){
  #   x_vals <- c(x_vals,splinefit$PC2[which.min(sqrt((pc_df$PC2[i]-splinefit$PC2)**2+(pc_df$PC1[i]-splinefit$PC1)**2))])
  # }
  # pc_df$nearest_spline_x_vals <- x_vals
  # pc_df <- pc_df[order(pc_df$nearest_spline_x_vals),]
  # # save the results
  # saveRDS(pc_df, paste0("pc_df_", background, "_background_", typeA, "_cellType.rds"))

  #######-version2 -- y=PC2-#######
  # plot the first two pcs (x=PC1)
  p <- ggplot(pc_df, aes(x=PC1, y=PC2, color=Tissue)) +
    geom_point(size=2) + scale_color_manual(values=color) + theme_ArchR()
  p
  #ggsave(plot=p,height=5,width=5, filename=paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_pca_on_diff_genes.pdf"), useDingbats=FALSE)

  # fit a spline (0.11 knot selected based on previous plot)
  require(splines)
  r <- summary(pc_df$PC1)
  r
  x <- seq(round(as.numeric(r[1]),2)-0.02,round(as.numeric(r[6]),2)+0.02, by = 0.001)
  spline_fit<-vector("list",length(x))
  names(spline_fit) = x
  A<-vector("list",length(x))
  names(A) = x
  for (i in 1:length(x)) {
    # set.seed(321)
    spline_fit[[i]] <- lm(PC2 ~ bs(PC1,knots = x[i]),data = pc_df )
    A[[i]] <-  AIC(spline_fit[[i]])
    # B <-  BIC(spline_fit[[i]])
    print(A[[i]]);
    #print(B)
  }
  unlist(A)
  a <- which.min(unlist(A))
  k <- as.numeric(names(a))
  k
  #spline_fit <- lm(PC2 ~ ns(PC1,3),data = pc_df)#
  spline_fit <- lm(PC2 ~ bs(PC1,knots = k),data = pc_df )
  fit      <-   lm(PC2 ~ bs(PC1,knots = c(0.1)),data = pc_df )
  AIC(fit, spline_fit)
  BIC(fit,spline_fit)
  scalef <- 1
  age.grid<-seq(from=(round(as.numeric(r[1]),2)-0.02)*10000*scalef, to = (round(as.numeric(r[6]),2)+0.02)*10000*scalef)/10000
  splinefit <- data.frame(age.grid,predict(spline_fit,newdata = list(PC1=age.grid)))
  colnames(splinefit) <- c("PC1", "PC2")

  # plot the pcs with the spline fit
  p <- ggplot(pc_df, aes(x=PC1, y=PC2, color=Tissue)) +
    geom_point(size=2) + scale_color_manual(values=color) +
    xlab(paste0("PC1 (",100*summary(pcs)$importance["Proportion of Variance","PC1"], "%)")) + ylab(paste0("PC2 (",100*summary(pcs)$importance["Proportion of Variance","PC2"], "%)")) + theme_ArchR()
  p <- p+ geom_line(data=splinefit, colour="#CC0000") + theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
    legend.title  = element_text(color = 'black', size  = 16),
    legend.text   = element_text(color = 'black', size   = 16),
    axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  )
  p
  ggsave(plot=p,height=5,width=6, filename=paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_pca_on_diff_genes.pdf"), useDingbats=FALSE)
  #get the nearest spline point for each sample
  # get x values of nearest points--this will order from "least to most cancerous"
  x_vals <- c()
  for (i in 1:length(pc_df[,1])){
    x_vals <- c(x_vals,splinefit$PC1[which.min(sqrt((pc_df$PC1[i]-splinefit$PC1)**2+(pc_df$PC2[i]-splinefit$PC2)**2))])
  }
  pc_df$nearest_spline_x_vals <- x_vals
  pc_df <- pc_df[order(pc_df$nearest_spline_x_vals),]

  # save the results
  saveRDS(pc_df, paste0("pc_df_", background, "_background_", typeA, "_cellType.rds"))
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 10) Cluster genes for heatmap
table(rownames(colon_full_diff_test_avg_logFC_significant)==rownames(colon_full_diff_test_p_val_adj_significant))
DEG <- cbind(colon_full_diff_test_avg_logFC_significant,colon_full_diff_test_p_val_adj_significant)
setdiff(new_sample_list,sig_sample)#"LP16STM"
n <- new_sample_list#sig_sample
df_empty <- NULL#data.frame()
df_deg <- data.frame(sample = n, deg = rep(0,length(n)))
df <- lapply(setNames(n, n), function(nameindex) {
  col.logfc <- paste0("avg_log2FC",nameindex)
  col.pval <- paste0('p_val_adj',nameindex)
  y <- DEG[,c(col.logfc,col.pval)]
  y <- y[abs(y[,1])>=0.5 & y[,2]<=0.05,]
  up_gene <- rownames(y)[y[,1]>0] #rownames(y)[y[,1]>=0.75 & y[,2]<=0.05]
  down_gene <- rownames(y)[y[,1]<0]#rownames(y)[y[col.pval]>0]
  # batch_name <<- c(batch_name, rep(x, ncol(y)-1))
  x <- data.frame(sample=nameindex,up_gene=length(up_gene),
                  down_gene=length(down_gene))
  df_empty <<- rbind(df_empty,x)
  return(df_empty)
})
df_deg <- left_join(df_deg,df_empty, by='sample')
df_deg$deg <- apply(df_deg, 1, function(x){
  deg <- sum(as.numeric(x["up_gene"]),as.numeric(x["down_gene"]))
  return(deg)
})
df_deg$sample <- gsub("STM",'',df_deg$sample)
df_deg$Tissue <- pc_df$Tissue[match(df_deg$sample,pc_df$SimplifiedSampleName)]
#colon_full$Tissue[match(df_deg$sample,colon_full$SimplifiedSampleName)]
dim(colon_full_diff_test_avg_logFC_significant)

Latt <- colon_full_diff_test_avg_logFC_significant
Latt$gene <- rownames(Latt)
Latt <- melt(Latt,id.vars = 'gene')

Patt <- colon_full_diff_test_p_val_adj_significant
Patt$gene <- rownames(Patt)
Patt <- melt(Patt,id.vars = 'gene')
colnames(Patt) <- c('Symbol','Replicates','Pvalue_Adj')
Patt$Replicates <- gsub('p_val_adj','avg_log2FC',Patt$Replicates)
Patt$Cell <- Patt$Replicates
table(Patt$Cell == Latt$variable) 
Patt$Log2FC <- Latt$value
table(Patt$Cell %in% rownames(pc_df))
Patt$Replicates <- pc_df$SimplifiedSampleName[match(Patt$Cell,rownames(pc_df))]
Patt$Species <- 'Human'
Patt$Organ <- 'Cervix'
Patt$Tissue <- pc_df$Tissue[match(Patt$Cell,rownames(pc_df))]
Patt$Malignancy <- pc_df$nearest_spline_x_vals[match(Patt$Cell,rownames(pc_df))]
Patt <- Patt[Patt$Pvalue_Adj<0.05,]
length(unique(Patt$Symbol))#2595
save(Patt, colon_full_diff_test_avg_logFC_significant,
    colon_full_diff_test_p_val_adj_significant,pc_df, file = 'Cervix_DEG.RData')
     
# -----Draw dot plot------
# View(pc_df[,c("PC1",'PC2','nearest_spline_x_vals',"Tissue")])
test.diff <- colon_full_diff_test_avg_logFC_significant[,colnames(colon_full_diff_test_avg_logFC_significant) %in% rownames(pc_df)]
dim(test.diff)
colnames(test.diff)
cohort_directory <- '/public/home/lorihan/lrh/database/Cervix/SCP/epithelial_results'
setwd(cohort_directory)
# Create directory
deg_folder <- 'deg_figures'
if (!dir.exists(paste0(deg_folder))){
  dir.create(paste0(deg_folder))
}
setwd(deg_folder)
getwd()
dim(test.diff)#2596   19
n <- rownames(test.diff)
p <- NULL
plot <-  lapply(setNames(n, n), function(nameindex) {
  gene <- nameindex
  gene_p <- pc_df[match(colnames(test.diff),rownames(pc_df)),]
  gene_p$gene <- t(test.diff[gene,])
  ylab <- 'Log2FC'
  p = ggplot(gene_p, aes(x=nearest_spline_x_vals, y=gene, color=Tissue)) +
    geom_point(size=2) + scale_color_manual(values=color)+ggtitle(gene) + 
    xlab("Malignancy Continuum")+ylab(ylab)+theme(
      plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
      plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
      plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
      axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
      axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
      axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
      legend.title  = element_text(color = 'black', size  = 16),
      legend.text   = element_text(color = 'black', size   = 16),
      axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
      axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
      #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
    )
  f <- paste0(gene,"_pca_on_diff_genes.png")
  # 0.2 保存为png等位图格式
  ggsave(plot=p,height=4,width=6, filename=f, dpi = 300, device = "png")
  return(p)
})

if (10 %in% execute_steps){
  # colon_full_diff_test <- readRDS("./colon_full_diff_test_UnaffectedSTM_background.rds")#"./colon_full_diff_test.rds")
  # background <- "Unaffected"
  # typeA <- "STM"
  #
  # colon_full_diff_test_avg_logFC <- colon_full_diff_test[,grepl("avg_log2FC", colnames(colon_full_diff_test))]
  # colon_full_diff_test_avg_logFC <- colon_full_diff_test_avg_logFC/0.69314718056
  # colon_full_diff_test_p_val_adj <- colon_full_diff_test[,grepl("p_val_adj", colnames(colon_full_diff_test))]
  # colon_full_diff_test_avg_logFC[is.na(colon_full_diff_test_avg_logFC)] <- 0
  # colon_full_diff_test_p_val_adj[is.na(colon_full_diff_test_p_val_adj)] <- 1
  # # 1/28/2020
  colon_full_diff_test_p_val_adj_significant <- colon_full_diff_test_p_val_adj_significant[,match(new_sample_list,sig_sample)]
  order_samples <- rownames(pc_df)[rownames(pc_df) %in% colnames(colon_full_diff_test_avg_logFC_significant)]
  
  #载入绘图包
  library(iheatmapr)
  library(datasets)
  library(reshape2)
  #devtools::install_github("ropensci/iheatmapr")
  test.diff <- as.matrix(colon_full_diff_test_avg_logFC_significant[,
                                                                    match(order_samples,colnames(colon_full_diff_test_avg_logFC_significant))])
  set.seed(1)
  gene_clusters <- kmeans(cbind(test.diff), 10, iter.max = 500, algorithm = "Lloyd")
}
clust_assign <- gene_clusters$cluster#kmeans(Indometh_matrix, 3)$cluster
main_heatmap(test.diff) %>%
  add_row_clusters(clust_assign) %>%
  add_col_clusters(clust_assign)

# # -----Draw dot plot------
test.diff <- colon_full_diff_test_avg_logFC_significant[,colnames(colon_full_diff_test_avg_logFC_significant) %in% rownames(pc_df)]
colnames(test.diff)
gene <- 'S100A2'
gene_p <- pc_df[match(colnames(test.diff),rownames(pc_df)),]
gene_p$gene <- t(test.diff[gene,])
p = ggplot(gene_p, aes(x=nearest_spline_x_vals, y=gene, color=Tissue)) +
  geom_point(size=2) + scale_color_manual(values=color)+
  xlab("Malignant continuum")+ylab(gene)+theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
    legend.title  = element_text(color = 'black', size  = 16),
    legend.text   = element_text(color = 'black', size   = 16),
    axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  )
p
ggsave(plot=p,height=5,width=5.5, filename=paste0(gene,"_pca_on_diff_genes.pdf"), useDingbats=FALSE)

gene = 'GNAS'
gene_p <- pc_df[match(colnames(test.diff),rownames(pc_df)),]
gene_p$gene <- t(test.diff[gene,])
p = ggplot(gene_p, aes(x=nearest_spline_x_vals, y=gene, color=Tissue)) +
  geom_point(size=2) + scale_color_manual(values=color)+
  xlab("Malignant continuum")+ylab(gene)+theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
    legend.title  = element_text(color = 'black', size  = 16),
    legend.text   = element_text(color = 'black', size   = 16),
    axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  )
p
ggsave(plot=p,height=5,width=5.5, filename=paste0(gene,"_pca_on_diff_genes.pdf"), useDingbats=FALSE)

# gene = 'GPX2'
# gene_p <- pc_df[match(colnames(test.diff),rownames(pc_df)),]
# gene_p$gene <- t(test.diff[gene,])
# p = ggplot(gene_p, aes(x=nearest_spline_x_vals, y=gene, color=Tissue)) +
#   geom_point(size=2) + scale_color_manual(values=color)+
#   xlab("Malignant continuum")+ylab(gene)+theme(
#     plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
#     plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
#     plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#     axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
#     axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
#     axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
#     axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
#     legend.title  = element_text(color = 'black', size  = 16),
#     legend.text   = element_text(color = 'black', size   = 16),
#     axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#     axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#     #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
#   )
# p
# ggsave(plot=p,height=5,width=5.5, filename=paste0(gene,"_pca_on_diff_genes.pdf"), useDingbats=FALSE)
# 
# gene = 'IGLC3'
# gene_p <- pc_df[match(colnames(test.diff),rownames(pc_df)),]
# gene_p$gene <- t(test.diff[gene,])
# p = ggplot(gene_p, aes(x=nearest_spline_x_vals, y=gene, color=Tissue)) +
#   geom_point(size=2) + scale_color_manual(values=color)+
#   xlab("Malignant continuum")+ylab(gene)+theme(
#     plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
#     plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
#     plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#     axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
#     axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
#     axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
#     axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
#     legend.title  = element_text(color = 'black', size  = 16),
#     legend.text   = element_text(color = 'black', size   = 16),
#     axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#     axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#     #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
#   )
# p
# ggsave(plot=p,height=5,width=5.5, filename=paste0(gene,"_pca_on_diff_genes.pdf"), useDingbats=FALSE)


###############################################################################################################################
# --------snRNA_immune_analysis.R--------
# save metadata table:
setwd(cohort_directory)
setwd('initial_clustering/')
Guo_i <- readRDS('/public/home/lorihan/lrh/database/Cervix/Guo/initial_clustering/Guo_immune_initial.rds')
table(Guo_i$batch)
Hua1_i <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua1/initial_clustering/Hua1_immune_initial.rds')
table(Hua1_i$batch)
table(Hua1_i$Tissue)
Hua2_i <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua2/initial_clustering/Hua2_immune_initial.rds')
table(Hua2_i$batch)
Hua3_i <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua3/initial_clustering/Hua3_immune_initial.rds')
table(Hua3_i$batch)
Hua4_i <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua4/initial_clustering/Hua4_immune_initial.rds')
table(Hua4_i$batch)
SCP_i <- readRDS('/public/home/lorihan/lrh/database/Cervix/SCP/SCP_immune_initial.rds')
table(SCP_i$batch)
table(SCP_i$Tissue)
# DefaultAssay(Hua1_i) <- 'RNA'
table(All_i$Celltype)
table(SCP_i$cell_type__ontology_label)
# SCP_i$Celltype <- gsub(',','',SCP_i$cell_type__ontology_label)
# sed -i 's/CD8-positive, alpha-beta T cell/CD8-positive alpha-beta T cell/g' metadata.csv
seurat_obj <- merge(All_i,c(Hua1_i,Hua2_i,Hua3_i,Hua4_i,Guo_i,SCP_i))
seurat_obj
colnames(seurat_obj@meta.data)
table(seurat_obj$batch)
table(seurat_obj$orig.ident[is.na(seurat_obj$batch)])
setwd("/public/home/lorihan/lrh/database/Cervix/SCP/")
if (!dir.exists(paste0('STAR-Methods'))){
  dir.create(paste0('STAR-Methods'))
}
dir.create(paste0('STAR-Methods/Imm'))
setwd('STAR-Methods/Imm/')
colnames(seurat_obj@meta.data)
names <- c( "orig.ident","nCount_RNA", "nFeature_RNA", 'batch','Tissue','sample_name','Celltype')
seurat_obj@meta.data <- seurat_obj@meta.data[, names]
seurat_obj$barcode <- colnames(seurat_obj)
# seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
# seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('./', 'counts.mtx'))
# write dimesnionality reduction matrix, in this example case pca matrix
# write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
# cell_embeddings<-Embeddings(seurat_obj, reduction = "umap")#使用Embeddings功能提取seurat UMAP或者TSNE
# write.csv(cell_embeddings, file = "cell_embeddings.csv")

# --------scRNA_stroma_analysis.R--------
# save metadata table:
setwd(cohort_directory)
Guo_s <- readRDS('/public/home/lorihan/lrh/database/Cervix/Guo/initial_clustering/Guo_stromal_initial.rds')
table(Guo_s$batch)
Hua1_s <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua1/initial_clustering/Hua1_stromal_initial.rds')
table(Hua1_s$batch)
table(Hua1_s$Tissue)
Hua2_s <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua2/initial_clustering/Hua2_stromal_initial.rds')
table(Hua2_s$batch)
Hua3_s <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua3/initial_clustering/Hua3_stromal_initial.rds')
table(Hua3_s$batch)
Hua4_s <- readRDS('/public/home/lorihan/lrh/database/Cervix/Hua4/initial_clustering/Hua4_stromal_initial.rds')
table(Hua4_s$batch)
SCP_s <- readRDS('/public/home/lorihan/lrh/database/Cervix/SCP/SCP_stromal_initial.rds')
table(SCP_s$batch)
table(SCP_s$Tissue)
table(All_s$Celltype)
table(SCP_s$cell_type__ontology_label)
SCP_s$Celltype <- SCP_s$cell_type__ontology_label
# DefaultAssay(Hua1_s) <- 'RNA'
# seurat_obj <- merge(All_s, y=c(Meng_s,Gao_s))
seurat_obj <- merge(All_s,c(Hua1_s,Hua2_s,Hua3_s,Hua4_s,Guo_s,SCP_s))
seurat_obj
colnames(seurat_obj@meta.data)
table(seurat_obj$batch)
table(seurat_obj$orig.ident[is.na(seurat_obj$batch)])

setwd("/public/home/lorihan/lrh/database/Cervix/SCP/")
if (!dir.exists(paste0('STAR-Methods'))){
  dir.create(paste0('STAR-Methods'))
}
dir.create(paste0('STAR-Methods/Str'))
setwd('STAR-Methods/Str/')
colnames(seurat_obj@meta.data)
names <- c( "orig.ident","nCount_RNA", "nFeature_RNA", 'batch','Tissue','sample_name','Celltype')
seurat_obj@meta.data <- seurat_obj@meta.data[, names]
seurat_obj$barcode <- colnames(seurat_obj)
# seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
# seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('./', 'counts.mtx'))
# write dimesnionality reduction matrix, in this example case pca matrix
# write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
# cell_embeddings<-Embeddings(seurat_obj, reduction = "umap")#使用Embeddings功能提取seurat UMAP或者TSNE
# write.csv(cell_embeddings, file = "cell_embeddings.csv")
# --------
setwd(cohort_directory)
setwd(analysis_parent_folder)
library(Seurat)
library(SeuratDisk)

# step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
seu = DietSeurat(
  colon_full,
  counts = TRUE, # so, raw counts save to adata.raw.X
  data = TRUE, # so, log1p counts save to adata.X
  scale.data = FALSE, # set to false, or else will save to adata.X
  features = rownames(colon_full), # export all genes, not just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)
# step 2: factor to character, or else your factor will be number in adata
i <- sapply(seu@meta.data, is.factor)
seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)

# step 3: convert
SaveH5Seurat(seu, filename = "srt.h5seurat", overwrite = TRUE)
Convert("srt.h5seurat", "cervix.epi.h5ad", assay="RNA", overwrite = TRUE)

