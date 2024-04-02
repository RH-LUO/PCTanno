#!/usr/local/bin/Rscript
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library('e1071') #必须加载，因为后面反卷积的算法是基于这个包里的数据

#####------Cibersort------#####
data_path <- '/data/rluo4/All/Output/Cluster/'
organ_all <- list.files(data_path)
indir <- '/data/rluo4/All/Output/TCGA/'
# for (organ in organ_all) {
#   res_file = paste0(indir,'Malignant-Transformation-',organ, '.rds')
#   if(organ == 'Chen'){
#     res_file = gsub('Chen/','CRC/',res_file)
#   }
#   load(file = res_file)
#   seu <- rdata_filter
#   # table(Idents(seu))
#   table(seu$Tissue)
#   seu$cell_class <- paste0("C",seu$leiden_res2)
#   if(organ =='Pancreas'){
#     seu$cell_class <- paste0("C",seu$leiden_res1.5)
#   }
#   outdir <- paste0(indir,organ)
#   outfile <- paste0(outdir,"/Output/Epi/rdata_subset_class.txt")
#   if ( file.exists(outfile) ) {
#     print(paste0('CIBERSORT for ', organ,  " is done, so skip !"))
#     next;
#   }
#   group.by <- "cell_class" 
#   mat <- as.data.frame(t(as.matrix(GetAssayData(seu, assay = "RNA", slot = "data")))) 
#   mat <- aggregate(mat, by=list(seu@meta.data[[group.by]]), FUN="mean") 
#   # -----
#   rownames(mat) <- mat$Group.1 
#   mat <- t(mat[,-1]) 
#   head(mat)
#   
#   if(organ == 'Chen'){
#     outfile = gsub('Chen/','CRC/',outfile)
#     outfile = gsub('rdata_subset','rdata_subset_Chen',outfile)
#   }
#   write.table(mat, outfile, sep = "\t",col.names = T,row.names = T)     #保存为sig.txt
#   length(unique(rownames(seu$RNA)))
#   dim(seu$RNA)
#   ex <- mat
#   qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#   LogC <- (qx[5] > 100) ||
#     (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#     (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
#   
#   if (LogC) { ex[which(ex <= 0)] <- NaN
#   exprSet <- log2(ex)
#   print("log2 transform finished")}else{print("log2 transform not needed")}
#   
# }


########------bulk RNAseq------#######
# Changed Chen into CRC (Becker)
cohort_directory = '/data/rluo4/All/Output/TCGA/'
setwd(cohort_directory)
load('TCGA-Mutation.RData')
files <- list.files()
files <- files[grep('result.Rdata', files)]
cohort_all <- files
cohort_all <- str_split(cohort_all, '_', simplify = T)[,1]
cohort_all <- cohort_all[! cohort_all %in% 'GBM']
my_ct <- c('BRCA','CESC','COAD','ESCA','HNSC','LIHC','LUAD','LUSC','PAAD','PRAD','READ','STAD','THCA','UCEC')
# my_tt <- c('Breast','Cervix','Chen','Esophagus','HNSCC','Liver','Lung','Lung','Pancreas','Prostate','Chen','GC','THCA','Endometrium')
my_tt <- c('Breast','Cervix','CRC','Esophagus','HNSCC','Liver','Lung','Lung','Pancreas','Prostate','CRC','GC','THCA','Endometrium')
index <- data.frame(TCGA = my_ct, PCT = my_tt)
library(data.table)
library(survival)
library(survminer)

# for (cohort in my_ct) { # already done in ts860
#   # cohort <- 'LUAD' # rsync -a -P -r -e 'ssh -p 42769' lorihan@ts860.hxdsjzx.com:/home/lorihan/lrh/All/Output/TCGA/cibersort/COAD ./
#   print(cohort)
#   f <- files[grep(cohort, files)]
#   load(paste0(cohort_directory,f))
#   outdir <- paste0(cohort_directory,"cibersort/",cohort)
#   if (!dir.exists(paste0(outdir))){
#     dir.create(paste0(outdir))
#   }
#   setwd(outdir)
#   outfile <- paste0(outdir,'/logtpm.txt')
#   dim(logtpm)#18767   423
#   exp2 = as.data.frame(logtpm[match(unique(rownames(logtpm)),
#                                     rownames(logtpm)),])
#   exp2 = rownames_to_column(exp2)
#   write.table(exp2,file = outfile, row.names = F,quote = F,sep = "\t")
#   
#   outfile <- paste0(outdir,'/tpms.txt')
#   dim(tpms)#18767   474
#   tpms=data.frame(tpms[match(unique(rownames(tpms)),
#                              rownames(tpms)),])
#   dim(tpms)#18268   474
#   exp2 = as.data.frame(tpms[match(unique(rownames(tpms)),
#                                   rownames(tpms)),])
#   exp2 = rownames_to_column(exp2)
#   write.table(exp2, file = outfile, row.names = F,quote = F,sep = "\t")
# }

source("/home/rluo4/Rcode/HPC/Rscripts_test/CIBERSORTx.R")  #激活function
for (cohort in my_ct[c(3,11)]) {
  print(cohort)
  organ = index$PCT[index$TCGA==cohort]
  # outdir <- paste0('/data/rluo4/database/',organ,'/Output/')
  # if(organ == 'Chen'){
  #   outdir = gsub('Chen/','CRC/',outdir)
  # }
  # setwd(outdir)
  
  indir <- paste0(cohort_directory,"cibersort/",cohort)
  if (!dir.exists(paste0(indir))){
    dir.create(paste0(indir))
  }
  setwd(indir)
  if ( file.exists(paste0(indir,"cibersort.logtpm.Rdata"))) {
    print(paste0('CIBERSORT for ', organ,  " is done, so skip !"))
    next;
  }
  # cd /data/rluo4/All/Output/TCGA/cibersort/READ
  # rsync -a -P -r -e 'ssh -p 42769' lorihan@ts860.hxdsjzx.com:/home/lorihan/lrh/database/CRC/Output/Epi/rdata_subset* ./
  outfile = paste0(indir, "/rdata_subset_class.txt")
  if(organ == 'Chen'){
    outfile = gsub('Chen/','CRC/',outfile)
    outfile = gsub('rdata_subset','rdata_subset_Chen',outfile)
  }
  if(T){
    TME.results = CIBERSORT(outfile,
                            paste0(indir, "/logtpm.txt"),
                            perm = 1000,
                            QN = FALSE)
    save(TME.results,file = paste0(indir,"/cibersort.logtpm.Rdata"))
  }
  
  #   setwd(outdir)
  #   if(T){
  #     TME.results = CIBERSORT(outfile,
  #                             paste0(indir, "/tpms.txt"),
  #                             perm = 1000,
  #                             QN = FALSE)
  #     save(TME.results,file = paste0(indir,"cibersort.tpms.Rdata"))
}