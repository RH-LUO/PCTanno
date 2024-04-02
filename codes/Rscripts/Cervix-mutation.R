#!/usr/local/bin/Rscript
# title: "Single-cell analysis of gastric pre-cancerous and cancer lesions reveals cell lineage diversity and intratumoral heterogeneity"
# author: "Ruihan Luo"
# date: "Aug 8th,2023"
# cn12
# rm(list=ls())
library(maftools, lib.loc = "/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1") ## 
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
load('~/lrh/database/All/organ12_Epi.RData')
cohort_directory <- '/public/home/lorihan/lrh/database/Cervix/Output/Epi'
setwd(cohort_directory)
Patt <- read.table('~/lrh/database/All/organ12-deg.txt',sep = '\t',fill = T,header = TRUE)
organ <- unique(Patt$Organ)
DEG <- Patt$Symbol[Patt$Organ==organ[1]]
for (i in organ[-1]) {
  j <- Patt$Symbol[Patt$Organ==i]
  DEG <- intersect(DEG,j) 
}
unique(Cervix_Epi$sample_name)
table(Cervix_Epi$Tissue,Cervix_Epi$sample_name)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
##################################################################################
# CHROM  Start   End     REF     ALT     FILTER  Cell_types      Up_context      Down_context    N_ALT   
# Dp      Nc   -Bc      Cc      VAF     CCF     BCp     CCp     Cell_types_min_BC       Cell_types_min_CC       Rest_BC Rest_CC Fisher_p      Cell_type_Filter        INFO  
NS <- unique(Cervix_Epi$orig.ident[Cervix_Epi$Tissue %in% c('AdultCervix','normal','NO_HPV')])
NS
Epi <- Cervix_Epi[Cervix_Epi$sample_name=='Guo',]
table(Epi$Tissue,Epi$orig.ident)
setwd("/public/home/lorihan/lrh/database/Cervix/Guo/anno.var/")
path = '/public/home/lorihan/lrh/database/Cervix/Guo/anno.var'
solution <- list.files(path,pattern='variants.avinput')
Nor <- paste(NS, 'variants.avinput', sep = '.')
solution <- solution[! solution %in% Nor]
solution
name <- str_split(solution,'[.]',simplify = T)[,1]
filePath <- sapply(solution, function(x){ 
  paste(path,x,sep='/')})  
solutions <- lapply(filePath, function(x){
  data.table::fread(file = x, stringsAsFactors = FALSE,sep = '\t')
})
names(solutions)
solutions.all <- solutions
# annovar_ops = list.files(pattern = "ann2.txt$") #all your annovar outputs per sample
solution <- list.files(path,pattern='hg38_multianno.txt')
Nor <- paste(NS, 'hg38_multianno.txt', sep = '.')
solution <- solution[! solution %in% Nor]
solution
name <- str_split(solution,'[.]',simplify = T)[,1]
filePath <- sapply(solution, function(x){ 
  paste(path,x,sep='/')})  
annovar_mafs = lapply(filePath, annovarToMaf) #convert to MAFs using annovarToMaf
names(annovar_mafs)
n <- names(annovar_mafs)
n
annovar_mafs <-  lapply(setNames(n, n), function(nameindex) {
  y <- as.data.frame(annovar_mafs[[nameindex]])
  z <- gsub('.hg38_multianno.txt','.variants.avinput',nameindex)
  w <- as.data.frame(solutions.all[[z]])
  VAF <- str_split(w$V6, '[-]', simplify = T)[,15]
  CCF <- str_split(w$V6, '[-]', simplify = T)[,16]
  y[,'VAF'] <- VAF
  y[,'CCF'] <- CCF
  return(y)  
})

dir = '/public/home/lorihan/lrh/database/Cervix/Guo/scMapping'
scMapping <- list.files(dir)
scMapping <- scMapping[! scMapping %in% NS]
Guo_anno <- vector("list",length(scMapping))
names(Guo_anno) <- scMapping
for (i in scMapping) {
  path <- paste(dir,i,sep='/')
  solution <- list.files(path,pattern='single_cell_genotype.tsv')
  # solution <- solution[!grepl('abn',solution)]
  name <- str_split(solution,'[.]',simplify = T)[,1]
  filePath <- sapply(solution, function(x){ 
    paste(path,x,sep='/')})  
  solutions <- lapply(filePath, function(x){
    data.table::fread(file = x, stringsAsFactors = FALSE)
  })
  names(solutions)
  n <- names(solutions)# sample='MSI_1-2020'
  solutions <-  lapply(setNames(n, n), function(nameindex) {
    a <- paste(i, '.hg38_multianno.txt', sep = '')
    b <- as.data.frame(annovar_mafs[[a]])
    c <- as.data.frame(solutions[[nameindex]])
    colnames(c)[1:2] <- colnames(b)[2:3]
    c$Start_Position <- as.character(c$Start_Position)
    intersect(colnames(b),colnames(c))
    y <- left_join(b, c[,-1],by='Start_Position')
    print(unique(y$Tumor_Sample_Barcode))
    print(unique(Epi$orig.ident))
    extract_last_element <- function(x) {
      split_string <- strsplit(x, "_")
      last_element <- sapply(split_string, function(y) tail(y, n = 1))
      return(last_element)
    }
    Epi$CB <- extract_last_element(str_split(rownames(Epi),'-1',simplify = T)[,1])
    print(table(y$CB %in% Epi$CB))
    Epi$CB <- paste(Epi$CB, Epi$orig.ident,sep = '-')#rownames(Epi)
    y$CB <- paste(y$CB, y$Tumor_Sample_Barcode,sep = '-')
    print(table(y$CB %in% Epi$CB))
    Epi$CB <- paste(Epi$CB, Epi$Cell_Type,sep = '-')#rownames(Epi)
    y$CB <- paste(y$CB, y$Cell_type_observed,sep = '-')
    print(table(y$CB %in% Epi$CB))
    # View(y[y$CB %ni% Epi$CB,])
    y <- y[y$CB %in% Epi$CB,]
    # filename <- sapply(rownames(Epi),shorten_names)
    return(y)
  })
  Guo_anno[[i]] = data.table::rbindlist(l = solutions, fill = TRUE) #Merge into single MAF
  }
scRNA.annovar <- NULL
allDat <-  lapply(Guo_anno,function(y){
  scRNA.annovar <<- rbind(scRNA.annovar, y)
  return(scRNA.annovar)
})
remove(allDat)
dim(scRNA.annovar)
colnames(scRNA.annovar)
unique(scRNA.annovar$Tumor_Sample_Barcode)#5

grep("FILTER",colnames(scRNA.annovar))
unique(scRNA.annovar$Tumor_Sample_Barcode)#62

table(na.omit(scRNA.annovar$X1000g2015aug_all)>0.05)#3-->1
filter.1000g <- which(scRNA.annovar$X1000g2015aug_all > 0.05)

filter.ExAC <- which(scRNA.annovar$ExAC_ALL > 0.05)
filter.gnomeAD <- which(scRNA.annovar$gnomAD_genome_ALL > 0.05)
filter.maf <- unique(c(filter.1000g,filter.ExAC,filter.gnomeAD))
# filter.maf <- unique(c(filter.1000g,filter.ExAC))
length(unique(filter.maf))#-->87
table(!is.na(scRNA.annovar$genomicSuperDups))#T:6552 
table(scRNA.annovar$rmsk)
table(!is.na(scRNA.annovar$rmsk))#33198
dim(scRNA.annovar)#761596     31-->40
scRNA.annovar <- scRNA.annovar[-filter.maf,]
dim(scRNA.annovar)#710112     33
# scRNA.annovar$Score <- round(scRNA.annovar$t_alt_count/scRNA.annovar$t_depth,3)
table(scRNA.annovar$VAF>0.02)
table(scRNA.annovar$VAF == '.')
table(scRNA.annovar$FATHMM_score>0.7)
# scRNA.annovar<-scRNA.annovar[scRNA.annovar$Score>=0.02,]#2%
maf.sc <- data.frame(scRNA.annovar[scRNA.annovar$VAF!='.', ])
# maf.sc <-  maf.sc[maf.sc$Score>=0.02,]#0.4%
unique(maf.sc$Tumor_Sample_Barcode)#35
table(maf.sc$FATHMM_score>0.7)
table(is.na(maf.sc$FATHMM_score))
# maf.sc <- maf.sc[maf.sc$FATHMM_score>0.7,]
dim(scRNA.annovar)#710112     33
dim(maf.sc)#2469894     33
table(scRNA.annovar$t_depth>2)
table(!is.na(maf.sc$genomicSuperDups))#21340
table(!is.na(maf.sc$rmsk))#198870
table(maf.sc$Variant_Type)
table(maf.sc$Variant_Type !="SNP")
table(maf.sc[!is.na(maf.sc$rmsk),]$Variant_Type)
table(maf.sc$Score>=0.004)
table(maf.sc[!is.na(maf.sc$genomicSuperDups),]$Score>0.04)
table(maf.sc[!is.na(maf.sc$rmsk) & maf.sc$Variant_Type !="SNP",]$Score<0.1)
table(maf.sc[!is.na(maf.sc$genomicSuperDups) & maf.sc$Variant_Type !="SNP",]$Score<0.1)
#F:458  T:43 
####
filter.rmsk <- which(!is.na(maf.sc$rmsk))
filter.dups <- which(!is.na(maf.sc$genomicSuperDups) )
filter.repeat <- which(!is.na(maf.sc$rmsk) & !is.na(maf.sc$genomicSuperDups))#union(filter.rmsk,filter.dups)
dim(maf.sc)
dim(maf.sc[-filter.repeat,])
################################################################################
maf.sc <- maf.sc[-filter.repeat,]
unique(maf.sc$Tumor_Sample_Barcode)#13
dim(maf.sc)#2469894     140
unique(maf.sc$Tumor_Sample_Barcode)
maf.Guo <- maf.sc

################################################################################
setwd("/public/home/lorihan/lrh/database/Cervix/Hua1/anno.var/")
path = '/public/home/lorihan/lrh/database/Cervix/Hua1/anno.var'
solution <- list.files(path,pattern='variants.avinput')
solution <- solution[!grepl('normal',solution)]
solution
name <- str_split(solution,'[.]',simplify = T)[,1]
filePath <- sapply(solution, function(x){ 
  paste(path,x,sep='/')})  
solutions <- lapply(filePath, function(x){
  data.table::fread(file = x, stringsAsFactors = FALSE)
})
names(solutions)
solutions.all <- solutions
# annovar_ops = list.files(pattern = "ann2.txt$") #all your annovar outputs per sample
solution <- list.files(path,pattern='hg38_multianno.txt')
solution <- solution[!grepl('normal',solution)]
name <- str_split(solution,'[.]',simplify = T)[,1]
filePath <- sapply(solution, function(x){ 
  paste(path,x,sep='/')})  
annovar_mafs = lapply(filePath, annovarToMaf) #convert to MAFs using annovarToMaf

names(annovar_mafs)
n <- names(annovar_mafs)
n
annovar_mafs <-  lapply(setNames(n, n), function(nameindex) {
  y <- as.data.frame(annovar_mafs[[nameindex]])
  z <- gsub('.hg38_multianno.txt','.variants.avinput',nameindex)
  w <- as.data.frame(solutions.all[[z]])
  VAF <- str_split(w$V6, '[-]', simplify = T)[,15]
  CCF <- str_split(w$V6, '[-]', simplify = T)[,16]
  y[,'VAF'] <- VAF
  y[,'CCF'] <- CCF
  return(y)  
})

dir = '/public/home/lorihan/lrh/database/Cervix/Hua1/scMapping'
scMapping <- list.files(dir)
scMapping <- scMapping[!grepl('normal',scMapping)]
Hua1_anno <- vector("list",length(scMapping))
names(Hua1_anno) <- scMapping

for (i in scMapping) {
  path <- paste(dir,i,sep='/')
  solution <- list.files(path,pattern='single_cell_genotype.tsv')
  # solution <- solution[!grepl('abn',solution)]
  name <- str_split(solution,'[.]',simplify = T)[,1]
  filePath <- sapply(solution, function(x){ 
    paste(path,x,sep='/')})  
  solutions <- lapply(filePath, function(x){
    data.table::fread(file = x, stringsAsFactors = FALSE)
  })
  names(solutions)
  n <- names(solutions)# sample='MSI_1-2020'
  solutions <-  lapply(setNames(n, n), function(nameindex) {
    a <- paste(i, '.hg38_multianno.txt', sep = '')
    b <- as.data.frame(annovar_mafs[[a]])
    c <- as.data.frame(solutions[[nameindex]])
    colnames(c)[1:2] <- colnames(b)[2:3]
    c$Start_Position <- as.character(c$Start_Position)
    intersect(colnames(b),colnames(c))
    y <- left_join(b, c[,-1],by='Start_Position')
    Epi <- Cervix_Epi[Cervix_Epi$sample_name=='Hua1',]
    rownames(Epi) <- gsub('Tumor','cervical',rownames(Epi))
    Epi$orig.ident <- gsub('Tumor','cervical',Epi$orig.ident)
    
    Epi$CB <- gsub("-1",'',str_split(rownames(Epi),'_',simplify = T)[,2])
    print(table(y$CB %in% Epi$CB))
    Epi$CB <- paste(Epi$CB, Epi$orig.ident,sep = '-')#rownames(Epi)
    y$CB <- paste(y$CB, y$Tumor_Sample_Barcode,sep = '-')
    print(table(y$CB %in% Epi$CB))
    Epi$CB <- paste(Epi$CB, Epi$Cell_Type,sep = '-')#rownames(Epi)
    y$CB <- paste(y$CB, y$Cell_type_observed,sep = '-')
    print(table(y$CB %in% Epi$CB))
    # View(y[y$CB %ni% Epi$CB,])
    y <- y[y$CB %in% Epi$CB,]
    # filename <- sapply(rownames(Epi),shorten_names)
    return(y)
  })
  Hua1_anno[[i]] = data.table::rbindlist(l = solutions, fill = TRUE) #Merge into single MAF
}
scRNA.annovar <- NULL
allDat <-  lapply(Hua1_anno,function(y){
  scRNA.annovar <<- rbind(scRNA.annovar, y)
  return(scRNA.annovar)
})
remove(allDat)
dim(scRNA.annovar)
colnames(scRNA.annovar)
unique(scRNA.annovar$Tumor_Sample_Barcode)#5

grep("FILTER",colnames(scRNA.annovar))
unique(scRNA.annovar$Tumor_Sample_Barcode)#62

table(na.omit(scRNA.annovar$X1000g2015aug_all)>0.05)#3-->1
filter.1000g <- which(scRNA.annovar$X1000g2015aug_all > 0.05)

filter.ExAC <- which(scRNA.annovar$ExAC_ALL > 0.05)
filter.gnomeAD <- which(scRNA.annovar$gnomAD_genome_ALL > 0.05)
filter.maf <- unique(c(filter.1000g,filter.ExAC,filter.gnomeAD))
# filter.maf <- unique(c(filter.1000g,filter.ExAC))
length(unique(filter.maf))#-->87
table(!is.na(scRNA.annovar$genomicSuperDups))#T:6552 
table(scRNA.annovar$rmsk)
table(!is.na(scRNA.annovar$rmsk))#33198
dim(scRNA.annovar)#761596     31-->40
scRNA.annovar <- scRNA.annovar[-filter.maf,]
dim(scRNA.annovar)#710112     33
# scRNA.annovar$Score <- round(scRNA.annovar$t_alt_count/scRNA.annovar$t_depth,3)
table(scRNA.annovar$VAF>0.02)
table(scRNA.annovar$VAF == '.')
table(scRNA.annovar$FATHMM_score>0.7)
# scRNA.annovar<-scRNA.annovar[scRNA.annovar$Score>=0.02,]#2%
maf.sc <- data.frame(scRNA.annovar[scRNA.annovar$VAF!='.', ])
# maf.sc <-  maf.sc[maf.sc$Score>=0.02,]#0.4%
unique(maf.sc$Tumor_Sample_Barcode)#35
table(maf.sc$FATHMM_score>0.7)
table(is.na(maf.sc$FATHMM_score))
# maf.sc <- maf.sc[maf.sc$FATHMM_score>0.7,]
dim(scRNA.annovar)#710112     33
dim(maf.sc)#2469894     33
table(scRNA.annovar$t_depth>2)
table(!is.na(maf.sc$genomicSuperDups))#21340
table(!is.na(maf.sc$rmsk))#198870
table(maf.sc$Variant_Type)
table(maf.sc$Variant_Type !="SNP")
table(maf.sc[!is.na(maf.sc$rmsk),]$Variant_Type)
table(maf.sc$Score>=0.004)
table(maf.sc[!is.na(maf.sc$genomicSuperDups),]$Score>0.04)
table(maf.sc[!is.na(maf.sc$rmsk) & maf.sc$Variant_Type !="SNP",]$Score<0.1)
table(maf.sc[!is.na(maf.sc$genomicSuperDups) & maf.sc$Variant_Type !="SNP",]$Score<0.1)
#F:458  T:43 
filter.rmsk <- which(!is.na(maf.sc$rmsk))
filter.dups <- which(!is.na(maf.sc$genomicSuperDups) )
filter.repeat <- which(!is.na(maf.sc$rmsk) & !is.na(maf.sc$genomicSuperDups))#union(filter.rmsk,filter.dups)
dim(maf.sc)
dim(maf.sc[-filter.repeat,])
################################################################################
maf.sc <- maf.sc[-filter.repeat,]
unique(maf.sc$Tumor_Sample_Barcode)#13
dim(maf.sc)#2469894     140
unique(maf.sc$Tumor_Sample_Barcode)
maf.Hua1 <- maf.sc


maf.Cervix <- rbind(maf.Hua1, maf.Guo)
################################################################################
# table(solutions$STM.single_cell_genotype.tsv$Start %in% y$Start_Position)
# table(solutions$STM.single_cell_genotype.tsv$Start %in% solutions.all$`MSI_1-2020.variants.avinput`$V2)
# table(maf.Cervix$ExonicFunc.refGene)
# dsdn=as.data.frame(table(maf.Cervix[-(grep("synonymous SNV",maf.Cervix$ExonicFunc.refGene)),]$Tumor_Sample_Barcode)/38)
# head(dsdn,20)
# str(dsdn)
# colnames(dsdn)<-c("sampleID","Freq")
# str(dsdn)
# dim(dsdn)#10
# dsdn=as.data.frame(table(y$ExonicFunc.refGene,y$Hugo_Symbol))
# dsdn <- dsdn[dsdn$Freq!=0,]
# length(unique(maf.Cervix$CB))
# table(y$ExonicFunc.refGene)
# stat <- aggregate(ExonicFunc.refGene~Hugo_Symbol,data=maf.Cervix ,FUN="table")
# unlist(stat$ExonicFunc.refGene)
# 
# drivers <- maf[,c(1:3,5,6)]
# unique = paste(drivers$Start_Position,drivers$Reference_Allele,drivers$Tumor_Seq_Allele2,sep ='-')
# length(unique(unique))
# colnames(drivers) <- c('sampleID','chr','pos','ref','mut')
# length(unique(drivers$pos))
# drivers$chr <- gsub("chr",'',drivers$chr)
# setwd('/public/home/lorihan/lrh/database/Cervix/Guo')
# write.table(drivers[!duplicated(unique),],"maf.Cervix.txt",quote = F,row.names = F)
# library(dndscv)
# # data("dataset_simbreast", package="dndscv")
# dndsout = dndscv(drivers[!duplicated(unique),])
# dndsout <- dndscv(mutations)
# dndsout <- dndscv(
#   drivers[!duplicated(unique),],
#   gene_list = NULL,
#   refdb = "hg38",
#   sm = "192r_3w",
#   kc = "cCervix81",
#   cv = "hg38",
#   max_muts_per_gene_per_sample = 3,
#   max_coding_muts_per_sample = 3000,
#   use_indel_sites = T,
#   min_indels = 5,
#   maxcovs = 20,
#   constrain_wnon_wspl = T,
#   outp = 3,
#   numcode = 1,
#   outmats = F,
#   mingenecovs = 500,
#   dc = NULL
# )
################################################################################
# n <- unique(maf.sc$Tumor_Sample_Barcode)
# scMapping <-  lapply(setNames(n, n), function(nameindex) {
#   y <- maf[maf$Tumor_Sample_Barcode == nameindex,]
#   STM <- y[y$Cell_type_observed=='STM',]
#   STM <- STM[! is.na(STM$FATHMM_score) & STM$FATHMM_score>0.7,]
#   unique(STM$Hugo_Symbol)
#   table(y$Cell_type_observed[y$Hugo_Symbol %in% unique(STM$Hugo_Symbol)[4]])
#   stat <- aggregate(ExonicFunc.refGene~Hugo_Symbol,data=y ,FUN="table")
#  
# })
## 准备和整理了两个文件maf，clincald data,保证有相同Tumor_Sample_Barcode
cohort_directory <- '/public/home/lorihan/lrh/database/Cervix/Output/Epi'
setwd(cohort_directory)
load('Cervix_DEG.RData')
table(Cervix_Epi$Tissue)
clin <- Cervix_Epi
clin$Tissue <- gsub('cervical cancer tissue','CC',clin$Tissue)
clin$Tissue <- gsub('CA_HPV','CC',clin$Tissue)
clin$Tissue <- gsub('neoplasm','CC',clin$Tissue)
clin$Tissue <- gsub('high','HSIL_HPV',clin$Tissue)
clin$Tissue <- gsub('metastatic','CC',clin$Tissue)
# clin$Tissue <- gsub('normal adjacent tissue','ADJ',clin$Tissue)
clin$Tissue <- gsub('normal','ADJ',clin$Tissue)
clin$Tissue <- gsub('AdultCervix','Healthy',clin$Tissue)
clin$Tissue <- gsub('NO_HPV','Healthy',clin$Tissue)
clin$orig.ident <- gsub('Tumor','cervical',clin$orig.ident)
table(clin$Tissue)
clin <- clin[! clin$Tissue %in% c('ADJ','Healthy'),]
clin$Malignancy <-  pc_df$nearest_spline_x_vals[match(clin$orig.ident,pc_df$SimplifiedSampleName)]

maf <- maf.Cervix[! is.na(maf.Cervix$FATHMM_score) & maf.Cervix$FATHMM_score>0.7,]
maf$Tumor_Sample_Barcode <- maf$CB
clin$Tumor_Sample_Barcode <- rownames(clin)
CB <- extract_last_element(str_split(rownames(clin),'-1',simplify = TRUE)[,1])
index <- clin$sample_name %in% c('Guo')
clin$Tumor_Sample_Barcode[index] <- CB[index]

CB <- gsub("-1",'',str_split(rownames(clin),'_',simplify = TRUE)[,2])
index <- clin$sample_name %in% c('Hua1')
clin$Tumor_Sample_Barcode[index] <- CB[index]

clin$Tumor_Sample_Barcode <- paste(clin$Tumor_Sample_Barcode,clin$orig.ident,clin$Cell_Type,sep = '-')
table(maf$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode)
setdiff(maf$Tumor_Sample_Barcode, clin$Tumor_Sample_Barcode)
CGC <- fread('~/lrh/database/All/Cosmic_CancerGeneCensus_v98_GRCh38.tsv.gz')
CGM <- fread('~/lrh/database/All/Cosmic_MutantCensus_v98_GRCh38.tsv.gz')

unique(maf$Hugo_Symbol[maf$Hugo_Symbol %in% CGM$GENE_SYMBOL])
unique(maf$Hugo_Symbol[maf$Hugo_Symbol %in% CGC$GENE_SYMBOL])
intersect(maf$Hugo_Symbol, Patt$Symbol)
intersect(maf$Hugo_Symbol, Patt$Symbol)
STM <- maf[maf$Cell_type_observed=='STM',]
unique(STM$Hugo_Symbol)
intersect(STM$Hugo_Symbol, DEG)
drivers <- intersect(maf$Hugo_Symbol, Patt$Symbol[abs(Patt$Log2FC)>=1])
drivers <- intersect(maf$Hugo_Symbol, Patt$Symbol[Patt$Log2FC>0.5])
drivers <- intersect(STM$Hugo_Symbol, Patt$Symbol)#[abs(Patt$Log2FC)>1])
drivers
# [1] "CLDN18"   "ATP5G3"   "ALDOA"    "GKN2"     "IGLL5"    "TM4SF20" 
# [7] "ANXA10"   "IFITM3"   "APOA1"    "TFF3"     "FTH1"     "TMPRSS15"
intersect(drivers, CGM$GENE_SYMBOL)
intersect(drivers, CGC$GENE_SYMBOL)
sort(intersect(drivers, DEG))
# [1] "ACTG1"    "ANXA2"    "ASH1L"    "CALM2"    "COX7A2"   "CXCL3"    "CYBA"    
# [8] "DBI"      "DNAJA1"   "EEF1A1"   "ELF3"     "FLNB"     "GAPDH"    "GNAS"    
# [15] "GPX4"     "HLA-B"    "HLA-C"    "HLA-DRA"  "HLA-DRB1" "HNRNPA1"  "HSP90AA1"
# [22] "HSPA5"    "HSPH1"    "IER2"     "IFI6"     "JUN"      "KRT19"    "LY6E"    
# [29] "MYL12A"   "MYL12B"   "NACA"     "NFIA"     "NPC2"     "NPM1"     "PABPC1"  
# [36] "PEBP1"    "PFN1"     "RPL10"    "RPL12"    "RPL13"    "RPL21"    "RPL23"   
# [43] "RPL5"     "RPL6"     "RPL8"     "RPS10"    "RPS13"    "RPS15A"   "RPS26"   
# [50] "RPS27A"   "RPS4X"    "RPS6"     "RPS9"     "RPSA"     "SLC25A6"  "SYNE2"   
# [57] "TIMP1"    "TMBIM6"   "TMSB4X"   "YWHAB"       
drivers <- STM[STM$Hugo_Symbol %in% drivers,]

mut <- clin#[clin$sample_name=='Becker',]
dim(mut)
table(mut$Tissue)
# Stat <- aggregate(AAChange.refGene~Tumor_Sample_Barcode,data=maf ,FUN="count")
stats <-  aggregate(AAChange.refGene~Tumor_Sample_Barcode,data=maf ,FUN="unique")
stats$Mut_count <- apply(stats,1,function(x){
  return(length(unlist(x["AAChange.refGene"])))
})
stats$Mut_gene <- apply(stats,1,function(x){
  y <- unlist(x["AAChange.refGene"])#[1]
  z <- str_split(y, ':', simplify = TRUE)[,1]         
  return(z)
})
stats$Mut_gene <- apply(stats,1,function(x){
  y <- unlist(x["AAChange.refGene"])#[1]
  z <- str_split(y, ':', simplify = TRUE)[,1]         
  z <- paste(z, collapse = ', ')
  return(z)
})
stats$AAChange <- apply(stats,1,function(x){
  y <- unlist(x["AAChange.refGene"])#[1]
  if (length(y)>1) {
    w <- paste(y, collapse = ', ')
  } else  {
    w <- y
  }
  return(w)
})
mut <- left_join(mut, stats[,-2],by='Tumor_Sample_Barcode')
rownames(mut) <- rownames(clin)#[clin$sample_name=='Becker',])
# data_as_strings <-  lapply(stats, function(column) sapply(column, toString))#tidyr::unnest(stats)#
write.table(mut[,c('Organ','Malignancy','Tumor_Sample_Barcode','Mut_count','Mut_gene','AAChange')], file = "Cervix_mut.txt", sep = "\t", quote = TRUE,row.names = TRUE)
mut <- read.table('Cervix_mut.txt')
table(mut$Tissue, mut$Cell_Type)
mut$Sample <- str_split(mut$Tumor_Sample_Barcode, '-', simplify = TRUE)[,2]
clin <- clin[clin$sample_name %in% c('Guo','Hua1'),]
table(mut$Sample %in% clin$orig.ident)
table(mut$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode)
# mut$Sample[! mut$Sample %in% clin$orig.ident ]
table(clin$sample_name, clin$Tissue)
mut <- mut[mut$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode,]
mut$Tissue <- clin$Tissue[match(mut$Tumor_Sample_Barcode, clin$Tumor_Sample_Barcode)]
# xlab <- expression(paste("PD-L1 expression[",log["2"],"(TPM+1)]"))
mut$Mut_count[is.na(mut$Mut_count)] <- 0
class(mut$Mut_count); class(mut$Malignancy)
table(is.na(mut$Malignancy))
ggscatter(mut[!is.na(mut$Malignancy),], y = "Mut_count", x = "Malignancy",
          color = "#CAB2D6", size = 1,
          add = "reg.line", add.params = list(color = "salmon", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
          conf.int = TRUE, # 添加回归线的置信区间
          cor.coef = TRUE, # 添加相关系数
          cor.coeff.args = list(method = "pearson", label.x = 0,label.y=7, label.sep =", ") #"\n")#选择Pearson相关
)+ xlab('Malignancy continumn')+ylab("No. of Mutants")#+ stat_cor( label.x = 3)
# 
# CRC_mut <- mut[mut$Tissue %in% c('CRC') & !is.na(mut$Mut_gene),]
# FAP_mut <- mut[mut$Tissue %in% c('FAP') & !is.na(mut$Mut_gene),]
# View(CRC_mut[,c('Malignancy','Tumor_Sample_Barcode','Mut_count','Mut_gene','AAChange')])
# CRC_genes <- unlist(str_split(CRC_mut$Mut_gene,', '))
# FAP_genes <- unlist(str_split(FAP_mut$Mut_gene,', '))
# intersect(unique(CRC_genes),unique(FAP_genes))
# setdiff(unique(CRC_genes),unique(FAP_genes))
# drivers <- intersect(setdiff(unique(CRC_genes),unique(FAP_genes)), Patt$Symbol)
# 
library(ggpubr)
table(mut$Tissue)
table(clin$sample_name, clin$Tissue)
compare_means(Mut_count~Tissue, data=mut,method="wilcox.test",p.adjust.method="BH")
compar<-list(c("CC","HSIL_HPV", "N_HPV"))
# pdf("Cervix_mut.pdf",width = 4,height = 4)
ylab<- "Number of mutations"#expression(paste("FCGR3A  ", "log"["2"], "(SCNA)"))
p1<-ggviolin(mut,"Tissue", "Mut_count", fill = "Tissue",
             palette =c("orange3","#1B9E77","salmon"#palette =c("salmon","#1B9E77"#"#E7298A","#E7B800"#"#E7B800"#)"#00AFBB" #"#8989FF","#00AFBB"
             ), bxp.errorbar = T, #是否添加error bar
             #bxp.errorbar.width = 0.2, #error bar的长度
             #palette = 'npg', #颜色风格
             add.params = list(fill = "white"),
             add = 'boxplot' # 是否添加boxplot上面的点点
)+xlab("")+ylab(ylab)+ggtitle("Cervix cohort")
#p1 + stat_compare_means(comparisons = compar,method = "wilcox.test", label = "p.signif")
p1+geom_signif(comparisons = compar, # 设置要对比的组
               y_position = c(34,36,38), #设置3个显著性标记的高度
               tip_length = c(0), #设置显著性那条横线两头向下的长度
               map_signif_level = F,#T, #设置是否标记显著性的*号，还是直接标记数值
               test =wilcox.test #设置显著性计算方式
) +#theme(plot.title = element_text(hjust = 0.45))+
  theme(legend.position='none')+
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
    #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  )
# dev.off()
p = ggstatsplot::ggbetweenstats(
  data = mut,
  x = Tissue,
  y = Mut_count,
  pairwise.display = 'all',
  p.adjust.method = 'bonferroni',#'BH',
  notch = TRUE, # show notched box plot
  mean.plotting = FALSE, # whether mean for each group is to be displayed
  mean.ci = TRUE, # whether to display confidence interval for means
  mean.label.size = 2.55, # size of the label for mean
  type = "np", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = FALSE, # whether outliers need to be tagged
  xlab = "Sample Type", # label for the x-axis variable
  ylab = "No. of Mutants", # label for the y-axis variable
  title = "Cervix cohort", # title text for the plot
  palette = "Dark2",#"Darjeeling1", # choosing a different color palette
  messages = FALSE
)+   theme(
  plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
  plot.subtitle = element_text(color = 'black', size   = 13,hjust = 0.5),
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
) + theme(legend.position='none')
getwd()
ggsave(plot = p,"Results/Cervix_mut.pdf",width = 8,height = 6)
f = paste0(cohort_directory, '/Results/Cervix_mut.png')
ggsave(plot=p,height=6.5,width=8.5, filename=f, dpi = 500, device = "png")

# d2 <- read.maf(maf=maf,clinicalData=clin)
# length(unique(d2@data$Tumor_Sample_Barcode))#36
# 
# levels(d2@data$Tumor_Sample_Barcode)#36
# table(d2@data$Hugo_Symbol %in% CGC$GENE_SYMBOL) #14
# table(d2@gene.summary$Hugo_Symbol %in% CGC$GENE_SYMBOL)
# d2@gene.summary$Hugo_Symbol[d2@gene.summary$Hugo_Symbol %in% CGC$GENE_SYMBOL]
# table(d2@gene.summary$Hugo_Symbol %in% CGM$GENE_SYMBOL)
# 
# d2@gene.summary$Hugo_Symbol[d2@gene.summary$Hugo_Symbol %in% CGM$GENE_SYMBOL]
# table(d2@data$FATHMM_score>0.7)
# str(d2@clinical.data)
# table(d2@data$Variant_Type)
# getClinicalData(d2) #查看临床信息
# getSampleSummary(d2) #查看每个样品发生突变的情况
# length(unique(d2@data$Tumor_Sample_Barcode))#1004
# dim(d2@data)
# min(!is.na(d2@data$Score))
# min(d2@data$Score)#0.02
# unique(d2@data$Tumor_Sample_Barcode)
# levels(d2@data$Tumor_Sample_Barcode)#369
# 
# table(as.character(d2@clinical.data$Histology))#368
# col = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
# names(col) = c( 'Nonsense_Mutation','Missense_Mutation','Frame_Shift_Del','In_Frame_Del','Frame_Shift_Ins','In_Frame_Ins'
#                 ,'Nonstop_Mutation','Splice_Site','Translation_Start_Site')
# #Color coding for FAB classification; try getAnnotations(x = laml) to see available annotations.
# fabcolors = RColorBrewer::brewer.pal(n = 9,name = 'Spectral')
# fabcolors = list(FAB_classification = fabcolors)
# mycol <- c("firebrick3","seagreen3","skyblue","lightgoldenrod","pink","slateblue","#1F78B4","orange3","magenta")
# vc_cols = mycol#RColorBrewer::brewer.pal(n = 9, name = 'Paired')
# names(vc_cols) = c(
#   'Frame_Shift_Ins',
#   'Frame_Shift_Del',
#   'Multi_Hit',
#   'In_Frame_Del',
#   'In_Frame_Ins',
#   'Splice_Site',
#   'Missense_Mutation',
#   'Nonsense_Mutation',
#   "Nonstop_Mutation"
# )
# #查看变异类型对应的颜色
# print(vc_cols)
# length(unique(d2@data$Hugo_Symbol))
# #gene.detected <- intersect(d2@data$Hugo_Symbol,as.character(str_split(freq_LUAD.SNV$Var1,"[ .]",simplify = T)[,1]))
# #length(unique(gene.detected))#14989
# fabcolors = c("orange4","seagreen4")#RColorBrewer::brewer.pal(n = 3,name = 'Spectral')
# names(fabcolors) = unique(d2@clinical.data$Histology)
# fabcolors = list(Histology = fabcolors)
# length(unique(d2@data$Tumor_Sample_Barcode))
# #png("LUAD_seq_genes.png",width = 940,height = 880)
# oncoplot(clinicalFeatures = c('Tumor_Sample_Barcode'),  #cex.main = 1,
#          fontSize = 1.15,top = 20,
#          maf = d2, #genes=gene.detected,
#          ## colors = col,  mutsig = laml.mutsig, mutsigQval = 0.01, clinicalFeatures = c('sex','pathologic_N'), 
#          annotationColor = fabcolors,
#          sortByAnnotation = TRUE, colors = vc_cols, drawRowBar = TRUE,sepwd_genes = 1,
#          showTumorSampleBarcodes = F,keepGeneOrder=F,removeNonMutated = TRUE,
#          titleFontSize = 0.1,legendFontSize = 2.2,gene_mar = 8,#barcode_mar = 0.1,
#          annotationFontSize = 2.2,additionalFeatureCex = 2.2,additionalFeaturePch = 20,
# )
# 
