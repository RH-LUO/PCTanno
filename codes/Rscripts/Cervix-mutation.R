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

################################################################################
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
################################################################################
# 1) load pdata of PCTanno database
################################################################################
# load('~/lrh/All/Output/organ13_Epi.RData')# from database-sum.R in ts860
# data_path = '/home/lorihan/lrh/All/Output'
data_path = '/data2/rluo4/All/Output'
setwd(data_path)
cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
tissue_summary <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
################################################################################
# NS <- tissue_summary$Replicates[tissue_summary$`Disease Stage` %in% c('Healthy','ADJ')]
# NS <- gsub('cd45+', 'cd45', NS)
# # writeLines(NS, "/home/lorihan/RPMfunc/PCTanno_pdata/Normal_Samples.txt")
# # scp -rp lorihan@10.91.221.59:/home/lorihan/RPMfunc/PCTanno_pdata/* /public/home/lorihan/PCTanno_pdata/
# setwd('/data2/rluo4/RPMfunc/PCTanno_pdata/')
# HTA11_DNA <- read.csv('HTA11_VCF_table.txt',sep = ';')
# HTA11_DNA$Filename <- sapply(HTA11_DNA$Filename, function(x) strsplit(x, '/')[[1]][3])
# # str_split(HTA11_DNA$Filename,'/',simplify = TRUE)[,3]
# HTA11_DNA$Atlas.Name <- sapply(HTA11_DNA$Filename, function(x) strsplit(x, '.rmdup.')[[1]][1])
# #str_split(HTA11_DNA$Filename,'.rmdup.',simplify = TRUE)[,1]
# HTA11_DNA$Atlas.Name <- sapply(HTA11_DNA$Filename, function(x) strsplit(x, '.001')[[1]][1])
# dis <- read.csv('Chen.dis_pdata.csv')
# val <- read.csv('Chen.val_pdata.csv')
# abn <- read.csv('Chen.abn_pdata.csv')
# Chen <- rbind(dis,val,abn)
# table(HTA11_DNA$Biospecimen %in% Chen$Samples)
# setdiff(HTA11_DNA$Biospecimen, Chen$Samples)
# nchar(HTA11_DNA$Biospecimen)
# nchar(Chen$Files)
# nchar(Chen$Samples)
# HTA11_DNA$short <- substr(HTA11_DNA$Biospecimen,1,16)
# Chen$short <- substr(Chen$Files,1,16)
# table(HTA11_DNA$short %in% Chen$short)
# # View(HTA11_DNA[HTA11_DNA$short %ni% Chen$short,])
# table(is.na(match(HTA11_DNA$short,Chen$short)))
# HTA11_DNA$Samples <- Chen$barcode[match(HTA11_DNA$short,Chen$short)]
# HTA11_DNA$Samples[is.na(HTA11_DNA$Samples)] <- HTA11_DNA$Atlas.Name[is.na(HTA11_DNA$Samples)]
# 
# # Chen$Tissue <- CRC_Epi$Tissue[match(Chen$Files,CRC_Epi$orig.ident)]
# # Chen$Malignancy <-  pc_df_Chen$nearest_spline_x_vals[match(Chen$Samples,pc_df_Chen$SimplifiedSampleName)]
# # On HPC:
# # MATCHED_FILES=/public/home/lorihan/PCTanno_pdata/tempfile
# # Delete the matched files or directories
# # while IFS= read -r FILE; do
# # rm -rf "$FILE"
# # done < "$MATCHED_FILES"
# # echo "Files removed."
# # (base) lorihan@cn06 ~/SCOMATIC$ rm */*/scMapping/*/*PASS.tsv
# # (base) lorihan@cn06 ~/SCOMATIC$ rm -rf */*/*/*ealthy* 
# # (base) lorihan@cn06 ~/SCOMATIC$ rm -rf */*/*/*ormal* 
# 
# # lorihan@cn06 ~/SCOMATIC/Cervix/Hua1/scMapping$ rm -r normal/
# # lorihan@cn06 ~/lrh/database/Liver/Losic$ rm -r scMapping
# # lorihan@cn06 ~/lrh/database/Liver/Losic$ cp -rp scMapping_new/ ~/SCOMATIC/Liver/Losic/scMapping
# 
# # lorihan@cn06 ~/lrh/database/Skin/Ji$ rm -r scMapping
# # lorihan@cn06 ~/lrh/database/Skin/Ji$ cp -rp scMapping_new/ ~/SCOMATIC/Skin/Ji/scMapping
# 
# # scp -rp  lorihan@10.91.221.30:/public/home/lorihan/SCOMATIC /home/lorihan/RPMfunc/
################################################################################
# 2) arrange the SNV pdata
################################################################################
# cohort_directory <- '/data2/rluo4/RPMfunc/SCOMATIC'
# setwd(cohort_directory)
# scMapping_dir <- list.dirs(cohort_directory); print(scMapping_dir)
# Diseased_Samples <- scMapping_dir[grepl('scMapping', scMapping_dir)]
# Diseased_Samples <- Diseased_Samples[!grepl("/scMapping$", Diseased_Samples)]
# Diseased_Samples <- data.frame(SNV_path = c(Diseased_Samples) )
# extract_last_element <- function(x) {
#   split_string <- strsplit(x, "/")
#   last_element <- sapply(split_string, function(y) tail(y, n = 1))
#   return(last_element)
# }
# extract_tissue_element <- function(x) {
#   split_string <- strsplit(x, "/")[[1]]
#   return(split_string[6])
# }
# extract_cohort_element <- function(x) {
#   split_string <- strsplit(x, "/")[[1]]
#   return(split_string[7])
# }
# Diseased_Samples$Tissue <- sapply(Diseased_Samples$SNV_path, extract_tissue_element)
# Diseased_Samples$Cohort <- sapply(Diseased_Samples$SNV_path, extract_cohort_element)
# Diseased_Samples$Replicates <- extract_last_element(Diseased_Samples$SNV_path)
# setdiff(Diseased_Samples$Replicates, tissue_summary$Replicates)
# intersect(Diseased_Samples$Replicates, Chen$Samples)
# Diseased_Samples$orig.ident <- Diseased_Samples$Replicates 
# # 1) change sample names of Breast-Wei and Skin-Ji
# Diseased_Samples$orig.ident <- sapply(Diseased_Samples$orig.ident, function(x) strsplit(x, '_scRNA')[[1]][1])
# # 2) change sample names of CRC-Chen
# table(Chen$Files %in% tissue_summary$Replicates)
# index = Diseased_Samples$Replicates %in% Chen$barcode
# Diseased_Samples$orig.ident[index] <- Chen$Files[match(Diseased_Samples$orig.ident[index], Chen$barcode)]
# setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)
# # 3) change sample names of Cervix-Hua1
# Diseased_Samples$orig.ident <- gsub('cervical','Tumor', Diseased_Samples$orig.ident)
# # 4) change sample names of GC-HY
# Diseased_Samples$orig.ident <- gsub('deep_cancer','Deep', Diseased_Samples$orig.ident)
# Diseased_Samples$orig.ident <- gsub('superficial_cancer','Superficial', Diseased_Samples$orig.ident)
# # 5) change sample names of Liver-Ramachandran
# Diseased_Samples$orig.ident <- gsub('cd45p','cd45+', Diseased_Samples$orig.ident)
# setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)
# # 6) change sample names of Prostate-Joseph
# remove_first_element <- function(char_vector) {
#   char_vector_parts <- strsplit(char_vector, "_")
#   char_vector_parts <- lapply(char_vector_parts, function(parts) {
#     if (length(parts) > 1) {
#       return(paste(parts[-1], collapse = "_"))
#     } else {
#       return(parts[1])
#     }
#   })
#   return(unlist(char_vector_parts))
# }
# Joseph <- tissue_summary[tissue_summary$Cohort=='Joseph',]
# index = Diseased_Samples$Replicates %in% remove_first_element(Joseph$Replicates)
# table(index)
# Diseased_Samples$orig.ident[index] <- Joseph$Replicates[match(Diseased_Samples$orig.ident[index], remove_first_element(Joseph$Replicates))]
# # 7) change sample names of Liver-Meng
# Meng <- tissue_summary[tissue_summary$Cohort=='Meng',]
# index = Diseased_Samples$Replicates %in% gsub('_Meng','', Meng$Replicates)
# table(index)
# index = Diseased_Samples$Cohort=='Meng'
# # Meng$Replicates[match(Diseased_Samples$orig.ident[index],  gsub('_Meng','', Meng$Replicates))]
# Diseased_Samples$orig.ident[index] <- paste0(Diseased_Samples$orig.ident[index] ,'_Meng')
# # 8) change sample names of Prostate-Dong
# Dong <- tissue_summary[tissue_summary$Cohort=='Dong',]
# index = Diseased_Samples$Replicates %in% remove_first_element(Dong$Replicates)
# table(index)
# Dong$Replicates[match(Diseased_Samples$orig.ident[index], remove_first_element(Dong$Replicates))]
# index = Diseased_Samples$Cohort=='Dong'
# Diseased_Samples$orig.ident[index] <- paste0('Dong_', Diseased_Samples$orig.ident[index] )
# setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)
# # UNC or NL: [1] "HTA11_8099_200000101113111" [2] "HTA11_392_300001101113111"[3] "HTA11_7956_200000101113211"
# # Remove: RNA-P23T2-P23T2-3
# remove_samples <- setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)# 4
# Diseased_Samples$GSE <- tissue_summary$Dataset[match(Diseased_Samples$Cohort, tissue_summary$Cohort)]
# write.table(Diseased_Samples, sep = "\t", file = '/data2/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples.txt')
# Diseased_Samples <- Diseased_Samples[Diseased_Samples$orig.ident %in% tissue_summary$Replicates,]
# table(Diseased_Samples$orig.ident %in% cell_summary$orig.ident) #703 samples from SNV calling
# unique(Diseased_Samples$Cohort)#32
################################################################################
# 3) arrange the SNV data
################################################################################
# CHROM  Start   End     REF     ALT     FILTER  Cell_types      Up_context      Down_context    N_ALT   
# Dp      Nc   -Bc      Cc      VAF     CCF     BCp     CCp     Cell_types_min_BC       Cell_types_min_CC       Rest_BC Rest_CC Fisher_p      Cell_type_Filter        INFO  
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
# Load Diseased Samples
Diseased_Samples <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples.txt')
# Load Normal Samples
normal_samples <- unique(cell_summary$orig.ident[cell_summary$Tissue %in% c('Healthy','ADJ')])
remove_samples <- setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)# 4
normal_samples <- c(normal_samples, remove_samples) #417
# Tissue <- cell_summary[cell_summary$Organ==tis,]
# unique(Tissue$sample_name)
# Helper function to load and filter files
load_and_filter_files <- function(dir, normal_samples, pattern) {
  solution <- list.files(dir, pattern = pattern)
  nor <- paste0(normal_samples, '.', pattern)  # Adjusted the separator
  solution <- solution[!solution %in% nor]
  file_paths <- file.path(dir, solution)
  # Use fread correctly within lapply
  solutions <- lapply(file_paths, function(fp) {
    fread(fp, sep = '\t', stringsAsFactors = FALSE)
  })
  names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.variants.avinput')[[1]][1])
  return(solutions)
}
# Helper function to load and filter Annovar files
load_and_filter_annovarToMaf <- function(dir, normal_samples, pattern) {
  solution <- list.files(dir, pattern = pattern)
  nor <- paste(normal_samples, pattern, sep = '.')
  solution <- solution[!solution %in% nor]
  file_paths <- file.path(dir, solution)
  # solutions <- lapply(file_paths, annovarToMaf)
  solutions <- lapply(file_paths, function(fp) {
    annovarToMaf(fp, refBuild = "hg38")
  })
  names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.hg38_multianno.txt')[[1]][1])
  return(solutions)
}
# Function to add VAF and CCF
add_vaf_ccf <- function(maf, solution, nameindex) {
  y <- as.data.frame(maf[[nameindex]])
  z <- gsub('.hg38_multianno.txt', '.variants.avinput', nameindex)
  w <- as.data.frame(solution[[z]])
  VAF <- sapply(w$V6, function(x) strsplit(x, '\\-')[[1]][15])
  CCF <- sapply(w$V6, function(x) strsplit(x, '\\-')[[1]][16])
  y <- mutate(y, VAF = VAF, CCF = CCF)
  return(y)
}
# Function to add VAF and CCF
add_vaf_ccf <- function(maf, solution, nameindex) {
  y <- as.data.frame(maf[[nameindex]])
  z <- gsub('.hg38_multianno.txt', '.variants.avinput', nameindex)
  w <- as.data.frame(solution[[z]])
  # Ensure w$V6 is split correctly and the lengths match
  split_v6 <- str_split(w$V6, '[-]', simplify = TRUE)
  if (nrow(split_v6) == nrow(y)) {
    VAF <- split_v6[, 15]
    CCF <- split_v6[, 16]
    y$VAF <- VAF
    y$CCF <- CCF
    return(y)
  } 
  # else {
  #   stop("Length of split V6 does not match the number of rows in y")
  # }
}
# Extract clean barcode
extract_clean_barcode <- function(barcodes) {
  patterns <- c(".*_([ACGT]+).\\d.*",
    ".*_([ACGT]+)-\\d.*", ".*_([ACGT]+)_\\d.*", ".*_([ACGT]+)-[A-Za-z0-9]+$", 
    ".*_([ACGT]+)$", ".*([ACGT]+)-[A-Za-z0-9]+$"
  )
  
  extract_barcode <- function(barcode) {
    for (pattern in patterns) {
      match <- regmatches(barcode, regexec(pattern, barcode))
      if (length(match[[1]]) > 1) {
        return(match[[1]][2])
      }
    }
    return(NA)
  }
  
  clean_barcodes <- sapply(barcodes, extract_barcode)
  return(clean_barcodes)
}

# Function to apply filters
apply_filters <- function(y) {
  filter_maf <- unique(c(
    which(y$X1000g2015aug_all > 0.05),
    which(y$ExAC_ALL > 0.05),
    which(y$gnomAD_genome_ALL > 0.05)
  ))
  y <- y[-filter_maf,]
  y <- y[y$VAF != '.' & y$VAF > 0.05, ]
  filter_repeat <- which(!is.na(y$rmsk) & !is.na(y$genomicSuperDups))
  y <- y[-filter_repeat, ]
  y <- y[!is.na(y$FATHMM_score) & y$FATHMM_score > 0.7, ]
  return(y)
}

# Process scMapping data
process_scMapping <- function(dir, annovar_mafs, cohort) {
  scMapping <- list.files(dir)
  scMapping <- scMapping[!scMapping %in% normal_samples]
  anno <- vector("list", length(scMapping))
  names(anno) <- scMapping
  
  for (i in scMapping) {
    path <- file.path(dir, i)
    solution <- list.files(path, pattern = 'single_cell_genotype.tsv')
    file_paths <- file.path(path, solution)
    solutions <- lapply(file_paths, fread, stringsAsFactors = FALSE)
    names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.')[[1]][1])
    
    solutions <- lapply(names(solutions), function(nameindex) {
      b <- as.data.frame(annovar_mafs[[i]])
      c <- as.data.frame(solutions[[nameindex]])
      colnames(c)[1:2] <- colnames(b)[2:3]
      c$Start_Position <- as.character(c$Start_Position)
      y <- left_join(b, c[,-1], by = 'Start_Position')
      y <- y[!is.na(y$Cell_type_observed),]
      y <- y[y$Tumor_Sample_Barcode %in% Diseased_Samples$Replicates,]
      index <- match(y$Tumor_Sample_Barcode, Diseased_Samples[Diseased_Samples$Cohort==cohort,]$Replicates)
      y$Tumor_Sample_Barcode <- Diseased_Samples[Diseased_Samples$Cohort==cohort,]$orig.ident[index]
      
      Epi <- cell_summary[cell_summary$sample_name == cohort,]
      # Epi$CB <- extract_clean_barcode(Epi$barcode)
      Epi$CB <- paste(Epi$CB, Epi$orig.ident, sep = '--')
      y$CB <- paste(y$CB, y$Tumor_Sample_Barcode, sep = '--')
      Epi$CB <- paste(Epi$CB, Epi$Cell_Type, sep = '--')
      y$CB <- paste(y$CB, y$Cell_type_observed, sep = '--')
      
      y <- y[y$CB %in% Epi$CB, ]
      y <- apply_filters(y)
      return(y)
    })
    
    anno[[i]] <- data.table::rbindlist(solutions, fill = TRUE)
  }
  
  scRNA_annovar <- NULL
  scRNA_annovar <- do.call(rbind, anno)
  
  return(scRNA_annovar)
}
# Process scMapping data for Cervix
# datasets <- unique(Diseased_Samples$Cohort)
# for(cohort in datasets){
#   # Define paths
#   base_dir <- '/data2/rluo4/RPMfunc/SCOMATIC'
#   start_time <- Sys.time()
#   tis <- Diseased_Samples$Tissue[match(cohort, Diseased_Samples$Cohort)]
#   print(paste("start parse:", tis, '--', cohort, ' at', start_time))
#   
#   annovar_dir <- file.path(base_dir, tis, cohort, "anno.var")
#   scMapping_dir <- file.path(base_dir, tis, cohort, "scMapping")
#   out_dir <-  '/data2/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
#   gse_acc <- Diseased_Samples$GSE[match(cohort, Diseased_Samples$Cohort)]
#   SNV_data <- file.path(out_dir, paste0(gse_acc, '_SCOMATIC_annovar.RData'))
#   if(! tis %in% c('Breast', 'Cervix', 'CRC', 'Endometrium', 'Esophagus', 'GC')){
#     print(paste0(tis, " run on ts860 !"))
#     next;
#   }
#   if( !file.exists(SNV_data)){
#     # Load and filter Guo files
#     test_solutions <- load_and_filter_files(annovar_dir, normal_samples, 'variants.avinput')
#     test_annovar_mafs <- load_and_filter_annovarToMaf(annovar_dir, normal_samples, 'hg38_multianno.txt')
#     test_annovar_mafs <- lapply(names(test_annovar_mafs), function(nameindex) {
#       add_vaf_ccf(test_annovar_mafs, test_solutions, nameindex)
#     })
#     names(test_annovar_mafs) <- names(test_solutions)
#     test_scRNA_annovar <- process_scMapping(scMapping_dir, test_annovar_mafs, cohort)
#     # Combine Guo and Hua1 data
#     # Save final data
#     save(test_scRNA_annovar, file = SNV_data)
#     end_time <- Sys.time()
#     print(end_time)
#     time_taken <- end_time - start_time
#     print(paste("Time taken to process", tis, '--', cohort, ':', time_taken, 'mins' ))
#     
#   }
# } #done in ts860 + 57
################################################################################
# 4) visualize the SNV data of Cervix
################################################################################
library(data.table)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggsignif)
library(ggstatsplot)
load('/data2/rluo4/All/Output/pheno_all.rds')#Cervix_Epi
# Base directory and output directory
base_dir <- '/data2/rluo4/RPMfunc/SCOMATIC'
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
# Get unique cohorts
combined_cohort <- unique(Diseased_Samples$Cohort)#[Diseased_Samples$Tissue == 'Prostate'])
# Initialize an empty data frame to store merged data
maf <- data.frame()
maf_results <-  vector("list", length(combined_cohort))
names(maf_results) <- combined_cohort
# Loop through each cohort to load and combine data
# for (cohort in combined_cohort) {
#   # Define paths
#   tis <- Diseased_Samples$Tissue[match(cohort, Diseased_Samples$Cohort)]
#   annovar_dir <- file.path(base_dir, tis, cohort, "anno.var")
#   scMapping_dir <- file.path(base_dir, tis, cohort, "scMapping")
#   gse_acc <- Diseased_Samples$GSE[match(cohort, Diseased_Samples$Cohort)]
#   SNV_data <- file.path(out_dir, paste0(gse_acc, '_SCOMATIC_annovar.RData'))
#   # Load data
#   if( file.exists(SNV_data)){
#     load(SNV_data)
#   maf_results[[cohort]] <- test_scRNA_annovar
#   SNV_data <- file.path(out_dir, paste0(gse_acc, '_SCOMATIC_annovar.rds'))
#   saveRDS(test_scRNA_annovar, file = SNV_data)
#   }
#   # # Convert to data frame and append to maf
#   # scRNA.annovar <- data.frame(test_scRNA_annovar)
#   # maf <- rbind(maf, scRNA.annovar)
# }
for (cohort in combined_cohort) {
  # Define paths
  tis <- Diseased_Samples$Tissue[match(cohort, Diseased_Samples$Cohort)]
  annovar_dir <- file.path(base_dir, tis, cohort, "anno.var")
  scMapping_dir <- file.path(base_dir, tis, cohort, "scMapping")
  gse_acc <- Diseased_Samples$GSE[match(cohort, Diseased_Samples$Cohort)]
  SNV_data <- file.path(out_dir, paste0(gse_acc, '_SCOMATIC_annovar.rds'))
  # Load data
  if( file.exists(SNV_data)){
    maf_results[[cohort]] <- readRDS(SNV_data)
    }

}
# Display unique sample barcodes
unique(maf$Tumor_Sample_Barcode)
# Table of mutations by tissue and cell type
table(mut$Tissue, mut$Cell_Type)
# Extract sample information
mut$Sample <- str_split(mut$CB, '--', simplify = TRUE)[, 2]

# Filter and prepare Epi data
Epi <- cell_summary %>%
  filter(sample_name %in% combined_cohort) %>%
  filter(!orig.ident %in% normal_samples) %>%
  mutate(CB = paste(extract_clean_barcode(barcode), orig.ident, sep = '--'),
         CB = paste(CB, Cell_Type, sep = '--'))
# Table of Epi CB matches
table(maf$CB %in% Epi$CB)
# Handle NA values in AAChange.refGene
maf$AAChange.refGene[is.na(maf$AAChange.refGene)] <- maf$Hugo_Symbol[is.na(maf$AAChange.refGene)]
# Aggregate mutation data
stats <- maf %>%
  group_by(CB) %>%
  summarise(AAChange.refGene = list(unique(AAChange.refGene))) %>%
  mutate(Mut_count = lengths(AAChange.refGene),
         Mut_gene = sapply(AAChange.refGene, function(x) paste(unique(str_split(x, ':', simplify = TRUE)[, 1]), collapse = ', ')),
         AAChange = sapply(AAChange.refGene, function(x) paste(x, collapse = ', ')))
# Join Epi and stats data
Epi <- left_join(Epi, stats, by = 'CB')
# Add malignancy information
tissue <- paste0(tis, '_Epi')
Epi$Malignancy <- pheno_all[[tissue]]$Malignancy[match(Epi$orig.ident, pheno_all[[tissue]]$orig.ident)]
# Prepare mutation data for analysis
mut <- Epi %>%
  select(Organ, Tissue, Malignancy, Cell_Type, CB, Mut_count, Mut_gene, AAChange) %>%
  mutate(Sample = str_split(Tumor_Sample_Barcode, '[--]', simplify = TRUE)[, 2],
         Mut_count = replace_na(Mut_count, 0))
# Statistical comparison
comparison <- list(c("CC", "HSIL_HPV", "N_HPV"))
ylab <- "Number of mutations"
# Plot violin plot with significance
p1 <- ggviolin(mut, "Tissue", "Mut_count", fill = "Tissue",
               palette = c("orange3", "#1B9E77", "salmon"),
               bxp.errorbar = TRUE,
               add = 'boxplot') +
  xlab("") + ylab(ylab) + ggtitle("Cervix cohort") +
  geom_signif(comparisons = comparison,
              y_position = c(34, 36, 38),
              tip_length = c(0),
              map_signif_level = FALSE,
              test = wilcox.test) +
  theme(legend.position = 'none',
        plot.title = element_text(color = 'black', size = 16, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size = 16, hjust = 0.5),
        plot.caption = element_text(color = 'black', size = 16, face = 'italic', hjust = 1),
        axis.text.x = element_text(color = 'black', size = 16),
        axis.text.y = element_text(color = 'black', size = 16),
        axis.title.x = element_text(color = 'black', size = 16),
        axis.title.y = element_text(color = 'black', size = 16, angle = 90),
        legend.title = element_text(color = 'black', size = 16),
        legend.text = element_text(color = 'black', size = 16),
        axis.line.y = element_line(color = 'black'),
        axis.line.x = element_line(color = 'black'))
# Save violin plot
ggsave(plot = p1, filename = "Results/Cervix_mut.pdf", width = 8, height = 6)

# Create and save box plot with ggbetweenstats
p <- ggstatsplot::ggbetweenstats(
  data = mut,
  x = Tissue,
  y = Mut_count,
  pairwise.display = 'all',
  p.adjust.method = 'bonferroni',
  notch = TRUE,
  mean.plotting = FALSE,
  mean.ci = TRUE,
  mean.label.size = 2.55,
  type = "np",
  k = 2,
  outlier.tagging = FALSE,
  xlab = "Sample Type",
  ylab = "No. of Mutants",
  title = "Cervix cohort",
  palette = "Dark2",
  messages = FALSE
) + theme(legend.position = 'none')
# Save box plot
ggsave(plot = p, filename = "Results/Cervix_mut.png", width = 8.5, height = 6.5, dpi = 500, device = "png")



base_dir <- '/data2/rluo4/RPMfunc/SCOMATIC'
combined_cohort <- unique(Diseased_Samples$Cohort[Diseased_Samples$Tissue=='Cervix'])
maf <- NULL
for (cohort in combined_cohort) {
  tis <- Diseased_Samples$Tissue[match(cohort, Diseased_Samples$Cohort)]
  annovar_dir <- file.path(base_dir, tis, cohort, "anno.var")
  scMapping_dir <- file.path(base_dir, tis, cohort, "scMapping")
  out_dir <-  '/data2/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
  gse_acc <- Diseased_Samples$GSE[match(cohort, Diseased_Samples$Cohort)]
  SNV_data <- file.path(out_dir, paste0(gse_acc, '_SCOMATIC_annovar.RData'))
  load(SNV_data)
  scRNA.annovar <- data.frame(test_scRNA_annovar)
  maf <<- rbind(maf, scRNA.annovar)
}
unique(maf$Tumor_Sample_Barcode)#13

table(mut$Tissue, mut$Cell_Type)
mut$Sample <- str_split(mut$Tumor_Sample_Barcode, '-', simplify = TRUE)[,2]
Epi <- cell_summary[cell_summary$sample_name %in% combined_cohort, ]
Epi <- Epi[! Epi$orig.ident %in% normal_samples,]
Epi$CB <- extract_clean_barcode(Epi$barcode)
Epi$CB <- paste(Epi$CB, Epi$orig.ident, sep = '--')
Epi$CB <- paste(Epi$CB, Epi$Cell_Type, sep = '--')
table(maf$CB %in% Epi$CB)

maf$AAChange.refGene[is.na(maf$AAChange.refGene)] <- maf$Hugo_Symbol[is.na(maf$AAChange.refGene)]
stats <-  aggregate(AAChange.refGene~CB,data=maf ,FUN="unique")
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
Epi <- left_join(Epi, stats[,-2],by='CB')
tissue <- paste0(tis, '_Epi')
Epi$Malignancy <- pheno_all[[tissue]]$Malignancy[match(Epi$orig.ident,pheno_all[[tissue]]$orig.ident)]
# rownames(Epi) <- rownames(clin)#[clin$sample_name=='Becker',])
mut <- Epi[,c('Organ', 'Tissue','Malignancy', 'Cell_Type','CB','Mut_count','Mut_gene','AAChange')]
# write.table(mut[,c('Organ','Malignancy','Tumor_Sample_Barcode','Mut_count','Mut_gene','AAChange')], file = "Cervix_mut.txt", sep = "\t", quote = TRUE,row.names = TRUE)
# mut <- read.table('Cervix_mut.txt')
table(mut$Tissue, mut$Cell_Type)
mut$Sample <- str_split(mut$Tumor_Sample_Barcode, '[--]', simplify = TRUE)[,2]

# xlab <- expression(paste("PD-L1 expression[",log["2"],"(TPM+1)]"))
mut$Mut_count[is.na(mut$Mut_count)] <- 0
library(ggpubr)
table(mut$Tissue)
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
