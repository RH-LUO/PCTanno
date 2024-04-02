##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) ########------ Define the key transitioning clusters based on clinical implications ----#######
library(ggsci)
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
Trace_Clusters =  vector("list", length(organ_all))
names(Trace_Clusters) <- organ_all
########################################################################
for (organ in organ_all) { # According to related ClinInfo_df of TCGA, we used  endpoint clusters in different paths (with two branches) or both endpoint and start-point clusters (only on branch) for comparisons of ssgsea scores  
  print(organ)
  if(organ == 'Breast'){# BRCA: Tumor Suppression - C14: (-) OS; C9: (+) DFS & OS 
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C4", precancer = "T1-C14") #c("T1-C4", "T1-C9")
  }
  if(organ == 'Cervix'){# CESC: Tumor Promotion ?
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C14", precancer = "T3-C6")#data.frame(cancer = "T1-C1", precancer = "T1-C5")#data.frame(cancer = "T3-C16", precancer = "T1-C5")#data.frame(cancer = "T2-C4", precancer = "T2-C10")#c("T1-C1", "T1-C5", "T2-C4", "T2-C10", "T3-C6", "T3-C5")
  }
  if(organ == 'Chen'){## COADREAD: Tumor Promotion - C1: (+) OS; 
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C1", precancer = "T3-C17")#c("T1-C1", "T1-C24", "T2-C26", "T2-C3", "T3-C17", "T3-C23")
  }
  if(organ == 'CRC'){# CRC: Tumor Promotion ? #which is Becker: 5
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C3", precancer = "T1-C16")#precancer = "T1-C19")#c("T1-C3", "T1-C19", "T1-C16", "T1-C17")
  }
  if(organ == 'Endometrium'){# UCEC: Tumor Promotion? C24: (+) Stage; C19: (-) Stage
    Trace_Clusters[[organ]] <- data.frame(cancer = "T2-C24", precancer = "T2-C4")# data.frame(cancer = "T1-C14", precancer = "T1-C12") #c("T1-C14", "T1-C12", "T2-C24", "T2-C19")
  }
  if(organ == 'Esophagus'){# ESCA: Tumor Promotion ?
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C10", precancer = "T1-C12")#c("T1-C10", "T1-C12")
  }
  if(organ == 'GC'){# STAD: Tumor Promotion ?
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C17", precancer = "T1-C2")#data.frame(cancer = "T2-C3", precancer = "T2-C7")#c("T1-C13", "T1-C2", "T2-C3", "T2-C12")
  }
  if(organ == 'HNSCC'){# HNSC: Tumor Suppression - C3: (+) DFS & OS
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C4", precancer = "T1-C16")# data.frame(cancer = "T1-C4", precancer = "T2-C3")# tumor suppression #data.frame(cancer = "T1-C4", precancer = "T1-C18") #tumor promotion #data.frame(cancer = "T2-C3", precancer = "T2-C19")#c("T1-C11", "T1-C18", "T2-C19", "T2-C3")#
  }
  if(organ == 'Liver'){# LIHC: Tumor Promotion - C0: (+) OS; C15: (-) DFS & OS
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C0", precancer = "T1-C1")#data.frame(cancer = "T1-C0", precancer = "T2-C7")
  }  
  if(organ == 'Lung'){# LUSC & LUAD: Tumor Promotion - C12, C14: (+) DFS & OS; C7: (-) DFS
    Trace_Clusters[[organ]] <-  data.frame(cancer = "T1-C14", precancer = "T1-C8")#c("T1-C12", "T1-C6", "T1-C4", "T1-C3")
  }
  if(organ == 'Pancreas'){# PAAD: Tumor Promotion - C4: (+) OS; C9: (-) DFS & OS; C3: (+) DFS
    Trace_Clusters[[organ]] <- data.frame(cancer = "T2-C4", precancer = "T2-C9")#c("T1-C8", "T1-C2", "T2-C4", "T2-C9")
  }
  if(organ == "Prostate"){# PRAD: Tumor Promotion - C34: (+) Outcome; C42: (-) Outcome; ? C16 + C34, C42-->C5; # C34: (-) OS; C42: (+) DFS & OS
    Trace_Clusters[[organ]] <-  data.frame(cancer = "T2-C34", precancer = "T2-C14")#data.frame(cancer = "T2-C34", precancer = "T2-C17")
  }
  if(organ == "Skin"){# Tumor Promotion  
    Trace_Clusters[[organ]] <- data.frame(cancer = "T2-C14", precancer = "T1-C4")#data.frame(cancer = "T1-C11", precancer = "T1-C9")#data.frame(cancer = "T1-C12", precancer = "T1-C11") #SEN-DN: data.frame(cancer = "T1-C12", precancer = "T1-C5")  #c("T1-C12", "T1-C9", "T2-C1", "T2-C11")
  } 
  if(organ == "THCA"){# THCA: Tumor Promotion  ?
    Trace_Clusters[[organ]] <- data.frame(cancer = "T1-C5", precancer = "T2-C2") #or  data.frame(cancer = "T1-C1", precancer = "T1-C15") #c("T1-C1", "T1-C15", "T2-C5", "T2-C2")
  }
  
}
load(file = paste0(data_path,'CF_cor_clin.rds'))
Survival_df <- surv_df#[surv_df$cancer_type %in% my_ct,]
Survival_df <- Survival_df[Survival_df$pval < 0.05, ]
ClinInfo_df <- categ_df#[categ_df$cancer_type %in% my_ct,]
ClinInfo_df <- ClinInfo_df[ClinInfo_df$pval < 0.05, ]

data_path = '/data/rluo4/All/Output/TCGA/'
ClinInfo <- Survival_df[Survival_df$cancer_type %in% c('BRCA', 'COAD', 'LIHC','LUSC','LUAD','PAAD', 'HNSC'),]#PRAD
clusters <- unlist(Trace_Clusters)
cohort_index$Clusters <- Trace_Clusters[match(cohort_index$PCT, names(Trace_Clusters))]
cohort_index$Clusters <- apply(cohort_index, 1, function(x){
  y <- str_split(paste0(unlist(x['Clusters'])), '-', simplify = T)[,2]
  z <- paste0(y, collapse = ', ')
  return(z)
})
cohort_index$Cc <- str_split(cohort_index$Clusters, ', ', simplify = T)[,1]
cohort_index$Cp <- str_split(cohort_index$Clusters, ', ', simplify = T)[,2]
index1 <- paste(ClinInfo$cancer_type, ClinInfo$cibersort_fraction)
cohort_index$Cc <- paste(cohort_index$TCGA, cohort_index$Cc)
cohort_index$Cp <- paste(cohort_index$TCGA, cohort_index$Cp)
index2 <- c('PAAD C1', 'PAAD C11')
ClinInfo <- ClinInfo[ ! index1 %in% index2,]#c(cohort_index$Cc, cohort_index$Cp) ,]
ClinInfo$cibersort_fraction <- paste(ClinInfo$cancer_type, ClinInfo$cibersort_fraction, sep = '-')

cor_clin <- matrix(data = 0, ncol = length(unique(ClinInfo$cibersort_fraction)), nrow = length(unique(ClinInfo$var_name)))
rownames(cor_clin) <- unique(ClinInfo$var_name)
colnames(cor_clin) <- unique(ClinInfo$cibersort_fraction)
for (i in rownames(cor_clin)) {
  for (j in colnames(cor_clin)) {
    index = ClinInfo$var_name==i & ClinInfo$cibersort_fraction==j
    if(length(unique(index))>1){
      cor_clin[i,j]<- as.numeric(ClinInfo[index,'eff'])# fill data
    }
  }
}

p_clin <- matrix(data = 1, ncol = length(unique(ClinInfo$cibersort_fraction)), nrow = length(unique(ClinInfo$var_name)))
rownames(p_clin) <- unique(ClinInfo$var_name)
colnames(p_clin) <- unique(ClinInfo$cibersort_fraction)
for (i in rownames(p_clin)) {
  for (j in colnames(p_clin)) {
    index = ClinInfo$var_name==i & ClinInfo$cibersort_fraction==j
    if(length(unique(index))>1){
      p_clin[i,j]<- as.numeric(ClinInfo[index,'pval'])# fill data
    }
  }
}

library(reshape2)
pmt.0 <- p_clin
cmt.0 <- cor_clin
if (!is.null(pmt.0)){
  sssmt <-  pmt.0< 0.001
  pmt.0[sssmt] <-'***'
  ssmt <- pmt.0 >=0.001 &pmt.0< 0.01
  pmt.0[ssmt] <-'**'
  smt <- pmt.0 >0.01 & pmt.0 <0.05
  pmt.0[smt] <- '*'
  pmt.0[!sssmt&!ssmt&!smt]<- ''
} else {
  pmt.0 <- F
}
class(cmt.0)
cmt <- as.matrix(cmt.0) # rbind(cmt.0,cmt.1)
pmt <- as.matrix(pmt.0) #rbind(pmt.0,pmt.1)
df <-melt(cmt,value.name="cor")
df$pvalue <-as.vector(pmt)
#df[df$pvalue<0.05 & abs(df$cor)>0.4,]
max(cmt)
min(cmt)
library(ComplexHeatmap)
library(grid)
library(circlize)
#?ComplexHeatmap::pheatmap()
mycol<-colorRampPalette(c("#1B9E77","lightyellow","firebrick"))(800)
mycol<-colorRamp2(c(-10, 0, 10),c("#1B9E77","white","firebrick"))#c("seagreen4", "lightyellow", "orange"))
p = ComplexHeatmap::pheatmap(cmt,scale = "none", cluster_row = F, cluster_col = F, name="logHR",
                             border=NA, fontsize_row=12, fontsize_col=12,
                             display_numbers = pmt, #display_numbers = round(cmt,2),
                             fontsize_number = 14, 
                             main="",fontsize = 10,
                             number_color = "white", cellwidth =21, cellheight =21,
                             color=mycol)
f = paste0(data_path,'res_fig/Cox_PCT_TCGA_heatmap.png')
png(height=2,width=5.5, filename=f, units = "in", res = 500)#
draw(p, heatmap_legend_side = "right", annotation_legend_side = "right", legend_grouping = "original")
dev.off()



library(ggsci)
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
Precancer_Clusters =  vector("list", length(organ_all))
names(Precancer_Clusters) <- organ_all
########################################################################
for (organ in organ_all) { # According to related ClinInfo_df of TCGA, we used  endpoint clusters in different paths (with two branches) or both endpoint and start-point clusters (only on branch) for comparisons of ssgsea scores
  print(organ)
  if(organ == 'Breast'){# BRCA: Tumor Suppression - C14: (-) OS; C9: (+) DFS & OS
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C14", precancer = "T1-C9") #c("T1-C4", "T1-C9")
  }
  if(organ == 'Cervix'){# CESC: Tumor Promotion ?
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C1", precancer = "T1-C5")#data.frame(cancer = "T1-C1", precancer = "T1-C5")#data.frame(cancer = "T3-C16", precancer = "T1-C5")#data.frame(cancer = "T2-C4", precancer = "T2-C10")#c("T1-C1", "T1-C5", "T2-C4", "T2-C10", "T3-C6", "T3-C5")
  }
  if(organ == 'Chen'){## COADREAD: Tumor Promotion - C1: (+) OS;
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C24", precancer = "T3-C17")#c("T1-C1", "T1-C24", "T2-C26", "T2-C3", "T3-C17", "T3-C23")
  }
  if(organ == 'CRC'){# CRC: Tumor Promotion ? #which is Becker: 5
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C17", precancer = "T1-C16")#precancer = "T1-C19")#c("T1-C3", "T1-C19", "T1-C16", "T1-C17")
  }
  if(organ == 'Endometrium'){# UCEC: Tumor Promotion? C24: (+) Stage; C19: (-) Stage
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T2-C4", precancer = "T2-C19")# data.frame(cancer = "T1-C14", precancer = "T1-C12") #c("T1-C14", "T1-C12", "T2-C24", "T2-C19")
  }
  if(organ == 'Esophagus'){# ESCA: Tumor Promotion ?
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C10", precancer = "T1-C12")#c("T1-C10", "T1-C12")
  }
  if(organ == 'GC'){# STAD: Tumor Promotion ?
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T2-C7", precancer = "T2-C12")#data.frame(cancer = "T2-C3", precancer = "T2-C7")#c("T1-C13", "T1-C2", "T2-C3", "T2-C12")
  }
  if(organ == 'HNSCC'){# HNSC: Tumor Suppression - C3: (+) DFS & OS
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C16", precancer = "T1-C18")# data.frame(cancer = "T1-C4", precancer = "T2-C3")# tumor suppression #data.frame(cancer = "T1-C4", precancer = "T1-C18") #tumor promotion #data.frame(cancer = "T2-C3", precancer = "T2-C19")#c("T1-C11", "T1-C18", "T2-C19", "T2-C3")#
  }
  if(organ == 'Liver'){# LIHC: Tumor Promotion - C0: (+) OS; C15: (-) DFS & OS
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T2-C15", precancer = "T2-C11")#data.frame(cancer = "T1-C0", precancer = "T2-C7")
  }
  if(organ == 'Lung'){# LUSC & LUAD: Tumor Promotion - C12, C14: (+) DFS & OS; C7: (-) DFS
    Precancer_Clusters[[organ]] <-  data.frame(cancer = "T1-C8", precancer = "T1-C6")#c("T1-C12", "T1-C6", "T1-C4", "T1-C3")
  }
  if(organ == 'Pancreas'){# PAAD: Tumor Promotion - C4: (+) OS; C9: (-) DFS & OS; C3: (+) DFS
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T2-C9", precancer = "T2-C6")#c("T1-C8", "T1-C2", "T2-C4", "T2-C9")
  }
  if(organ == "Prostate"){# PRAD: Tumor Promotion - C34: (+) Outcome; C42: (-) Outcome; ? C16 + C34, C42-->C5; # C34: (-) OS; C42: (+) DFS & OS
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C42", precancer = "T2-C14")#c("T1-C6", "T1-C9", "T2-C34", "T2-C14")
  }
  if(organ == "Skin"){# Tumor Promotion
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C11", precancer = "T1-C4")#data.frame(cancer = "T1-C11", precancer = "T1-C9")#data.frame(cancer = "T1-C12", precancer = "T1-C11") #SEN-DN: data.frame(cancer = "T1-C12", precancer = "T1-C5")  #c("T1-C12", "T1-C9", "T2-C1", "T2-C11")
  }
  if(organ == "THCA"){# THCA: Tumor Promotion  ?
    Precancer_Clusters[[organ]] <- data.frame(cancer = "T1-C15", precancer = "T2-C2") #or  data.frame(cancer = "T1-C1", precancer = "T1-C15") #c("T1-C1", "T1-C15", "T2-C5", "T2-C2")
  }
  
}

library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(ggpubr)
load('/data/rluo4/All/Output/pheno_all.rds')
load('/data/rluo4/All/Output/Imm_all.rds')
load('/data/rluo4/All/Output/Str_all.rds')

my_genesets <- read.csv('/data/rluo4/All/my_genesets.csv')
my_genesets <- my_genesets[-grep("KEGG", my_genesets$gs_name),-1]
# my_genesets <- rbind(Senescense, Metaplasia, H_genesets, cellcycle.gmt, reactome.gmt)
colnames(my_genesets) <- c("term","gene")
my_genesets$term <- as.factor(my_genesets$term)
unique(my_genesets$term)
my.gmt <- split(my_genesets$gene, my_genesets$term)
unique(my_genesets$term[grep('MHC', my_genesets$term)])
unique(my_genesets$term[grep('ANTIGEN', my_genesets$term)])


data_path = '/data/rluo4/All/Output'
setwd(data_path)
organ_all <- list.files('Cluster/')
# organ_all <- organ_all[organ_all!='Chen']
organ_all
cell_all <- c(pheno_all, Imm_all, Str_all)
pheno <- paste(organ_all, 'Epi', sep = '_')
pheno <- c(pheno, 'Chen_Epi')
df.enrich_ADJ <- NULL
cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
for (organ in organ_all) {
  # organ = 'Breast'
  epi <- paste0(organ, '_Epi')
  print(organ)
  epi <- gsub('Chen','CRC',paste(organ, 'Epi', sep = "_"))
  Epi <- cell_all[[epi]]
  Epi$barcode = rownames(Epi)
  Epi$cell_type <- Epi$Cell_Type
  table(Epi$Organ, Epi$Tissue)
  Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
  Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                              gsub('Tumor','PRAD', Epi$Tissue)))
  Epi$Major_type <- 'Epithelial'
  print(table(Epi$Organ))
  print(table(Epi$Tissue))
  
  load(paste0(cohort_directory, organ,'_ADJ_ssgsea.rds'))
  dim(ssgsea.res)
  ssgsea.res_ADJ <- ssgsea.res
  
  load(paste0(cohort_directory, organ,'_PCT_ssgsea.rds'))
  dim(ssgsea.res)
  
  terms <- intersect(rownames(ssgsea.res), rownames(ssgsea.res_ADJ))
  ssgsea.res <- cbind(ssgsea.res_ADJ[terms, ], ssgsea.res[terms,])
  
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds')
  load(file)
  hub_cluster = Precancer_Clusters[[organ]]
  hub_cluster
  
  
  set.seed(123)
  group_ABN <- rdata_filter@meta.data
  table(group_ABN$cell_class)
  table(group_ABN$Sample_Type)
  group_ABN$cell_subtype <- as.character(group_ABN$cell_class)
  table(group_ABN$cell_subtype)
  group_ABN <- group_ABN[group_ABN$cell_subtype %in% hub_cluster,]
  group_ADJ <- Epi[Epi$Tissue %in% c('Healthy','ADJ') & Epi$Cell_Type == 'STM',]
  if(organ == 'Chen'){
    group_ADJ <- group_ADJ[group_ADJ$sample_name=='Chen',]
  }
  if(organ == 'CRC'){
    group_ADJ <- group_ADJ[group_ADJ$sample_name=='Becker',]
  }
  group_ADJ$cell_subtype <- group_ADJ$Tissue
  colN <- c('orig.ident', 'cell_subtype')
  group_all <- rbind(group_ABN[, colN], group_ADJ[, colN])
  table(group_all$cell_subtype)
  
  df.organ <- NULL
  for (cancer in hub_cluster) {
    # group_ABN <- group_ABN[group_ABN$cell_subtype %in% c(cancer, precancer),]
    # cancer <- hub_cluster
    group_ABN <- group_all[group_all$cell_subtype %in% c(cancer, c('Healthy','ADJ')),]
    library(edgeR)
    library(limma)
    print(table(group_ABN$cell_subtype))
    # print(table(group_ABN$cell_subtype %in% c(cancer, precancer)))
    # # group_ABN$cell_subtype[group_ABN$cell_subtype=='T1-C14'] = 'T1-C1'
    # precancer <- hub_cluster$precancer
    group <- as.factor(ifelse(group_ABN$cell_subtype %in% cancer,
                              'Bad','Good'))
    # subtype <- data.frame(subset = subtype)
    # subtype <- subtype[subtype$subset.Freq !=0,]
    table(group)
    desigN <- model.matrix(~group  + 0)
    head(desigN)
    # Activation of the UPR can enhance the cell's ability to cope with nutrient deprivation, hypoxia, and other stressors commonly encountered in solid tumors.
    rownames(desigN) <- rownames(group_ABN)
    head(desigN)
    comparE <- makeContrasts(groupBad - groupGood, levels=desigN)
    enrich.res <- ssgsea.res[, colnames(ssgsea.res) %in% rownames(group_ABN)]
    # enrich.res<-gsva.res[,colnames(gsva.res) %in% rownames(group_ABN)]
    fiT <- lmFit(enrich.res, desigN)
    fiT2 <- contrasts.fit(fiT, comparE)
    fiT3 <- eBayes(fiT2)
    veloGood <- topTable(fiT3, coef=1, number=200)
    head(veloGood, n=3)
    max(abs(veloGood$logFC))
    
    dat.enrich=enrich.res
    df.res=veloGood
    df.res=df.res[df.res$P.Value<0.05 ,]
    df.res$Pathways <- rownames(df.res)
    df.res$organ <- organ
    df.res$cluster[df.res$t>0] <- cancer#cancer #'Precancer'
    df.res$cluster[df.res$t<0] <- paste0('ADJ v.s ', cancer)
    
    df.organ <<- rbind(df.organ, df.res)
  }
  
  df.enrich_ADJ <<- rbind(df.enrich_ADJ, df.organ)
  
}
save(df.enrich_ADJ, file = paste0('/data/rluo4/All/Output/Enrichment_STM_ADJ.rds'))
df.enrich_ADJ$organ <- gsub('HNSCC','OralCalvity', gsub('THCA','Thyroid',df.enrich_ADJ$organ))
write.csv(df.enrich_ADJ, file ="/data/rluo4/All/Output/df.enrich_ADJ.csv", row.names = F) # Table S5
write.xlsx(df.enrich_ADJ, file = paste0('/data/rluo4/All/Output/table_xlsx/TableS5_1.xlsx'), rowNames = FALSE)#,

dup.pathways <- as.data.frame(table(df.enrich_ADJ$Pathways))
ADJ_path <- dup.pathways$Var1[dup.pathways$Freq>=8]
include_path <- c(
  as.character(unique(my_genesets$term[grep('JAK', my_genesets$term)])),
  as.character(unique(my_genesets$term[grep('TNF', my_genesets$term)])),
  as.character(unique(my_genesets$term[grep('NFK', my_genesets$term)])),
  as.character(unique(my_genesets$term[grep('IFN', my_genesets$term)])),
  as.character(unique(my_genesets$term[grep('TLR', my_genesets$term)])),
  as.character(unique(my_genesets$term[grep('MAPK', my_genesets$term)]))
)
include_path <- include_path[include_path %ni% c('REACTOME_NEGATIVE_REGULATION_OF_MAPK_PATHWAY')]
ADJ_path <- unique(c(as.character(ADJ_path), include_path))
ADJ_path

df_plot <- df.enrich_ADJ[df.enrich_ADJ$Pathways %in% include_path,]
exclude_compare <- ! grepl('ADJ',df_plot$cluster)
df_plot <- df_plot[exclude_compare, ]
dim(df_plot)
# Fig.4
library(stringr)
df_plot$id <- df_plot$Pathways
df_plot$id <- str_replace(df_plot$id , "HALLMARK_","")
df_plot$id <- str_replace(df_plot$id , "REACTOME_","")
# df_plot$id <- gsub('EpiSen','EPITHELIAL SENESCENCE',df_plot$id)
# df_plot$id <- gsub('Metaplasia','METAPLASIA AND DAMAGE RESPONSE',df_plot$id)
# df_plot$id <- gsub('TGF_BETA_RECEPTOR_SIGNALING_IN_EMT_EPITHELIAL_TO_MESENCHYMAL_TRANSITION','TGF_BETA_RECEPTOR_SIGNALING_IN_EMT',df_plot$id)
unique(df_plot$id)
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}

df_plot$id = sapply(df_plot$id ,shorten_names)
df_plot$trace <- str_split(df_plot$cluster, '-', simplify = T)[,2]
df_plot$trace <- paste(df_plot$organ, df_plot$trace, sep = '-')
df_plot$trace <- gsub('CRC','Colon',gsub('Chen','Colon', df_plot$trace))
p_dotplot <- ggplot(df_plot, aes(x = trace, y = id, size = -log10(adj.P.Val), fill = t)) +
  geom_point(data = df_plot,#df[!grepl("Cancer|HALLMARK", df$gene),], 
             pch = 21, color = "white") +
  scale_fill_distiller("GSVA\nt value",#"Fold\nchange\n(log2)",
                       palette = "PiYG", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = max(abs(as.numeric(df_plot$t))) * c(-1, 1),#max(abs(as.numeric(df$avg_log2FC[!grepl("HALLMARK", df$gene)]))) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top")) +
  geom_point(data = df_plot,#df[!grepl("Cancer|HALLMARK", df$gene),], 
             pch = 1, color = ifelse(df_plot$adj.P.Val < 0.05, "grey50", "white")) +#ifelse(df[!grepl("Cancer|HALLMARK", df$gene),]$p_val_adj < 0.05, "grey50", "white")) +
  scale_size("FDR\n(-log10)", range = c(2, 4.5), guide = guide_legend(title.position = "top",
                                                                      override.aes = list(pch = 16))) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "bottom") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_text(color = 'black',size = 11,face = "italic"),
        axis.text.x   = element_text(color = 'black', family="Arial",# face = "italic",#color = 'black',  family="Arial",# colour="red"
                                     size = 10, angle = 335),
        axis.text.y = element_text(color = 'black',size = 10, family="Arial"), #face = "italic",
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  ylab("") +
  xlab("")
p_dotplot
ggsave("/data/rluo4/All/Output/res_fig/df.enrich_ADJ.png", width = 10, height = 8, dpi = 500)
ggsave("/data/rluo4/All/Output/res_fig/df.enrich_ADJ.pdf", width = 10, height = 8, dpi = 500)
# ########################################################################
# library(ComplexHeatmap)
# library(circlize)
# ann_colors = list( Group = c("High"="#D95F02","Low"="#1B9E77"))#"#0072B5FF", "white", "#BC3C29FF"
# # ha <- HeatmapAnnotation(df.res=ac, col = ann_colors,show_legend = TRUE, annotation_legend_param = list(Group = list(direction = "horizontal"))
# # )
# head( df.res$Pathways, 100)
# 
# dat_high<-dat.enrich[pathways,group=="Dediff"]
# dat_low<-dat.enrich[pathways,group=="Diff"]
# low_mean<-apply(dat_low,1,mean)
# high_mean<-apply(dat_high,1,mean)
# 
# df1<-data.frame(group=rep("High",length(pathways)),Pathway=pathways,exp_mean=high_mean,x=rep("Low",length(pathways)),xend=rep("High",length(pathways)),
#                 y=low_mean,yend=high_mean)
# df2<-data.frame(group=rep("Low",length(pathways)),Pathway=pathways,exp_mean=low_mean,x=rep("Low",length(pathways)),xend=rep("Low",length(pathways)),
#                 y=low_mean,yend=low_mean)
# 
# Df<-rbind(df1,df2)
# Df$Pathway <- factor(Df$Pathway,levels=unique(Df$Pathway))
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer")
#   library(RColorBrewer)
# }
# library(RColorBrewer)
# library(ggplot2)
# Df$Pathway <- gsub('HALLMARK_','',Df$Pathway)
# Df$Pathway <- gsub('REACTOME_','',Df$Pathway)
# Df$Pathway <- gsub('_',' ', Df$Pathway)
# Df$Pathway  <- gsub('EpiSen','EPITHELIAL SENESCENCE', Df$Pathway)
# Df$Pathway  <- gsub('Metaplasia','METAPLASIA AND DAMAGE RESPONSE',Df$Pathway)
# 
# p <- ggplot(Df,aes(x=group,y=exp_mean)) +
#   geom_violin(fill="#C7B3E5",color="#C7B3E5",trim=F,alpha=0.2,width=0.8) +
#   geom_segment(aes(x=x,xend=xend,y=y,yend=yend,color=Pathway),size=0.5) +
#   geom_point(aes(fill=Pathway),shape=21,size=2,color="black") +
#   ylab("GSVA Enrichment Score") +xlab("Mutation number in stem-like cells") 
# 
# p+ theme(
#   legend.title=element_text(#face="italic", family="Times",# colour="red",
#     size=12),
#   legend.text=element_text(#face="italic", family="Times", #colour="red",
#     size=9.5),
#   axis.ticks=element_line(color='black'),
#   axis.text=element_text(colour="black",size=14),
#   panel.background = element_rect(fill='transparent'),
#   axis.line=element_line(color='black'),
#   axis.title=element_text(color='black',size=14),
#   plot.title = element_text(size=13,hjust = 0.5))


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) GSVA for STM cell clusters in rdata_filter objects
library(Rcpp)
#Import GSVA results
df.enrich_all <- NULL
cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
for (organ in organ_all) {
  
  load(paste0(cohort_directory, organ,'_PCT_ssgsea.rds'))
  dim(ssgsea.res)
  
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds')
  load(file)
  hub_cluster = Trace_Clusters[[organ]]
  hub_cluster
  # cancer = "T1-C1"; precancer = "T1-C15"
  cancer = hub_cluster$cancer
  precancer = hub_cluster$precancer
  
  set.seed(123)
  group_ABN <- rdata_filter@meta.data
  table(group_ABN$cell_class)
  table(group_ABN$Sample_Type)
  group_ABN$cell_subtype <- as.character(group_ABN$cell_class)
  table(group_ABN$cell_subtype)
  library(edgeR)
  library(limma)
  print(table(group_ABN$cell_subtype %in% c(cancer, precancer)))
  # group_ABN$cell_subtype[group_ABN$cell_subtype=='T1-C14'] = 'T1-C1'
  group_ABN <- group_ABN[group_ABN$cell_subtype %in% c(cancer, precancer),]
  group <- as.factor(ifelse(group_ABN$cell_subtype %in% cancer, 
                            'Bad','Good'))
  # subtype <- data.frame(subset = subtype)
  # subtype <- subtype[subtype$subset.Freq !=0,]
  table(group)
  desigN <- model.matrix(~group  + 0)
  head(desigN)
  # Activation of the UPR can enhance the cell's ability to cope with nutrient deprivation, hypoxia, and other stressors commonly encountered in solid tumors.
  rownames(desigN) <- rownames(group_ABN)
  head(desigN)
  comparE <- makeContrasts(groupBad - groupGood, levels=desigN)
  enrich.res<-ssgsea.res[,colnames(ssgsea.res) %in% rownames(group_ABN)]
  # enrich.res<-gsva.res[,colnames(gsva.res) %in% rownames(group_ABN)]
  fiT <- lmFit(enrich.res, desigN)
  fiT2 <- contrasts.fit(fiT, comparE)
  fiT3 <- eBayes(fiT2)
  veloGood <- topTable(fiT3, coef=1, number=200)
  head(veloGood, n=3)
  max(abs(veloGood$logFC))
  
  dat.enrich=enrich.res
  df.res=veloGood
  df.res=df.res[df.res$P.Value<0.05 ,] 
  df.res$Pathways <- rownames(df.res)
  df.res$organ <- organ
  df.res$cluster[df.res$t>0] <- cancer
  df.res$cluster[df.res$t<0] <- precancer
  
  df.enrich_all <<- rbind(df.enrich_all, df.res)
  
  ## barplot
  df.res=df.res[df.res$adj.P.Val<0.05 ,]  # sort according to p-value
  dat_plot <- data.frame(id = row.names(df.res),
                         t = df.res$t)
  dat_plot<-dat_plot[order(dat_plot$t,decreasing = T),]#[1:15,]
  top_up <- head(dat_plot$id,15)
  top_down <- tail(dat_plot$id,15)
  
  include_path <- c('EpiSen','SenMayo', 'FRIDMAN_SENESCENCE_UP','CancerG0Arrest','HALLMARK_HYPOXIA', 
                    'REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE','REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP',
                    'REACTOME_CELLULAR_SENESCENCE','REACTOME_ONCOGENE_INDUCED_SENESCENCE',
                    'REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE', 
                    as.character(unique(my_genesets$term[grep('MHC', my_genesets$term)])),
                    as.character(unique(my_genesets$term[grep('AUTOPHAGY', my_genesets$term)])),
                    as.character(unique(my_genesets$term[grep('TGF_BETA', my_genesets$term)]))
  )
  dat_plot <- dat_plot[dat_plot$id %in% c(top_up,top_down,include_path),]
  
  # 去掉"HALLMARK_"
  library(stringr)
  dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
  dat_plot$id <- str_replace(dat_plot$id , "REACTOME_","")
  dat_plot$id <- gsub('EpiSen','EPITHELIAL SENESCENCE',dat_plot$id)
  dat_plot$id <- gsub('Metaplasia','METAPLASIA AND DAMAGE RESPONSE',dat_plot$id)
  
  shorten_names <- function(x, n_word=4){
    if (length(strsplit(x, "_")[[1]]) > n_word )
      x <- paste(paste(strsplit(x, "_")[[1]][1:n_word], collapse="_"),
                 paste(strsplit(x, "_")[[1]][(n_word+1):length(strsplit(x,"_")[[1]])],
                       collapse="_"), sep="\n_")
    return(x)
  } 
  shorten_names <- function(x, n_word=4, n_char=50){
    if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 50))
    {
      if (nchar(x) > 50) x <- substr(x, 1, 50)
      x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")
      return(x)
    }
    else
    {
      return(x)
    }
  }
  dat_plot$term <- dat_plot$id
  dat_plot$id <- dat_plot$term
  dat_plot$id = sapply(dat_plot$id ,shorten_names)
  # dat_plot$id <- shorten_names(dat_plot$id)
  # 新增一列 根据t阈值分类
  dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
  # 排序
  dat_plot <- dat_plot %>% arrange(t)
  # 变成因子类型
  dat_plot$id <- gsub('_',' ',dat_plot$id)
  dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
  # 绘制
  library(ggplot2)
  library(ggthemes)
  # install.packages("ggprism")
  library(ggprism)
  if(organ %in% c('Breast','Chen')){
    color = c('Up'= '#f89e81','NoSignifi'='#cccccc','Down'='#7bcd7b')
  }else{
    color = c('Up'= '#C17188','NoSignifi'='#cccccc','Down'='#748EC2')
  }
  ylab = paste0('GSVA: t value in ', organ, '(', cancer, ' v.s ', precancer, ')')#'(HCC:STM:19 v.s HCC:STM:7)'
  p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
    geom_col()+
    coord_flip() +
    scale_fill_manual(values = color ) +
    # scale_fill_manual(values = c('Up'= '#C17188','NoSignifi'='#cccccc','Down'='#748EC2')) +
    geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
    xlab('') + 
    ylab(ylab) + #注意坐标轴旋转了
    guides(fill=F)+ # 不显示图例
    theme_prism(border = T) +
    theme( 
      plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
      plot.subtitle = element_text(color = 'black', size   = 13,hjust = 0.5),
      plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title  = element_text(color = 'black', size  = 16),
      legend.text   = element_text(color = 'black', size   = 16),
      # axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
      axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
      panel.border = element_rect(linetype = 'dotted', size = 1.2,fill = NA) # 图四周框起来
    )
  p
  # 添加标签 # 此处参考了：https://mp.weixin.qq.com/s/eCMwWCnjTyQvNX2wNaDYXg
  # 小于-2的数量
  low1 <- dat_plot %>% filter(t < -2) %>% nrow()
  # 小于0总数量
  low0 <- dat_plot %>% filter( t < 0) %>% nrow()
  # 小于2总数量
  high0 <- dat_plot %>% filter(t < 2) %>% nrow()
  # 总的柱子数量
  high1 <- nrow(dat_plot)
  # 依次从下到上添加标签
  p1 <- p + 
    geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'black') + # 小于-1的为黑色标签
    # geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
    #           hjust = 0,color = 'grey') + # 灰色标签
    # geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
    #           hjust = 1,color = 'grey') + # 灰色标签
    geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'black') # 大于1的为黑色标签
  p1
  ggsave(p1, file = paste0(cohort_directory, organ, "_res_fig/", cancer, '_', precancer, "_ssgsea_bar.png"),dpi = 500,
         width = 10,height  = 7)
  
}
save(df.enrich_all, file = paste0('/data/rluo4/All/Output/Enrichment_STM_all.rds'))
load(paste0('/data/rluo4/All/Output/Enrichment_STM_all.rds'))
write.csv(df.enrich_all, file ="/data/rluo4/All/Output/df.enrich_all.csv", row.names = F) # Table S5
write.xlsx(df.enrich_all, file = paste0('/data/rluo4/All/Output/table_xlsx/TableS5_2.xlsx'), rowNames = FALSE)#,
# Fig.4
load(paste0('/data/rluo4/All/Output/Enrichment_STM_all.rds'))
df.enrich_select <- df.enrich_all
df.enrich_select$organ <- gsub('HNSCC','OralCalvity', gsub('THCA','Thyroid',df.enrich_select$organ))
dup.pathways <- as.data.frame(table(df.enrich_select$Pathways))
select_path <- dup.pathways$Var1[dup.pathways$Freq>=8]
include_path <- c('EpiSen','SenMayo', 'FRIDMAN_SENESCENCE_UP','HALLMARK_HYPOXIA', #'CancerG0Arrest',
                  'REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE','REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP',
                  'REACTOME_CELLULAR_SENESCENCE','REACTOME_ONCOGENE_INDUCED_SENESCENCE',
                  'REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE', 
                  as.character(unique(my_genesets$term[grep('MHC', my_genesets$term)])),
                  as.character(unique(my_genesets$term[grep('MESENCHYMAL', my_genesets$term)])),
                  as.character(unique(my_genesets$term[grep('AUTOPHAGY', my_genesets$term)])),
                  as.character(unique(my_genesets$term[grep('TGF_BETA', my_genesets$term)]))
)
select_path <- unique(c(as.character(select_path), include_path))
select_path
df.enrich_select <- df.enrich_select[df.enrich_select$Pathways %in% include_path,]
df.enrich_select <- df.enrich_select[df.enrich_select$Pathways %ni% c('DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING','REACTOME_ANTIGEN_PRESENTATION:_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC'),]
# Fig.4
library(stringr)
df.enrich_select$id <- df.enrich_select$Pathways
df.enrich_select$id <- str_replace(df.enrich_select$id , "HALLMARK_","")
df.enrich_select$id <- str_replace(df.enrich_select$id , "REACTOME_","")
df.enrich_select$id <- gsub('EpiSen','EPITHELIAL SENESCENCE',df.enrich_select$id)
df.enrich_select$id <- gsub('Metaplasia','METAPLASIA AND DAMAGE RESPONSE',df.enrich_select$id)
df.enrich_select$id <- gsub('TGF_BETA_RECEPTOR_SIGNALING_IN_EMT_EPITHELIAL_TO_MESENCHYMAL_TRANSITION','TGF_BETA_RECEPTOR_SIGNALING_IN_EMT',df.enrich_select$id)
unique(df.enrich_select$id)
df.enrich_select$trace <- str_split(df.enrich_select$cluster, '-', simplify = T)[,2]
df.enrich_select$trace <- paste(df.enrich_select$organ, df.enrich_select$trace, sep = '-')
df.enrich_select$trace <- gsub('CRC','Colon',gsub('Chen','Colon', df.enrich_select$trace))

p_dotplot <- ggplot(df.enrich_select, aes(x = trace, y = id, size = -log10(adj.P.Val), fill = t)) +
  geom_point(data = df.enrich_select,#df[!grepl("Cancer|HALLMARK", df$gene),], 
             pch = 21, color = "white") +
  scale_fill_distiller("GSVA\nt value",#"Fold\nchange\n(log2)",
                       palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = max(abs(as.numeric(df.enrich_select$t))) * c(-1, 1),#max(abs(as.numeric(df$avg_log2FC[!grepl("HALLMARK", df$gene)]))) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top")) +
  geom_point(data = df.enrich_select,#df[!grepl("Cancer|HALLMARK", df$gene),], 
             pch = 1, color = ifelse(df.enrich_select$adj.P.Val < 0.05, "grey50", "white")) +#ifelse(df[!grepl("Cancer|HALLMARK", df$gene),]$p_val_adj < 0.05, "grey50", "white")) +
  scale_size("FDR\n(-log10)", range = c(2, 4.5), guide = guide_legend(title.position = "top",
                                                                      override.aes = list(pch = 16))) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "bottom") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x   = element_text(color = 'black', family="Arial",# face = "italic",#color = 'black',  family="Arial",# colour="red"
                                     size = 10, angle = 335),
        # axis.text.x = element_blank(),
        axis.text.y = element_text(color = 'black',size = 10, family="Arial"), #face = "italic",
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  ylab("") +
  xlab("")
p_dotplot
# ggsave("/data/rluo4/All/Output/res_fig/df.enrich_select.pdf", width = 10, height = 6)
ggsave("/data/rluo4/All/Output/res_fig/df.enrich_select.png", width = 10, height = 5)

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) # load the epithelial phenotypes and the MPs
load('/data/rluo4/All/Output/pheno_all.rds')
library(Seurat)
data_path <- '/data/rluo4/All/Output/Cluster/'
load(paste0(data_path,'../MP_new.rds'))
setwd(data_path)
organ_all <- list.files(data_path)
module_GO = vector("list", length(organ_all))
names(module_GO) <- organ_all
########################################################################
for (organ in organ_all) {
  print(organ)
  if(organ == 'Breast'){#: 11
    module_GO[[organ]] <- c("Stress","MES","Secreted","Unfolded-protein-response",
                            "Interferon/MHC-II","Proteasomal-degradation","Cellcycle-G2M","Epithelial-senescence",#7:"Cellcycle-G2M/G0"
                            "Unknown1","Hypoxia", "Respiration")
  }
  if(organ == 'Cervix'){#: 11
    module_GO[[organ]] <- c("Stress","Immune-response/Metabolism","MES","EMT-related",
                            "Interferon/MHC-II","Epithelial-senescence","Cellcycle-G2M","Respiration",##7:"Cellcycle-G2M/G0"
                            "Hypoxia","Immune-activation/Secreted","Hemato-related")
  }
  if(organ == 'Chen'){#: 9
    module_GO[[organ]] <- c("MYC",'Epithelial-senescence',"Metal-response","Secreted",#1: "Translation-initiation",
                            "WNT","Cellcycle-G2M","Translation-initiation", "Unknown1",#5: CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION
                            "Respiration")
  }
  if(organ == 'CRC'){# which is Becker: 5
    module_GO[[organ]] <- c("MES",'Hypoxia',"Respiration","MYC",#1:WNT
                            "WNT")#Proteasomal-degradation
  }
  if(organ == 'Endometrium'){#12
    module_GO[[organ]] <- c("Proteasomal-degradation","Interferon/MHC-II","Cilia1","Cilia2",#cilium
                            "Immune-activation/Secreted","Respiration","Cellcycle-G2M","Stress",
                            "MYC", "Hypoxia","Secreted","EMT-related")#2: EpiSen 3: TUBULIN_BINDING 11:HORMONE
  }
  if(organ == 'Esophagus'){#: 14
    module_GO[[organ]] <- c('MES',"Unknown","MYC","EMT-related",#1: ECM/Angiogenesis 3: EpiSen
                            "Hypoxia","Stress",'Epithelial-senescence','Metabolism',#5: Interferon/MHC-II
                            "Respiration","Immune-activation/Secreted", "Cellcycle-G2M","Virus-defense-response",#10: Interferon 11: Protein-maturation
                            "Immune-response/Metabolism","Interferon/MHC-I") #12: Interferon; MP_12_CXCL10_ETV7_GBP4_GBP5_IDO1_LAP3_TAP1_UBD_WARS
  }
  if(organ == 'GC'){#: 7
    module_GO[[organ]] <- c("Interferon/MHC-II","Epithelial-senescence","Stress","Cellcycle-G2M", #4:"Cellcycle-G2M/G0"
                            "MYC","Translation-initiation","Metal-response")
  }
  if(organ == 'HNSCC'){#: 17
    module_GO[[organ]] <- c("Immune-response/Metabolism","EMT-related1", "EMT-related2","Cellcycle-G1S",
                            "Secreted","Immune-activation/Secreted","Cellcycle-G2M","Hypoxia",#5: B-cell 6: T/B-cell 7: "Cellcycle-G2M/G0"
                            "MES", "Interferon/MHC-I",  "Epithelial-senescence", "Stress",
                            "Virus-defense-response","EMT-related3",'Respiration',"Interferon/MHC-II",'MYC')#16; MHC-II
  }
  if(organ=='Liver'){#: 10
    module_GO[[organ]] <- c("Platelet-activation","EMT-related","Interferon/MHC-II","Metabolism",
                            "Cellcycle-G2M", 'Immune-activation/Secreted','Stress',"MYC",#5: "Cellcycle-G2M/G0" 8:"Translation-initiation"
                            "Respiration","Virus-defense-response")
  }  
  if(organ=='Lung'){#: 5
    module_GO[[organ]] <- c( "MYC", 'Respiration', "Proteasomal-degradation",'Cellcycle-G1S',
                             "Stress",'Secreted')
  }
  if(organ == 'Pancreas'){#: 5
    module_GO[[organ]] <- c( "Pancreas-classical", 'MYC', "Proteasomal-degradation","Hemato-related")#2: "Translation-initiation"
  }
  if(organ=="Prostate"){#: 10
    module_GO[[organ]] <- c("Cellcycle-G0Arrest","Stress","EMT-related1","EMT-related2",#1: "Respiration"
                            "Epithelial-senescence", "Unfolded-protein-response", "Neuronal", "MES",
                            "Interferon/MHC-II", "MYC") #10: "Translation-initiation"
  }
  if(organ=="Skin"){#: 7
    module_GO[[organ]] <- c("Epithelial-senescence","MES","Stress","Cellcycle-G2M",#2: "Wnt-signaling" 4:"CancerG0Arrest"
                            "Metal-response","Interferon/MHC-II","Respiration")
  } 
  
  if(organ=="THCA"){#: 8
    module_GO[[organ]] <- c("Cellcycle-G2M","Epithelial-senescence","Secreted","Cellcycle-G0Arrest",#HEMATOPOIETIC_STEM_CELL_PROLIFERATION, NEGATIVE_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS
                            "Translation-initiation","Interferon/MHC-II", "Stress", "Metabolism")
  }
  
}
MPs <- unique(unlist(module_GO))
MPs <- MPs[-grep("Unknown",MPs)]
MPs#31
MPs <- unique(gsub('EMT-related1','EMT-related',MPs))
MPs#30
load(paste0(data_path,"../Cor_MP_Malignancy_new.rds"))
Cor_MP_Malignancy[["HNSCC"]]$mRNAs <- gsub('EMT-related1','EMT-related',Cor_MP_Malignancy[["HNSCC"]]$mRNAs)
Cor_MP_Malignancy[["Prostate"]]$mRNAs <- gsub('EMT-related1','EMT-related',Cor_MP_Malignancy[["Prostate"]]$mRNAs)
MP_Tissue <- data.frame(MP=c(MPs,'CytoTRACE'))
n <- organ_all
p <- lapply(setNames(n, n), function(nameindex){
  Cor_table <-  Cor_MP_Malignancy[[nameindex]]
  # Cor_table <- Cor_table[Cor_table$p.value <0.05, ]
  Tissue<- Cor_table$cor[match(MP_Tissue$MP,Cor_table$mRNAs)]
  # MP_Tissue[,nameindex] <- Tissue
  return(Tissue)
})
Cor_table <- data.frame(p)
rownames(Cor_table) <- MP_Tissue$MP
# write.csv(Cor_table,file = paste0(data_path,'../Cor_MP_Malignancy_new.csv'),sep = ',',quote = F,row.names = T,col.names = T)
Cor_table_Malignancy <- Cor_table
load(paste0(data_path,"../Cor_MP_CytoTRACE_new.rds"))
Cor_MP_CytoTRACE[["HNSCC"]]$mRNAs <- gsub('EMT-related1','EMT-related',Cor_MP_CytoTRACE[["HNSCC"]]$mRNAs)
Cor_MP_CytoTRACE[["Prostate"]]$mRNAs <- gsub('EMT-related1','EMT-related',Cor_MP_CytoTRACE[["Prostate"]]$mRNAs)
MP_Tissue <- data.frame(MP=c(MPs,'Malignancy'))
n <- organ_all
p <- lapply(setNames(n, n), function(nameindex){
  Cor_table <-  Cor_MP_CytoTRACE[[nameindex]]
  # Cor_table <- Cor_table[Cor_table$p.value <0.05, ]
  Tissue<- Cor_table$cor[match(MP_Tissue$MP,Cor_table$mRNAs)]
  # MP_Tissue[,nameindex] <- Tissue
  return(Tissue)
})
Cor_table <- data.frame(p)
rownames(Cor_table) <- MP_Tissue$MP
# write.csv(Cor_table,file = paste0(data_path,'../Cor_MP_CytoTRACE_new.csv'),sep = ',',quote = F,row.names = T,col.names = T)
Cor_table_CytoTRACE <- Cor_table

Cor_table <- Cor_table_CytoTRACE#[-nrow(Cor_table_CytoTRACE),]*Cor_table_Malignancy[-nrow(Cor_table_Malignancy),]
#Cor_table <- Cor_table[,-2]
#Cor_table <- Cor_table_Malignancy
library(ggplot2)
library(RColorBrewer)
coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

for(organ in organ_all){
  # organ = 'Liver'
  Lab<-data.frame(
    MP=rownames(Cor_table_CytoTRACE)[1:30],#c("biology" , "english" ,"math" ,  "music" , "R-coding" ),
    id=c(1:30) ,
    Stemness=Cor_table_CytoTRACE[[organ]][1:30],
    Malignancy=Cor_table_Malignancy[[organ]][1:30]
  )
  index <- is.na(Lab$Stemness) & is.na(Lab$Malignancy)
  label_data <- Lab[! index,]
  label_data$id <- 1:nrow(label_data)
  AddRow<-c(NA,nrow(label_data)+1,label_data[1,ncol(label_data)-1],label_data[1,ncol(label_data)])
  mydata<-rbind(label_data, AddRow)
  
  myAngle<- 360- 360 * (label_data$id-1) /nrow(label_data)  
  
  mydata<-melt(mydata,id=c("MP", "id"))
  
  ggplot(data=mydata,aes(x=id, y=value,group=variable,fill=variable)) + 
    geom_polygon(colour="black",alpha=0.1)+
    geom_point(size=4,shape=21,color = 'black')+
    coord_radar()+
    #coord_polar() +
    scale_x_continuous(breaks =label_data$id,labels=label_data$MP)+
    theme_bw() +
    ylim(0,1)+
    theme(axis.text.x=element_text(size = 11,colour="black",angle = myAngle),
          axis.title=element_text(size=15,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          panel.grid.major = element_line(color="grey80"),
          axis.line = element_line(color="black"),
          axis.ticks =  element_line(color="black"))
  
  library(tidyverse, lib.loc = "/usr/local/lib/R/site-library")
  library(ggradar)   #要用remotes::install_github("ricardo-bion/ggradar")安装
  
  data_radar <- as.data.frame( t(label_data[,-(1:2)]) )
  colnames(data_radar) <- label_data$MP
  #data_radar$group <- rownames(data_radar)
  data_radar <- abs(data_radar)
  data_radar <- data_radar[, sort(colnames(data_radar), decreasing = F)]
  # ggradar(
  #   data_radar[1, ], 
  #   values.radar = c("0", "0.3", "1"),
  #   grid.min = 0, grid.mid = 0.3, grid.max = 1,
  #   # Polygons
  #   group.line.width = 1, 
  #   group.point.size = 3,
  #   group.colours = "#00AFBB",
  #   # Background and grid lines
  #   background.circle.colour = "white",
  #   gridline.mid.colour = "grey"
  # )
  df <- data_radar %>% rownames_to_column("group")
  # ggradar(
  #   df, plot.extent.x.sf = 1.4,
  #   values.radar = c("0", "0.3", "1"),
  #   grid.min = 0, grid.mid = 0.3, grid.max = 1, legend.title = "",#"Cell feature",
  #   # Polygons
  #   group.line.width = 1, 
  #   group.point.size = 3,
  #   group.colours = c("#00AFBB", "#E7B800"),#, "#FC4E07"),
  #   # Background and grid lines
  #   background.circle.colour = "white",
  #   gridline.mid.colour = "#FC4E07",#"grey",
  #   legend.position = "none",#"bottom",
  #   legend.text.size = 15
  # )+ 
  #   # theme(
  #   #   # legend.position = c(1, 0),
  #   #   # legend.justification = c(1, 0),
  #   #   legend.text = element_text(size = 12),
  #   #   legend.key = element_rect(fill = NA, color = NA),
  #   #   legend.background = element_blank()
  #   # )+
  #   # labs(title = "Correlation of MPs with liver cell features") + 
  #   theme( axis.text=element_text(size = 11,colour="black",angle = myAngle),
  #          # plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
  #          # panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
  #          # plot.title.position = "plot", # slightly different from default
  #          plot.title = element_text(
  #            size = 20,
  #            face = "bold",
  #            #color = "#2a475e"
  #          )
  #   )
  ggradar(df,base.size = 12,
          plot.extent.x.sf = 2,
          plot.extent.y.sf = 1.2,
          x.centre.range = 0.02 * (grid.max - centre.y),
          label.centre.y = FALSE,
          values.radar = c("0", "0.3", "1"),
          grid.min = 0, grid.mid = 0.3, grid.max = 1,
          legend.title = "Cell feature", 
          legend.text.size = 12,
          legend.position = "none") + theme(
            legend.title=element_text(#face="italic", family="Times",# colour="red",
              size=12),
            legend.text=element_text(#face="italic", family="Times", #colour="red",
              size=15),
            # axis.ticks=element_line(color='black'),
            # axis.text=element_text(colour="black",size=14),
            panel.background = element_rect(fill='transparent'),
            # axis.line=element_line(color='black'),
            # axis.title=element_text(color='black',size=14),
            plot.title = element_text(size=16,hjust = 0.5))
  outdir = '/data/rluo4/All/Output/res_fig/radar_fig/'
  
  f <- paste0(outdir, organ, '_Cor_table_radar_plot.png')
  # ggsave(filename = f, height = 5, width = 7.2, dpi = 500, device = "png")#, useDingbats=FALSE)
  ggsave(filename = f, height = 5, width = 8.8, dpi = 500, device = "png")#, useDingbats=FALSE)
  
}

#############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) Load the results of intersected genes from ts860's final steps & regenerate DEGs from the pheno-related clusters
Patt <- read.table('/data/rluo4/All/Output/organ13-deg.txt',sep = '\t',fill = TRUE,header = TRUE)
data_path <- '/data/rluo4/All/Output/Cluster/'
setwd(data_path)
load(paste0(data_path,'../MP_new.rds'))
organ_all <- list.files(data_path)
# print(t)
for (organ in organ_all) {
  print(organ)
  # save(Trans_G,DEG_list, Sig_G, Dat_sig, rdata_filter, file =paste0(cohort_directory, 'Epi/Results/','Malignant-Transformation-',organ, '.rds' ))
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds')
  load(file)
  ##############################################################################################################################
  #............................................................................................................................#
  ##############################################################################################################################
  # 0)---load results of AUCell_MP/new.R---
  analysis_parent_folder <- '/data/rluo4/All/Output/cor_auc_MPs/'
  if (!dir.exists(paste0(analysis_parent_folder))){
    dir.create(paste0(analysis_parent_folder))
  }
  setwd(analysis_parent_folder)
  getwd()
  range = 4:12 # 迭代次数
  gmin = 5 #最低关联数量
  ncores = 5 #运行核数量
  load(paste0(organ,'_cor_auc_MPs_new.rds'))#save(EpithelialCell_obj, AUC, file=paste0(organ,'_cor_auc_MPs_new.rds'))
  
  table(colnames(EpithelialCell_obj) %in% colnames(AUC))
  table(colnames(AUC) %in% colnames(rdata_filter))
  table(colnames(EpithelialCell_obj) %in% colnames(rdata_filter))
  
  AUC <- AUC[, colnames(AUC) %in% colnames(rdata_filter)]
  library(ggpubr)
  library(ggsignif)
  cell.data <- as.data.frame(t(AUC))
  modules <- rownames(AUC)
  print(modules)
  ##############################################################################################################################
  #............................................................................................................................#
  ##############################################################################################################################
  # 1) Load the results of DEGs against healthy tissues
  # data_path <- '/data/rluo4/All/Output/Cluster/'
  setwd(data_path)
  # load(paste0(data_path,'../MP_new.rds'))
  Dat_sig <- vector("list", length(modules))
  names(Dat_sig) <- modules
  
  mt_genes <- grep("^MT-", Patt$Symbol, value = TRUE)
  rps_genes <- grep("^RPS", Patt$Symbol, value = TRUE)
  rpl_genes <- grep("^RPL", Patt$Symbol, value = TRUE)#grep("^RP[SL]",rownames(sce))
  # hla_genes <- grep("^HLA-", rownames(colon), value = TRUE)
  blacklist <- c(mt_genes, rps_genes, rpl_genes)#, hla_genes)
  Patt_all <- Patt
  blacklist <- unique(blacklist)
  
  for (i in modules){
    # i = modules[1]
    data_path <- '/data/rluo4/All/Output/Cluster/'
    setwd(data_path)
    organ_all <- list.files(data_path)
    Dat_MP <-  vector("list", length(organ_all))
    names(Dat_MP) <- organ_all
    # for(organ in organ_all){
    # organ = 'Breast'
    indir <- '/data/rluo4/All/Output/sdata_ABN/'
    pc_file = paste0(indir, organ,'_DEG.RData')
    if(organ !='Chen'){
      load(pc_file)
      if(organ=='CRC'){
        pc_df <-  pc_df_Becker
        Patt <- Patt_Becker
      }
    } else{
      pc_file <- gsub('Chen','CRC',pc_file)
      load(pc_file)
      pc_df <- pc_df_Chen
      Patt <- Patt_Chen
    }
    unique(pc_df$Sample)
    pc_sample <- pc_df$SimplifiedSampleName
    table(pc_df$Tissue)
    head(Patt)
    
    analysis_parent_folder <- paste0(data_path,organ)
    if (!dir.exists(paste0(analysis_parent_folder))){
      dir.create(paste0(analysis_parent_folder))
    }
    setwd(analysis_parent_folder)
    # load(paste0(organ,'_cor_auc_MPs.rds'))
    dd_auc <- MP_all[[i]][[organ]]
    Gx <- dd_auc[dd_auc$p.value <0.05 & ! is.na(dd_auc$p.value),]# dd_auc$cor >0 & 
    Gx <- unique(Gx$mRNAs)
    Gy <- Patt[Patt$Pvalue_Adj < 0.05 & Patt$Log2FC >0.25,]
    Gy <- unique(Gy$Symbol)
    Gy <- Gy[Gy %in% Patt_all$Symbol]
    G_MP <- intersect(Gx, Gy)
    G_MP <- G_MP[! G_MP %in% blacklist]
    print(paste0(i,"-related gene number for ",organ,' is ',length(G_MP)))    
    Dat_MP[[organ]] <- dd_auc[dd_auc$mRNAs %in% G_MP,]
    rownames(Dat_MP[[organ]]) <- Dat_MP[[organ]]$mRNAs
    print(Cor_table[i,])
    
    # just check the Cor with Cytotrace Score> 0:
    names(Dat_MP)
    index = which(as.numeric(Cor_table_CytoTRACE[i,])>0)
    print(paste0(i," in ",organ, ' first positively corr with Cytotrace'))
    Dat <- Dat_MP[[organ]][,-1]
    rownames(Dat) <- Dat[,1]
    #Dat <- Dat[,-1]
    print(head(Dat))
    Dat_sig[[i]] <- Dat #save common genes between MP and DEGs as the Dat_sig list
  }
  
  ##############################################################################################################################
  #............................................................................................................................#
  ##############################################################################################################################
  # 2) # DEGS between STM subclusters
  print(table(rdata_filter$cell_class, rdata_filter$Trajectory))
  Trajectory <- unique(rdata_filter$Trajectory)
  if(organ %in% c('Breast','CRC')){ # original: T1 only
    Trajectory <- c('T1','T2')
  }
  if(organ %in% c('HNSCC','Lung')){ # origincal: T2, T1 only
    Trajectory <- c('T1','T2','T3')
  }
  print(paste0(organ, ': ', Trajectory))
  ########
  cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
  DEG_list <-  vector("list", length(Trajectory))
  names(DEG_list) <- Trajectory
  for(t in Trajectory){
    if(organ == 'Breast'){#?
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C14','T1-C9'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C4','T1-C9'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds') #paste0(cohort_directory, 'Epi/Results/', organ, '_',t, '_deg.rds')
      if(t=='T2'){
        f <- gsub('T2','T1',f)
      }
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C14"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T1-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C14"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C9"]
        annotation_col=data.frame(Cluster=c(rep("T1-C14",length(a)),rep("T1-C9",length(b))))
        ann_colors = list(Type= c('T1-C14'="darkred",'T1-C9'='lightblue'))
        col = list(State = c("T1-C14" = "#E41A1C","T1-C9" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C4"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C9"]
        annotation_col=data.frame(Cluster=c(rep("T1-C4",length(a)),rep("T1-C9",length(b))))
        ann_colors = list(Type= c('T1-C4'="darkred",'T1-C9'='lightblue'))
        col = list(State = c("T1-C4" = "#E41A1C","T1-C9" = "#377EB8"))
      }
    }
    
    if(organ == 'Cervix'){
      if(t == 'T3'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        rdata_T$cell_class <- gsub('T1','T3', rdata_T$cell_class)
        
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else{
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T3-C8"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T3-C5"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C14"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        } else{# # Trajectory 2:
          marker_C = marker%>%filter(cluster == c("T2-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T2-C10"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        }
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      if(t == 'T3'){# Trajectory 3: 
        a = rownames(dat@meta.data)[dat$cell_class=="T3-C8"]
        b = rownames(dat@meta.data)[dat$cell_class=="T3-C5"]
        annotation_col=data.frame(Cluster=c(rep("T3-C8",length(a)),rep("T3-C5",length(b))))
        ann_colors = list(Type= c('T3-C8'="darkred",'T3-C5'='lightblue'))
        col = list(State = c("T3-C8" = "#E41A1C","T3-C5" = "#377EB8"))
      } else{
        if(t == 'T1'){# Trajectory 1: 
          a = rownames(dat@meta.data)[dat$cell_class=="T1-C14"]
          b = rownames(dat@meta.data)[dat$cell_class=="T1-C1"]
          annotation_col=data.frame(Cluster=c(rep("T1-C14",length(a)),rep("T1-C1",length(b))))
          ann_colors = list(Type= c('T1-C14'="darkred",'T1-C1'='lightblue'))
          col = list(State = c("T1-C14" = "#E41A1C","T1-C1" = "#377EB8"))
        } else{# Trajectory 2: 
          a = rownames(dat@meta.data)[dat$cell_class=="T2-C4"]
          b = rownames(dat@meta.data)[dat$cell_class=="T2-C10"]
          annotation_col=data.frame(Cluster=c(rep("T2-C4",length(a)),rep("T2-C10",length(b))))
          ann_colors = list(Type= c('T2-C4'="darkred",'T2-C10'='lightblue'))
          col = list(State = c("T2-C4" = "#E41A1C","T2-C10" = "#377EB8"))
        }
      }
    } 
    
    if(organ == 'Chen'){
      
      if(t == 'T3'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else {
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T3-C17"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T3-C23"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T1-C21"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        } else{
          marker_C = marker%>%filter(cluster == c("T2-C26"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T2-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        }
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      if(t == 'T3'){# Trajectory 3: 
        a = rownames(dat@meta.data)[dat$cell_class=="T3-C17"]
        b = rownames(dat@meta.data)[dat$cell_class=="T3-C23"]
        annotation_col=data.frame(Cluster=c(rep("T3-C17",length(a)),rep("T3-C23",length(b))))
        ann_colors = list(Type= c('T3-C17'="darkred",'T3-C23'='lightblue'))
        col = list(State = c("T3-C17" = "#E41A1C","T3-C23" = "#377EB8"))
      } else{
        if(t == 'T1'){# Trajectory 1: 
          a = rownames(dat@meta.data)[dat$cell_class=="T1-C1"]
          b = rownames(dat@meta.data)[dat$cell_class=="T1-C21"]
          annotation_col=data.frame(Cluster=c(rep("T1-C1",length(a)),rep("T1-C21",length(b))))
          ann_colors = list(Type= c('T1-C1'="darkred",'T1-C21'='lightblue'))
          col = list(State = c("T1-C1" = "#E41A1C","T1-C21" = "#377EB8"))
        } else{# Trajectory 2: 
          a = rownames(dat@meta.data)[dat$cell_class=="T2-C26"]
          b = rownames(dat@meta.data)[dat$cell_class=="T2-C3"]
          annotation_col=data.frame(Cluster=c(rep("T2-C26",length(a)),rep("T2-C3",length(b))))
          ann_colors = list(Type= c('T2-C26'="darkred",'T2-C3'='lightblue'))
          col = list(State = c("T2-C26" = "#E41A1C","T2-C3" = "#377EB8"))
        }
      }
    }
    
    if(organ == 'CRC'){
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C3','T1-C19'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C16','T1-C17'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      if(t=='T2'){
        f <- gsub('T2','T1',f)
      }
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C19"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T1-C16"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C17"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C3"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C19"]
        annotation_col=data.frame(Cluster=c(rep("T1-C3",length(a)),rep("T1-C19",length(b))))
        ann_colors = list(Type= c('T1-C3'="darkred",'T1-C19'='lightblue'))
        col = list(State = c("T1-C3" = "#E41A1C","T1-C19" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C16"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C17"]
        annotation_col=data.frame(Cluster=c(rep("T1-C16",length(a)),rep("T1-C17",length(b))))
        ann_colors = list(Type= c('T1-C16'="darkred",'T1-C17'='lightblue'))
        col = list(State = c("T1-C16" = "#E41A1C","T1-C17" = "#377EB8"))
      }
      
    }
    
    if(organ == 'Endometrium'){
      
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C14"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C12"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C24"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T2-C19"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C14"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C12"]
        annotation_col=data.frame(Cluster=c(rep("T1-C14",length(a)),rep("T1-C12",length(b))))
        ann_colors = list(Type= c('T1-C14'="darkred",'T1-C12'='lightblue'))
        col = list(State = c("T1-C14" = "#E41A1C","T1-C12" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T2-C24"]
        b = rownames(dat@meta.data)[dat$cell_class=="T2-C19"]
        annotation_col=data.frame(Cluster=c(rep("T2-C24",length(a)),rep("T2-C19",length(b))))
        ann_colors = list(Type= c('T2-C24'="darkred",'T2-C19'='lightblue'))
        col = list(State = c("T2-C24" = "#E41A1C","T2-C19" = "#377EB8"))
      }
      
    }
    
    if(organ == 'Esophagus'){
      
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, RNA_snn_res.1 %in% rdata_filter$RNA_snn_res.1[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, RNA_snn_res.1 %in% rdata_filter$RNA_snn_res.1[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C10"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C12"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } 
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C10"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C12"]
        annotation_col=data.frame(Cluster=c(rep("T1-C10",length(a)),rep("T1-C12",length(b))))
        ann_colors = list(Type= c('T1-C10'="darkred",'T1-C12'='lightblue'))
        col = list(State = c("T1-C10" = "#E41A1C","T1-C12" = "#377EB8"))
      } 
      
    }
    
    if(organ == 'GC'){
      
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C17"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C2"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T2-C7"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      # a = rownames(dat@meta.data)[dat$cell_class=="C16"]
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C17"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C2"]
        annotation_col=data.frame(Cluster=c(rep("T1-C17",length(a)),rep("T1-C2",length(b))))
        ann_colors = list(Type= c('T1-C17'="darkred",'T1-C2'='lightblue'))
        col = list(State = c("T1-C17" = "#E41A1C","T1-C2" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T2-C3"]
        b = rownames(dat@meta.data)[dat$cell_class=="T2-C7"]
        annotation_col=data.frame(Cluster=c(rep("T2-C3",length(a)),rep("T2-C7",length(b))))
        ann_colors = list(Type= c('T2-C3'="darkred",'T2-C7'='lightblue'))
        col = list(State = c("T2-C3" = "#E41A1C","T2-C7" = "#377EB8"))
      }
      
    }
    
    if(organ == 'HNSCC'){
      if(t == 'T3'){
        rdata_filter$Trajectory[rdata_filter$cell_class %in% c('T1-C18','T1-C16')] <- 'T3'
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else {
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      if(t=='T3'){
        f <- gsub('T3','T1',f)
      }
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T1-C16"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C18"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T1-C11"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        } else{
          marker_C = marker%>%filter(cluster == c("T2-C19"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T2-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        }
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      if(t == 'T3'){# Trajectory 3: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C16"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C18"]
        annotation_col=data.frame(Cluster=c(rep("T1-C16",length(a)),rep("T1-C18",length(b))))
        ann_colors = list(Type= c('T1-C16'="darkred",'T1-C18'='lightblue'))
        col = list(State = c("T1-C16" = "#E41A1C","T1-C18" = "#377EB8"))
      } else{
        if(t == 'T1'){# Trajectory 1: 
          a = rownames(dat@meta.data)[dat$cell_class=="T1-C4"]
          b = rownames(dat@meta.data)[dat$cell_class=="T1-C11"]
          annotation_col=data.frame(Cluster=c(rep("T1-C4",length(a)),rep("T1-C11",length(b))))
          ann_colors = list(Type= c('T1-C4'="darkred",'T1-C11'='lightblue'))
          col = list(State = c("T1-C4" = "#E41A1C","T1-C11" = "#377EB8"))
        } else{# Trajectory 2: 
          a = rownames(dat@meta.data)[dat$cell_class=="T2-C19"]
          b = rownames(dat@meta.data)[dat$cell_class=="T2-C3"]
          annotation_col=data.frame(Cluster=c(rep("T2-C19",length(a)),rep("T2-C3",length(b))))
          ann_colors = list(Type= c('T2-C19'="darkred",'T2-C3'='lightblue'))
          col = list(State = c("T2-C19" = "#E41A1C","T2-C3" = "#377EB8"))
        }
      }
    }
    
    if(organ == 'Liver'){
      
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds') #paste0(cohort_directory, 'Epi/Results/', organ, '_',t, '_deg.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C0"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C7"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T2-C15"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C0"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C1"]
        annotation_col=data.frame(Cluster=c(rep("T1-C0",length(a)),rep("T1-C1",length(b))))
        ann_colors = list(Type= c('T1-C0'="darkred",'T1-C1'='lightblue'))
        col = list(State = c("T1-C0" = "#E41A1C","T1-C1" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T2-C7"]
        b = rownames(dat@meta.data)[dat$cell_class=="T2-C15"]
        annotation_col=data.frame(Cluster=c(rep("T2-C7",length(a)),rep("T2-C15",length(b))))
        ann_colors = list(Type= c('T2-C7'="darkred",'T2-C15'='lightblue'))
        col = list(State = c("T2-C7" = "#E41A1C","T2-C15" = "#377EB8"))
      }
      
    }
    
    if(organ == 'Lung'){
      if(t == 'T3'){
        rdata_filter$Trajectory[rdata_filter$cell_class %in% c('T1-C8','T1-C9')] <- 'T3'
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C14','T1-C6'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else{
          rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C7','T1-C3'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
      
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      if(t=='T3'){
        f <- gsub('T3','T1',f)
      }
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T1-C8"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C14"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T1-C6"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        } else{
          marker_C = marker%>%filter(cluster == c("T1-C7"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
          marker_P = marker%>%filter(cluster==c("T1-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        }
      }  
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T3'){# Trajectory 3: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C8"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C9"]
        annotation_col=data.frame(Cluster=c(rep("T1-C8",length(a)),rep("T1-C9",length(b))))
        ann_colors = list(Type= c('T1-C8'="darkred",'T1-C9'='lightblue'))
        col = list(State = c("T1-C8" = "#E41A1C","T1-C9" = "#377EB8"))
      } else{
        if(t == 'T1'){# Trajectory 1: 
          a = rownames(dat@meta.data)[dat$cell_class=="T1-C14"]
          b = rownames(dat@meta.data)[dat$cell_class=="T1-C6"]
          annotation_col=data.frame(Cluster=c(rep("T1-C14",length(a)),rep("T1-C6",length(b))))
          ann_colors = list(Type= c('T1-C14'="darkred",'T1-C6'='lightblue'))
          col = list(State = c("T1-C14" = "#E41A1C","T1-C6" = "#377EB8"))
        } else{# Trajectory 2: 
          a = rownames(dat@meta.data)[dat$cell_class=="T1-C7"]
          b = rownames(dat@meta.data)[dat$cell_class=="T1-C3"]
          annotation_col=data.frame(Cluster=c(rep("T1-C7",length(a)),rep("T1-C3",length(b))))
          ann_colors = list(Type= c('T1-C7'="darkred",'T1-C3'='lightblue'))
          col = list(State = c("T1-C7" = "#E41A1C","T1-C3" = "#377EB8"))
        }
      }
    }
    
    if(organ == 'Pancreas'){
      
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res1.5 %in% rdata_filter$leiden_res1.5[rdata_filter$Trajectory==t])#
      } else{
        rdata_T <- subset(rdata_filter, leiden_res1.5 %in% rdata_filter$leiden_res1.5[rdata_filter$Trajectory==t])#
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C8"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C2"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T2-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C8"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C2"]
        annotation_col=data.frame(Cluster=c(rep("T1-C2",length(a)),rep("T1-C8",length(b))))
        ann_colors = list(Type= c('T1-C2'="darkred",'T1-C8'='lightblue'))
        col = list(State = c("T1-C2" = "#E41A1C","T1-C8" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T2-C4"]
        b = rownames(dat@meta.data)[dat$cell_class=="T2-C9"]
        annotation_col=data.frame(Cluster=c(rep("T2-C4",length(a)),rep("T2-C9",length(b))))
        ann_colors = list(Type= c('T2-C4'="darkred",'T2-C9'='lightblue'))
        col = list(State = c("T2-C4" = "#E41A1C","T2-C9" = "#377EB8"))
      }
      
    }
    
    if(organ == 'Prostate'){
      
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T2'){# Trajectory 2:
        marker_C = marker%>%filter(cluster==c("T2-C34"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T2-C16"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T1-C5"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C42"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T2'){# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T2-C34"]
        b = rownames(dat@meta.data)[dat$cell_class=="T2-C16"]
        annotation_col=data.frame(Cluster=c(rep("T2-C34",length(a)),rep("T2-C16",length(b))))
        ann_colors = list(Type= c('T2-C34'="darkred",'T2-C16'='lightblue'))
        col = list(State = c("T2-C34" = "#E41A1C","T2-C16" = "#377EB8"))
      } else{# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C5"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C42"]
        annotation_col=data.frame(Cluster=c(rep("T1-C5",length(a)),rep("T1-C42",length(b))))
        ann_colors = list(Type= c('T1-C5'="darkred",'T1-C42'='lightblue'))
        col = list(State = c("T1-C5" = "#E41A1C","T1-C42" = "#377EB8"))
      }
      
    }
    
    if(organ == 'Skin'){
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, cell_class %in% c('T2-C1','T1-C11'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C9','T1-C11'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T2'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T2-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C11"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C11"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T2'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T2-C1"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C11"]
        annotation_col=data.frame(Cluster=c(rep("T2-C1",length(a)),rep("T1-C11",length(b))))
        ann_colors = list(Type= c('T2-C1'="darkred",'T1-C11'='lightblue'))
        col = list(State = c("T2-C1" = "#E41A1C","T1-C11" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C9"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C11"]
        annotation_col=data.frame(Cluster=c(rep("T1-C9",length(a)),rep("T1-C11",length(b))))
        ann_colors = list(Type= c('T1-C9'="darkred",'T1-C11'='lightblue'))
        col = list(State = c("T1-C9" = "#E41A1C","T1-C11" = "#377EB8"))
      }
    }
    
    if(organ == 'THCA'){
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
      print(table(rdata_T$cell_class, rdata_T$Trajectory))
      options(stringsAsFactors = F)
      dat = rdata_T
      DefaultAssay(dat)="RNA"
      table(dat$cell_class)
      dat$cell_class <- gsub('T1-C5','T2-C5',dat$cell_class) #
      Idents(dat) = factor(dat$cell_class, levels = unique(dat$cell_class))
      f <- paste0(cohort_directory, organ, '_',t, '_deg.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T1-C15"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C5"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
        marker_P = marker%>%filter(cluster==c("T2-C2"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.5)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
      if(t == 'T1'){# Trajectory 1: 
        a = rownames(dat@meta.data)[dat$cell_class=="T1-C1"]
        b = rownames(dat@meta.data)[dat$cell_class=="T1-C15"]
        annotation_col=data.frame(Cluster=c(rep("T1-C1",length(a)),rep("T1-C15",length(b))))
        ann_colors = list(Type= c('T1-C1'="darkred",'T1-C15'='lightblue'))
        col = list(State = c("T1-C1" = "#E41A1C","T1-C15" = "#377EB8"))
      } else{# Trajectory 2: 
        a = rownames(dat@meta.data)[dat$cell_class=="T2-C5"]
        b = rownames(dat@meta.data)[dat$cell_class=="T2-C2"]
        annotation_col=data.frame(Cluster=c(rep("T2-C5",length(a)),rep("T2-C2",length(b))))
        ann_colors = list(Type= c('T2-C5'="darkred",'T2-C2'='lightblue'))
        col = list(State = c("T2-C5" = "#E41A1C","T2-C2" = "#377EB8"))
      }
      
    }
    
    DEG_list[[t]] <- rbind(marker_C, marker_P)
    # save up- and down-regulated genes in the DEG_list
  }
  
  ##############################################################################################################################
  #............................................................................................................................#
  ##############################################################################################################################
  # 3) Identifying transition genes
  library(Rcpp)
  library(harmony)
  library(future)#, lib.loc = "/home/lorihan/miniconda3/lib/R/library")
  library(Matrix)
  options(stringsAsFactors = FALSE)
  options(future.globals.maxSize = 120000 * 1024^2)
  library(future)
  # plan(multisession, workers=100)
  availableCores() #12 #查看几个核可用
  nbrOfWorkers() #4 当前可用的核有多少个
  pid <- Sys.getpid()
  pid #2109512
  library(future.apply)
  cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
  
  batch_cor <- function(gene){
    y <- as.numeric(Dat.all[gene,])
    rownames <- rownames(Dat.all)
    do.call(rbind,future_lapply(rownames, function(x){
      dd  <- cor.test(as.numeric(Dat.all[x,]),y,type='spearman')
      data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
    }))
  }
  for(t in Trajectory){
    print(t)
    if( ! organ %in% c('Breast','Cervix','Chen','CRC','Esophagus','Pancreas','HNSCC','Skin','Lung') ){
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])#
      } else{
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])#
      }
    }
    
    if(organ == 'Breast'){#?
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C14','T1-C9'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C4','T1-C9'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
    }
    if(organ == 'Cervix'){
      if(t == 'T3'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        rdata_T$cell_class <- gsub('T1','T3', rdata_T$cell_class)
        
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else{
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
    }
    if(organ == 'Chen'){
      if(t == 'T3'){
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else {
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
    }
    if(organ == 'CRC'){
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C3','T1-C19'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C16','T1-C17'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
    }
    if(organ == 'Esophagus'){
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, RNA_snn_res.1 %in% rdata_filter$RNA_snn_res.1[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, RNA_snn_res.1 %in% rdata_filter$RNA_snn_res.1[rdata_filter$Trajectory==t])
      }
    }
    if(organ == 'Lung'){#?
      if(t == 'T3'){
        rdata_filter$Trajectory[rdata_filter$cell_class %in% c('T1-C8','T1-C9')] <- 'T3'
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C14','T1-C6'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else{
          rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C7','T1-C3'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
    }
    if(organ == 'Pancreas'){
      if(t == 'T1'){
        rdata_T <- subset(rdata_filter, leiden_res1.5 %in% rdata_filter$leiden_res1.5[rdata_filter$Trajectory==t])#
      } else{
        rdata_T <- subset(rdata_filter, leiden_res1.5 %in% rdata_filter$leiden_res1.5[rdata_filter$Trajectory==t])#
      }
    }
    if(organ == 'HNSCC'){
      if(t == 'T3'){
        rdata_filter$Trajectory[rdata_filter$cell_class %in% c('T1-C18','T1-C16')] <- 'T3'
        rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        if(t == 'T1'){
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        } else {
          rdata_T <- subset(rdata_filter, leiden_res2 %in% rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
        }
      }
    }
    if(organ == 'Skin'){
      if(t == 'T2'){
        rdata_T <- subset(rdata_filter, cell_class %in% c('T2-C1','T1-C11'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      } else{
        rdata_T <- subset(rdata_filter, cell_class %in% c('T1-C9','T1-C11'))#rdata_filter$leiden_res2[rdata_filter$Trajectory==t])
      }
    }  
    print(table(rdata_T$cell_class, rdata_T$Trajectory))
    expression <- GetAssayData(object = rdata_T, slot = "data")#EpithelialCell_obj@assays$RNA@data
    expression[1:5,1:5]
    table(rowSums(expression)>0)
    expression <- as.matrix(expression)
    dim(expression)
    #   if(length(unique(rownames(AUC)==i))>1)  {
    f <-  paste0(cohort_directory, organ, '_cor_',t,'_pseudotime.rds')
    print(f)
    # f <- gsub('/','_',f)
    if (!file.exists(f)){
      print(paste0("add Cor for pseudotime: ", organ))
      pseudo <- rdata_T$latent_time
      Dat.all<-as.matrix(rbind(pseudo,expression))
      str(Dat.all)
      getwd()
      Dat.all[1:5,1:5]
      set.seed(211)
      system.time(dd_auc <- batch_cor('pseudo'))
      save(dd_auc, file=f)
    }
  }
  
  # Transition genes:
  Trans <- Cor_table_CytoTRACE[organ]
  str(Trans)
  Trans <- na.omit(Trans)
  Trans_G <-   vector("list", length(Trajectory))
  names(Trans_G) <- Trajectory
  Sig_G <- vector("list", length(Trajectory))
  names(Sig_G) <- Trajectory
  for(t in Trajectory){
    print(t)
    Trans_G[[t]] <- Dat_sig  #Assign MP-genes to a Trans_G list
    cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
    f <-  paste0(cohort_directory, organ, '_cor_',t,'_pseudotime.rds')
    load(f)
    table(is.na(dd_auc$p.value))
    sig_T <- dd_auc[dd_auc$p.value <0.05 & ! is.na(dd_auc$p.value),]
    sig_T <- sig_T[sig_T$mRNAs !='pseudo',]
    sig_T <- sig_T[! sig_T$mRNAs %in% blacklist,]
    sig_T <- sig_T[order(abs(sig_T$cor),decreasing = T),]
    dim(sig_T)
    Num <- round(nrow(sig_T)*0.1)
    sig_T <- head(sig_T, Num)
    inter <- intersect(sig_T$mRNAs, DEG_list[[t]]$gene) #common genes for Transition-genes and DEGs between precancer and cancers
    print(table(modules %in% rownames(Trans)))
    sig_GL <- NULL
    for (i in modules){
      sect <- Dat_sig[[i]]$mRNAs %in% inter
      print(paste0(i,"-related gene number in ",organ,' for ',t,' is ', table(sect)[2]))
      j = paste0('Num_Genes_',t)
      Trans[i, j] <- table(sect)[2]
      Trans_G[[t]][[i]] <- Dat_sig[[i]][sect,]
      gL <- Dat_sig[[i]][sect,1] #final genes from 3 parts
      sig_GL <<- c(sig_GL, unlist(gL))
    }
    Sig_G[[t]] <- unique(sig_GL)
  }
  Trans <- Trans[-nrow(Trans),]
  sig_glist <- unique(unlist(Sig_G))#511
  PCT_path <- paste0('/data/rluo4/All/Output/PCT/')
  save(Trans_G, DEG_list, Sig_G, Dat_sig, rdata_filter, file = paste0(PCT_path, 'PCT-',organ, '.rds' ))
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) # Choose the final 100 hub genes for each MP 
#######################################################################
library(purrr)
Trans_MP <- vector("list", length(MPs))
names(Trans_MP) <- MPs
for(i in MPs){
  # i = MPs[1]
  print(paste0("MP: ",i))
  data_path <- '/data/rluo4/All/Output/Cluster/'
  setwd(data_path)
  organ_all <- list.files(data_path)
  Dat_MP <-  vector("list", length(organ_all))
  names(Dat_MP) <- organ_all
  
  for(organ in organ_all){
    # organ="Prostate"
    analysis_parent_folder <- paste0(data_path,organ)
    if (!dir.exists(paste0(analysis_parent_folder))){
      dir.create(paste0(analysis_parent_folder))
    }
    setwd(analysis_parent_folder)
    getwd()
    path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
    f = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds') # old version
    f = paste0(path_PCT, 'PCT-',organ,'.rds') # new version
    
    if (file.exists(f)){
      load(f)
      sig_glist <- unique(unlist(Sig_G))
      #  print(f)
      j = 'EMT-related1'
      index = grepl(j, names(Dat_sig))
      if( i == 'EMT-related' & length(unique(index)) >1 ){
        print(paste0("MP: ",j," in ",organ))
        sect <- Dat_sig[[j]]$mRNAs %in% sig_glist
        print(paste0(j,"-related gene number in ",organ, ' is ', table(sect)[2]))
        Dat_MP[[organ]] = Dat_sig[[j]][sect,]
      } else{
        print(paste0("MP: ",i," in ",organ))
        sect <- Dat_sig[[i]]$mRNAs %in% sig_glist
        print(paste0(i,"-related gene number in ",organ, ' is ', table(sect)[2]))
        Dat_MP[[organ]] = Dat_sig[[i]][sect,]
      }
    }
  }
  Trans_MP[[i]] <- discard(Dat_MP, is.null)
}
# save(Trans_MP, file=paste0(data_path,'../Transition_MP.rds'))
# load(paste0(data_path,'../Transition_MP_old.rds'))

library(dplyr)
setwd(data_path)
load(paste0(data_path,'../MP_new.rds'))
Dat_MP <- vector("list", length(MPs))
names(Dat_MP) <- MPs
for (i in MPs){
  # i = MPs[1]
  data_path <- '/data/rluo4/All/Output/Cluster/'
  setwd(data_path)
  organ_all <- list.files(data_path)
  # organ = 'Breast'
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  print(Cor_table[i,])
  # 1.1) Cor > 0:
  # index = which(as.numeric(Cor_table[i,])>0)
  index = which(as.numeric(Cor_table[i,])!=0)
  print(index)
  n <- colnames(Cor_table)[index] 
  organ <- n[1]
  if(length(index)>0){
    print(paste0(i," in ",organ, ' first positively corr with CytoTRACE'))
  } else{
    print(paste0(i," in ",organ, ' first negatively corr with CytoTRACE'))
  }
  Dat <- Trans_MP[[i]][[organ]]#[,-1]
  colnames(Dat) <- gsub('cor',paste0('R.',n[1]),colnames(Dat))
  colnames(Dat) <- gsub('p.value',paste0('P.',n[1]),colnames(Dat))
  if(length(n)>1){
    n <- n[-1]
    for(j in 1:length(n)){
      organ = n[j]
      print(paste0(i," in ",organ, ' significantly corr with CytoTRACE'))
      Dat <- merge(Dat, Trans_MP[[i]][[organ]], by=1, all=T)#
      # Dat <- bind_rows(Dat, Trans_MP[[i]][[organ]]) %>% distinct()#left_join(Dat, Trans_MP[[i]][[organ]],  by='mRNAs')
      colnames(Dat) <- gsub('cor',paste0('R.',n[j]),colnames(Dat))
      colnames(Dat) <- gsub('p.value',paste0('P.',n[j]),colnames(Dat))
      index <- !is.na(Dat[,ncol(Dat)-1]) & !is.na(Dat[,ncol(Dat)])
      print(table(index))
      # Dat <- Dat[index,]
    }   
  }
  rownames(Dat) <- Dat[,1]
  Dat <- Dat[,-1]
  head(Dat)
  library(psych)
  Dat_abs <- apply(Dat[seq(1, ncol(Dat), by = 2),drop=FALSE], 2, abs)#apply(Dat[grep('R.',colnames(Dat)),drop=FALSE], 2, abs)
  print(head(sort(apply(Dat_abs, 1, geometric.mean),decreasing = T)))
  len <- na.omit(as.numeric(Cor_table[i,]))
  if(length(len) >1){
    Dat_MP[[i]] <- Dat[seq(1, ncol(Dat), by = 2),drop=FALSE]#list(Dat_p, Dat_n)
  } else if(length(len)==1){
    print(paste0(i, ' only in ', organ))
    Dat_MP[[i]] <- Dat[seq(1, ncol(Dat), by = 2),drop=FALSE]#list(Dat)
  } else{
    Dat_MP[[i]] <- NULL
  }
  Dat_MP[[i]]['MP_gene'] <- rownames(Dat_MP[[i]])
}
# There are >=28 genes in each MP (Dat_MP[[i]])
# unique(names(DEG_Trans))
# DEG_Trans <- DEG_Trans[unique(names(DEG_Trans))]

# Top 100 genes for each MP
library(psych)
Trans_genes <- NULL
Dat_final <- vector("list", length(MPs))
names(Dat_final) <- MPs
for (i in MPs){
  Dat <- Dat_MP[[i]]
  # Second Choice: select genes in more than 2 tissues #save as Transition_MP_m2ts.rds
  # if(ncol(Dat)>=2){
  #   n_organ <- apply(Dat, 1, function(x){
  #     N <- na.omit(x)
  #     return(length(N))
  #   })
  #   # Dat$n_organ <- n_organ
  #   Dat <- Dat[n_organ>=2,]
  # }
  # Dat_abs <- apply(Dat[seq(1, ncol(Dat), by = 2),drop=FALSE], 2, abs)
  Dat_abs <- apply(Dat[-ncol(Dat),drop=FALSE], 2, abs)
  order_val <- sort(apply(Dat_abs, 1, geometric.mean), decreasing = T)
  sect = head(names(order_val), 100)
  print( sect)
  order_data <- Dat[sect, , drop = FALSE]
  # order_data$MP_Gene <- rownames(order_data)
  order_data$MP <- i
  order_data$index <- 1:nrow(order_data)
  Dat_final[[i]] <- order_data
  MP_Gene <- order_data$MP
  names(MP_Gene) <- order_data$MP_gene#order_data$MP_Gene
  Trans_genes <<- c(Trans_genes,MP_Gene)
}

Trans_genes <- data.frame(MPs = Trans_genes, MP_Gene = names(Trans_genes))#unique(Trans_genes)
length(unique(Trans_genes$MP_Gene))#1365-->1083 (top 100), 712 (top 50) -->1039

# load(paste0(data_path,'../Transition_MP.rds'))
## load(paste0(data_path,'../Transition_MP_m2ts.rds'))

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) # load the cell clusters of PCT
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}

DEG_Trans <- NULL
Mut_Trans <- NULL
Mut_genes <- NULL
Mut_all <- NULL
data_path <- '/data/rluo4/All/Output/Cluster/'
setwd(data_path)
organ_all <- list.files(data_path)
for (organ in organ_all) {
  # organ = 'Breast'
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  # file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds') # old version
  file = paste0(path_PCT, 'PCT-',organ,'.rds') # new version
  
  load(file)
  # save the DEGs and Mutations in PCT
  names(DEG_list) <- paste0(organ,'_', names(DEG_list))
  DEG_Trans <- c(DEG_Trans, DEG_list)
  Mutants <- list(rdata_filter@meta.data[,c('barcode','orig.ident','SimplifiedSampleName','Sample_Type','cnv_score','Mut_count','Mut_gene','AAChange')])
  names(Mutants) <- paste0(organ,'_Mutants')
  Mut_Trans <- c(Mut_Trans, Mutants)
  Mut_G <- na.omit(rdata_filter$Mut_gene)
  Mut_genes <- c(Mut_genes, paste0(Mut_G[Mut_G!='']))
  Mut_df <-rdata_filter@meta.data[,c('barcode','orig.ident','SimplifiedSampleName','Sample_Type','cnv_score','Mut_count','Mut_gene','AAChange')]
  Mut_df$Tissue <- organ
  Mut_all <- rbind(Mut_all, Mut_df)
}
# Trans_genes
Mut_G <- paste0(unlist(strsplit(Mut_genes,', ')))
unique(Mut_G)
Mut_df <- as.data.frame(table(Mut_G))
final_G <- intersect(Mut_G, Trans_genes$MP_Gene)
# View(Trans_genes[Trans_genes$MP_Gene %in% final_G,])
final_G <- Mut_df[Mut_df$Mut_G %in% final_G,]
# View(final_G[final_G$Mut_G %in% Ligand_Receptor$Ligand_Receptor,])
# load(file=paste0(data_path,'../Transition_MP.rds'))

cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
DEG_list_all <- NULL
data_path <- '/data/rluo4/All/Output/Cluster/'
setwd(data_path)
organ_all <- list.files(data_path)
for (organ in organ_all) {
  # organ = 'Breast'
  print(organ)
  # save(Trans_G,DEG_list, Sig_G, Dat_sig, rdata_filter, file =paste0(cohort_directory, 'Epi/Results/','Malignant-Transformation-',organ, '.rds' ))
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds')
  load(file)
  print(table(rdata_filter$cell_class, rdata_filter$Trajectory))
  Trajectory <- unique(rdata_filter$Trajectory)
  if(organ %in% c('Breast','CRC')){ # original: T1 only
    Trajectory <- c('T1','T2')
  }
  if(organ %in% c('HNSCC','Lung')){ # origincal: T2, T1 only
    Trajectory <- c('T1','T2','T3')
  }
  
  print(paste0(organ, ': ', Trajectory))
  ########
  cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
  DEG_list <-  vector("list", length(Trajectory))
  names(DEG_list) <- Trajectory
  for(t in Trajectory){
    if(organ == 'Breast'){#?
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C14"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T1-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
    }
    
    if(organ == 'Cervix'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T3-C8"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C5"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T1-C5"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        } else{
          marker_C = marker%>%filter(cluster == c("T2-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T2-C10"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        }
      }
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
    } 
    
    if(organ == 'Chen'){
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T3-C17"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T3-C23"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T1-C21"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        } else{
          marker_C = marker%>%filter(cluster == c("T2-C26"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T2-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        }
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
    }
    
    if(organ == 'CRC'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C19"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T1-C16"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C17"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'Endometrium'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C14"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C12"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C24"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T2-C19"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'Esophagus'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C10"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C12"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } 
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'GC'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C17"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C2"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T2-C7"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'HNSCC'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T1-C16"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C18"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T1-C11"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        } else{
          marker_C = marker%>%filter(cluster == c("T2-C19"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T2-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        }
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'Liver'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds') #paste0(cohort_directory, 'Epi/Results/', organ, '_',t, '_deg.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C0"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C7"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T2-C15"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'Lung'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T3'){
        # Trajectory 3:
        marker_C = marker%>%filter(cluster==c("T1-C8"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        if(t == 'T1'){
          # Trajectory 1:
          marker_C = marker%>%filter(cluster==c("T1-C14"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T1-C6"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        } else{
          marker_C = marker%>%filter(cluster == c("T1-C7"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
          marker_P = marker%>%filter(cluster==c("T1-C3"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        }
      }  
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'Pancreas'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C8"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C2"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C4"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T2-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'Prostate'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T2'){# Trajectory 2:
        marker_C = marker%>%filter(cluster==c("T2-C34"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T2-C16"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T1-C5"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C42"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'Skin'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C9"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C11"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C11"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    if(organ == 'THCA'){
      
      f <-  paste0(cohort_directory, organ, '_deg_all.rds')
      load(file = f)
      print(t)
      marker = dat.markers
      marker$avg_log2FC = round(marker$avg_log2FC,2)
      marker$pct.ratio = marker$pct.1/marker$pct.2
      table(marker$cluster)
      if(t == 'T1'){
        # Trajectory 1:
        marker_C = marker%>%filter(cluster==c("T1-C1"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T1-C15"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      } else{
        marker_C = marker%>%filter(cluster == c("T2-C5"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
        marker_P = marker%>%filter(cluster==c("T2-C2"),avg_log2FC>=0.25,pct.1>=0.25,pct.ratio>=1.025)
      }
      
      del = intersect(marker_C$gene,marker_P$gene)
      print(del)
      
    }
    
    DEG_list[[t]] <- rbind(marker_C,#[ ! marker_C$gene %in% del,],
                           marker_P[ ! marker_P$gene %in% del,])
  }
  # save the DEGs and Mutations in PCT
  DEG_list$all <- marker
  names(DEG_list) <- paste0(organ,'_', names(DEG_list))
  DEG_list_all <- c(DEG_list_all, DEG_list)
}

save(Trans_MP,Trans_genes, Mut_all,  DEG_Trans, DEG_list_all, file=paste0(data_path,'../Transition_MP.rds'))

############################################################################################################
uniprot <- read.csv('/data/rluo4/All/uniprot-hs.tsv',fill = T,header = T,sep = '\t')
library(org.Hs.eg.db)
eg2trans <- toTable(org.Hs.egENSEMBLTRANS2EG)
eg2symbol=toTable(org.Hs.egSYMBOL)
eg2uniprot=toTable(org.Hs.egUNIPROT)
length(unique(eg2uniprot$uniprot_id))
length(unique(eg2uniprot$uniprot_id))
GeneList=mappedLkeys(org.Hs.egSYMBOL)
if( GeneList[1] %in% eg2symbol$symbol ){
  symbols=GeneList
  geneIds=eg2symbol[match(symbols,eg2symbol$symbol),'gene_id']
}else{
  geneIds=GeneList
  symbols=eg2symbol[match(geneIds,eg2symbol$gene_id),'symbol']
}
geneIds <- eg2symbol$gene_id
geneUniprot=eg2uniprot[match(eg2symbol$gene_id,eg2uniprot$gene_id),'uniprot_id']
head(geneUniprot)
table(is.na(match(eg2symbol$gene_id,eg2trans$gene_id)))
geneTrans=eg2trans[match(eg2symbol$gene_id,eg2trans$gene_id),'trans_id']
head(geneUniprot)

gene_info <- data.frame(   symbols=eg2symbol$symbol,
                           geneTrans=geneTrans,
                           geneUniprots=geneUniprot,
                           geneIds = geneIds,
                           stringsAsFactors = F
) 
x <- DEG_list_all$Breast_all[DEG_list_all$Breast_all$cluster=='T1-C4',]
x <- DEG_list_all$Breast_all[DEG_list_all$Breast_all$cluster=='T1-C14',]
x <- DEG_list_all$Liver_all[DEG_list_all$Liver_all$cluster=='T1-C0',]
x <- DEG_list_all$Liver_all[DEG_list_all$Liver_all$cluster=='T1-C1',]

x <- DEG_Trans$Liver_T1[DEG_Trans$Liver_T1$cluster=='T1-C0',]
x <- DEG_Trans$Liver_T1[DEG_Trans$Liver_T1$cluster=='T1-C1',]
# x <- DEG_Trans$Esophagus_T1[DEG_Trans$Esophagus_T1$cluster=='T1-C10',]
# x <- DEG_Trans$Prostate_T2[DEG_Trans$Prostate_T2$cluster=='T2-C34',]
# x <- DEG_Trans$GC_T1[DEG_Trans$GC_T1$cluster=='T1-C17',]
# x <- DEG_list_all$GC_all[DEG_list_all$GC_all$cluster=='T1-C17',]
x <- DEG_Trans$Breast_T2[DEG_Trans$Breast_T2$cluster=='T1-C4',]
x <- DEG_Trans$Breast_T1[DEG_Trans$Breast_T1$cluster=='T1-C14',]

table(x$gene %in% gene_info$symbols)
table(x$gene %in% Patt$Symbol)
x$uniprot <- gene_info$geneUniprots[match(x$gene,gene_info$symbols)]
table(x$uniprot %in% uniprot$Entry)
table(x$gene %in% uniprot$Gene.Names)
x$localization <- uniprot$Subcellular.location..CC.[match(x$uniprot,uniprot$Entry)]
table(x$localization=='')#3288
x$localization <- gsub('SUBCELLULAR LOCATION: ','',x$localization)
x$Nucleus <- ifelse(grepl('ucleus',x$localization),1,0)
x$Secreted <- ifelse(grepl('ecreted',x$localization),1,0)
x$Cytoplasm <- ifelse(grepl('ytoplasm',x$localization),1,0)
x$Surface <- ifelse(grepl('ell membrane',x$localization),1,0)
multi.loc <- apply(x[, c('Nucleus','Secreted','Cytoplasm','Surface')],1,function(x){
  s <- sum(x)
  return(s)
})
dim(x);length(multi.loc)
x$multi.loc <- multi.loc  
table(x$multi.loc)
table(x$multi.loc %in% c(0,1))
x <- x[x$multi.loc %in% c(0,1),]
x$localization[is.na(x$localization)] = ''
table(x$localization=='')#3638
# x <- x[x$localization!='',]
length(strsplit(x$localization[100], "[;]")[[1]])
x$loc[x$Nucleus==1] <- 'Nucleus'
x$loc[x$Secreted==1] <- 'Secreted'
x$loc[x$Cytoplasm==1] <- 'Cytoplasm'
x$loc[x$Surface==1] <- 'Surface'
# table(x$loc)
x$loc <- str_split(x$loc,'[;]',simplify = T)[,1]
x$loc <- str_split(x$loc,'[.]',simplify = T)[,1]
x$loc <- str_split(x$loc,'[{]',simplify = T)[,1]
x$loc <- str_split(x$loc,'[,]',simplify = T)[,1]
table(x$loc)
table(x$loc[x$loc %ni% c('Nucleus','Secreted','Cytoplasm','Surface')])
x$loc[x$loc %ni% c('Nucleus','Secreted','Cytoplasm','Surface')] <- 'Rest'
table(x$loc)
loc_SEN_up <- data.frame(table(x$loc))
loc_SEN_up$pct <- loc_SEN_up$Freq/sum(loc_SEN_up$Freq)
loc_SEN_up
loc_SEN_down <- data.frame(table(x$loc))
loc_SEN_down$pct <- loc_SEN_down$Freq/sum(loc_SEN_down$Freq)
loc_SEN_down


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 9) # load CellChat object
library(CellChat)
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
Cell_Clusters =  vector("list", length(organ_all))
names(Cell_Clusters) <- organ_all
load(paste0(data_path,'../Transition_MP.rds'))
########################################################################
for (organ in organ_all) {
  # Define Cell_Clusters
  ########
  print(organ)
  if(organ == 'Breast'){#
    Cell_Clusters[[organ]] <- c("T1-C14", "T1-C9",'T1-C4')#c("T1-C4", "T1-C9")
  }
  if(organ == 'Cervix'){#
    Cell_Clusters[[organ]] <- c("T1-C1",  'T1-C14', 'T3-C8', "T3-C5", "T2-C4", "T2-C10")
  }
  if(organ == 'Chen'){# 
    Cell_Clusters[[organ]] <- c("T1-C1", "T1-C21", "T1-C24",   "T2-C26", "T2-C3",   "T3-C17", "T3-C23")#c("T1-C1", "T1-C24", "T2-C26", "T2-C3", "T3-C17", "T3-C23")
  }
  if(organ == 'CRC'){# which is Becker
    Cell_Clusters[[organ]] <- c("T1-C3", "T1-C19",   "T1-C16", "T1-C17")
  }
  if(organ == 'Endometrium'){#
    Cell_Clusters[[organ]] <- c("T2-C24", "T2-C4", 'T2-C19',   "T1-C14", "T1-C12")#c("T1-C14", "T1-C12", "T2-C24", "T2-C19")
  }
  if(organ == 'Esophagus'){#
    Cell_Clusters[[organ]] <- c("T1-C10", "T1-C12") 
  }
  if(organ == 'GC'){#
    Cell_Clusters[[organ]] <- c("T1-C17", "T1-C2", "T2-C3", "T2-C7")#c("T1-C13", "T1-C2", "T2-C3", "T2-C12")
  }
  if(organ == 'HNSCC'){#
    Cell_Clusters[[organ]] <- c("T1-C4", 'T1-C11', "T2-C19", "T2-C3", "T1-C18", 'T1-C16')# ,c("T1-C11", "T1-C18", "T2-C19", "T2-C3")
  }
  if(organ=='Liver'){#
    Cell_Clusters[[organ]] <- c("T1-C0", "T1-C1", "T2-C7", "T2-C15")#c("T1-C0", "T1-C1", "T2-C7", "T2-C11")
  }  
  if(organ=='Lung'){#
    Cell_Clusters[[organ]] <- c("T1-C14",  "T1-C6",   "T1-C7", "T1-C3", "T1-C8",  "T1-C9")#c("T1-C12", "T1-C6", "T1-C4", "T1-C3")
  }
  if(organ == 'Pancreas'){#
    Cell_Clusters[[organ]] <- c("T1-C8", "T1-C2",  "T2-C4", "T2-C9")#c("T1-C8", "T1-C2", "T2-C4", "T2-C9")
  }
  if(organ == "Prostate"){#
    Cell_Clusters[[organ]] <- c("T1-C42", "T1-C5", "T2-C34", "T2-C14")#c("T1-C6", "T1-C9", "T2-C34", "T2-C14")
  }
  if(organ == "Skin"){#
    Cell_Clusters[[organ]] <- c( "T2-C1", "T1-C11", 'T1-C9') #, "T2-C14"# c("T1-C12", "T1-C9", "T2-C1", "T2-C11")
  } 
  if(organ=="THCA"){#
    Cell_Clusters[[organ]] <- c("T1-C1", "T1-C15", "T2-C5", "T2-C2")
  }
  
}
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
setwd(data_path)
organ_all <- list.files(data_path)
CCI_summary <- vector("list",length(organ_all))
names(CCI_summary) <- organ_all
for (organ in organ_all) {
  # organ = 'Breast'
  CellChat <- NULL
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  # file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds') # old version
  file = paste0(path_PCT, 'PCT-',organ,'.rds') # new version
  load(file)
  # print(table(rdata_filter$SimplifiedSampleName))
  rdata_filter$cell_subtype <- extract_last_element(rdata_filter$cell_class)
  subtype <- table(rdata_filter$orig.ident, rdata_filter$cell_subtype)
  subtype <- data.frame(subset = subtype)
  subtype <- subtype[subtype$subset.Freq !=0,]
  subtype <- subtype[subtype$subset.Freq >=10,]
  if(organ=='Liver'){
    subtype <- table(rdata_filter$SimplifiedSampleName, rdata_filter$cell_subtype)
    subtype <- data.frame(subset = subtype)
    subtype <- subtype[subtype$subset.Freq !=0,]
    subtype <- subtype[subtype$subset.Freq >=10,]
  }
  if(organ=='Esophagus'){
    subtype <- table(rdata_filter$orig.ident, rdata_filter$cell_subtype)
    subtype <- data.frame(subset = subtype)
    subtype <- subtype[subtype$subset.Freq !=0,]
    # subtype <- subtype[subtype$subset.Freq >=10,]
  }
  print(Cell_Clusters[[organ]])
  
  for (t in Cell_Clusters[[organ]]) {
    print(t)
    subcluster <- extract_last_element(t)
    patient = unique(subtype$subset.Var1[subtype$subset.Var2==subcluster])
    TissueType <- paste0(patient, collapse = '_')
    
    if(organ == 'Esophagus' & t == 'T1-C12'){
      TissueType = 'LZE11D_LZE20T_LZE21D1_LZE21T_LZE22D1_LZE22D3_LZE22T_LZE24D1_LZE24T_LZE2D_LZE2T_LZE3D_LZE4T_LZE5T_LZE6T_LZE7T_LZE8T_P107T-E_P130T-E_P16T-E_P75T-E_P76T-E'
    }
    if(organ == 'Esophagus' & t == 'T1-C10'){
      TissueType = 'P104T-E_P127T-E_P128T-E_P15T-E_P16T-E_P1T-E_P21T-E_P23T-E_P2T-E_P31T-E_P39T-E_P47T-E_P54T-E_P57T-E_P61T-E_P65T-E_P74T-E_P75T-E_P76T-E_P79T-E_P82T-E_P84T-E_P8T-E'
    }
    
    Sample = paste0(TissueType,'_cluster',subcluster,'_CCI.RData')
    print(Sample)
    file = paste0(data_path,'../CellChat/',Sample)
    
    if(file.exists(file)){
      load(file)
    } else{
      print(paste0(organ, ' ', Sample, " doesn't exist!"))
    }
    
    print(table(cellchat@idents))
    signaling.name <- cellchat@netP$pathways
    # print(table(pathways.show.all %in% pathways.all))
    # print(pathways.all[!pathways.all %in% pathways.show.all])
    if (is.null(signaling.name)) {
      signaling.name <- signaling
    }
    # CCI <- vector("list",length(signaling.name))
    # names(CCI) <- signaling.name
    data.cci=NULL
    thresh <- 0.05
    net <- cellchat@net
    pairLR.use.name <- dimnames(net$prob)[[3]]
    for (i in signaling.name) {
      pairLR <- searchPair(signaling = i, pairLR.use = cellchat@LR$LRsig,
                           key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
      pairLR <- pairLR[pairLR.name, ]
      prob <- net$prob
      pval <- net$pval
      prob[pval > thresh] <- 0
      if (length(pairLR.name) > 1) {
        pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 3, sum) != 0]
      }else{
        pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 0]
      }
      print(setdiff(pairLR.name, pairLR.name.use))
      if (length(pairLR.name.use) == 0) {
        stop(paste0("There is no significant communication of ",
                    signaling.name))
      }else {
        pairLR <- pairLR[pairLR.name.use, ]
      }
      # nRow <- length(pairLR.name.use)
      # prob <- prob[, , pairLR.name.use]
      # pval <- pval[, , pairLR.name.use]
      # if (length(pairLR.name.use) > 1){
      #   CCI <- apply(prob,3,function(x){
      #     c(t(x))   ## 逐行转换为向量# c(dat) ## 逐列转换为向量
      #   })
      # }else {
      #   CCI <- as.matrix(c(t(prob)))
      #   colnames(CCI) <- pairLR.name.use
      # }
      data.cci <<- rbind(data.cci,pairLR)
    }
    data.cci <- data.cci[,c(3,4,1,2)]
    colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
    data.cci$Tissue <- organ
    data.cci$DiseaseStage <- paste0(TissueType,'_cluster',subcluster)
    data.cci$Trace <- t
    CellChat <<- rbind(CellChat, data.cci)
  }
  CCI_summary[[organ]] <- CellChat
}
# setdiff(CellChat_GC$LRpair[CellChat$DiseaseStage=='GC'],
#         CellChat_GC$LRpair[CellChat$DiseaseStage=='Precancer'])
# 
# setdiff(CellChat_GC$LRpair[CellChat$DiseaseStage=='Precancer'],
#         CellChat_GC$LRpair[CellChat$DiseaseStage=='Healthy'])
df_empty <- NULL
CCI_summary_All <-  lapply(CCI_summary,function(y){
  df_empty <<- rbind(df_empty, y)
  return(df_empty)
})
CCI_summary_All <- df_empty
CCI_summary_All$Index <- 1:nrow(CCI_summary_All)
# write.table(CCI_summary_All[,c(7,1:6)], '/data/rluo4/All/Output/Index_LR_All.txt',#paste0(TissueType,'_Index_CellChat.txt'),
#             sep='\t',quote=F,row.name=F,col.name=T)
library(stringr)
LRpairs <- CCI_summary_All#read.csv('Index_LR_All.txt', sep = '\t')
L <- strsplit(LRpairs$Ligand,split = '')
L <- unlist(L)
non_alnum_chars <- grepl("[^[:alnum:]]", L)
table(L[non_alnum_chars])
L1 <- str_split(LRpairs$Ligand, "[_]", simplify = T)[,1]
L1 <- L1[!grepl('_',L1)]
L2 <- str_split(LRpairs$Ligand, "[_]", simplify = T)[,2]
L2 <- unique(L2)
Ligand <- c(L1, L2)
Ligand <- unique(Ligand)

R <- strsplit(LRpairs$Receptor,split = '')
R <- unlist(R)
non_alnum_chars <- grepl("[^[:alnum:]]", R)
table(R[non_alnum_chars])
R1 <- str_split(LRpairs$Receptor, "[:]", simplify = T)[,1]
R1 <- R1[!grepl('_',R1)]
R2 <- str_split(LRpairs$Receptor, "[:]", simplify = T)[,2]
R2 <- unique(R2)
R1.1 <- str_split(LRpairs$Receptor, "[_]", simplify = T)[,1]
R1.1 <- R1.1[!grepl(':',R1.1)]
R2.2 <- str_split(LRpairs$Receptor, "[_]", simplify = T)[,2]
R2.2 <- unique(R2.2)
Receptor <- c(R1,R2,R1.1,R2.2)
Receptor <- unique(Receptor)

Ligand_Receptor <- data.frame(Index= 1:length(c(Ligand, Receptor)),
                              Ligand_Receptor=c(Ligand, Receptor))
# write.table(Ligand_Receptor, 'Ligand_Receptor.txt',sep = '\t',quote = F)
table(Ligand_Receptor$Ligand_Receptor %in% Trans_genes$MP_Gene)
unique(intersect(Ligand_Receptor$Ligand_Receptor, Trans_genes$MP_Gene))# 70 genes --> 49 genes
# View(Trans_genes[Trans_genes$MP_Gene %in% Ligand_Receptor$Ligand_Receptor,])

table(Ligand_Receptor$Ligand_Receptor=='')

intersect(Ligand_Receptor$Ligand_Receptor, Trans_genes$MP_Gene)
MP_LR <- Trans_genes[Trans_genes$MP_Gene %in% Ligand_Receptor$Ligand_Receptor,]
MP_LR$MP_Gene[MP_LR$MPs=='Epithelial-senescence']
unique(MP_LR$MPs) # except for the G2M-Cellcycle

Index_LR <- NULL
summary <- apply(LRpairs,1,function(x){
  if(grepl('_',x['Ligand'])){
    L1 <- str_split(x['Ligand'], "[_]", simplify = T)[,1]
    # L1 <- L1[!grepl('_',L1)]
    L2 <- str_split(x['Ligand'], "[_]", simplify = T)[,2]
    #L2 <- unique(L2)
  } else{
    L1 <- L2 <- x['Ligand']
  }
  if(grepl(':',x['Receptor'])){
    R1 <- str_split(x['Receptor'], "[:]", simplify = T)[,1]
    #R1 <- R1[!grepl('_',R1)]
    R2 <- str_split(x['Receptor'], "[:]", simplify = T)[,2]
    #R2 <- unique(R2)
  } else{
    R1 <- R2 <- x['Receptor']
  }
  if(grepl('_',x['Receptor'])){
    R1.1 <- str_split(x['Receptor'], "[_]", simplify = T)[,1]
    #R1.1 <- R1.1[!grepl(':',R1.1)]
    R2.2 <- str_split(x['Receptor'], "[_]", simplify = T)[,2]
    # R2.2 <- unique(R2.2)
  } else{
    R1.1 <- R2.2 <- x['Receptor']
  }
  if( L1 %in% Trans_genes$MP_Gene | L2 %in% Trans_genes$MP_Gene |
      R1 %in% Trans_genes$MP_Gene | R2 %in% Trans_genes$MP_Gene |
      R1.1 %in% Trans_genes$MP_Gene | R2.2 %in% Trans_genes$MP_Gene ){
    Index_LR <<- rbind(Index_LR, x)
  }
})
Index_LR <- as.data.frame(Index_LR)
Index_LR$Index <- 1:nrow(Index_LR)
write.table(Index_LR[, c(8, 1:7)], 'Index_LR.txt', sep='\t',quote=F,row.name=F,col.name=T)
sort(table(Trans_genes[Trans_genes$MP_Gene %in% Ligand_Receptor$Ligand_Receptor,]$MP_Gene))
# CCL2      CCL3      CD74      CD99    COL4A2     CXCL8       GRN  HLA-DPA1  HLA-DQB1     ITGA6 
# 1         1         1         1         1         1         1         1         1         1 
# NECTIN2     PRSS3      SDC2     VEGFA      AREG     CCL20     CXCL2     CXCL3       FN1     HLA-A 
# 1         1         1         1         2         2         2         2         2         2 
# IL1R2       MDK       MIF       NCL TNFRSF12A   TNFSF10       APP     GDF15     HLA-C  HLA-DRB1 
# 2         2         2         2         2         2         3         3         3         3 
# IL6     IL6ST     ITGB4     NTRK2   RARRES2      SAA1    SEMA3C      DSC2     EFNA5      NRG1 
# 3         3         3         3         3         3         3         4         4         4 
# SPP1     ANXA1   CEACAM5  HLA-DPB1   HLA-DRA     LAMB3      DSC3     PTPRF      SDC1 
# 4         5         5         5         5         6         7         8        10 
save(CCI_summary,Cell_Clusters, CCI_summary_All, Ligand_Receptor, Index_LR, file = paste0(data_path,'../CCI_summary_MP.rds'))

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) save CCI_all object
library(CellChat)
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
CCI_Clusters =  vector("list", length(organ_all))
names(CCI_Clusters) <- organ_all
########################################################################
for (organ in organ_all) {
  # Define CCI_Clusters
  ########
  print(organ)
  if(organ == 'Breast'){#
    CCI_Clusters[[organ]] <- c("T1-C14", "T1-C9",'T1-C4')
  }
  if(organ == 'Cervix'){#
    CCI_Clusters[[organ]] <- c("T3-C5", "T1-C1",  'T1-C14', 'T3-C8')# "T2-C4", "T2-C10")
  }
  if(organ == 'Chen'){# 
    CCI_Clusters[[organ]] <- c("T1-C1", "T1-C21", "T1-C24",   "T2-C26", "T2-C3",   "T3-C17", "T3-C23")
  }
  if(organ == 'CRC'){# which is Becker
    CCI_Clusters[[organ]] <- c("T1-C16", "T1-C17", "T1-C19", "T1-C3" )
  }
  if(organ == 'Endometrium'){#
    CCI_Clusters[[organ]] <- c('T1-C12','T1-C14', 'T2-C19', "T2-C4",  "T2-C24")
  }
  if(organ == 'Esophagus'){#
    CCI_Clusters[[organ]] <- c("T1-C10", "T1-C12") 
  }
  if(organ == 'GC'){#
    CCI_Clusters[[organ]] <- c('T2-C12', "T2-C7", "T1-C2", "T1-C17",  "T2-C3")
  }
  if(organ == 'HNSCC'){#
    CCI_Clusters[[organ]] <- c("T1-C4", 'T1-C11', "T2-C19", "T2-C3", "T1-C18", 'T1-C16')
  }
  if(organ=='Liver'){#
    CCI_Clusters[[organ]] <- c("T1-C0", "T1-C1", 'T2-C11', "T2-C7", "T2-C15")
  }  
  if(organ=='Lung'){#
    CCI_Clusters[[organ]] <- c("T1-C6", "T1-C9", "T1-C8", "T1-C14" , "T1-C3","T1-C7")
  }
  if(organ == 'Pancreas'){#
    CCI_Clusters[[organ]] <- c( "T2-C9", "T2-C4",  "T1-C8", "T1-C2")#c("T2-C6",'T2-C5', 'T2-C0', 'T2-C1', "T2-C9", "T2-C4",  "T1-C8", "T1-C2")
  }
  if(organ == "Prostate"){#
    CCI_Clusters[[organ]] <- c("T1-C42", "T1-C5", "T2-C34", "T2-C14")
  }
  if(organ == "Skin"){#
    CCI_Clusters[[organ]] <- c( "T2-C1", 'T2-C14', "T1-C11", 'T1-C9') 
  } 
  if(organ=="THCA"){#
    CCI_Clusters[[organ]] <- c("T2-C2",  "T1-C15", "T1-C5", "T1-C1")#c("T2-C2", "T2-C16", "T1-C15", "T1-C5", "T1-C1")
  }
  
}

library(CellChat)
library(stringr)
indir <-'/data/rluo4/database/'
setwd(data_path)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}

load(file = paste0(data_path,'../CCI_summary_MP.rds'))
# save(CCI_STM_Myeloid, CCI_all, Cell_Clusters, file = paste0(data_path,'../CCI_ANXA1.rds'))
# load(file =  paste0(data_path,'../CCI_ANXA1.rds') )

CCI_summary <- vector("list",length(organ_all))
names(CCI_summary) <- organ_all
CellChat_PCT <- NULL
for (organ in organ_all) {
  # organ = 'Breast'
  CellChat <- NULL
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  # file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds') # old version
  file = paste0(path_PCT, 'PCT-',organ,'.rds') # new version
  load(file)
  # print(table(rdata_filter$SimplifiedSampleName))
  rdata_filter$cell_subtype <- extract_last_element(rdata_filter$cell_class)
  subtype <- table(rdata_filter$orig.ident, rdata_filter$cell_subtype)
  subtype <- data.frame(subset = subtype)
  subtype <- subtype[subtype$subset.Freq !=0,]
  subtype <- subtype[subtype$subset.Freq >=10,]
  if(organ=='Liver'){
    subtype <- table(rdata_filter$SimplifiedSampleName, rdata_filter$cell_subtype)
    subtype <- data.frame(subset = subtype)
    subtype <- subtype[subtype$subset.Freq !=0,]
    subtype <- subtype[subtype$subset.Freq >=10,]
  }
  if(organ=='Esophagus'){
    subtype <- table(rdata_filter$orig.ident, rdata_filter$cell_subtype)
    subtype <- data.frame(subset = subtype)
    subtype <- subtype[subtype$subset.Freq !=0,]
    # subtype <- subtype[subtype$subset.Freq >=10,]
  }
  print(CCI_Clusters[[organ]])
  
  # for (t in Cell_Clusters[[organ]]) {
  for (t in CCI_Clusters[[organ]]) {
    print(t)
    subcluster <- extract_last_element(t)
    patient = unique(subtype$subset.Var1[subtype$subset.Var2==subcluster])
    TissueType <- paste0(patient, collapse = '_')
    
    if(organ == 'Esophagus' & t == 'T1-C12'){
      TissueType = 'LZE11D_LZE20T_LZE21D1_LZE21T_LZE22D1_LZE22D3_LZE22T_LZE24D1_LZE24T_LZE2D_LZE2T_LZE3D_LZE4T_LZE5T_LZE6T_LZE7T_LZE8T_P107T-E_P130T-E_P16T-E_P75T-E_P76T-E'
    }
    if(organ == 'Esophagus' & t == 'T1-C10'){
      TissueType = 'P104T-E_P127T-E_P128T-E_P15T-E_P16T-E_P1T-E_P21T-E_P23T-E_P2T-E_P31T-E_P39T-E_P47T-E_P54T-E_P57T-E_P61T-E_P65T-E_P74T-E_P75T-E_P76T-E_P79T-E_P82T-E_P84T-E_P8T-E'
    }
    
    Sample = paste0(TissueType,'_cluster',subcluster,'_CCI.RData')
    print(Sample)
    file = paste0(data_path,'../CellChat/',Sample)
    
    if(file.exists(file)){
      load(file)
    } else{
      print(paste0(organ, ' ', Sample, " doesn't exist!"))
    }
    
    print(table(cellchat@idents))
    signaling.name <- cellchat@netP$pathways
    # print(table(pathways.show.all %in% pathways.all))
    # print(pathways.all[!pathways.all %in% pathways.show.all])
    if (is.null(signaling.name)) {
      signaling.name <- signaling
    }
    # CCI <- vector("list",length(signaling.name))
    # names(CCI) <- signaling.name
    print(signaling.name)
    data.cci=NULL
    data.cci <- subsetCommunication(cellchat)
    data.cci <- data.cci[, c(1:7, 9)]#data.cci[,c(3,4,1,2)]
    # colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
    data.cci$Tissue <- organ
    data.cci$DiseaseStage <- paste0(TissueType,'_cluster',subcluster)
    data.cci$Trace <- t
    
    data.cci$TissueType <- paste(data.cci$Tissue, extract_last_element(data.cci$Trace), sep = '-')
    CellChat <<- rbind(CellChat, data.cci)
    cellchat_obj <-  vector("list", 1)
    # names(cellchat_obj) <- c('obj')
    names(cellchat_obj) <- unique(data.cci$TissueType)
    cellchat_obj[[1]] <- cellchat
    CellChat_PCT <- c(CellChat_PCT, cellchat_obj)
    
  }
  CCI_summary[[organ]] <- CellChat
}


df_empty <- NULL
CCI_all <-  lapply(CCI_summary,function(y){
  df_empty <<- rbind(df_empty, y)
  return(df_empty)
})
CCI_all <- df_empty
CCI_all$Index <- 1:nrow(CCI_all)
unique(CCI_all$target[CCI_all$ligand=='ANXA1'])
CCI_all$target <- gsub('MDSCs','MDSC',CCI_all$target)
save(CellChat_PCT,  CCI_Clusters, CCI_all, file = paste0(data_path,'../CellChat_PCT.rds'))
write.csv(CCI_all, file ="/data/rluo4/All/Output/CCI_all.csv", row.names = F, quote = F) #Table S6
write.xlsx(CCI_all[, -13], file = paste0('/data/rluo4/All/Output/table_xlsx/TableS6.xlsx'), rowNames = FALSE)#,
#############################################################################################################
# i) STM -> Myeloid 
load(file = paste0(data_path,'../CellChat_PCT.rds'))

CCI_all_ANXA1 <- CCI_all[CCI_all$ligand=='ANXA1',]
unique(CCI_all_ANXA1$target)# 22: Myeloid + TREG, TFH, MSC.ADIPO
Lineage <- read.csv('/data/rluo4/All/Output/organ13-celltypes.csv')
CCI_all_ANXA1 <- CCI_all_ANXA1[CCI_all_ANXA1$target %in% Lineage$Minor_type[Lineage$Major_type=='Myeloid'],]
CCI_all_ANXA1 <- CCI_all_ANXA1[CCI_all_ANXA1$source=='STM',]
unique(CCI_all_ANXA1$Tissue)# 11
unique(CCI_all_ANXA1$target)# 22
# CCI_STM_Myeloid <- CCI_all[CCI_all$source=='STM' & CCI_all$target %in% unique(CCI_all_ANXA1$target),]
# unique(CCI_STM_Myeloid$target)
CCI_STM_Myeloid <- NULL
for (organ in organ_all) {
  dataPlot <- CCI_all[CCI_all$source=='STM' &  CCI_all$Tissue==organ,]
  dataPlot <-  dataPlot[dataPlot$target %in%  unique(CCI_all_ANXA1$target[CCI_all_ANXA1$Tissue==organ]), ] 
  CCI_STM_Myeloid <<- rbind(CCI_STM_Myeloid, dataPlot)
}
# CCI_STM_Myeloid <- CCI_all[CCI_all$source=='STM' & CCI_all$target %in% Lineage$Minor_type[Lineage$Major_type=='Myeloid'],]
unique(CCI_STM_Myeloid$target)
unique(CCI_STM_Myeloid$ligand)
table(CCI_STM_Myeloid$pathway_name)
pathway_sub <- sort(table(CCI_STM_Myeloid$pathway_name),decreasing = T)
pathway_sub
CCI_STM_Myeloid <- CCI_STM_Myeloid[CCI_STM_Myeloid$pathway_name %in% c(head(names(pathway_sub),15), 'ANXA1'),]
table(CCI_STM_Myeloid$Tissue)
unique(CCI_STM_Myeloid$pathway_name)
CCI_STM_Myeloid <- CCI_STM_Myeloid[! CCI_STM_Myeloid$pathway_name %in% c("MHC-I", "MHC-II", 'FN1','COLLAGEN','ICAM'),]
unique(CCI_STM_Myeloid$pathway_name)

CCI_STM_Myeloid$trace <- str_split(CCI_STM_Myeloid$Trace, '-', simplify = T)[,2]
CCI_STM_Myeloid$Tissue <- gsub('CRC','Colon',gsub('Chen','Colon', CCI_STM_Myeloid$Tissue))
CCI_STM_Myeloid$Tissue <- gsub('GC','Stomach',gsub('HNSCC','OralCalvity', CCI_STM_Myeloid$Tissue))
CCI_STM_Myeloid$Tissue <- gsub('THCA','Thyroid', CCI_STM_Myeloid$Tissue)
# CCI_STM_Myeloid$trace <- paste(CCI_STM_Myeloid$Tissue, CCI_STM_Myeloid$trace, sep = '-')
unique(CCI_STM_Myeloid$trace)
CCI_STM_Myeloid$interaction_name <- as.character(CCI_STM_Myeloid$interaction_name)
CCI_STM_Myeloid$source.target <- paste(CCI_STM_Myeloid$trace, CCI_STM_Myeloid$target, 
                                       sep = " -> ")
unique(CCI_STM_Myeloid$source.target )
# pathway_sub <- as.data.frame(table(CCI_STM_Myeloid$pathway_name, CCI_STM_Myeloid$interaction_name))
# pathway_sub <- pathway_sub[pathway_sub$Freq!=0,]
# Fig.5
# Dot plot of STM-Myeloid from 11 tissues
# load libraries
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(data.table)
library(ggsci)
library(ggrepel)
library(viridis)
# library(fgsea)
library(dplyr)
library(tibble)
library(ggpointdensity)
library(RColorBrewer)
library(ggnewscale)
library(wesanderson)
unique(plotdata$pathway_name)
plotdata <- CCI_STM_Myeloid[order(CCI_STM_Myeloid$interaction_name,decreasing = F),]
plotdata <- CCI_STM_Myeloid[order(CCI_STM_Myeloid$pathway_name,decreasing = F),]
plotdata <- plotdata[plotdata$Tissue %in% c('Breast', 'Liver', 'Cervix'),]
source.targets <- unique(as.character(plotdata$source.target))
Tissues <- unique(as.character(plotdata$Tissue))
interaction_names <- unique(as.character(plotdata$interaction_name))
plotdata <- plotdata %>%
  # filter(gene %in% genelist) %>%
  filter(pval < 0.05) %>%
  mutate(source.target = factor(source.target, levels = source.targets)) %>%
  # mutate(lof_gof = ifelse(grepl("GOF", cell_line), "GOF", "LOF")) %>% 
  mutate(Tissue = factor(Tissue, levels = Tissues)) %>% 
  # mutate(cell_line = gsub("_", " ", gsub("_GOF", "", cell_line))) %>% 
  mutate(interaction_name = factor(interaction_name, levels = interaction_names))

ggplot(plotdata, aes(x = source.target, y = interaction_name, size = -log10(pmax(pval, 1e-10)) )) + #-log10(pval))) +
  geom_point(aes(fill = prob), pch = 21) +
  scale_fill_distiller("Commun. Prob.", palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = c(quantile(plotdata$prob[plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]],
                                                         0, na.rm = T), quantile(plotdata$prob[plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]], 1, na.rm = T)), guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +#, limits = LIMITS) +
  ggnewscale::new_scale_fill() +
  geom_point(data = plotdata[!plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)],], aes(fill = prob), pch = 21) +
  scale_fill_distiller("Commun. Prob.", palette = "PiYG", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = c(quantile(plotdata$prob[!plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]],
                                                         0, na.rm = T), quantile(plotdata$prob[!plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]], 1, na.rm = T)), direction = 1, guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +#, limits = LIMITS) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  theme_bw() + 
  
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',#"right", #
        legend.key.width = unit(0.3, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        plot.title = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed")) +
  labs(size = "P value (-log10)") +
  guides(size = guide_legend(override.aes = list(pch = 16), title.position = "top", title.hjust = 0.5)) +
  ylab("") +
  xlab("") +
  facet_grid(. ~ Tissue, scales = "free_x", space = "free_x")
ggsave("/data/rluo4/All/Output/res_fig/CCI_STM_Myeloid_dotplot.png", height = 14.75, width = 19.75, dpi = 500, device = "png")
ggsave("/data/rluo4/All/Output/res_fig/CCI_STM_Myeloid_dotplot_3tissues.png", height = 7, width = 7, dpi = 500, device = "png")

ggsave("/data/rluo4/All/Output/res_fig/CCI_STM_Myeloid_dotplot.pdf", height = 12.75, width = 23.75)

#################################################################################################################
# ii) CAF -> Myeloid
CCI_all_ANXA1 <- CCI_all[CCI_all$ligand=='ANXA1',]
unique(CCI_all_ANXA1$target)# 22: Myeloid + TREG, TFH, MSC.ADIPO
CCI_all_ANXA1 <- CCI_all_ANXA1[CCI_all_ANXA1$source %in% Lineage$Minor_type[Lineage$Major_type=='Mesenchymal'],]
CCI_all_ANXA1 <- CCI_all_ANXA1[CCI_all_ANXA1$target %in% Lineage$Minor_type[Lineage$Major_type=='Myeloid'],]
unique(CCI_all_ANXA1$Tissue)# 11
unique(CCI_all_ANXA1$target)# 10
CCI_CAF_Myeloid <- NULL
for (organ in organ_all) {
  dataPlot <- CCI_all[CCI_all$source %in% Lineage$Minor_type[Lineage$Major_type=='Mesenchymal'] &  CCI_all$Tissue==organ,]
  dataPlot <-  dataPlot[dataPlot$target %in%  unique(CCI_all_ANXA1$target[CCI_all_ANXA1$Tissue==organ]), ] 
  CCI_CAF_Myeloid <<- rbind(CCI_CAF_Myeloid, dataPlot)
}

unique(CCI_CAF_Myeloid$target)
unique(CCI_CAF_Myeloid$ligand)
table(CCI_CAF_Myeloid$pathway_name)
pathway_sub <- sort(table(CCI_CAF_Myeloid$pathway_name),decreasing = T)
pathway_sub
CCI_CAF_Myeloid <- CCI_CAF_Myeloid[CCI_CAF_Myeloid$pathway_name %in% c(head(names(pathway_sub),12), 'ANXA1','CXCL'),]
table(CCI_CAF_Myeloid$Tissue)
unique(CCI_CAF_Myeloid$pathway_name)
CCI_CAF_Myeloid <- CCI_CAF_Myeloid[! CCI_CAF_Myeloid$pathway_name %in% c("MHC-I", "MHC-II", 'APP', 'CD99'),]
CCI_CAF_Myeloid$trace <- str_split(CCI_CAF_Myeloid$Trace, '-', simplify = T)[,2]
CCI_CAF_Myeloid$Tissue <- gsub('CRC','Colon',gsub('Chen','Colon', CCI_CAF_Myeloid$Tissue))
CCI_CAF_Myeloid$Tissue <- gsub('GC','Stomach',gsub('HNSCC','OralCalvity', CCI_CAF_Myeloid$Tissue))
CCI_CAF_Myeloid$Tissue <- gsub('THCA','Thyroid', CCI_CAF_Myeloid$Tissue)
# CCI_CAF_Myeloid$trace <- paste(CCI_CAF_Myeloid$Tissue, CCI_CAF_Myeloid$trace, sep = '-')
unique(CCI_CAF_Myeloid$trace)
CCI_CAF_Myeloid$interaction_name <- as.character(CCI_CAF_Myeloid$interaction_name)
CCI_CAF_Myeloid$source.target <- paste(CCI_CAF_Myeloid$source, CCI_CAF_Myeloid$target, 
                                       sep = " -> ")
unique(CCI_CAF_Myeloid$source.target )
CCI_CAF_Myeloid$index <- paste(CCI_CAF_Myeloid$TissueType, CCI_CAF_Myeloid$source.target, sep=': ')
CCI_CAF_Myeloid$source.target <- paste(CCI_CAF_Myeloid$trace, CCI_CAF_Myeloid$source.target, sep=': ')

LR_annexin <- CCI_CAF_Myeloid$index[CCI_CAF_Myeloid$pathway_name=='ANNEXIN']
index <- paste(CCI_CAF_Myeloid$Tissue, CCI_CAF_Myeloid$source.target, sep = '-')
table(CCI_CAF_Myeloid$index %in% LR_annexin)
CCI_CAF_Myeloid <- CCI_CAF_Myeloid[CCI_CAF_Myeloid$index %in% unique(LR_annexin), ]

# Fig.5
plotdata <- CCI_CAF_Myeloid[order(CCI_CAF_Myeloid$interaction_name,decreasing = F),]
plotdata <- CCI_CAF_Myeloid[order(CCI_CAF_Myeloid$pathway_name,decreasing = F),]
source.targets <- unique(as.character(plotdata$source.target))
Tissues <- unique(as.character(plotdata$Tissue))
interaction_names <- unique(as.character(plotdata$interaction_name))
table(plotdata$interaction_name, plotdata$pathway_name)
plotdata <- plotdata %>%
  # filter(gene %in% genelist) %>%
  filter(pval < 0.05) %>%
  mutate(source.target = factor(source.target, levels = source.targets)) %>%
  # mutate(lof_gof = ifelse(grepl("GOF", cell_line), "GOF", "LOF")) %>% 
  mutate(Tissue = factor(Tissue, levels = Tissues)) %>% 
  # mutate(cell_line = gsub("_", " ", gsub("_GOF", "", cell_line))) %>% 
  mutate(interaction_name = factor(interaction_name, levels = interaction_names))
sort(table(plotdata$ligand))
ligand.df <- data.frame(table(plotdata$ligand))
plotdata <- plotdata[plotdata$ligand %in% ligand.df$Var1[ligand.df$Freq>100],]
ggplot(plotdata, aes(x = source.target, y = interaction_name, size = -log10(pmax(pval, 1e-10)) )) + #-log10(pval))) +
  geom_point(aes(fill = prob), pch = 21) +
  scale_fill_distiller("Commun. Prob.", palette = 1, values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = c(quantile(plotdata$prob[plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]],
                                                         0, na.rm = T), quantile(plotdata$prob[plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]], 1, na.rm = T)), guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +#, limits = LIMITS) +
  ggnewscale::new_scale_fill() +
  geom_point(data = plotdata[!plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)],], aes(fill = prob), pch = 21) +
  scale_fill_distiller("Commun. Prob.", palette = 3, values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = c(quantile(plotdata$prob[!plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]],
                                                         0, na.rm = T), quantile(plotdata$prob[!plotdata$Tissue %in% Tissues[c(1,3,5,7,9,11)]], 1, na.rm = T)), direction = 1, guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +#, limits = LIMITS) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  theme_bw() + 
  
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',#"bottom", 
        legend.key.width = unit(0.3, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        plot.title = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed")) +
  labs(size = "P value (-log10)") +
  guides(size = guide_legend(override.aes = list(pch = 16), title.position = "top", title.hjust = 0.5)) +
  ylab("") +
  xlab("") +
  facet_grid(. ~ Tissue, scales = "free_x", space = "free_x")
ggsave("/data/rluo4/All/Output/res_fig/CCI_CAF_Myeloid_dotplot.png", height = 19.75, width = 44, dpi = 500, device = "png")
ggsave("/data/rluo4/All/Output/res_fig/CCI_CAF_Myeloid_dotplot.pdf", height = 20.75, width = 42)

save(CCI_all, CCI_STM_Myeloid, CCI_CAF_Myeloid, file = paste0(data_path,'../CCI_ANXA1.rds'))

CCI_all_ANXA1 <- CCI_all#[CCI_all$ligand=='ANXA1',]
unique(CCI_all_ANXA1$target)# 22: Myeloid + TREG, TFH, MSC.ADIPO
CCI_all_ANXA1 <- CCI_all_ANXA1[CCI_all_ANXA1$source %in% Lineage$Minor_type[Lineage$Major_type=='Mesenchymal'],]
CCI_all_ANXA1 <- CCI_all_ANXA1[CCI_all_ANXA1$target %in% Lineage$Minor_type[Lineage$Major_type=='Myeloid'],]
unique(CCI_all_ANXA1$Tissue)# 11
unique(CCI_all_ANXA1$target)# 10
CCI_CAF_Myeloid <- NULL
for (organ in organ_all) {
  dataPlot <- CCI_all_ANXA1[CCI_all_ANXA1$Tissue==organ,]
  dataPlot <- dataPlot[dataPlot$Trace %in% Trace_Clusters[[organ]]$precancer,]
  CCI_CAF_Myeloid <<- rbind(CCI_CAF_Myeloid, dataPlot)
}
CCI_CAF_Myeloid <- CCI_all[CCI_all$source %in% Lineage$Minor_type[Lineage$Major_type=='Mesenchymal'] & CCI_all$target %in% unique(CCI_all_ANXA1$target),]

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) more specific CCI-figures
library(CellChat)
data_path = '/data/rluo4/All/Output/res_fig/'
# typeCCI_use <- c("CA_HPV_1_CA_HPV_2_CA_HPV_3_HSIL_HPV_1_clusterC8", "CRC1_8810_HTA11_10711_200000101113111_HTA11_99999974143_8462013111_clusterC1","C09_C43_clusterC11", "C08_C21_C43_C46_C57_EOLP-2_clusterC3", )
load(paste0(data_path,'../CCI_ANXA1.rds'))                        
##############################################################################################################################
# i)  ANNEXIN_CCI
##############################################################################################################################
# CCI_all_ANNEXIN <- CCI_all[CCI_all$ligand=='ANNEXIN',]
CCI_all_ANNEXIN <- CCI_all[CCI_all$pathway_name == 'ANNEXIN',]
unique(CCI_all_ANNEXIN$target)# 22: Myeloid + TREG, TFH, MSC.ADIPO
unique(CCI_all_ANNEXIN$interaction_name) # ANXA1_FPR1, ANXA1_FPR2
Lineage <- read.csv('/data/rluo4/All/Output/organ13-celltypes.csv')
# CCI_all_ANNEXIN <- CCI_all_ANNEXIN[CCI_all_ANNEXIN$target %in% Lineage$Minor_type[Lineage$Major_type=='Myeloid'],]
# CCI_all_ANNEXIN <- CCI_all_ANNEXIN[CCI_all_ANNEXIN$source=='STM',]
unique(CCI_all_ANNEXIN$Tissue)
typeCCI_use <- unique(CCI_all_ANNEXIN$DiseaseStage) # 32 TissueType
typeCCI_use_df <- CCI_all_ANNEXIN[!duplicated(CCI_all_ANNEXIN$DiseaseStage), ]
typeCCI_use_df$Index <- 1:nrow(typeCCI_use_df)

extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
typeCCI_use_df$typeCCI_use <- paste(typeCCI_use_df$Tissue, extract_last_element(typeCCI_use_df$Trace), sep = '-')

CCI_dir <- "/data/rluo4/All/Output/res_fig/ANNEXIN_CCI/"
if (!dir.exists(paste0(CCI_dir))){
  dir.create(paste0(CCI_dir))
}
setwd(CCI_dir)
getwd()

CellChat_ANNEXIN <- NULL
CellChat_ANXA1 <- vector("list", length(typeCCI_use_df$DiseaseStage))
names(CellChat_ANXA1) <- typeCCI_use_df$typeCCI_use

for (i in 1:nrow(typeCCI_use_df)) {  
  TissueType = typeCCI_use_df$typeCCI_use[i]
  DiseaseStage = typeCCI_use_df$DiseaseStage[i]
  # TissueType = 'CA_HPV_1_CA_HPV_2_CA_HPV_3_HSIL_HPV_1_HSIL_HPV_2_clusterC14'
  print(TissueType)
  setwd(CCI_dir)
  
  # cohort_directory <- paste0(data_path, '../CellChat/',TissueType, '_CCI')
  print(paste0(TissueType,' from : ', DiseaseStage))
  load(paste0(data_path, '../CellChat/', DiseaseStage, '_CCI.RData'))
  #展示每个亚群作为source的信号传递
  print(table(cellchat@idents))
  
  library(ArchR)
  # 我们也可以可视化由单个配体-受体对介导的细胞间通信。我们提供了一个函数 extractEnrichedLR 来提取给定信号通路的所有显著相互作用（L-R对）和相关的信号基因。
  library(CellChat)                      
  pathways.show = 'ANNEXIN'
  # pdf(paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.pdf'), height = 5 , width = length(unique(cellchat@idents))*0.35)                     
  p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
  # p1
  # dev.off()  
  ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.pdf'), height = 5, width = length(unique(cellchat@idents))*0.35, dpi = 500)
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  # show one ligand-receptor pair 'ANNEXIN-FPR1/2' # Hierarchy plot
  vertex.receiver = seq(1,4) # a numeric vector
  LR.show = 'ANXA1_FPR1'  #pairLR.CXCL[3,]    
  for(LR.show in pairLR.CXCL$interaction_name){
    pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 6.5)                    
    netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
    dev.off()
  }
  # LR.show <- 'ANXA1_FPR2'
  # if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>=1){
  #   pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 6.5)
  #   netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  #   dev.off()
  # }                     
  data.cci <- subsetCommunication(cellchat)
  data.cci <- data.cci[, c(1:7, 9)]#data.cci[,c(3,4,1,2)]
  # colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
  data.cci$DiseaseStage <- DiseaseStage
  data.cci$TissueType <- TissueType
  data.cci$Organ <- str_split(data.cci$TissueType, '-', simplify = T)[,1]
  
  CellChat_ANNEXIN <<- rbind(CellChat_ANNEXIN, data.cci)
  CellChat_ANXA1[[TissueType]] <- cellchat
  
}

##############################################################################################################################
# ii)  TGFb_CCI
##############################################################################################################################
# CCI_all_TGFb <- CCI_all[CCI_all$ligand=='TGFb',]
CCI_all_TGFb <- CCI_all[CCI_all$pathway_name == 'TGFb',]
unique(CCI_all_TGFb$target)# 22: Myeloid + TREG, TFH, MSC.ADIPO
unique(CCI_all_TGFb$interaction_name) # ANXA1_FPR1, ANXA1_FPR2
Lineage <- read.csv('/data/rluo4/All/Output/organ13-celltypes.csv')
# CCI_all_TGFb <- CCI_all_TGFb[CCI_all_TGFb$target %in% Lineage$Minor_type[Lineage$Major_type=='Myeloid'],]
# CCI_all_TGFb <- CCI_all_TGFb[CCI_all_TGFb$source=='STM',]
unique(CCI_all_TGFb$Tissue)
typeCCI_use <- unique(CCI_all_TGFb$DiseaseStage) # 32 TissueType
typeCCI_use_df <- CCI_all_TGFb[!duplicated(CCI_all_TGFb$DiseaseStage), ]
typeCCI_use_df$Index <- 1:nrow(typeCCI_use_df)

extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
typeCCI_use_df$typeCCI_use <- paste(typeCCI_use_df$Tissue, extract_last_element(typeCCI_use_df$Trace), sep = '-')

CCI_dir <- "/data/rluo4/All/Output/res_fig/TGFb_CCI/"
if (!dir.exists(paste0(CCI_dir))){
  dir.create(paste0(CCI_dir))
}
setwd(CCI_dir)
getwd()

CellChat_TGFb <- NULL
CellChat_TGFB <- vector("list", length(typeCCI_use_df$DiseaseStage))
names(CellChat_TGFB) <- typeCCI_use_df$typeCCI_use

for (i in 1:nrow(typeCCI_use_df)) {  
  TissueType = typeCCI_use_df$typeCCI_use[i]
  DiseaseStage = typeCCI_use_df$DiseaseStage[i]
  # TissueType = 'CA_HPV_1_CA_HPV_2_CA_HPV_3_HSIL_HPV_1_HSIL_HPV_2_clusterC14'
  print(TissueType)
  setwd(CCI_dir)
  
  # cohort_directory <- paste0(data_path, '../CellChat/',TissueType, '_CCI')
  print(paste0(TissueType,' from : ', DiseaseStage))
  load(paste0(data_path, '../CellChat/', DiseaseStage, '_CCI.RData'))
  #展示每个亚群作为source的信号传递
  print(table(cellchat@idents))
  
  library(ArchR)
  # 我们也可以可视化由单个配体-受体对介导的细胞间通信。我们提供了一个函数 extractEnrichedLR 来提取给定信号通路的所有显著相互作用（L-R对）和相关的信号基因。
  library(CellChat)                      
  pathways.show = 'TGFb'
  # pdf(paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.pdf'), height = 5 , width = length(unique(cellchat@idents))*0.35)                     
  p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
  # p1
  # dev.off()  
  ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.pdf'), height = 5, width = length(unique(cellchat@idents))*0.35, dpi = 500)
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  # show one ligand-receptor pair 'TGFb-ITGAV+ITGB1' # Hierarchy plot
  vertex.receiver = seq(1,4) # a numeric vector
  # LR.show = 'ANXA1_FPR1'  #pairLR.CXCL[3,]    
  for(LR.show in pairLR.CXCL$interaction_name){
    pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 6.5)                    
    netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
    dev.off()
  }
  # LR.show <- 'ANXA1_FPR2'
  # if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>=1){
  #   pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 6.5)
  #   netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  #   dev.off()
  # }                     
  data.cci <- subsetCommunication(cellchat)
  data.cci <- data.cci[, c(1:7, 9)]#data.cci[,c(3,4,1,2)]
  # colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
  data.cci$DiseaseStage <- DiseaseStage
  data.cci$TissueType <- TissueType
  data.cci$Organ <- str_split(data.cci$TissueType, '-', simplify = T)[,1]
  
  CellChat_TGFb <<- rbind(CellChat_TGFb, data.cci)
  CellChat_TGFB[[TissueType]] <- cellchat
  
}
save(CellChat_ANNEXIN, CellChat_ANXA1, CellChat_TGFb, CellChat_TGFB, file =  paste0(data_path,'../CCI_TGFb.rds'))
load(paste0(data_path,'../CCI_TGFb.rds'))

##############################################################################################################################
# iii)  FN1 - CCI
##############################################################################################################################
typeCCI_use <- unique(CCI_all$DiseaseStage[CCI_all$source == 'STM' & CCI_all$ligand == 'FN1' & 
                                             CCI_all$receptor %in% c('ITGAV_ITGB1', 'SDC1')]) #cirrhotic1_cirrhotic2_cirrhotic3_clusterC15 has SDC1 receptor interaction, so we will skip this one to show the difference between benign and malignant clusters ,#c('SDC1','ITGAV_ITGB1')])                         
print(typeCCI_use)
CCI_dir <- "/data/rluo4/All/Output/res_fig/FN1_CCI/"
if (!dir.exists(paste0(CCI_dir))){
  dir.create(paste0(CCI_dir))
}
setwd(CCI_dir)
getwd()
typeCCI_use_df <- CCI_all %>% filter(source == 'STM' & ligand == 'FN1' & receptor %in% c('ITGAV_ITGB1', 'SDC1'))
typeCCI_use_df$Index <- 1:nrow(typeCCI_use_df)

extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
typeCCI_use_df$typeCCI_use <- paste(typeCCI_use_df$Tissue, extract_last_element(typeCCI_use_df$Trace), sep = '-')

CellChat_FN1 <- NULL
for (i in unique(typeCCI_use_df$typeCCI_use)) {  
  TissueType = i#typeCCI_use_df$typeCCI_use[i]
  DiseaseStage = unique(typeCCI_use_df$DiseaseStage[typeCCI_use_df$typeCCI_use==i])
  print(TissueType)
  setwd(CCI_dir)
  # cohort_directory <- paste0(data_path, '../CellChat/',TissueType, '_CCI')
  print(paste0(TissueType,' from : ', DiseaseStage))
  load(paste0(data_path, '../CellChat/', DiseaseStage, '_CCI.RData'))
  #展示每个亚群作为source的信号传递
  print(table(cellchat@idents))
  library(ArchR)
  #     我们也可以可视化由单个配体-受体对介导的细胞间通信。我们提供了一个函数 extractEnrichedLR 来提取给定信号通路的所有显著相互作用（L-R对）和相关的信号基因。
  library(CellChat)                      
  pathways.show = 'FN1'
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  print(pairLR.CXCL)
  # show one ligand-receptor pair 'FN1-ITGAV+ITGB1' # Hierarchy plot
  levels(cellchat@idents)
  sources.use = na.omit( match( unique(typeCCI_use_df$source[typeCCI_use_df$typeCCI_use==TissueType]), levels(cellchat@idents))  )
  vertex.receiver = na.omit( match( unique(typeCCI_use_df$target[typeCCI_use_df$typeCCI_use==TissueType]), levels(cellchat@idents))  )
  
  p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
  # dev.off()  
  # ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.pdf'), height = nrow(pairLR.CXCL), width = length(unique(cellchat@idents))*0.35, dpi = 500)
  ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.png'), height = nrow(pairLR.CXCL)*0.8, width = length(unique(cellchat@idents))*0.28, dpi = 500)
  
  LR.show = 'FN1_ITGAV_ITGB1'  #pairLR.CXCL[3,]    
  if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>1 ){
    # pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 7)
    png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = 5.5, width = 5.5, units = 'in',  res = 500)
    par(mfrow = c(1,1))
    netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 0.5, #vertex.weight = 15, edge.weight.max = 10,
                         layout = "circle" , vertex.receiver = vertex.receiver)
    dev.off()
  }   
  LR.show <- 'FN1_SDC1'#pairLR.CXCL[5,]
  if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>1){
    # pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 7)
    png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = 5.5, width = 5.5, units = 'in',  res = 500)
    par(mfrow = c(1,1))
    netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 0.5, #vertex.weight = 15, edge.weight.max = 10,
                         layout = "circle" , vertex.receiver = vertex.receiver)
    dev.off()
  }   
  
  # pairLR.use = pairLR.CXCL[pairLR.CXCL$interaction_name == LR.show, , drop = FALSE]
  pairLR.use = pairLR.CXCL[pairLR.CXCL$interaction_name %in% c('FN1_ITGAV_ITGB1', 'FN1_SDC1'), , drop = FALSE]
  png(paste0(CCI_dir, TissueType, '_', pathways.show, '.png'), width = 5, height = 5, units = 'in', res = 500)
  # par(mfrow = c(1,1))
  a <- netVisual_chord_gene(cellchat, targets.use = vertex.receiver, #sources.use =sources.use,
                            pairLR.use = pairLR.use, legend.pos.x = 8, legend.pos.y = 40#signaling = pathways.show,  
  )
  print(a)
  dev.off()
  png(paste0(CCI_dir, TissueType, '_STM_', pathways.show, '.png'), width = 4, height = 4, units = 'in', res = 500)
  par(mfrow = c(1,1))
  a <- netVisual_chord_gene(cellchat, targets.use = vertex.receiver, sources.use =sources.use,
                            pairLR.use = pairLR.use, legend.pos.x = 10, legend.pos.y = -4#signaling = pathways.show,  
  )
  print(a)
  dev.off()
  data.cci <- subsetCommunication(cellchat)
  data.cci <- data.cci[, c(1:7, 9)]#data.cci[,c(3,4,1,2)]
  # colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
  data.cci$DiseaseStage <- DiseaseStage
  data.cci$TissueType <- TissueType
  data.cci$Organ <- str_split(data.cci$TissueType, '-', simplify = T)[,1]
  CellChat_FN1 <<- rbind(CellChat_FN1, data.cci)
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# iv)  NECTIN - CCI
##############################################################################################################################
TissueType = 'brca10_brca2_brca3_DCIS2_M1_NCCBC14_NCCBC3_NCCBC5_P1_P2_clusterC4'
print(TissueType)
names(TissueType)
cohort_directory <- paste0(data_path,'/', TissueType, '_CCI')
print(paste0(TissueType, ' from dir: ', cohort_directory))
load(paste0(data_path,'/../CellChat/',TissueType, '_CCI.RData'))
#展示每个亚群作为source的信号传递
print(table(cellchat@idents))
# subtype <- data.frame(subset =  table(cellchat@idents))
# subtype <- subtype[subtype$subset.Freq !=0,]
# subtype <- subtype[subtype$subset.Freq >=10,]
# table(cellchat@idents)
# cellchat <-  subsetCellChat(cellchat, idents.use = subtype$subset.Var1)
library(ArchR)
#     我们也可以可视化由单个配体-受体对介导的细胞间通信。我们提供了一个函数 extractEnrichedLR 来提取给定信号通路的所有显著相互作用（L-R对）和相关的信号基因。
library(CellChat)                      
pathways.show = 'NECTIN'
CCI_dir <- paste0("/data/rluo4/All/Output/res_fig/Breast_CCI_clusterC4/")
if (!dir.exists(paste0('',CCI_dir))){
  dir.create(paste0(CCI_dir))
}
setwd(CCI_dir)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
# pdf(paste0(CCI_dir, TissueType, '_', pathways.show, '.pdf'), height = nrow(pairLR.CXCL)*1.35, width = length(unique(cellchat@idents))*0.35)                     
p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
# p1
# dev.off() 
# ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.pdf'), height = nrow(pairLR.CXCL), width = length(unique(cellchat@idents))*0.35, dpi = 500)
ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.png'), height = nrow(pairLR.CXCL)*1.2, width = length(unique(cellchat@idents))*0.28, dpi = 500)

# show one ligand-receptor pair 'FN1-ITGAV+ITGB1' # Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
pairLR.CXCL
LR.show = 'NECTIN3_TIGIT'  #pairLR.CXCL[3,]    
pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 8.5, width = 8.5)  
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()

LR.show <- 'NECTIN2_TIGIT'#pairLR.CXCL[5,]
pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 8.5, width = 8.5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()

table(cellchat@idents)   
levels(cellchat@idents)   
vertex.receiver = c(4,6,15) 
sources.use = c(16, 23, 26, 33)
netVisual_aggregate(cellchat, signaling = pathways.show,# sources.use = sources.use,
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = sources.use, targets.use = vertex.receiver,  slot.name = "netP", legend.pos.x = 10)
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# group.cellType <- levels(cellchat@idents)[c(vertex.receiver, sources.use)]
# names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, sources.use = sources.use,  #group = group.cellType, 
                     targets.use = vertex.receiver, title.name = paste0(pathways.show, " signaling network"))
png(filename=paste0(CCI_dir, pathways.show, '_', LR.show, '.png'), width = 5.5, height = 4.5, units = 'in', res = 500)
netVisual_chord_gene(cellchat, sources.use = sources.use, targets.use = vertex.receiver, signaling = pathways.show, legend.pos.x = 6, legend.pos.y = 40)
dev.off()
# ggsave(filename=paste0(CCI_dir, pathways.show, '_', LR.show, '.png'), width = 4, height = 3, units = 'in', dpi = 500)

cellchat@netP$pathways
pathways.show <- c("TGFb")  
library(CellChat)                      
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
pairLR.CXCL
p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
# p1
# dev.off() 
ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.png'), height = nrow(pairLR.CXCL)*1.2, width = length(unique(cellchat@idents))*0.28, dpi = 500)

# show one ligand-receptor pair 'FN1-ITGAV+ITGB1' # Hierarchy plot
# vertex.receiver = seq(1,4) # a numeric vector
pairLR.CXCL
LR.show = 'TGFB1_TGFBR1_TGFBR2'  #pairLR.CXCL[3,]    
# LR.show = 'TGFB3_TGFBR1_TGFBR2'  #pairLR.CXCL[3,]    
# pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 8.5, width = 8.5)  
png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = 5.5, width = 5.5, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 0.4, #vertex.weight = 15, edge.weight.max = 10,
                     layout = "chord" )
dev.off()

LR.show <- 'TGFB1_ACVR1B_TGFBR2'#pairLR.CXCL[5,]
# pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 8.5, width = 8.5)
png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = 5.5, width = 5.5, units = 'in',  res = 500)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

# levels(cellchat@idents)   
# sources.use = which(levels(cellchat@idents) %in% c('PVA','FIB','INMON','TH17','CD8TRM',
#                                                'DC','CD8TEREX','CD8TEXP','pDC','CD8TEREXINT','CD8TCM',
#                                                'M1MAC' , 'CD8TEX')) 
# vertex.receiver = which(levels(cellchat@idents) %in% c('ECM'))
# netVisual_aggregate(cellchat, signaling = pathways.show,# sources.use = sources.use,
#                     vertex.receiver = vertex.receiver,layout = "hierarchy")
# # show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# netVisual_chord_gene(cellchat, sources.use = sources.use, targets.use = vertex.receiver,  slot.name = "netP", legend.pos.x = 10)
# # Chord diagram
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# # group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# # group.cellType <- levels(cellchat@idents)[c(vertex.receiver, sources.use)]
# # names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, sources.use = sources.use,  #group = group.cellType, 
#                      targets.use = vertex.receiver, title.name = paste0(pathways.show, " signaling network"))
# png(filename=paste0(CCI_dir, pathways.show, '_', LR.show, '.png'), width = 5.5, height = 4.5, units = 'in', res = 500)
# netVisual_chord_gene(cellchat, sources.use = sources.use, targets.use = vertex.receiver, signaling = pathways.show, legend.pos.x = 6, legend.pos.y = 40)
# dev.off()

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# v)  TNF-TNFRSF1A, CCL2-ACKR1, CXCL2-ACKR1, CXCL12_ACKR3
##############################################################################################################################
data_path = '/data/rluo4/All/Output/res_fig/'
typeCCI_use <- unique(CCI_all$DiseaseStage[CCI_all$source == 'STM' & CCI_all$ligand %in% c('TNF', 'CCL2', 'CXCL2', 'CXCL12') & 
                                             CCI_all$receptor %in% c('TNFRSF1A', 'ACKR1', 'ACKR3')]) #cirrhotic1_cirrhotic2_cirrhotic3_clusterC15 has SDC1 receptor interaction, so we will skip this one to show the difference between benign and malignant clusters ,#c('SDC1','ITGAV_ITGB1')])                         
print(typeCCI_use)
CCI_dir <- "/data/rluo4/All/Output/res_fig/Inflam_CCI/"
if (!dir.exists(paste0(CCI_dir))){
  dir.create(paste0(CCI_dir))
}
setwd(CCI_dir)
getwd()
typeCCI_use_df <- CCI_all %>% filter(# source == 'STM' & 
  CCI_all$ligand %in% c('TNF', 'CCL2', 'CXCL2', 'CXCL12') & 
    CCI_all$receptor %in% c('TNFRSF1A', 'ACKR1', 'CXCR4', 'ACKR3'))
typeCCI_use_df$Index <- 1:nrow(typeCCI_use_df)

extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
typeCCI_use_df$typeCCI_use <- paste(typeCCI_use_df$Tissue, extract_last_element(typeCCI_use_df$Trace), sep = '-')

CellChat_inflam <- NULL
for (i in unique(typeCCI_use_df$typeCCI_use)) {  
  TissueType = i#typeCCI_use_df$typeCCI_use[i]
  DiseaseStage = unique(typeCCI_use_df$DiseaseStage[typeCCI_use_df$typeCCI_use==i])
  print(TissueType)
  setwd(CCI_dir)
  # cohort_directory <- paste0(data_path, '../CellChat/',TissueType, '_CCI')
  print(paste0(TissueType,' from : ', DiseaseStage))
  load(paste0(data_path, '../CellChat/', DiseaseStage, '_CCI.RData'))
  #展示每个亚群作为source的信号传递
  print(table(cellchat@idents))
  library(ArchR)
  #     我们也可以可视化由单个配体-受体对介导的细胞间通信。我们提供了一个函数 extractEnrichedLR 来提取给定信号通路的所有显著相互作用（L-R对）和相关的信号基因。
  library(CellChat)                      
  for(pathways.show in c('TNF', 'CCL', 'CXCL')){
    
    paths <- unique(typeCCI_use_df$pathway_name[typeCCI_use_df$typeCCI_use==TissueType])
    if(pathways.show %in% paths){
      print(pathways.show)
      pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
      print(pairLR.CXCL)
      # show one ligand-receptor pair 'FN1-ITGAV+ITGB1' # Hierarchy plot
      levels(cellchat@idents)
      sources.use = na.omit( match( unique(typeCCI_use_df$source[typeCCI_use_df$typeCCI_use==TissueType]), levels(cellchat@idents))  )
      vertex.receiver = na.omit( match( unique(typeCCI_use_df$target[typeCCI_use_df$typeCCI_use==TissueType]), levels(cellchat@idents))  )
      if(nrow(pairLR.CXCL)>=1 ){
        p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
        hei <- ifelse(nrow(pairLR.CXCL) <=3, nrow(pairLR.CXCL)*6,  nrow(pairLR.CXCL)*0.6)
        
        ggsave(plot = p1, file = paste0(CCI_dir, TissueType, '_', pathways.show, '_expression.png'), height = hei, width = length(unique(cellchat@idents))*0.28, dpi = 500)
        
        hei = wid = nlevels(cellchat@idents)*0.28
        
        if(pathways.show == 'TNF'){
          LR.show = 'TNF_TNFRSF1A'  #pairLR.CXCL[3,]    
          if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>=1){
            # pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 7)
            png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = hei, width = wid, units = 'in',  res = 500)
            par(mfrow = c(1,1))
            netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 1, edge.width.max = 1, point.size = 2.5,
                                 vertex.size.max = 1, vertex.label.cex = 1, layout = "circle" , vertex.receiver = vertex.receiver)
            dev.off()
          }  
        }
        if(pathways.show == 'CCL'){
          
          LR.show = 'CCL2_ACKR1'  #pairLR.CXCL[3,]    
          if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>1){
            # pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 7)
            png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = hei, width = wid, units = 'in',  res = 500)
            par(mfrow = c(1,1))
            netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 1, edge.width.max = 1, point.size = 2.5,
                                 vertex.size.max = 1, vertex.label.cex = 1, layout = "circle" , vertex.receiver = vertex.receiver)
            dev.off()
          }  
        }
        
        if(pathways.show == 'CXCL'){
          
          LR.show = 'CXCL2_ACKR1'  #pairLR.CXCL[3,]    
          if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>1 ){
            # pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 7)
            png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = hei, width = wid, units = 'in',  res = 500)
            par(mfrow = c(1,1))
            netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 1, edge.width.max = 1, point.size = 2.5,
                                 vertex.size.max = 1, vertex.label.cex = 1, layout = "circle" , vertex.receiver = vertex.receiver)
            dev.off()
          }   
          LR.show <- 'CXCL12_CXCR4'#pairLR.CXCL[5,]
          if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>1){
            # pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 7)
            png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = hei, width = wid, units = 'in',  res = 500)
            par(mfrow = c(1,1))
            netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 1, edge.width.max = 1, point.size = 2.5,
                                 vertex.size.max = 1, vertex.label.cex = 1, layout = "circle" , vertex.receiver = vertex.receiver)
            dev.off()
          }  
          LR.show <- 'CXCL12_ACKR3'#pairLR.CXCL[5,]
          if(length(unique(grepl(LR.show, pairLR.CXCL[,1])))>1){
            # pdf(paste0(CCI_dir, TissueType, '_', LR.show, '.pdf'), height = 7, width = 7)
            png(paste0(CCI_dir, TissueType, '_', LR.show, '.png'), height = hei, width = wid, units = 'in',  res = 500)
            par(mfrow = c(1,1))
            netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, arrow.size = 1, edge.width.max = 1, point.size = 2.5,
                                 vertex.size.max = 1, vertex.label.cex = 1, layout = "circle" , vertex.receiver = vertex.receiver)
            dev.off()
          }
        }
      }
    }
    
  }
  data.cci <- subsetCommunication(cellchat)
  data.cci <- data.cci[, c(1:7, 9)]#data.cci[,c(3,4,1,2)]
  # colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
  data.cci$DiseaseStage <- DiseaseStage
  data.cci$TissueType <- TissueType
  data.cci$Organ <- str_split(data.cci$TissueType, '-', simplify = T)[,1]
  CellChat_inflam <<- rbind(CellChat_inflam, data.cci)
}
load(file = paste0('/data/rluo4/All/Output/CellChat_inflam.rds'))

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 13) Check the gene expression of ANNEXIN in ST data
# load( file = paste0('/data/rluo4/All/Output/ANXA1_exp_df.rds'))
markers_dir <- "/data/rluo4/All/Output/res_fig/ST_markers_fig/"

Lineage <- read.csv('/data/rluo4/All/Output/organ13-celltypes.csv')
markers <- c('ANXA1','FPR1','FPR2','FPR3', 
             'CTGF','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','ACVR1B', 'ACVR1',
             'FAP','MMP1','ACTA2','MYH11','COL6A3')#MP_LR$MP_Gene[MP_LR$MPs==mp]
data_path = '/data/rluo4/All/Output/res_fig/'
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}

organ = 'Breast'
marker_Breast <- as.data.frame(matrix(NA, nrow = 0, ncol = (length(markers) + 2)))
colnames(marker_Breast) <- c(markers, 'Cell_Type', 'TissueType')
for (sample in samples_Breast) {
  min.cells <- 10
  n <- paste0(sample, '_obj')
  seu <- SPOTlight_obj_Breast[[n]]#PCT_seu[[organ]]
  table(Idents(seu))
  genes <- intersect(markers, rownames(seu@assays$SCT$data))
  genes
  
  PCT_data <- as.data.frame(t(seu@assays$SCT$data[genes, ]))
  # PCT_data$cell_type <- seu$labels
  PCT_data$Cell_Type <- seu$labels
  # SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj.rds'))
  # PCT_data$spotCoordinates <- rownames(SpaCET_obj@input$spotCoordinates)
  # 
  subtype <- data.frame(subset =  table(PCT_data$Cell_Type))
  subtype <- subtype[subtype$subset.Freq !=0,]
  subtype <- subtype[subtype$subset.Freq >=5,]
  PCT_data$Cell_Type <- as.character(PCT_data$Cell_Type)
  table(PCT_data$Cell_Type)
  extract_last_element <- function(x) {
    split_string <- strsplit(x, "-")
    last_element <- sapply(split_string, function(y) tail(y, n = 1))
    return(last_element)
  }
  celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
  celltype_choose <- c(extract_last_element(PCT_Clusters[[organ]]),celltype_choose )
  celltype_choose
  PCT_data <- PCT_data  %>% 
    filter(Cell_Type %in% subtype$subset.Var1) %>% #& MP %in% c('Epithelial-senescence')) %>%
    filter(Cell_Type %in%  celltype_choose)
  table( PCT_data$Cell_Type)
  level_compare <- celltype_choose[celltype_choose %in% PCT_data$Cell_Type]
  PCT_data$Cell_Type <- factor(PCT_data$Cell_Type, levels = level_compare )
  PCT_data$TissueType = sample
  # p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
  # View(paletteer::palettes_d_names)
  coln <- setdiff(colnames(marker_Breast), colnames(PCT_data))
  if(length(coln)>=1){
    PCT_data[coln] <- NA
  }
  marker_Breast <<- rbind(marker_Breast, PCT_data[, colnames(marker_Breast)])
  # marker_Breast <<-  marker_Breast %>%
  #   inner_join(PCT_data, by = coln)
  
  # markers_dir <- "/data/rluo4/All/Output/res_fig/ST_markers_fig/"
  # for (gene in genes) {
  #   # gene = 'FPR1'#'ANXA1'#'FPR1'
  #   # any_zero <- any(PCT_data[[gene]] == 0)
  #   # Filter out groups with all zero values
  #   # Convert 'gene' column to numeric if it's stored as character
  #   PCT_data[[gene]] <- as.numeric(PCT_data[[gene]])
  #   # Define the formula string for the t-test
  #   formula_str <- reformulate("Cell_Type", response = gene)
  #   stat.test <- PCT_data %>%
  #     # group_by(Cell_Type) %>%
  #     t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
  #     adjust_pvalue(method = "bonferroni") %>%
  #     add_significance()
  #   # stat.test <- PCT_data %>%
  #   #   group_by(Cell_Type) %>%
  #   #   t_test(TGFB1 ~ Cell_Type, ref.group = "Breast-C14")
  #   stat.test <- stat.test %>% add_y_position(fun = "mean_se")
  #   bp <- ggbarplot(
  #     PCT_data, x = "Cell_Type", y = gene, fill = 'Cell_Type',  package = "palettetown", #View(paletteer::palettes_d_names)
  #     palette = "quilava",
  #     add = c('mean_se'), #facet.by = "Cell_Type"
  #   ) + xlab("") + #ylab("") +
  #     stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  #     scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
  #       text = element_text(size = 12),
  #       panel.grid.major = element_line(colour = "grey90", size=0.2),
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_blank(),
  #       axis.line = element_line(colour = "black"),
  #       # axis.text.x = element_blank(),
  #       axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
  #       axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
  #       legend.position="none",
  #       legend.title=element_text(size=12)) 
  #   bp
  #   ggsave(plot = bp, file = paste0(markers_dir, sample, '_', gene, '_expression.png'), height =  4, width = length(unique(PCT_data$Cell_Type))*0.6, dpi = 500)
  # }
}

unique(marker_Breast$TissueType)
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
celltype_choose <- c(extract_last_element(PCT_Clusters[[organ]]),celltype_choose )
level_compare <- celltype_choose[celltype_choose %in% marker_Breast$Cell_Type]
# marker_Breast$Cell_Type <- factor(marker_Breast$Cell_Type, levels = level_compare )
marker_Breast$Cell_Type <- as.character(marker_Breast$Cell_Type)
unique(marker_Breast$Cell_Type)
marker_Breast$Lineage <-   Lineage$Major_type[match(marker_Breast$Cell_Type, Lineage$Minor_type)]
marker_Breast$Lineage[is.na(marker_Breast$Lineage)] <- 'STM'
# marker_Breast$Cell_Type[marker_Breast$Lineage=='STM'] <- factor(marker_Breast$Cell_Type[marker_Breast$Lineage=='STM'], levels = extract_last_element(PCT_Clusters[[organ]]) )
for (gene in genes) {
  # gene = 'ANXA1'#'FPR1'
  marker_Breast[[gene]][is.na(marker_Breast[[gene]])] = 0
  marker_Breast[[gene]] <- as.numeric(marker_Breast[[gene]])
  
  # for (majortype in unique(marker_Breast$Lineage)) {
  ggplotdata <- marker_Breast #%>% filter(Lineage==majortype)
  levels <- c(extract_last_element(PCT_Clusters[[organ]]), 
              unique(marker_Breast$Cell_Type[marker_Breast$Lineage=='Myeloid']),
              unique(marker_Breast$Cell_Type[marker_Breast$Lineage=='Mesenchymal']))
  ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = levels)
  # Define the formula string for the t-test
  formula_str <- reformulate("Cell_Type", response = gene)
  stat.test <- ggplotdata %>%
    # group_by(Lineage) %>%
    t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
  # stat.test <- PCT_data %>%
  #   group_by(Cell_Type) %>%
  #   t_test(TGFB1 ~ Cell_Type, ref.group = "Breast-C14")
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  
  stat.test <- stat.test %>% add_y_position(fun = "mean_se")
  bp <- ggbarplot(
    ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "rickandmorty",#"palettetown", #View(paletteer::palettes_d_names)
    #palette = "quilava",
    add = c('mean_se'), #facet.by = "Lineage"
  ) + xlab("") + #ylab("") +
    # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
      text = element_text(size = 12),
      panel.grid.major = element_line(colour = "grey90", size=0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # axis.text.x = element_blank(),
      axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
      legend.position="none",
      legend.title=element_text(size=12))
  bp
  ggsave(plot = bp, file = paste0(markers_dir, organ, '_ST_',  gene, '_expression.png'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
}

for (gene in genes) {
  # gene = 'FPR1'#'ANXA1'#'FPR1'
  # any_zero <- any(PCT_data[[gene]] == 0)
  # Filter out groups with all zero values
  # Convert 'gene' column to numeric if it's stored as character
  marker_Breast[[gene]][is.na(marker_Breast[[gene]])] = 0
  marker_Breast[[gene]] <- as.numeric(marker_Breast[[gene]])
  
  for (majortype in unique(marker_Breast$Lineage)) {
    ggplotdata <- marker_Breast %>% filter(Lineage==majortype)
    if(majortype=='STM'){
      ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = extract_last_element(PCT_Clusters[[organ]]) )
    }
    
    # Define the formula string for the t-test
    formula_str <- reformulate("Cell_Type", response = gene)
    stat.test <- ggplotdata %>%
      # group_by(Lineage) %>%
      t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    # stat.test <- PCT_data %>%
    #   group_by(Cell_Type) %>%
    #   t_test(TGFB1 ~ Cell_Type, ref.group = "Breast-C14")
    color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
    
    stat.test <- stat.test %>% add_y_position(fun = "mean_se")
    bp <- ggbarplot(
      ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "jco",#"palettetown", #View(paletteer::palettes_d_names)
      #palette = "quilava",
      add = c('mean_se'), #facet.by = "Lineage"
    ) + xlab("") + #ylab("") +
      stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
        text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_blank(),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="none",
        legend.title=element_text(size=12))
    bp
    ggsave(plot = bp, file = paste0(markers_dir, organ, '_', majortype, '_', gene, '_expression.png'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
  }
}
organ = 'Breast'
sample = 'Visium_BC'
n <- paste0(sample, '_obj')
seu <- SPOTlight_obj_Breast[[n]]#PCT_seu[[organ]]
table(Idents(seu))
outputDir = '/data/rluo4/All/Output/Epi_Results'
outfile <- paste0(outputDir,"/",organ,'_ST_ssgsea.rds')
# if (! file.exists(outfile)) {
#   save(ssgsea.res, file = outfile)
# }
load(outfile)
dim(ssgsea.res)
ssgsea.res[1:5,1:5]
set.seed(123)
group_ABN <- seu@meta.data
table(group_ABN$labels)
table(colnames(ssgsea.res) == rownames(group_ABN))

library(edgeR)
library(limma)
group_ABN <- group_ABN[group_ABN$labels %in% c('C4', 'C9'),]
ssgsea.res <- ssgsea.res[,colnames(ssgsea.res) %in% rownames(group_ABN)]
group <- as.factor(ifelse(group_ABN$labels %in% 'C9', 
                          'Bad','Good'))
# subtype <- data.frame(subset = subtype)
# subtype <- subtype[subtype$subset.Freq !=0,]
table(group)
desigN <- model.matrix(~group  + 0)
head(desigN)
# Activation of the UPR can enhance the cell's ability to cope with nutrient deprivation, hypoxia, and other stressors commonly encountered in solid tumors.
rownames(desigN) <- rownames(group_ABN)
head(desigN)
comparE <- makeContrasts(groupBad - groupGood, levels=desigN)
enrich.res<-ssgsea.res[,colnames(ssgsea.res) %in% rownames(group_ABN)]
# enrich.res<-gsva.res[,colnames(gsva.res) %in% rownames(group_ABN)]
fiT <- lmFit(enrich.res, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
veloGood <- topTable(fiT3, coef=1, number=200)
head(veloGood, n=3)
max(abs(veloGood$logFC))
dat.enrich=enrich.res
df.res=veloGood
df.res=df.res[df.res$P.Value<0.05 ,] 
df.res$Pathways <- rownames(df.res)
df.res$organ <- organ
df.res$cluster[df.res$t>0] <- 'C9'
df.res$cluster[df.res$t<0] <- 'C4'
## barplot
# df.res=df.res[df.res$adj.P.Val<0.05 ,]  # sort according to p-value
dat_plot <- data.frame(id = row.names(df.res),
                       t = df.res$t)
dat_plot<-dat_plot[order(dat_plot$t,decreasing = T),]#[1:15,]
top_up <- head(dat_plot$id,15)
top_down <- tail(dat_plot$id,20)[-c(18:20)]

my_genesets <- read.csv('/data/rluo4/All/my_genesets.csv')
my_genesets <- my_genesets[-grep("KEGG", my_genesets$gs_name),-1]
# my_genesets <- rbind(Senescense, Metaplasia, H_genesets, cellcycle.gmt, reactome.gmt)
colnames(my_genesets) <- c("term","gene")
my_genesets$term <- as.factor(my_genesets$term)
unique(my_genesets$term)
my.gmt <- split(my_genesets$gene, my_genesets$term)

include_path <- c('EpiSen','SenMayo', 'FRIDMAN_SENESCENCE_UP','CancerG0Arrest','HALLMARK_HYPOXIA', 
                  'REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE','REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP',
                  'REACTOME_CELLULAR_SENESCENCE','REACTOME_ONCOGENE_INDUCED_SENESCENCE',
                  'REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE', 
                  as.character(unique(my_genesets$term[grep('EGFR', my_genesets$term)])),
                  as.character(unique(my_genesets$term[grep('AUTOPHAGY', my_genesets$term)])),
                  as.character(unique(my_genesets$term[grep('TGF_BETA', my_genesets$term)]))
)
dat_plot <- dat_plot[dat_plot$id %in% c(top_up,top_down,include_path),]
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
dat_plot$id <- str_replace(dat_plot$id , "REACTOME_","")
dat_plot$id <- gsub('EpiSen','EPITHELIAL SENESCENCE',dat_plot$id)
dat_plot$id <- gsub('Metaplasia','METAPLASIA AND DAMAGE RESPONSE',dat_plot$id)

shorten_names <- function(x, n_word=4){
  if (length(strsplit(x, "_")[[1]]) > n_word )
    x <- paste(paste(strsplit(x, "_")[[1]][1:n_word], collapse="_"),
               paste(strsplit(x, "_")[[1]][(n_word+1):length(strsplit(x,"_")[[1]])],
                     collapse="_"), sep="\n_")
  return(x)
} 
shorten_names <- function(x, n_word=4, n_char=50){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 50))
  {
    if (nchar(x) > 50) x <- substr(x, 1, 50)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}
dat_plot$term <- dat_plot$id
dat_plot$id <- dat_plot$term
dat_plot$id = sapply(dat_plot$id ,shorten_names)
dat_plot <- dat_plot[! dat_plot$id %in% c('PD_1_SIGNALING', 'ADAPTIVE_IMMUNE_SYSTEM',#'FORMATION_OF_THE_CORNIFIED_ENVELOPE',
                                          'KERATINIZATION','DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING'),]

# dat_plot$id <- shorten_names(dat_plot$id)
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- gsub('_',' ',dat_plot$id)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
color = c('Up'= '#f89e81','NoSignifi'='#cccccc','Down'='#7bcd7b')

ylab = paste0('GSVA: t value in Breast ST data (', 'C9', ' v.s ', 'C4', ')')#'(HCC:STM:19 v.s HCC:STM:7)'
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = color ) +
  # scale_fill_manual(values = c('Up'= '#C17188','NoSignifi'='#cccccc','Down'='#748EC2')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme( 
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 13,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = 'black', size   = 11,face = 'italic', hjust = 1),
    axis.ticks.y = element_blank(),
    legend.title  = element_text(color = 'black', size  = 16),
    legend.text   = element_text(color = 'black', size   = 16),
    # axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    panel.border = element_blank()#element_rect(linetype = 'dotted', size = 1.2,fill = NA) # 图四周框起来
  )
p
ggsave(p, file = paste0(outdir, 'C4_C9_ssgsea_bar_NA.png'),dpi = 500,
       width = 5,height  = 4)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = color ) +
  # scale_fill_manual(values = c('Up'= '#C17188','NoSignifi'='#cccccc','Down'='#748EC2')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab(ylab) + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme( 
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 13,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title  = element_text(color = 'black', size  = 16),
    legend.text   = element_text(color = 'black', size   = 16),
    # axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    panel.border = element_blank()#panel.border = element_rect(linetype = 'dotted', size = 1.2,fill = NA) # 图四周框起来
  )
p# 添加标签 # 此处参考了：https://mp.weixin.qq.com/s/eCMwWCnjTyQvNX2wNaDYXg
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)
# 依次从下到上添加标签
p1 <- p + 
  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black', size   = 4.5, face = 'arial') + # 小于-1的为黑色标签
  # geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
  #           hjust = 0,color = 'grey') + # 灰色标签
  # geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
  #           hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black', size   = 4.5, face = 'arial') # 大于1的为黑色标签
p1
outdir <- "/data/rluo4/database/Breast/Barkley/vis_seu/"
ggsave(p1, file = paste0(outdir, 'C4_C9_ssgsea_bar.png'),dpi = 500,
       width = 12.5,height  = 8)


organ = 'Cervix'
marker_Cervix <- as.data.frame(matrix(NA, nrow = 0, ncol = (length(markers) + 2)))
colnames(marker_Cervix) <- c(markers, 'Cell_Type', 'TissueType')
for (sample in pdata.ST$geo_accession) {
  min.cells <- 10
  n <- paste0(sample, '_obj')
  seu <- SPOTlight_obj_Cervix[[n]]#PCT_seu[[organ]]
  table(Idents(seu))
  genes <- intersect(markers, rownames(seu@assays$SCT$data))
  genes
  
  PCT_data <- as.data.frame(t(seu@assays$SCT$data[genes, ]))
  # PCT_data$cell_type <- seu$labels
  PCT_data$Cell_Type <- seu$labels
  # SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj.rds'))
  # PCT_data$spotCoordinates <- rownames(SpaCET_obj@input$spotCoordinates)
  # 
  subtype <- data.frame(subset =  table(PCT_data$Cell_Type))
  subtype <- subtype[subtype$subset.Freq !=0,]
  subtype <- subtype[subtype$subset.Freq >=5,]
  PCT_data$Cell_Type <- as.character(PCT_data$Cell_Type)
  table(PCT_data$Cell_Type)
  extract_last_element <- function(x) {
    split_string <- strsplit(x, "-")
    last_element <- sapply(split_string, function(y) tail(y, n = 1))
    return(last_element)
  }
  celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
  celltype_choose <- c(extract_last_element(PCT_Clusters[[organ]]),celltype_choose ,'Epithelia','Mesenchymal')
  celltype_choose
  PCT_data <- PCT_data  %>% 
    filter(Cell_Type %in% subtype$subset.Var1) %>% #& MP %in% c('Epithelial-senescence')) %>%
    filter(Cell_Type %in%  celltype_choose)
  table( PCT_data$Cell_Type)
  level_compare <- celltype_choose[celltype_choose %in% PCT_data$Cell_Type]
  PCT_data$Cell_Type <- factor(PCT_data$Cell_Type, levels = level_compare )
  PCT_data$TissueType = sample
  # p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
  # View(paletteer::palettes_d_names)
  coln <- setdiff(colnames(marker_Cervix), colnames(PCT_data))
  if(length(coln)>=1){
    PCT_data[coln] <- NA
  }
  marker_Cervix <<- rbind(marker_Cervix, PCT_data[, colnames(marker_Cervix)])
  # marker_Cervix <<-  marker_Cervix %>%
  #   inner_join(PCT_data, by = coln)
}
table(marker_Cervix$Cell_Type)
marker_Cervix$Cell_Type <- gsub('Epithelia', 'Healthy-STM', marker_Cervix$Cell_Type)
marker_Cervix$Cell_Type <- gsub('Mesenchymal', 'FIB', marker_Cervix$Cell_Type)
unique(marker_Cervix$Cell_Type)

unique(marker_Cervix$TissueType)
marker_Cervix$TissueType <- gsub('GSM6360689', 'Healthy', gsub('GSM6360690','HPV_HSIL', gsub('GSM6360691', 'CC', gsub('GSM6360692', 'CC',
                                                                                                                      marker_Cervix$TissueType ))))
unique(marker_Cervix$TissueType)
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
celltype_choose <- c('Healthy-STM',extract_last_element(PCT_Clusters[[organ]]),celltype_choose )
level_compare <- celltype_choose[celltype_choose %in% marker_Cervix$Cell_Type]
# marker_Cervix$Cell_Type <- factor(marker_Cervix$Cell_Type, levels = level_compare )
marker_Cervix$Cell_Type <- as.character(marker_Cervix$Cell_Type)
unique(marker_Cervix$Cell_Type)

marker_Cervix$Lineage <-   Lineage$Major_type[match(marker_Cervix$Cell_Type, Lineage$Minor_type)]
marker_Cervix$Lineage[is.na(marker_Cervix$Lineage)] <- 'STM'
df.ABN <- marker_Cervix %>% group_by(Cell_Type, TissueType) %>%  summarise( sd= sd(ANXA1), Med_exp = median(ANXA1), Avg_exp= mean(ANXA1), max = max(ANXA1), min = min(ANXA1))
markers_dir <- "/data/rluo4/All/Output/res_fig/ST_markers_fig/"


for (gene in genes) {
  # gene = 'FPR1'#'ANXA1'#'FPR1'
  # any_zero <- any(PCT_data[[gene]] == 0)
  # Filter out groups with all zero values
  # Convert 'gene' column to numeric if it's stored as character
  marker_Cervix[[gene]][is.na(marker_Cervix[[gene]])] = 0
  marker_Cervix[[gene]] <- as.numeric(marker_Cervix[[gene]])
  
  for (majortype in unique(marker_Cervix$Lineage)) {
    ggplotdata <- marker_Cervix %>% filter(Lineage==majortype)
    if(majortype=='STM'){
      ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = c('Healthy-STM',extract_last_element(PCT_Clusters[[organ]])) )
    }
    levels = c('Healthy','HPV_HSIL','CC')
    levels = levels[levels %in% unique(ggplotdata$TissueType)]
    ggplotdata$TissueType<- factor(ggplotdata$TissueType, levels = levels)
    
    if(majortype !='STM'){
      ggplotdata$Cell_Type <- paste(ggplotdata$TissueType, ggplotdata$Cell_Type, sep = '-')
      # factor(ggplotdata$Cell_Type, levels = extract_last_element(PCT_Clusters[[organ]]) )
    } 
    print(table(ggplotdata$Cell_Type))
    if(length(unique(ggplotdata$Cell_Type)) <=1){
      print(paste0('There is no not enough cell type to compare', ' in ', majortype,  " , so skip !"))
      next;
    }
    # Define the formula string for the t-test
    formula_str <- reformulate("Cell_Type", response = gene)
    stat.test <- ggplotdata %>%
      # group_by(TissueType) %>%
      t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    # stat.test <- PCT_data %>%
    #   group_by(Cell_Type) %>%
    #   t_test(TGFB1 ~ Cell_Type, ref.group = "Breast-C14")
    color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
    
    stat.test <- stat.test %>% add_y_position(fun = "mean_se")
    bp <- ggbarplot(
      ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "jco",#"palettetown", #View(paletteer::palettes_d_names)
      #palette = "quilava",
      add = c('mean_se'), #facet.by = "TissueType"
    ) + xlab("") + #ylab("") +
      stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
        text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_blank(),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="none",
        legend.title=element_text(size=12))
    bp
    wid = ifelse(length(unique(ggplotdata$Cell_Type)) < 4, length(unique(ggplotdata$Cell_Type))*0.58, length(unique(ggplotdata$Cell_Type))*0.5)
    ggsave(plot = bp, file = paste0(markers_dir, organ, '_', majortype, '_', gene, '_expression.png'), height =  3, width = wid, dpi = 500)
  }
}

for (gene in genes) {
  # gene = 'ANXA1'#'FPR1'
  marker_Cervix[[gene]][is.na(marker_Cervix[[gene]])] = 0
  marker_Cervix[[gene]] <- as.numeric(marker_Cervix[[gene]])
  
  # for (majortype in unique(marker_Cervix$Lineage)) {
  ggplotdata <- marker_Cervix #%>% filter(Lineage==majortype)
  unique(ggplotdata$Cell_Type)
  
  levels <- c('Healthy-STM',extract_last_element(PCT_Clusters[[organ]]), 
              unique(marker_Cervix$Cell_Type[marker_Cervix$Lineage=='Myeloid']),
              unique(marker_Cervix$Cell_Type[marker_Cervix$Lineage=='Mesenchymal']))
  levels <- levels[levels %in% unique(ggplotdata$Cell_Type)]
  ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = levels)
  # Define the formula string for the t-test
  formula_str <- reformulate("Cell_Type", response = gene)
  stat.test <- ggplotdata %>%
    # group_by(Lineage) %>%
    t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
  # stat.test <- PCT_data %>%
  #   group_by(Cell_Type) %>%
  #   t_test(TGFB1 ~ Cell_Type, ref.group = "Cervix-C14")
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  
  stat.test <- stat.test %>% add_y_position(fun = "mean_se")
  bp <- ggbarplot(
    ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "rickandmorty",#"palettetown", #View(paletteer::palettes_d_names)
    #palette = "quilava",
    add = c('mean_se'), #facet.by = "Lineage"
  ) + xlab("") + #ylab("") +
    # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
      text = element_text(size = 12),
      panel.grid.major = element_line(colour = "grey90", size=0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # axis.text.x = element_blank(),
      axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
      legend.position="none",
      legend.title=element_text(size=12))
  bp
  ggsave(plot = bp, file = paste0(markers_dir, organ, '_ST_',  gene, '_expression.png'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
}


organ = 'Liver'
marker_Liver <- as.data.frame(matrix(NA, nrow = 0, ncol = (length(markers) + 2)))
colnames(marker_Liver) <- c(markers, 'Cell_Type', 'TissueType')
indir = paste0('/data/rluo4/database/', organ, '/Barkley/')
outdir = paste0('/data/rluo4/database/', organ, '/Barkley/vis_seu/')
samples = list.files(indir)
sample = samples[grep('NYU',samples)]
print(sample)#'NYU_LIHC1'
# for (sample in pdata.ST$geo_accession) {
min.cells <- 10
n <- paste0(sample, '_obj')
seu <- SPOTlight_obj_Liver[[n]]#PCT_seu[[organ]]
table(Idents(seu))
genes <- intersect(markers, rownames(seu@assays$SCT$data))
genes

PCT_data <- as.data.frame(t(seu@assays$SCT$data[genes, ]))
# PCT_data$cell_type <- seu$labels
PCT_data$Cell_Type <- seu$labels
# SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj.rds'))
# PCT_data$spotCoordinates <- rownames(SpaCET_obj@input$spotCoordinates)
# 
subtype <- data.frame(subset =  table(PCT_data$Cell_Type))
subtype <- subtype[subtype$subset.Freq !=0,]
subtype <- subtype[subtype$subset.Freq >=5,]
PCT_data$Cell_Type <- as.character(PCT_data$Cell_Type)
table(PCT_data$Cell_Type)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid', 'Endothelia')])
celltype_choose <- c('C1','C0','C11','C15','C7',celltype_choose,'Myeloid' )
celltype_choose
PCT_data <- PCT_data  %>% 
  filter(Cell_Type %in% subtype$subset.Var1) %>% #& MP %in% c('Epithelial-senescence')) %>%
  filter(Cell_Type %in%  celltype_choose)
table( PCT_data$Cell_Type)
level_compare <- celltype_choose[celltype_choose %in% PCT_data$Cell_Type]
PCT_data$Cell_Type <- factor(PCT_data$Cell_Type, levels = level_compare )
PCT_data$TissueType = sample
# p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
# View(paletteer::palettes_d_names)
coln <- setdiff(colnames(marker_Liver), colnames(PCT_data))
if(length(coln)>=1){
  PCT_data[coln] <- NA
}
marker_Liver <<- rbind(marker_Liver, PCT_data[, colnames(marker_Liver)])
# marker_Liver <<-  marker_Liver %>%
#   inner_join(PCT_data, by = coln)
# }
table(marker_Liver$Cell_Type)

celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
celltype_choose <- c('C1','C0','C11','C15','C7',celltype_choose,'Myeloid' )
level_compare <- celltype_choose[celltype_choose %in% marker_Liver$Cell_Type]
# marker_Liver$Cell_Type <- factor(marker_Liver$Cell_Type, levels = level_compare )
marker_Liver$Cell_Type <- as.character(marker_Liver$Cell_Type)
unique(marker_Liver$Cell_Type)
marker_Liver$Lineage <-   Lineage$Major_type[match(marker_Liver$Cell_Type, Lineage$Minor_type)]
marker_Liver$Lineage[marker_Liver$Cell_Type=='Myeloid'] <- 'Myeloid'
marker_Liver$Lineage[is.na(marker_Liver$Lineage)] <- 'STM'
for (gene in genes) {
  # gene = 'ANXA1'#'FPR1'
  marker_Liver[[gene]][is.na(marker_Liver[[gene]])] = 0
  marker_Liver[[gene]] <- as.numeric(marker_Liver[[gene]])
  
  # for (majortype in unique(marker_Liver$Lineage)) {
  ggplotdata <- marker_Liver #%>% filter(Lineage==majortype)
  levels <- c(extract_last_element(PCT_Clusters[[organ]])[extract_last_element(PCT_Clusters[[organ]]) %in%  unique(marker_Liver$Cell_Type)], 
              unique(marker_Liver$Cell_Type[marker_Liver$Lineage=='Myeloid']),
              unique(marker_Liver$Cell_Type[marker_Liver$Lineage=='Endothelia']),
              unique(marker_Liver$Cell_Type[marker_Liver$Lineage=='Mesenchymal']))
  # levels <- c(levels, 'Myeloid')
  ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = levels)
  unique(ggplotdata$Cell_Type)
  # Define the formula string for the t-test
  formula_str <- reformulate("Cell_Type", response = gene)
  # stat.test <- ggplotdata %>%
  #   # group_by(Lineage) %>%
  #   t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
  #   adjust_pvalue(method = "bonferroni") %>%
  #   add_significance()
  # stat.test <- PCT_data %>%
  #   group_by(Cell_Type) %>%
  #   t_test(TGFB1 ~ Cell_Type, ref.group = "Liver-C14")
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  
  stat.test <- stat.test %>% add_y_position(fun = "mean_se")
  bp <- ggbarplot(
    ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "rickandmorty",#"palettetown", #View(paletteer::palettes_d_names)
    #palette = "quilava",
    add = c('mean_se'), #facet.by = "Lineage"
  ) + xlab("") + #ylab("") +
    # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
      text = element_text(size = 12),
      panel.grid.major = element_line(colour = "grey90", size=0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # axis.text.x = element_blank(),
      axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
      legend.position="none",
      legend.title=element_text(size=12))
  bp
  ggsave(plot = bp, file = paste0(markers_dir, organ, '_ST_',  gene, '_expression.png'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
}
for (gene in genes) {
  marker_Liver[[gene]][is.na(marker_Liver[[gene]])] = 0
  marker_Liver[[gene]] <- as.numeric(marker_Liver[[gene]])
  
  for (majortype in unique(marker_Liver$Lineage)) {
    ggplotdata <- marker_Liver %>% filter(Lineage==majortype)
    level_compare <- c('C1','C0','C11','C15','C7')
    level_compare <- level_compare[level_compare %in% unique(ggplotdata$Cell_Type)]
    if(majortype=='Myeloid'){
      level_compare <- c(level_compare, 'Myeloid')
    }
    if(majortype=='STM'){
      ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = level_compare )
    }
    print(unique(ggplotdata$Cell_Type))
    if(length(unique(ggplotdata$Cell_Type)) <=1){
      print(paste0('There is no not enough cell type to compare', ' in ', majortype,  " , so skip !"))
      next;
    }
    
    # Define the formula string for the t-test
    formula_str <- reformulate("Cell_Type", response = gene)
    stat.test <- ggplotdata %>%
      # group_by(Lineage) %>%
      t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    # stat.test <- PCT_data %>%
    #   group_by(Cell_Type) %>%
    #   t_test(TGFB1 ~ Cell_Type, ref.group = "Liver-C14")
    color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
    
    stat.test <- stat.test %>% add_y_position(fun = "mean_se")
    bp <- ggbarplot(
      ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "jco",#"palettetown", #View(paletteer::palettes_d_names)
      #palette = "quilava",
      add = c('mean_se'), #facet.by = "Lineage"
    ) + xlab("") + #ylab("") +
      stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
        text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_blank(),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="none",
        legend.title=element_text(size=12))
    bp
    ggsave(plot = bp, file = paste0(markers_dir, organ, '_', majortype, '_', gene, '_expression.png'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
  }
}


organ = 'Endometrium'
marker_Endometrium <- as.data.frame(matrix(NA, nrow = 0, ncol = (length(markers) + 2)))
colnames(marker_Endometrium) <- c(markers, 'Cell_Type', 'TissueType')
indir = paste0('/data/rluo4/database/', organ, '/Barkley/')
outdir = paste0('/data/rluo4/database/', organ, '/Barkley/vis_seu/')
samples = list.files(indir)
sample = samples[grep('NYU',samples)]
print(sample)#'NYU_LIHC1'
# for (sample in pdata.ST$geo_accession) {
min.cells <- 10
n <- paste0(sample, '_obj')
seu <- SPOTlight_obj_Endometrium[[n]]#PCT_seu[[organ]]
table(Idents(seu))
genes <- intersect(markers, rownames(seu@assays$SCT$data))
genes

PCT_data <- as.data.frame(t(seu@assays$SCT$data[genes, ]))
# PCT_data$cell_type <- seu$labels
PCT_data$Cell_Type <- seu$labels
# SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj.rds'))
# PCT_data$spotCoordinates <- rownames(SpaCET_obj@input$spotCoordinates)
# 
subtype <- data.frame(subset =  table(PCT_data$Cell_Type))
subtype <- subtype[subtype$subset.Freq !=0,]
subtype <- subtype[subtype$subset.Freq >=5,]
PCT_data$Cell_Type <- as.character(PCT_data$Cell_Type)
table(PCT_data$Cell_Type)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
celltype_choose <- c(extract_last_element(PCT_Clusters[[organ]]),celltype_choose )
celltype_choose
PCT_data <- PCT_data  %>% 
  filter(Cell_Type %in% subtype$subset.Var1) %>% #& MP %in% c('Epithelial-senescence')) %>%
  filter(Cell_Type %in%  celltype_choose)
table( PCT_data$Cell_Type)
level_compare <- celltype_choose[celltype_choose %in% PCT_data$Cell_Type]
PCT_data$Cell_Type <- factor(PCT_data$Cell_Type, levels = level_compare )
PCT_data$TissueType = sample
# p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
# View(paletteer::palettes_d_names)
coln <- setdiff(colnames(marker_Endometrium), colnames(PCT_data))
if(length(coln)>=1){
  PCT_data[coln] <- NA
}
marker_Endometrium <<- rbind(marker_Endometrium, PCT_data[, colnames(marker_Endometrium)])
# marker_Endometrium <<-  marker_Endometrium %>%
#   inner_join(PCT_data, by = coln)
# }
unique(marker_Endometrium$TissueType)

table(marker_Endometrium$Cell_Type)
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
celltype_choose <- c(extract_last_element(PCT_Clusters[[organ]]),celltype_choose )
level_compare <- celltype_choose[celltype_choose %in% marker_Endometrium$Cell_Type]
# marker_Endometrium$Cell_Type <- factor(marker_Endometrium$Cell_Type, levels = level_compare )
marker_Endometrium$Cell_Type <- as.character(marker_Endometrium$Cell_Type)
unique(marker_Endometrium$Cell_Type)
marker_Endometrium$Lineage <- Lineage$Major_type[match(marker_Endometrium$Cell_Type, Lineage$Minor_type)]
marker_Endometrium$Lineage[is.na(marker_Endometrium$Lineage)] <- 'STM'
for (gene in genes) {
  # gene = 'ANXA1'#'FPR1'
  marker_Endometrium[[gene]][is.na(marker_Endometrium[[gene]])] = 0
  marker_Endometrium[[gene]] <- as.numeric(marker_Endometrium[[gene]])
  
  # for (majortype in unique(marker_Endometrium$Lineage)) {
  ggplotdata <- marker_Endometrium #%>% filter(Lineage==majortype)
  levels <- c(extract_last_element(PCT_Clusters[[organ]]), 
              unique(marker_Endometrium$Cell_Type[marker_Endometrium$Lineage=='Myeloid']),
              unique(marker_Endometrium$Cell_Type[marker_Endometrium$Lineage=='Mesenchymal']))
  ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = levels)
  # Define the formula string for the t-test
  formula_str <- reformulate("Cell_Type", response = gene)
  # stat.test <- ggplotdata %>%
  #   # group_by(Lineage) %>%
  #   t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
  #   adjust_pvalue(method = "bonferroni") %>%
  #   add_significance()
  # # stat.test <- PCT_data %>%
  #   group_by(Cell_Type) %>%
  #   t_test(TGFB1 ~ Cell_Type, ref.group = "Endometrium-C14")
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  
  stat.test <- stat.test %>% add_y_position(fun = "mean_se")
  bp <- ggbarplot(
    ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "rickandmorty",#"palettetown", #View(paletteer::palettes_d_names)
    #palette = "quilava",
    add = c('mean_se'), #facet.by = "Lineage"
  ) + xlab("") + #ylab("") +
    # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
      text = element_text(size = 12),
      panel.grid.major = element_line(colour = "grey90", size=0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # axis.text.x = element_blank(),
      axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
      legend.position="none",
      legend.title=element_text(size=12))
  bp
  ggsave(plot = bp, file = paste0(markers_dir, organ, '_ST_',  gene, '_expression.pdf'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
}

for (gene in genes) {
  marker_Endometrium[[gene]][is.na(marker_Endometrium[[gene]])] = 0
  marker_Endometrium[[gene]] <- as.numeric(marker_Endometrium[[gene]])
  
  for (majortype in unique(marker_Endometrium$Lineage)) {
    ggplotdata <- marker_Endometrium %>% filter(Lineage==majortype)
    level_compare <- extract_last_element(PCT_Clusters[[organ]])
    level_compare <- level_compare[level_compare %in% unique(ggplotdata$Cell_Type)]
    if(majortype=='STM'){
      ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = level_compare )
    }
    print(table(ggplotdata$Cell_Type))
    if(length(unique(ggplotdata$Cell_Type)) <=1){
      print(paste0('There is no not enough cell type to compare', ' in ', majortype,  " , so skip !"))
      next;
    }
    # Define the formula string for the t-test
    formula_str <- reformulate("Cell_Type", response = gene)
    stat.test <- ggplotdata %>%
      # group_by(Lineage) %>%
      t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    # stat.test <- PCT_data %>%
    #   group_by(Cell_Type) %>%
    #   t_test(TGFB1 ~ Cell_Type, ref.group = "Endometrium-C14")
    color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
    
    stat.test <- stat.test %>% add_y_position(fun = "mean_se")
    bp <- ggbarplot(
      ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "jco",#"palettetown", #View(paletteer::palettes_d_names)
      #palette = "quilava",
      add = c('mean_se'), #facet.by = "Lineage"
    ) + xlab("") + #ylab("") +
      stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
        text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_blank(),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="none",
        legend.title=element_text(size=12))
    bp
    ggsave(plot = bp, file = paste0(markers_dir, organ, '_', majortype, '_', gene, '_expression.pdf'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
  }
}

organ = 'Pancreas'
marker_Pancreas <- as.data.frame(matrix(NA, nrow = 0, ncol = (length(markers) + 2)))
colnames(marker_Pancreas) <- c(markers, 'Cell_Type', 'TissueType')
indir = paste0('/data/rluo4/database/', organ, '/Barkley/')
outdir = paste0('/data/rluo4/database/', organ, '/Barkley/vis_seu/')
samples = list.files(indir)
sample = samples[grep('NYU',samples)]
print(sample)#'NYU_LIHC1'
# for (sample in pdata.ST$geo_accession) {
min.cells <- 10
n <- paste0(sample, '_obj')
seu <- SPOTlight_obj_Pancreas[[n]]#PCT_seu[[organ]]
table(Idents(seu))
genes <- intersect(markers, rownames(seu@assays$SCT$data))
genes

PCT_data <- as.data.frame(t(seu@assays$SCT$data[genes, ]))
# PCT_data$cell_type <- seu$labels
PCT_data$Cell_Type <- seu$labels
# SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj.rds'))
# PCT_data$spotCoordinates <- rownames(SpaCET_obj@input$spotCoordinates)
# 
subtype <- data.frame(subset =  table(PCT_data$Cell_Type))
subtype <- subtype[subtype$subset.Freq !=0,]
subtype <- subtype[subtype$subset.Freq >=5,]
PCT_data$Cell_Type <- as.character(PCT_data$Cell_Type)
table(PCT_data$Cell_Type)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
celltype_choose <- c(extract_last_element(PCT_Clusters[[organ]]),celltype_choose )
celltype_choose
PCT_data <- PCT_data  %>% 
  filter(Cell_Type %in% subtype$subset.Var1) %>% #& MP %in% c('Epithelial-senescence')) %>%
  filter(Cell_Type %in%  celltype_choose)
table( PCT_data$Cell_Type)
level_compare <- celltype_choose[celltype_choose %in% PCT_data$Cell_Type]
PCT_data$Cell_Type <- factor(PCT_data$Cell_Type, levels = level_compare )
PCT_data$TissueType = sample
# p1 = plotGeneExpression(cellchat, signaling = pathways.show)#"TGFb")
# View(paletteer::palettes_d_names)
coln <- setdiff(colnames(marker_Pancreas), colnames(PCT_data))
if(length(coln)>=1){
  PCT_data[coln] <- NA
}
marker_Pancreas <<- rbind(marker_Pancreas, PCT_data[, colnames(marker_Pancreas)])
# marker_Pancreas <<-  marker_Pancreas %>%
#   inner_join(PCT_data, by = coln)
# }
table(marker_Pancreas$Cell_Type)
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Mesenchymal', 'Myeloid')])
celltype_choose <- c(extract_last_element(PCT_Clusters[[organ]]),celltype_choose )
level_compare <- celltype_choose[celltype_choose %in% marker_Pancreas$Cell_Type]
# marker_Pancreas$Cell_Type <- factor(marker_Pancreas$Cell_Type, levels = level_compare )
marker_Pancreas$Cell_Type <- as.character(marker_Pancreas$Cell_Type)
unique(marker_Pancreas$Cell_Type)
marker_Pancreas$Lineage <-   Lineage$Major_type[match(marker_Pancreas$Cell_Type, Lineage$Minor_type)]
marker_Pancreas$Lineage[is.na(marker_Pancreas$Lineage)] <- 'STM'
for (gene in genes) {
  # gene = 'ANXA1'#'FPR1'
  marker_Pancreas[[gene]][is.na(marker_Pancreas[[gene]])] = 0
  marker_Pancreas[[gene]] <- as.numeric(marker_Pancreas[[gene]])
  
  # for (majortype in unique(marker_Pancreas$Lineage)) {
  ggplotdata <- marker_Pancreas #%>% filter(Lineage==majortype)
  levels <- c(extract_last_element(PCT_Clusters[[organ]]), 
              unique(marker_Pancreas$Cell_Type[marker_Pancreas$Lineage=='Myeloid']),
              unique(marker_Pancreas$Cell_Type[marker_Pancreas$Lineage=='Mesenchymal']))
  ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = levels)
  # Define the formula string for the t-test
  formula_str <- reformulate("Cell_Type", response = gene)
  # stat.test <- ggplotdata %>%
  #   # group_by(Lineage) %>%
  #   t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
  #   adjust_pvalue(method = "bonferroni") %>%
  #   add_significance()
  # stat.test <- PCT_data %>%
  #   group_by(Cell_Type) %>%
  #   t_test(TGFB1 ~ Cell_Type, ref.group = "Pancreas-C14")
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  
  stat.test <- stat.test %>% add_y_position(fun = "mean_se")
  bp <- ggbarplot(
    ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "rickandmorty",#"palettetown", #View(paletteer::palettes_d_names)
    #palette = "quilava",
    add = c('mean_se'), #facet.by = "Lineage"
  ) + xlab("") + #ylab("") +
    # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
      text = element_text(size = 12),
      panel.grid.major = element_line(colour = "grey90", size=0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # axis.text.x = element_blank(),
      axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
      legend.position="none",
      legend.title=element_text(size=12))
  bp
  ggsave(plot = bp, file = paste0(markers_dir, organ, '_ST_',  gene, '_expression.pdf'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
}
for (gene in genes) {
  marker_Pancreas[[gene]][is.na(marker_Pancreas[[gene]])] = 0
  marker_Pancreas[[gene]] <- as.numeric(marker_Pancreas[[gene]])
  
  for (majortype in unique(marker_Pancreas$Lineage)) {
    ggplotdata <- marker_Pancreas %>% filter(Lineage==majortype)
    level_compare <- extract_last_element(PCT_Clusters[[organ]])
    level_compare <- level_compare[level_compare %in% unique(ggplotdata$Cell_Type)]
    if(majortype=='STM'){
      ggplotdata$Cell_Type<- factor(ggplotdata$Cell_Type, levels = level_compare )
    }
    print(table(ggplotdata$Cell_Type))
    if(length(unique(ggplotdata$Cell_Type)) <=1){
      print(paste0('There is no not enough cell type to compare', ' in ', majortype,  " , so skip !"))
      next;
    }
    
    # Define the formula string for the t-test
    formula_str <- reformulate("Cell_Type", response = gene)
    stat.test <- ggplotdata %>%
      # group_by(Lineage) %>%
      t_test(formula_str) %>%  # t_test(ACVR1B ~ Cell_Type) %>% # t_test(!!sym(gene) ~ Cell_Type) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    # stat.test <- PCT_data %>%
    #   group_by(Cell_Type) %>%
    #   t_test(TGFB1 ~ Cell_Type, ref.group = "Pancreas-C14")
    color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
    
    stat.test <- stat.test %>% add_y_position(fun = "mean_se")
    bp <- ggbarplot(
      ggplotdata, x = "Cell_Type", y = gene, fill = 'Cell_Type',  palette = "jco",#"palettetown", #View(paletteer::palettes_d_names)
      #palette = "quilava",
      add = c('mean_se'), #facet.by = "Lineage"
    ) + xlab("") + #ylab("") +
      stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
        text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_blank(),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="none",
        legend.title=element_text(size=12))
    bp
    ggsave(plot = bp, file = paste0(markers_dir, organ, '_', majortype, '_', gene, '_expression.pdf'), height =  3, width = length(unique(ggplotdata$Cell_Type))*0.6, dpi = 500)
  }
}




##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 12) Check the gene expression of ANNEXIN and TGFb
library(CellChat)
library(stringr)
library(Seurat)
library(ArchR)
library(dplyr)
########################################################################
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
PCT_Clusters =  vector("list", length(organ_all))
names(PCT_Clusters) <- organ_all
########################################################################
for (organ in organ_all) {
  # Define PCT_Clusters
  ########
  print(organ)
  if(organ == 'Breast'){#
    PCT_Clusters[[organ]] <- c( "T1-C9","T1-C14",'T1-C4')
  }
  if(organ == 'Cervix'){#
    PCT_Clusters[[organ]] <- c("T1-C5", 'T3-C8', "T1-C1",  'T1-C14')# "T2-C4", "T2-C10",  
  }
  if(organ == 'Chen'){# 
    PCT_Clusters[[organ]] <- c('T1-C24',  "T1-C21", "T1-C1")
  }
  if(organ == 'CRC'){# which is Becker
    PCT_Clusters[[organ]] <- c("T1-C17", "T1-C19", "T1-C3")
  }
  if(organ == 'Endometrium'){#
    PCT_Clusters[[organ]] <- c('T1-C12','T1-C14', 'T2-C19', "T2-C4",  "T2-C24")#c( 'T2-C19', "T2-C4",  "T2-C24")
  }
  if(organ == 'Esophagus'){#
    PCT_Clusters[[organ]] <- c("T1-C12", "T1-C10")  # 'T1-C7'
  }
  if(organ == 'GC'){#
    PCT_Clusters[[organ]] <- c("T2-C12", "T2-C7", "T1-C2", "T1-C17")#c("T1-C13", "T1-C2", "T2-C3", "T2-C12")
  }
  if(organ == 'HNSCC'){#
    PCT_Clusters[[organ]] <- c("T1-C18", 'T1-C16', "T1-C11", 'T1-C4')# ,c("T1-C11", "T1-C18", "T2-C19", "T2-C3")
  }
  if(organ=='Liver'){#
    PCT_Clusters[[organ]] <-  c( "T2-C11", 'T2-C15', 'T2-C7', "T1-C1", "T1-C0")#c( "T2-C11", "T1-C1", "T1-C0")
  }  
  if(organ=='Lung'){#
    PCT_Clusters[[organ]] <- c('T1-C6',"T1-C9",  "T1-C8",  "T1-C14")#c("T1-C12", "T1-C6", "T1-C4", "T1-C3")
  }
  if(organ == 'Pancreas'){#
    PCT_Clusters[[organ]] <- c("T1-C2", "T1-C8", "T2-C9", "T2-C4") #c("T2-C6", "T2-C5",  "T2-C0", "T2-C1", "T2-C9", 'T2-C4')#
  }
  if(organ == "Prostate"){#
    PCT_Clusters[[organ]] <- c( "T2-C14", 'T2-C16',  "T2-C34")#c("T1-C6", "T1-C9", "T2-C34", "T2-C14")
  }
  if(organ == "Skin"){#
    PCT_Clusters[[organ]] <- c( 'T1-C9', "T1-C11", "T2-C1",  "T2-C14" ) #,# c("T1-C12", "T1-C9", "T2-C1", "T2-C11")
  } 
  if(organ=="THCA"){#
    PCT_Clusters[[organ]] <- c( "T2-C2",  "T2-C16", "T1-C15", "T1-C5", "T1-C1")
  }
  
}
######################################################################################################
Lineage <- read.csv('/data/rluo4/All/Output/organ13-celltypes.csv')
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
remove_last_element <- function(char_vector) {
  char_vector_parts <- strsplit(char_vector, "_")
  char_vector_parts <- lapply(char_vector_parts, function(parts) {
    if (length(parts) > 1) {
      return(paste(parts[-length(parts)], collapse = "_"))
    } else {
      return(parts[1])
    }
  })
  return(unlist(char_vector_parts))
}

markers <- c('ANXA1','FPR1','FPR2','FPR3', 
             'CTGF','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','ACVR1B', 'ACVR1',
             'FAP','MMP1','ACTA2','MYH11','COL6A3')#MP_LR$MP_Gene[MP_LR$MPs==mp]
H_dir <- '/data/rluo4/All/Output/sdata_H/'
data_path = '/data/rluo4/All/Output/res_fig/'
# load(file = paste0(data_path,'../PCT_seu_CCI.rds'))
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}

load( file = paste0(data_path,'../PCT_SCT_CCI.rds'))
marker_df <- NULL
for (organ in organ_all) {
  # organ = 'Breast'
  min.cells <- 10
  # load(file = paste0('/data/rluo4/All/Output/PCT/PCT-', organ, '_SCT_CCI.rds'))
  # save(seu, file = paste0('/data/rluo4/All/Output/PCT/PCT-', organ, '_SCT_CCI.rds'))
  seu <- PCT_seu[[organ]]
  table(Idents(seu))
  genes <- intersect(markers, rownames(seu@assays$SCT$data))
  p <- DotPlot(seu, features = rev(genes))
  TME_pdata <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
  colnames(TME_pdata) <- c('id','features.plot','Pct.exp','Avg.exp.scaled')
  TME_pdata$id <- as.character(TME_pdata$id)
  TME_pdata$Cell_Type <- extract_last_element(TME_pdata$id)
  TME_pdata$TissueType <- remove_last_element(TME_pdata$id)
  TME_pdata$Lineage <- seu$Lineage[match(TME_pdata$TissueType, seu$TissueType)]
  TME_pdata$Organ <- organ
  TME_pdata$MP <- ifelse(TME_pdata$features.plot %in%  c('ANXA1','FPR1','FPR2','FPR3'), 'ANNEXIN', 'TGFb')
  marker_df <<- rbind(marker_df, TME_pdata)
  
  PCT_data <- as.data.frame(t(seu@assays$SCT$data[genes, ]))
  PCT_data$cell_type <- seu$cell_type
  PCT_data$Cell_Type <- seu$Cell_Type
  PCT_data$cell_class <- seu$cell_class
  PCT_data$TissueType <- seu$TissueType
  subtype <- data.frame(subset =  table(PCT_data$cell_class))
  subtype <- subtype[subtype$subset.Freq !=0,]
  subtype <- subtype[subtype$subset.Freq >=10,]
  PCT_data <- PCT_data[PCT_data$cell_class %in% subtype$subset.Var1,]
  extract_last_element <- function(x) {
    split_string <- strsplit(x, "-")
    last_element <- sapply(split_string, function(y) tail(y, n = 1))
    return(last_element)
  }
  # level_compare <- paste(organ, c(extract_last_element(PCT_Clusters[[organ]])), sep = '-')
  level_compare <- paste(organ, c('Healthy',extract_last_element(PCT_Clusters[[organ]])), sep = '-')
  print(table(PCT_data$TissueType %in% level_compare))
  PCT_data <- PCT_data[PCT_data$TissueType %in% level_compare,]
  CT_TT <- data.frame(table(PCT_data$Cell_Type, PCT_data$TissueType) )
  CT_TT <- CT_TT[CT_TT$Freq!=0,]
  CT_remove <- names(table(CT_TT$Var1))[which(table(CT_TT$Var1)<2)]
  print( paste0('remove cell types for grouping: ', paste(CT_remove, collapse = ', ')))
  CT_choose <- names(table(CT_TT$Var1))[which(table(CT_TT$Var1)>=2)]
  PCT_data <- PCT_data[PCT_data$Cell_Type %in% CT_choose,]
  print(table(PCT_data$Cell_Type))
  
  PCT_data$TissueType <- factor(PCT_data$TissueType, levels = level_compare )
  library(rstatix)
  library(rlang)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
  markers_dir <- "/data/rluo4/All/Output/res_fig/markers_fig/"
  for (gene in genes) {
    # gene = 'FPR1'#'ANXA1'#'FPR1'
    # any_zero <- any(PCT_data[[gene]] == 0)
    # Filter out groups with all zero values
    # Convert 'gene' column to numeric if it's stored as character
    PCT_data[[gene]] <- as.numeric(PCT_data[[gene]])
    # Define the formula string for the t-test
    formula_str <- reformulate("TissueType", response = gene)
    stat.test <- PCT_data %>%
      group_by(Cell_Type) %>%
      t_test(formula_str) %>%  # t_test(ACVR1B ~ TissueType) %>% # t_test(!!sym(gene) ~ TissueType) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    # stat.test <- PCT_data %>%
    #   group_by(Cell_Type) %>%
    #   t_test(TGFB1 ~ TissueType, ref.group = "Breast-C14")
    stat.test <- stat.test %>% add_y_position(fun = "mean_se")
    color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
    
    bp <- ggbarplot(
      PCT_data, x = "TissueType", y = gene, fill = 'TissueType', palette = color_cluster,#"jco",#"#FC4E07",
      add = c('mean_se'), facet.by = "Cell_Type"
    ) + xlab("") + #ylab("") +
      stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
        text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        # axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="top",
        legend.title=element_text(size=12)) 
    ggsave(plot = bp, file = paste0(markers_dir, organ, '_', gene, '_expression.pdf'), height =  length(unique(PCT_data$Cell_Type))*0.35, width = length(unique(PCT_data$Cell_Type))*0.35, dpi = 500)
  }
  
}
save(marker_df, file = paste0('/data/rluo4/All/Output/ANXA1_exp_df.rds'))
load( file = paste0('/data/rluo4/All/Output/ANXA1_exp_df.rds'))


########################################################################
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
Lineage <- read.csv('/data/rluo4/All/Output/organ13-celltypes.csv')
celltype_choose <- c(Lineage$Minor_type[Lineage$Major_type %in% c('Epithelia', 'Mesenchymal', 'Myeloid')])
class(celltype_choose)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}

# load( file = paste0(data_path,'../PCT_SCT_CCI.rds'))
markers <- c('ANXA1','FPR1','FPR2','FPR3', 
             'CTGF','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','ACVR1B', 'ACVR1',
             'FAP','MMP1','ACTA2','MYH11','COL6A3')#MP_LR$MP_Gene[MP_LR$MPs==mp]
remove_last_element <- function(char_vector) {
  char_vector_parts <- strsplit(char_vector, "_")
  char_vector_parts <- lapply(char_vector_parts, function(parts) {
    if (length(parts) > 1) {
      return(paste(parts[-length(parts)], collapse = "_"))
    } else {
      return(parts[1])
    }
  })
  return(unlist(char_vector_parts))
}
library(stringr)
marker_df <- NULL
for (organ in organ_all) {
  # organ = 'Breast'
  min.cells <- 10
  # load(file = paste0('/data/rluo4/All/Output/PCT/PCT-', organ, '_SCT_CCI.rds'))
  # save(seu, file = paste0('/data/rluo4/All/Output/PCT/PCT-', organ, '_SCT_CCI.rds'))
  seu <- PCT_seu[[organ]]
  table(Idents(seu))
  genes <- intersect(markers, rownames(seu@assays$SCT$data))
  seu$Major_Type <- Lineage$Major_type[match(seu$Cell_Type, Lineage$Minor_type)]
  seu$Major_Type <- paste(seu$Major_Type, str_split(Idents(seu), '_', simplify = T)[,1], sep = ': ')#paste(seu$Major_Type, remove_last_element(Idents(seu)), sep = ': ')
  Idents(seu) <- seu$Major_Type
  print(  table(Idents(seu)) )
  p <- DotPlot(seu, features = rev(genes))
  TME_pdata <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
  colnames(TME_pdata) <- c('id','features.plot','Pct.exp','Avg.exp.scaled')
  TME_pdata$id <- as.character(TME_pdata$id)
  TME_pdata$Cell_Type <- str_split(TME_pdata$id, ': ', simplify = T)[, 1]
  TME_pdata$TissueType <- str_split(TME_pdata$id, ': ', simplify = T)[, 2]
  # TME_pdata$Lineage <- seu$Lineage[match(TME_pdata$TissueType, seu$TissueType)]
  TME_pdata$Organ <- organ
  TME_pdata$MP <- ifelse(TME_pdata$features.plot %in%  c('ANXA1','FPR1','FPR2','FPR3'), 'ANNEXIN', 'TGFb')
  marker_df <<- rbind(marker_df, TME_pdata)
}
# save(marker_df, file = paste0('/data/rluo4/All/Output/ANXA1_exp_df_Lineage.rds'))
color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
load(file = paste0('/data/rluo4/All/Output/ANXA1_exp_df_Lineage.rds'))
for (organ in organ_all) {
  # organ = 'Breast'
  TME_pdata <- marker_df  %>% 
    filter(Organ %in% organ) %>% #& MP %in% c('Epithelial-senescence')) %>%
    # filter(Cell_Type %in%  celltype_choose)
    filter(Cell_Type %in%  c('Epithelia','Endothelia', 'Mesenchymal', 'Myeloid'))
  
  # TME_pdata$id <- extract_last_element(TME_pdata$id)
  # TME_pdata$Major_Type <- Lineage$Major_type[match(TME_pdata$Cell_Type, Lineage$Minor_type)]
  # TME_pdata$Major_Type <- paste(TME_pdata$Major_Type, remove_last_element(TME_pdata$Lineage), sep = ': ')
  table(TME_pdata$TissueType)
  extract_last_element <- function(x) {
    split_string <- strsplit(x, "-")
    last_element <- sapply(split_string, function(y) tail(y, n = 1))
    return(last_element)
  }
  TME_pdata$TissueType <- extract_last_element(TME_pdata$TissueType)
  table(TME_pdata$TissueType)
  # TME_pdata <- TME_pdata[order(TME_pdata$Cell_Type), ]
  level_compare <- c('Healthy',extract_last_element(PCT_Clusters[[organ]]) )
  TME_pdata <- TME_pdata[TME_pdata$TissueType %in% level_compare,]
  TME_pdata$TissueType <- factor(TME_pdata$TissueType, levels = rev(level_compare ))
  # TME_pdata$TissueType <- factor(TME_pdata$TissueType, levels = level_compare )
  
  library(rstatix)
  library(rlang)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
  markers_dir <- "/data/rluo4/All/Output/res_fig/markers_fig/"
  colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
  plotx <- ggplot(TME_pdata, aes(x = TissueType, y = features.plot)) +        ## global aes
    geom_point(aes(fill = Avg.exp.scaled, size = Pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.005) +
    facet_grid(.~Cell_Type) +
    # scale_y_discrete(breaks=0:11, labels=paste0("Treg_c", 0:11), limits = rev) +
    scale_x_discrete(limits = rev) +
    xlab("") + ylab("") +
    scale_fill_gradientn(
      colors = c( "#0077c1", "#6DCCFF", "lightyellow", colorsForDataType[3],"#ff1620")) +
    
    # colors = c("#84c7b3","#99a9cc", "lightyellow", "#EDCAE0" , "#F484AE"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    theme( text = element_text(size = 12),
           panel.grid.major = element_line(colour = "grey90", size=0.2),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black"),
           # axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
           axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1, color = 'black', size = 10),
           axis.text.x   = element_text(color = 'black', size = 12.5, angle = 40,
                                        hjust = 1,vjust = 1),
           axis.title.x  = element_text(color = 'black', size = 14, angle = 0),
           axis.title.y  = element_text(color = 'black', size = 12.5, angle = 90),
           axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
           axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
           legend.text   = element_text(color = 'black', size   = 10),
           legend.position="top", legend.key.height = unit(0.25, "cm"),
           legend.title=element_text(size=12)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               ncol = 3,
                               byrow = T,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
  ggsave(plot = plotx, file = paste0(markers_dir, organ, '_ANNEXIN_TGFb_expression.pdf'), height =  4, width = length(unique(TME_pdata$id))*0.4, dpi = 500)
  
  
  TME_pdata <- TME_pdata[TME_pdata$features.plot %in% c('FPR1','FPR2','FPR3'),]
  TME_pdata$Cell_Type <- factor(TME_pdata$Cell_Type, levels = c('Myeloid','Mesenchymal','Epithelia','Endothelia'))
  plotx <- ggplot(TME_pdata, aes(x = TissueType, y = features.plot)) +        ## global aes
    geom_point(aes(fill = Avg.exp.scaled, size = Pct.exp),
               color = 'black',
               shape = 21,
               stroke = 0.005) +
    facet_grid(.~Cell_Type) +
    # scale_y_discrete(breaks=0:11, labels=paste0("Treg_c", 0:11), limits = rev) +
    scale_x_discrete(limits = rev) +
    xlab("") + ylab("") +
    scale_fill_gradientn(
      colors = c( "#0077c1", "#6DCCFF", "lightyellow", colorsForDataType[3],"#ff1620")) +
    
    # colors = c("#84c7b3","#99a9cc", "lightyellow", "#EDCAE0" , "#F484AE"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    theme( text = element_text(size = 12),
           panel.grid.major = element_line(colour = "grey90", size=0.2),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black"),
           # axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
           axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1, color = 'black', size = 10),
           axis.text.x   = element_text(color = 'black', size = 12.5, angle = 40,
                                        hjust = 1,vjust = 1),
           axis.title.x  = element_text(color = 'black', size = 14, angle = 0),
           axis.title.y  = element_text(color = 'black', size = 12.5, angle = 90),
           axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
           axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
           legend.text   = element_text(color = 'black', size   = 10),
           legend.position="top", legend.key.height = unit(0.25, "cm"),
           legend.title=element_text(size=12)) +
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               ncol = 3,
                               byrow = T,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
  wid = ifelse(length(unique(TME_pdata$id))<=3,  length(unique(TME_pdata$id))*0.38, length(unique(TME_pdata$id))*0.32)
  ggsave(plot = plotx, file = paste0(markers_dir, organ, '_ANNEXIN_receptor_expression.png'), height =  3, width = wid, dpi = 500)
  
}

load(paste0('/data/rluo4/All/Output/ANNEXIN_exp_df.rds'))#save(ANNEXIN_exp, file = paste0('/data/rluo4/All/Output/ANNEXIN_exp_df.rds'))
for (organ in organ_all) {
  print(organ)
  PCT_data <- ANNEXIN_exp[[organ]]
  subtype <- data.frame(subset =  table(PCT_data$cell_class))
  subtype <- subtype[subtype$subset.Freq !=0,]
  subtype <- subtype[subtype$subset.Freq >=10,]
  PCT_data <- PCT_data[PCT_data$cell_class %in% subtype$subset.Var1,]
  print(table(PCT_data$Cell_Type))
  extract_last_element <- function(x) {
    split_string <- strsplit(x, "-")
    last_element <- sapply(split_string, function(y) tail(y, n = 1))
    return(last_element)
  }
  # level_compare <- paste(organ, c(extract_last_element(PCT_Clusters[[organ]])), sep = '-')
  PCT_data$TissueType <- gsub('GC', 'Stomach', PCT_data$TissueType)
  PCT_data$TissueType <- gsub('HNSCC', 'Oral cavity', PCT_data$TissueType)
  level_compare <- paste(organ, c('Healthy',extract_last_element(PCT_Clusters[[organ]])), sep = '-')
  level_compare <-  gsub('GC', 'Stomach', level_compare)
  level_compare <-  gsub('HNSCC', 'Oral cavity', level_compare)
  print(table(PCT_data$TissueType %in% level_compare))
  PCT_data <- PCT_data[PCT_data$TissueType %in% level_compare,]
  CT_TT <- data.frame(table(PCT_data$Cell_Type, PCT_data$TissueType) )
  CT_TT <- CT_TT[CT_TT$Freq!=0,]
  CT_remove <- names(table(CT_TT$Var1))[which(table(CT_TT$Var1)<2)]
  print( paste0('remove cell types for grouping: ', paste(CT_remove, collapse = ', ')))
  CT_choose <- names(table(CT_TT$Var1))[which(table(CT_TT$Var1)>=2)]
  PCT_data <- PCT_data[PCT_data$Cell_Type %in% CT_choose,]
  print(table(PCT_data$Cell_Type))
  
  PCT_data$TissueType <- factor(PCT_data$TissueType, levels = level_compare )
  
  markers_dir <- "/data/rluo4/All/Output/res_fig/markers_fig/"
  genes <- intersect(c('FPR1','FPR2','FPR3'), colnames(PCT_data))
  for (gene in genes) {
    # gene = 'FPR1'#'ANXA1'#'FPR1'
    # any_zero <- any(PCT_data[[gene]] == 0)
    # Filter out groups with all zero values
    # Convert 'gene' column to numeric if it's stored as character
    PCT_data[[gene]] <- as.numeric(PCT_data[[gene]])
    # Define the formula string for the t-test
    formula_str <- reformulate("TissueType", response = gene)
    ggplotdata <- PCT_data[PCT_data$Cell_Type %in% unique(Lineage$Minor_type[Lineage$Major_type=='Myeloid']),]
    stat.test <- ggplotdata %>%
      group_by(Cell_Type) %>%
      t_test(formula_str) %>%  # t_test(ACVR1B ~ TissueType) %>% # t_test(!!sym(gene) ~ TissueType) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    # stat.test <- PCT_data %>%
    #   group_by(Cell_Type) %>%
    #   t_test(TGFB1 ~ TissueType, ref.group = "Breast-C14")
    stat.test <- stat.test %>% add_y_position(fun = "mean_se")
    color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
    
    bp <- ggbarplot(
      ggplotdata, x = "TissueType", y = gene, fill = 'TissueType', palette = color_cluster,#"jco",#"#FC4E07",
      add = c('mean_se'), facet.by = "Cell_Type"
    ) + xlab("") + #ylab("") +
      stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
        text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        # axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="top",
        legend.title=element_text(size=12)) 
    bp
    wid = ifelse(length(unique(ggplotdata$Cell_Type))<=3, length(unique(ggplotdata$Cell_Type))*3.28, length(unique(ggplotdata$Cell_Type))*0.66)
    ggsave(plot = bp, file = paste0(markers_dir, organ, '_myeloid_', gene, '_expression.png'), height = length(unique(ggplotdata$Cell_Type))*0.8, width = wid, dpi = 500)
  }
}


#Load ANNEXIN_exp results
load(paste0('/data/rluo4/All/Output/ANNEXIN_exp_df.rds'))#save(ANNEXIN_exp, file = paste0('/data/rluo4/All/Output/ANNEXIN_exp_df.rds'))

#Import GSVA results
cor.TGFb <- NULL
df.ANXA1.TGFb <- NULL
cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
for (organ in organ_all) {
  # organ = 'Breast'
  load(paste0(cohort_directory, organ,'_PCT_ssgsea.rds'))
  print(paste('GSVA score for', organ, " STM cell number is: ", dim(ssgsea.res))[2])
  dim(ANNEXIN_exp[[organ]])
  
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds')
  load(file)
  
  set.seed(123)
  group_ABN <- rdata_filter@meta.data
  print(table(group_ABN$cell_class))
  table(group_ABN$Sample_Type)
  group_ABN$cell_subtype <- as.character(group_ABN$cell_class)
  print(table(group_ABN$cell_subtype))
  PCT_data <- ANNEXIN_exp[[organ]]
  PCT_data$cluster <- extract_last_element(PCT_data$TissueType)
  table(PCT_data$cluster)
  table(str_split(rownames(PCT_data), ':', simplify = T)[,2] %in% colnames(ssgsea.res))
  PCT_data <- PCT_data[PCT_data$cluster %in%  unique(extract_last_element(group_ABN$cell_subtype)),]
  PCT_data <- PCT_data[PCT_data$Cell_Type %in% 'STM',]
  dim(PCT_data)
  # setdiff(colnames(ssgsea.res), str_split(rownames(PCT_data), ':', simplify = T)[,2] )
  include_path <- c(as.character(unique(my_genesets$term[grep('TGF_BETA', my_genesets$term)])))
  include_path <- include_path[include_path %ni% c('DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING')]
  enrich.res <- ssgsea.res[rownames(ssgsea.res) %in% include_path, ]
  TGFb.res <- apply(enrich.res,2,mean)
  index <- match(str_split(rownames(PCT_data), ':', simplify = T)[,2], names(TGFb.res))
  print(table((str_split(rownames(PCT_data), ':', simplify = T)[,2] %in% names(TGFb.res))))
  PCT_data$TGFb.response <- TGFb.res[index] 
  PCT_data$organ <- organ
  # df.ANXA1.TGFb <- PCT_data[, c('ANXA1','TGFb.response','organ')]
  dd  <- cor.test(as.numeric(PCT_data$ANXA1), TGFb.res[index] ,type='spearman')
  dd <- data.frame(gene='ANXA1', pathway='TGFb', organ = organ, cor=dd$estimate,p.value=dd$p.value )
  cor.TGFb <<- rbind(cor.TGFb, dd)
  df.ANXA1.TGFb <<- rbind(df.ANXA1.TGFb, PCT_data[, c('ANXA1','TGFb.response','organ')])
}

# hei = ifelse(length(unique(test$Lineage))<=4, length(unique(test$Lineage))*20,length(unique(test$Lineage))*15)
# ggsave(file.path(outfile), plotx, width = 120, height = hei, units = "mm")
xlab = 'log-normalized expression of ANXA1'
ylab = paste0('GSVA score of TGF-β pathway')
outdir <- '/data/rluo4/All/Output/res_fig/ANXA1_TGFb/'
for (organ in organ_all) {
  # organ = 'Cervix'
  PCT_data <- df.ANXA1.TGFb[df.ANXA1.TGFb$organ == organ,]
  # summary(df.ANXA1.TGFb[df.ANXA1.TGFb$organ == organ,'TGFb.response'])
  label.y <- max(df.ANXA1.TGFb[df.ANXA1.TGFb$organ == organ,'TGFb.response'])
  ggscatter(PCT_data, x = "ANXA1", y = 'TGFb.response',
            color = "#CAB2D6", size = 2,
            add = "reg.line", add.params = list(color = "salmon", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
            conf.int = TRUE, # 添加回归线的置信区间
            cor.coef = TRUE, # 添加相关系数
            cor.coeff.args = list(method = "pearson", label.x = 0,label.y=label.y, label.sep =", ") #"\n")#选择Pearson相关
  )+ xlab('')+ylab('') + theme(legend.title  = element_text(color = 'black', size  = 10),
                               legend.text   = element_text(color = 'black', size   = 9.5),
                               panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) 
  )
  
  fn <- paste0(outdir, organ, '_ANXA1_cor_TGFb.png')
  ggsave(file=fn, width = 3,height = 3, dpi=500)
}

ggscatterhist(
  df.ANXA1.TGFb[df.ANXA1.TGFb$organ == organ,], x = "ANXA1", y = 'TGFb.response',
  color = "#A6CEE3", size = 1, #color = "site", palette = "nature",#杂志nature的配色
  add = "reg.line", add.params = list(color = "salmon", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
  conf.int = TRUE, # 添加回归线的置信区间
  cor.coef = TRUE, # 添加相关系数
  cor.coeff.args = list(method = "pearson", label.x = 0.45,label.y=0.4, label.sep =", "),
  #  palette = c("#00AFBB", "#E7B800", "#FC4E07"),#分组填色，如果不用杂志的颜色可以自行设置
  margin.params = list(fill ="#CAB2D6", color = "black", size = 0.2),
  # xlab(xlab) + ylab(ylab)#title="TCGA-SKCM cohort",
)# + stat_cor(aes(color = site), label.y = 0.5)
# ggsave(file="cor_ANXA1.png",width = 4,height = 3)
# ggscatter(df.ANXA1.TGFb[df.ANXA1.TGFb$organ == organ,], x = "ANXA1", y = 'TGFb.response',
#           color = "#CAB2D6", 
#           # palette = c("#00BA38", "#BBBBBB", "#F8766D"),
#           size = 2,
#           # label = ,
#           # font.label = 10, 
#           repel = T,
#           xlab = xlab, 
#           ylab = ylab) + 
#   theme_base() + xlab() + ylab(ylab)+ 
#   theme(legend.position = c(0.133,0.85),
#         legend.title  = element_text(color = 'black', size  = 10),
#         legend.text   = element_text(color = 'black', size   = 9.5))+#legend.position = "none",
#   geom_hline(yintercept = 1, linetype="dashed") +
#   geom_vline(xintercept = c(-0.15,0.15), linetype="dashed")


library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(ggpubr)
load('/data/rluo4/All/Output/pheno_all.rds')
load('/data/rluo4/All/Output/Imm_all.rds')
load('/data/rluo4/All/Output/Str_all.rds')

my_genesets <- read.csv('/data/rluo4/All/my_genesets.csv')
my_genesets <- my_genesets[-grep("KEGG", my_genesets$gs_name),-1]
# my_genesets <- rbind(Senescense, Metaplasia, H_genesets, cellcycle.gmt, reactome.gmt)
colnames(my_genesets) <- c("term","gene")
my_genesets$term <- as.factor(my_genesets$term)
unique(my_genesets$term)
my.gmt <- split(my_genesets$gene, my_genesets$term)
unique(my_genesets$term[grep('MHC', my_genesets$term)])
unique(my_genesets$term[grep('ANTIGEN', my_genesets$term)])

data_path = '/data/rluo4/All/Output'
setwd(data_path)
organ_all <- list.files('Cluster/')
# organ_all <- organ_all[organ_all!='Chen']
organ_all
cell_all <- c(pheno_all, Imm_all, Str_all)

pheno <- paste(organ_all, 'Epi', sep = '_')
pheno <- c(pheno, 'Chen_Epi')
cohort_directory = '/data/rluo4/All/Output/Epi_Results/'
cor.inflam <- NULL
df.ANXA1.inflam <- NULL
for (organ in organ_all) {
  # organ = 'Breast'
  epi <- paste0(organ, '_Epi')
  print(organ)
  epi <- gsub('Chen','CRC',paste(organ, 'Epi', sep = "_"))
  Epi <- cell_all[[epi]]
  Epi$barcode = rownames(Epi)
  Epi$cell_type <- Epi$Cell_Type
  table(Epi$Organ, Epi$Tissue)
  Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
  Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                              gsub('Tumor','PRAD', Epi$Tissue)))
  Epi$Major_type <- 'Epithelial'
  print(table(Epi$Organ))
  print(table(Epi$Tissue))
  
  load(paste0(cohort_directory, organ,'_ADJ_ssgsea.rds'))
  dim(ssgsea.res)
  ssgsea.res_ADJ <- ssgsea.res
  
  load(paste0(cohort_directory, organ,'_PCT_ssgsea.rds'))
  dim(ssgsea.res)
  
  terms <- intersect(rownames(ssgsea.res), rownames(ssgsea.res_ADJ))
  ssgsea.res <- cbind(ssgsea.res_ADJ[terms, ], ssgsea.res[terms,])
  
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'Malignant-Transformation-',organ,'.rds')
  load(file)
  
  set.seed(123)
  group_ABN <- rdata_filter@meta.data
  print(table(group_ABN$cell_class))
  table(group_ABN$Sample_Type)
  group_ABN$cell_subtype <- as.character(group_ABN$cell_class)
  print(table(group_ABN$cell_subtype))
  PCT_data <- ANNEXIN_exp[[organ]]
  PCT_data$cluster <- extract_last_element(PCT_data$TissueType)
  table(PCT_data$cluster)
  table(str_split(rownames(PCT_data), ':', simplify = T)[,2] %in% colnames(ssgsea.res))
  PCT_data <- PCT_data[PCT_data$cluster %in%  c('Healthy',unique(extract_last_element(group_ABN$cell_subtype))),]
  PCT_data <- PCT_data[PCT_data$Cell_Type %in% 'STM',]
  dim(PCT_data)
  # setdiff(colnames(ssgsea.res), str_split(rownames(PCT_data), ':', simplify = T)[,2] )
  # include_path <- c(as.character(unique(my_genesets$term[grep('TGF_BETA', my_genesets$term)])))
  include_path <- c(
    as.character(unique(my_genesets$term[grep('JAK', my_genesets$term)])),
    as.character(unique(my_genesets$term[grep('TNF', my_genesets$term)])),
    as.character(unique(my_genesets$term[grep('NFK', my_genesets$term)])),
    as.character(unique(my_genesets$term[grep('IFN', my_genesets$term)])),
    as.character(unique(my_genesets$term[grep('TLR', my_genesets$term)])),
    as.character(unique(my_genesets$term[grep('MAPK', my_genesets$term)]))
  )
  # include_path <- include_path[include_path %ni% c('DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING')]
  enrich.res <- ssgsea.res[rownames(ssgsea.res) %in% include_path, ]
  inflam.res <- apply(enrich.res,2,median)
  index <- match(str_split(rownames(PCT_data), ':', simplify = T)[,2], names(inflam.res))
  print(table((str_split(rownames(PCT_data), ':', simplify = T)[,2] %in% names(inflam.res))))
  PCT_data$inflam.response <- inflam.res[index] 
  PCT_data$organ <- organ
  # df.ANXA1.inflam <- PCT_data[, c('ANXA1','inflam.response','organ')]
  dd  <- cor.test(as.numeric(PCT_data$ANXA1), inflam.res[index] ,type='spearman')
  dd <- data.frame(gene='ANXA1', pathway='inflam', organ = organ, cor=dd$estimate,p.value=dd$p.value )
  cor.inflam <<- rbind(cor.inflam, dd)
  df.ANXA1.inflam <<- rbind(df.ANXA1.inflam, PCT_data[, c('ANXA1','inflam.response','organ','cluster')])
}
save(df.ANXA1.inflam, file = paste('/data/rluo4/All/Output/df.ANXA1.inflam.rds'))

load( paste('/data/rluo4/All/Output/df.ANXA1.inflam.rds'))
xlab = 'log-normalized expression of ANXA1'
ylab = paste0('GSVA score of inflammatory response')
outdir <- '/data/rluo4/All/Output/res_fig/ANXA1_TGFb/'
for (organ in organ_all) {
  # organ = 'Cervix'
  PCT_data <- df.ANXA1.inflam[df.ANXA1.inflam$organ == organ,]
  hub_cluster = Precancer_Clusters[[organ]]
  precancer <- extract_last_element(unlist(hub_cluster))
  PCT_data$Tissue <- PCT_data$cluster
  PCT_data$Tissue[PCT_data$cluster %in% precancer] <- 'Precancer'
  PCT_data$Tissue[! PCT_data$cluster %in% c('Healthy',precancer)] <- 'Cancer'
  print(table(PCT_data$Tissue))
  PCT_data <- PCT_data[PCT_data$Tissue %in% c('Healthy','Precancer'), ]
  level_compare <- c(extract_last_element(unlist(hub_cluster)), 'Healthy')#extract_last_element(PCT_Clusters[[organ]]) )
  level_compare <- level_compare[level_compare %in% PCT_data$cluster]
  PCT_data$cluster <- factor(PCT_data$cluster, levels = rev(level_compare ))
  
  dd  <- cor.test(as.numeric(PCT_data$ANXA1), as.numeric(PCT_data$inflam.response) ,type='spearman')
  dd <- data.frame(gene='ANXA1', pathway='inflam.response', organ = organ, cor=dd$estimate,p.value=dd$p.value )
  print(dd)
  # summary(df.ANXA1.TGFb[df.ANXA1.TGFb$organ == organ,'TGFb.response'])
  # label.y <- max(PCT_data[,'TGFb.response'])
  # ggscatter(PCT_data, x = "ANXA1", y = 'inflam.response',
  #           color = "#CAB2D6", size = 1,
  #           add = "reg.line", add.params = list(color = "salmon", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
  #           conf.int = TRUE, # 添加回归线的置信区间
  #           cor.coef = TRUE, # 添加相关系数
  #           cor.coeff.args = list(method = "spearman", label.x = 0,label.y=label.y, label.sep =", ") #"\n")#选择Pearson相关
  # )+ xlab(xlab)+ylab(ylab)
  
  library(Rmisc)
  ANXA1.HP <- summarySE(PCT_data, measurevar = 'ANXA1', groupvars = 'cluster')
  ANXA1.HP$var <- 'ANXA1'
  Inflam.HP <- summarySE(PCT_data, measurevar = 'inflam.response', groupvars = 'cluster')
  Inflam.HP$var <- 'Inflammation'
  colnames(Inflam.HP)[3] <- 'ANXA1'
  ggplot()+
    geom_line(data=ANXA1.HP, aes(x=cluster,y=ANXA1, color = var,  group=var), color = "#dd9bc5", linewidth = 1)+
    geom_point(data=ANXA1.HP, aes(x=cluster,y=ANXA1), fill = "#dc5035", shape = 21, size = 3)+
    geom_errorbar(data = ANXA1.HP, aes(x = cluster, ymin=ANXA1-se, ymax=ANXA1+se),  width = .1, color = "#dd9bc5")  +
    geom_line(data=Inflam.HP, aes(x = cluster, y = ANXA1*5, color = var, group=var), color = "#99a9cc", linewidth = 1)+
    geom_point(data=Inflam.HP, aes(x = cluster, y = ANXA1*5), fill =  "#76daff", shape = 21, size = 3)+
    geom_errorbar(data = Inflam.HP, aes(x = cluster, ymin=(ANXA1-se)*5, ymax=(ANXA1+se)*5),  width = .1, color = "#99a9cc") + 
    
    # scale_y_continuous(name = "log-normalized expression of ANXA1",#y1特征
    #                    sec.axis = sec_axis( trans=~./5, name = "GSVA score of inflammatory pathway"))+ #+#y2特征+
    scale_y_continuous(name = "",#y1特征
                       sec.axis = sec_axis( trans=~./5, name = ""))+ #+#y2特征+
    
    theme(
      text = element_text(size = 14),
      panel.grid.major = element_line(colour = "grey90", size=0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      # axis.line = element_line(colour = "black"),
      # axis.text.x = element_blank(),
      axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1, size=14, colour = "black"),
      axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1, size=12),
      # legend.position="top",
      legend.justification=c(1,0),# 这一项很关键,如果没有这个参数,图例会偏移,读者可以试一试
      legend.position=c(1,0),
      legend.title=element_text(size=12))
  fn <- paste0(outdir, organ, '_ANXA1_cor_Inflam.png')
  ggsave(file=fn, width = length(unique(PCT_data$cluster))*1, height = 3.5, dpi=500)
}
# ggscatterhist(
#   PCT_data, x = "ANXA1", y = 'inflam.response',
#   color = "Tissue", palette = "nature", size = 1, #杂志nature的配色
#   add = "reg.line", add.params = list(color = "salmon", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
#   conf.int = TRUE, # 添加回归线的置信区间
#   cor.coef = TRUE, # 添加相关系数
#   cor.coeff.args = list(method = "pearson", label.x = 0.45,label.y=0.4, label.sep =", "),
#   #  palette = c("#00AFBB", "#E7B800", "#FC4E07"),#分组填色，如果不用杂志的颜色可以自行设置
#   margin.params = list(fill ="#CAB2D6", color = "black", size = 0.2),
#   # xlab(xlab) + ylab(ylab)#title="TCGA-SKCM cohort",
# )# + stat_cor(aes(color = site), label.y = 0.5)
# # ggsave(file="cor_ANXA1.png",width = 4,height = 3)
# 
# 


library(Rmisc)
# tgc <- summarySE(tg, measurevar="len", groupvars=c("supp","dose"))
# tgc
ANXA1.HP <- summarySE(PCT_data, measurevar = 'ANXA1', groupvars = 'cluster')
ANXA1.HP$var <- 'ANXA1'
Inflam.HP <- summarySE(PCT_data, measurevar = 'inflam.response', groupvars = 'cluster')
Inflam.HP$var <- 'Inflammation'
colnames(Inflam.HP)[3] <- 'ANXA1'
# df.HP <- rbind(ANXA1.HP, Inflam.HP)
pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(ANXA1.HP, aes(x=cluster, y=ANXA1, colour=var, group=var)) +
  geom_errorbar(aes(ymin=ANXA1-se, ymax=ANXA1+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white")
ggplot()+
  geom_line(data=Inflam.HP, aes(x=cluster,y=ANXA1, color = var,  group=var), color = "#dd9bc5", linewidth = 1)+
  geom_point(data=Inflam.HP, aes(x=cluster,y=ANXA1), fill = "#dc5035", shape = 21, size = 3)+
  geom_errorbar(data = Inflam.HP, aes(x = cluster, ymin=ANXA1-se, ymax=ANXA1+se),  width = .1, color = "#dd9bc5")  +
  geom_line(data=ANXA1.HP, aes(x = cluster, y = ANXA1/5, color = var, group=var), color = "#99a9cc", linewidth = 1)+
  geom_point(data=ANXA1.HP, aes(x = cluster, y = ANXA1/5), fill =  "#76daff", shape = 21, size = 3)+
  geom_errorbar(data = ANXA1.HP, aes(x = cluster, ymin=(ANXA1-se)/5, ymax=(ANXA1+se)/5),  width = .1, color = "#99a9cc")

ggplot()+
  geom_line(data=Inflam.HP, aes(x=cluster,y=ANXA1, color = var,  group=var), color = "#dd9bc5", linewidth = 1)+
  geom_point(data=Inflam.HP, aes(x=cluster,y=ANXA1), fill = "#dc5035", shape = 21, size = 3)+
  geom_errorbar(data = Inflam.HP, aes(x = cluster, ymin=ANXA1-se, ymax=ANXA1+se),  width = .1, color = "#dd9bc5")  +
  geom_line(data=ANXA1.HP, aes(x = cluster, y = ANXA1/5, color = var, group=var), color = "#99a9cc", linewidth = 1)+
  geom_point(data=ANXA1.HP, aes(x = cluster, y = ANXA1/5), fill =  "#76daff", shape = 21, size = 3)+
  geom_errorbar(data = ANXA1.HP, aes(x = cluster, ymin=(ANXA1-se)/5, ymax=(ANXA1+se)/5),  width = .1, color = "#99a9cc") + 
  
  scale_y_continuous(name = "GSVA score of inflammatory pathway",#y1特征
                     sec.axis = sec_axis( trans=~.*5, name="log-normalized expression of ANXA1"))+ #+#y2特征+
  theme_classic() 

ggplot()+
  geom_line(data=ANXA1.HP, aes(x=cluster,y=ANXA1, color = var,  group=var), color = "#dd9bc5", linewidth = 1)+
  geom_point(data=ANXA1.HP, aes(x=cluster,y=ANXA1), fill = "#dc5035", shape = 21, size = 3)+
  geom_errorbar(data = ANXA1.HP, aes(x = cluster, ymin=ANXA1-se, ymax=ANXA1+se),  width = .1, color = "#dd9bc5")  +
  geom_line(data=Inflam.HP, aes(x = cluster, y = ANXA1*5, color = var, group=var), color = "#99a9cc", linewidth = 1)+
  geom_point(data=Inflam.HP, aes(x = cluster, y = ANXA1*5), fill =  "#76daff", shape = 21, size = 3)+
  geom_errorbar(data = Inflam.HP, aes(x = cluster, ymin=(ANXA1-se)*5, ymax=(ANXA1+se)*5),  width = .1, color = "#99a9cc") + 
  
  scale_y_continuous(name = "log-normalized expression of ANXA1",#y1特征
                     sec.axis = sec_axis( trans=~./5, name = "GSVA score of inflammatory pathway"))+ #+#y2特征+
  
  theme(
    text = element_text(size = 14),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    # axis.text.x = element_blank(),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1, size=15, colour = "black"),
    axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1, size=12),
    # legend.position="top",
    legend.justification=c(1,0),# 这一项很关键,如果没有这个参数,图例会偏移,读者可以试一试
    legend.position=c(1,0),
    legend.title=element_text(size=12))
#   theme_classic() +
#   # guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
#   theme(legend.position = "bottom",
#         # legend.title = element_blank(),
#         legend.text = element_text(size = 9),
#         plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
#         axis.text = element_text(size = 10),
#         # axis.title.y = element_blank()
#   ) -> f3splot
# print(f3splot)
# scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
#                  breaks=c("OJ", "VC"),
#                  labels=c("Orange juice", "Ascorbic acid"),
#                  l=40) +                    # Use darker colors, lightness=40
# ggtitle("The Effect of Vitamin C on\nTooth Growth in Guinea Pigs") +
# expand_limits(y=0) +                        # Expand y range
# scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
theme_bw() +
  theme(legend.justification=c(1,0),# 这一项很关键,如果没有这个参数,图例会偏移,读者可以试一试
        legend.position=c(1,0))               # Position legend in bottom right
# scale_x_continuous(name = "This is x1-axis",#x1特征
#                    sec.axis = sec_axis( trans=~.*50, name="This is x2-axis"))#x2特征

##############################################################################################################################
library(ggplot2)
data_path = '/data/rluo4/All/Output'
setwd(data_path)
organ_all <- list.files('Cluster/')
organ_all <- organ_all[organ_all!='Chen']
organ_all
cell_all <- c(pheno_all, Imm_all, Str_all)
pheno <- paste(organ_all, 'Epi', sep = '_')
pheno <- c(pheno, 'Chen_Epi')
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
remove_last_element <- function(char_vector) {
  char_vector_parts <- strsplit(char_vector, "_")
  char_vector_parts <- lapply(char_vector_parts, function(parts) {
    if (length(parts) > 1) {
      return(paste(parts[-length(parts)], collapse = "_"))
    } else {
      return(parts[1])
    }
  })
  return(unlist(char_vector_parts))
}

load(paste0(data_path,'/Transition_MP.rds'))
load(file = paste0(data_path,'/CCI_summary_MP.rds'))
# write.csv(Trans_genes, file = paste0(data_path, '/Transition_MP.csv'), quote = F, row.names = F)

MP_LR <- Trans_genes[Trans_genes$MP_Gene %in% Ligand_Receptor$Ligand_Receptor,]
# saveRDS(scRNA_sub, file = paste0(outputDir, organ ,'_TME_count.rds'))
indir = '/data/rluo4/All/Output/sdata_TME/'
exp_df <- NULL
for (obs in pheno){
  # obs <- 'Esophagus_Epi'
  organ <- str_split(obs,'_', simplify = T)[,1]
  print(organ)
  # organ = 'Liver'#'Esophagus'#'CRC'#'Breast'#'HNSCC'#'Cervix'
  path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
  file = paste0(path_PCT, 'PCT-',organ,'.rds') # new version
  load(file)
  scRNA_sub <- readRDS(file = paste0(indir, organ ,'_TME_count.rds'))
  table(scRNA_sub$cell_class)
  # scRNA_sub$cell_class <- gsub('ADJ','Healthy', scRNA_sub$cell_class)
  Idents(scRNA_sub) <- scRNA_sub$cell_class
  table(Idents(scRNA_sub))
  subtype <- data.frame(subset =  table(scRNA_sub$cell_class))
  subtype <- subtype[subtype$subset.Freq !=0,]
  subtype <- subtype[subtype$subset.Freq >=10,]
  scRNA_sub <- scRNA_sub[, scRNA_sub$cell_class %in% subtype$subset.Var1]
  table(scRNA_sub$cell_class)
  # summary(scRNA_sub@assays$RNA@data['ANXA1',scRNA_sub$Sample_Type=='HSIL_HPV'])
  # summary(scRNA_sub@assays$RNA@data['ANXA1',scRNA_sub$cell_class=='HSIL_HPV_INCAF'])
  # summary(scRNA_sub@assays$RNA@data['ANXA1',scRNA_sub$cell_class=='Healthy_INCAF'])
  # summary(scRNA_sub@assays$RNA@data['ANXA1',scRNA_sub$cell_class=='CC_INCAF'])
  table(scRNA_sub$Cell_Type %in% Lineage$Minor_type)
  scRNA_sub$Lineage <- Lineage$Major_type[match(scRNA_sub$Cell_Type, Lineage$Minor_type)]
  table(scRNA_sub$Lineage)
  index = scRNA_sub$Cell_Type =='STM'
  scRNA_sub$Lineage[index] = scRNA_sub$cell_class[index]
  table(scRNA_sub$Lineage)
  
  # test <- subset(scRNA_sub, Lineage %ni% c('Tcell','Others','Bcell') )
  test <- subset(scRNA_sub, Lineage %in% c('Tcell','Others','Bcell', 'Endothelia','Mesenchymal','Myeloid') )
  table(test$Lineage)
  index = test$Cell_Type !='STM'
  test$Lineage[index] = paste(test$Sample_Type, test$Lineage, sep = '_')[index]
  table(test$Lineage)
  test$Lineage <- gsub('ADJ','Healthy', test$Lineage)
  table(Idents(test))
  Idents(test) = test$Lineage
  test <- SCTransform(test, ncells = 3000, verbose = FALSE) 
  # This function calls sctransform::vst. The sctransform package is available at https://github.com/satijalab/sctransform. Use this function as an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow. Results are saved in a new assay (named SCT by default) with counts being (corrected) counts, data being log1p(counts), scale.data being pearson residuals; sctransform::vst intermediate results are saved in misc slot of new assay.
  # ##############################################################################################################################
  
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='Healthy_Mesenchymal'])
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='BRCA1-mut_Mesenchymal'])
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='DCIS_Mesenchymal'])
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='IDC_Mesenchymal'])
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='Healthy_STM'])
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='T1-C14'])
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='T1-C9'])
  #   summary(test@assays$RNA@data['NECTIN2',test$Lineage=='T1-C4'])
  
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='Healthy_Mesenchymal'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='HSIL_HPV_Mesenchymal'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='CC_Mesenchymal'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='Healthy_STM'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='T1-C5'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='T1-C1'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='T1-C14'])
  
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_STM'])#Endometrium
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='ADJ_STM'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T2-C4'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T2-C19'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T2-C24'])
  
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_STM'])#CRC
  #   # summary(test@assays$RNA@data['ANXA1',test$Lineage=='ADJ_STM'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C21'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C1'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T3-C17'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T3-C23'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_Myeloid'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='AD_Myeloid'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='CRC_Myeloid'])
  
  
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_Mesenchymal'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='EOLP_Mesenchymal'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='OSCC_Mesenchymal'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_STM'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C18'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C11'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C4'])
  
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_STM'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C12'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C10'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C7'])
  #   summary(test@assays$RNA@data['FPR1',test$Lineage=='Healthy_Myeloid'])
  #   # summary(test@assays$RNA@data['FPR1',test$Lineage=='ESCC_Myeloid'])
  
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_Mesenchymal'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Cirrhotic_Mesenchymal'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='HCC_Mesenchymal'])
  #   summary(test@assays$RNA@data['ANXA1',test$cell_class=='Healthy_PFIB'])
  #   summary(test@assays$RNA@data['ANXA1',test$cell_class=='Cirrhotic_PFIB'])
  #   summary(test@assays$RNA@data['ANXA1',test$cell_class=='HCC_PFIB'])
  
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='Healthy_STM'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C1'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T1-C0'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T2-C15'])
  #   summary(test@assays$RNA@data['ANXA1',test$Lineage=='T2-C7'])
  
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='Healthy_Myeloid'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='Cirrhotic_Myeloid'])
  #   summary(test@assays$RNA@data['FPR2',test$Lineage=='HCC_Myeloid'])
  ##############################################################################################################################
  for (mp in unique(MP_LR$MPs)) {
    # mp <-'Epithelial-senescence'
    print(mp)
    markers <- MP_LR$MP_Gene[MP_LR$MPs==mp]
    colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
    colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
    gene <- intersect(markers, rownames(test))
    ggplotdata <- as.data.frame(t(test@assays$RNA@counts[gene, ]))
    str(ggplotdata)
    
    p <- DotPlot(test, features = rev(gene))
    TME_pdata <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
    colnames(TME_pdata) <- c('id','features.plot','Pct.exp','Avg.exp.scaled')
    TME_pdata$id <- as.character(TME_pdata$id)
    TME_pdata$Cell_Type <- extract_last_element(TME_pdata$id)
    TME_pdata$Tissue <- remove_last_element(TME_pdata$id)
    TME_pdata$Organ <- organ
    TME_pdata$MP <- mp
    exp_df <<- rbind(exp_df, TME_pdata)
    ## data$id <- factor(data$id, levels = rev(dataTypeLevel))
    # data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
    ## data$id <- factor(data$id, levels = rev(dataTypeLevel))
    levels <- sort(unique(TME_pdata$id))
    TME_pdata$id <- factor(TME_pdata$id, levels = levels)
    plotx <- ggplot(TME_pdata, aes(y = id, x = features.plot)) +        ## global aes
      geom_point(aes(fill = Avg.exp.scaled, size = Pct.exp),
                 color = 'black',
                 shape = 21,
                 stroke = 0.005) +
      # scale_y_discrete(breaks=0:11, labels=paste0("Treg_c", 0:11), limits = rev) +
      scale_x_discrete(limits = rev) +
      xlab("") + ylab("") +
      scale_fill_gradientn(
        colors = c("#5DBCFF", 'lightblue', "white", "#EDCAE0" , "#F484AE"))+
      scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
      ## scale_x_reverse() + scale_y_reverse() +
      theme(
        text = element_text(size = 13),
        panel.grid.major = element_line(colour = "grey90", size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1),
        legend.position="bottom",
        legend.title=element_text(size=12)) +
      guides(size = guide_legend(title.position="top",
                                 title.hjust = 0.5,
                                 ncol = 3,
                                 byrow = T,
                                 override.aes = list(stroke = 0.4)),
             fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
    plotx
    
    outdir <- paste0(indir, organ, "_MP_LR/")
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    
    mp = gsub('/','_',mp)
    # outfile = paste0(outdir, mp, "_TME_bubbleplot.pdf")
    outfile = paste0(outdir, mp, "_STM_bubbleplot.pdf")
    
    hei = ifelse(length(unique(test$Lineage))<=4, length(unique(test$Lineage))*20,length(unique(test$Lineage))*15)
    ggsave(file.path(outfile), plotx, width = 120, height = hei, units = "mm")
    
  }
  
}
# save(exp_df, file = paste0('/data/rluo4/All/Output/MP_exp_df.rds'))
save(exp_df, file = paste0('/data/rluo4/All/Output/MP_STM_exp_df.rds'))

########################################################################
data_path <- '/data/rluo4/All/Output/Cluster/'
indir <-'/data/rluo4/database/'
setwd(data_path)
organ_all <- list.files(data_path)
PCT_Clusters =  vector("list", length(organ_all))
names(PCT_Clusters) <- organ_all
########################################################################
for (organ in organ_all) {
  # Define PCT_Clusters
  ########
  print(organ)
  if(organ == 'Breast'){#
    PCT_Clusters[[organ]] <- c( "T1-C9","T1-C14",'T1-C4')
  }
  if(organ == 'Cervix'){#
    PCT_Clusters[[organ]] <- c("T1-C5", 'T3-C8', "T1-C1",  'T1-C14')# "T2-C4", "T2-C10",  
  }
  if(organ == 'Chen'){# 
    PCT_Clusters[[organ]] <- c('T1-C24',  "T1-C21", "T1-C1")
  }
  if(organ == 'CRC'){# which is Becker
    PCT_Clusters[[organ]] <- c("T1-C17", "T1-C19", "T1-C3")
  }
  if(organ == 'Endometrium'){#
    PCT_Clusters[[organ]] <- c('T1-C12','T1-C14', 'T2-C19', "T2-C4",  "T2-C24")#c( 'T2-C19', "T2-C4",  "T2-C24")
  }
  if(organ == 'Esophagus'){#
    PCT_Clusters[[organ]] <- c("T1-C12", "T1-C10")  # 'T1-C7'
  }
  if(organ == 'GC'){#
    PCT_Clusters[[organ]] <- c("T2-C12", "T2-C7", "T1-C2", "T1-C17")#c("T1-C13", "T1-C2", "T2-C3", "T2-C12")
  }
  if(organ == 'HNSCC'){#
    PCT_Clusters[[organ]] <- c("T1-C18", 'T1-C16', "T1-C11", 'T1-C4')# ,c("T1-C11", "T1-C18", "T2-C19", "T2-C3")
  }
  if(organ=='Liver'){#
    PCT_Clusters[[organ]] <-  c( "T2-C11", 'T2-C15', 'T2-C7', "T1-C1", "T1-C0")#c( "T2-C11", "T1-C1", "T1-C0")
  }  
  if(organ=='Lung'){#
    PCT_Clusters[[organ]] <- c('T1-C6',"T1-C9",  "T1-C8",  "T1-C14")#c("T1-C12", "T1-C6", "T1-C4", "T1-C3")
  }
  if(organ == 'Pancreas'){#
    PCT_Clusters[[organ]] <- c("T2-C6", "T2-C5",  "T2-C0", "T2-C1", "T2-C9", 'T2-C4')#c("T1-C2", "T1-C8", "T2-C9", "T2-C4")
  }
  if(organ == "Prostate"){#
    PCT_Clusters[[organ]] <- c("T2-C14", "T2-C34", 'T1-C42','T1-C5')#c("T1-C6", "T1-C9", "T2-C34", "T2-C14")
  }
  if(organ == "Skin"){#
    PCT_Clusters[[organ]] <- c( 'T1-C9', "T1-C11", "T2-C1",  "T2-C14" ) #,# c("T1-C12", "T1-C9", "T2-C1", "T2-C11")
  } 
  if(organ=="THCA"){#
    PCT_Clusters[[organ]] <- c( "T2-C2",  "T2-C16", "T1-C15", "T1-C5", "T1-C1")
  }
  
}
########################################################################
load(file = paste0('/data/rluo4/All/Output/MP_exp_df.rds'))
load(file = paste0('/data/rluo4/All/Output/MP_STM_exp_df.rds'))
for (organ in organ_all) {
  plotdata <- exp_df  %>% 
    filter(Organ %in% organ & features.plot=='ANXA1' & MP %in% c('Epithelial-senescence')) %>%
    filter(Cell_Type=='STM' | Cell_Type %in% PCT_Clusters[[organ]])
  # index = TME_pdata$Cell_Type=='Mesenchymal'
  # index = TME_pdata$Cell_Type=='STM' | TME_pdata$Cell_Type %in% PCT_Clusters[[organ]]
  # plotdata <- TME_pdata[index,]
  # plotdata <- plotdata[plotdata$features.plot=='ANXA1' & plotdata$MP=='Epithelial-senescence',]
  # plotdata <- plotdata[ plotdata$MP=='Epithelial-senescence',]
  extract_last_element <- function(x) {
    split_string <- strsplit(x, "-")
    last_element <- sapply(split_string, function(y) tail(y, n = 1))
    return(last_element)
  }
  plotdata$Tissue <- extract_last_element(plotdata$Tissue)
  # plotdata$Tissue <- factor(plotdata$Tissue, levels = c('Healthy','C14','C9','C4'))
  plotdata$Tissue <- factor(plotdata$Tissue, levels = c('Healthy', extract_last_element(PCT_Clusters[[organ]])
  ))
  # df<- dt %>% group_by(dose,sex) %>%  summarise( sd= sd(len),  Avg_len= mean(len)) 
  #绘制点线图；
  p1<- ggplot(data=plotdata, aes(x = Tissue, y = Pct.exp, group = features.plot)) +
    geom_line(aes(color=features.plot),size=0.5,show.legend=T) +
    geom_point(aes(color=features.plot, size = Pct.exp),shape=21,fill= "white",
               size=2,show.legend=F)
  # mycolor<- c( "#0077c1", "#00a99e", "#6bc72b", "#ff5a20", "#ff1620", "#752995")
  # p2<- p1 + scale_colour_manual(values=alpha(mycolor[c(3,4)],1))
  # #添加误差线；
  # p3<- p2+geom_errorbar(aes(ymin = Avg_len-sd/3,
  #                           ymax= Avg_len+sd/3,
  #                           color=sex),
  #                       width=0.1,show.legend=F)
  # plotdata$Tissue[plotdata$Tissue!='Healthy'] <- paste(plotdata$Organ, plotdata$Tissue, collapse = "-")[plotdata$Tissue!='Healthy']
  p = ggplot(plotdata, aes(x = Tissue,  group = features.plot)) +
    # geom_point(aes(y = Avg.exp.scaled, color = "Avg.exp.scaled"), size = 3) +  # Dot plot for Avg.exp
    geom_line(aes(y = Pct.exp), size = 0.5, color = 'gray50') +  # Line plot for Avg.exp
    # geom_point(aes(y = Avg.exp.scaled, size = Pct.exp, color = "Pct.exp")) +  # Line plot with proportion as size
    geom_point(aes(y = Pct.exp, fill = Avg.exp.scaled, size = Pct.exp),
               color = 'black',
               shape = 25,
               stroke = 0.005) +
    # scale_x_discrete(limits = rev) +
    xlab("") + ylab("Percentage of cells (%)") + #Pct. of cells expressing each gene
    scale_fill_gradientn(
      colors = c( "#0077c1", "#6DCCFF", "white", colorsForDataType[3],"#ff1620"),
      limits = max(abs(plotdata$Avg.exp.scaled)) * c(-1, 1),
      guide = guide_colorbar(title.position = "top"),#oob = squish
    )+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    ## scale_x_reverse() + scale_y_reverse() +
    theme( text = element_text(size = 12),
           panel.grid.major = element_line(colour = "grey90", size=0.2),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black"),
           # axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
           axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1, color = 'black', size = 10),
           axis.text.x   = element_text(color = 'black', size = 13.5, angle = 40,
                                        hjust = 1,vjust = 1),
           axis.title.x  = element_text(color = 'black', size = 14, angle = 0),
           axis.title.y  = element_text(color = 'black', size = 12.5, angle = 90),
           axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
           axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
           legend.text   = element_text(color = 'black', size   = 10),
           legend.position="top", legend.key.height = unit(0.25, "cm"),
           legend.title=element_text(size=12)) +
    # theme(panel.background = element_rect(fill = NA),
    #       plot.margin = margin(t=10,r=10,b=5,l=5,unit = "mm"),
    #       # axis.ticks.y = element_blank,
    #       axis.ticks.x = element_line(colour = "grey40",size = 0.5),
    #       axis.line = element_line(colour = "grey40",size = 0.5),
    #       panel.grid.major.y = element_line(colour = NA,size = 0.5),
    #       panel.grid.major.x = element_blank()) + 
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               ncol = 3,
                               byrow = T,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
  p
  outdir = paste0('/data/rluo4/All/Output/res_fig/ANXA1_fig/')
  wid = length(unique(plotdata$Cell_Type))
  wid = ifelse(wid <=3, wid*1.38, ifelse(wid>4, wid*0.78, wid*0.98))
  ggsave(plot=p,height = 3.5, width=wid,
         filename=paste0(outdir, organ, "_STM.png"), dpi = 500, device = "png")
  
}

precancer <- c('AAH','AD','AEH','AK',#'AIS',
               'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
               'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
               'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
               'SER','SIM','WIM')

load(file = paste0('/data/rluo4/All/Output/MP_exp_df.rds'))

for (organ in organ_all) {
  plotdata <- exp_df  %>% 
    filter(Organ %in% organ & features.plot=='ANXA1' & MP %in% c('Epithelial-senescence')) %>%
    filter(Cell_Type=='Mesenchymal')
  
  # extract_last_element <- function(x) {
  #   split_string <- strsplit(x, "-")
  #   last_element <- sapply(split_string, function(y) tail(y, n = 1))
  #   return(last_element)
  # }
  # plotdata$Tissue <- extract_last_element(plotdata$Tissue)
  # plotdata$Tissue <- factor(plotdata$Tissue, levels = c('Healthy','C14','C9','C4'))
  plotdata$Tissue <- factor(plotdata$Tissue, levels = c('Healthy', plotdata$Tissue[plotdata$Tissue %in% precancer],
                                                        plotdata$Tissue[! plotdata$Tissue %in% c('Healthy',precancer)]) )
  
  # df<- dt %>% group_by(dose,sex) %>%  summarise( sd= sd(len),  Avg_len= mean(len)) 
  #绘制点线图；
  p1<- ggplot(data=plotdata, aes(x = Tissue, y = Pct.exp, group = features.plot)) +
    geom_line(aes(color=features.plot),size=0.5,show.legend=T) +
    geom_point(aes(color=features.plot, size = Pct.exp),shape=21,fill= "white",
               size=2,show.legend=F)
  # mycolor<- c( "#0077c1", "#00a99e", "#6bc72b", "#ff5a20", "#ff1620", "#752995")
  # p2<- p1 + scale_colour_manual(values=alpha(mycolor[c(3,4)],1))
  # #添加误差线；
  # p3<- p2+geom_errorbar(aes(ymin = Avg_len-sd/3,
  #                           ymax= Avg_len+sd/3,
  #                           color=sex),
  #                       width=0.1,show.legend=F)
  p = ggplot(plotdata, aes(x = Tissue,  group = features.plot)) +
    # geom_point(aes(y = Avg.exp.scaled, color = "Avg.exp.scaled"), size = 3) +  # Dot plot for Avg.exp
    geom_line(aes(y = Pct.exp), size = 0.5, color = 'gray50') +  # Line plot for Avg.exp
    # geom_point(aes(y = Avg.exp.scaled, size = Pct.exp, color = "Pct.exp")) +  # Line plot with proportion as size
    geom_point(aes(y = Pct.exp, fill = Avg.exp.scaled, size = Pct.exp),
               color = 'black',
               shape = 25,
               stroke = 0.005) +
    # scale_x_discrete(limits = rev) +
    xlab("") + ylab("Percentage of cells (%)") + #Pct. of cells expressing each gene
    scale_fill_gradientn(
      colors = c( "#0077c1", "#6DCCFF", "white", colorsForDataType[3],"#ff1620"))+
    scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
    ## scale_x_reverse() + scale_y_reverse() +
    theme( text = element_text(size = 12),
           panel.grid.major = element_line(colour = "grey90", size=0.2),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black"),
           # axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
           axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1, color = 'black', size = 10),
           axis.text.x   = element_text(color = 'black', size = 13.5, angle = 40,
                                        hjust = 1,vjust = 1),
           axis.title.x  = element_text(color = 'black', size = 14, angle = 0),
           axis.title.y  = element_text(color = 'black', size = 12.5, angle = 90),
           axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
           axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
           legend.text   = element_text(color = 'black', size   = 10),
           legend.position="top", legend.key.height = unit(0.25, "cm"),
           legend.title=element_text(size=12)) +
    # theme(panel.background = element_rect(fill = NA),
    #       plot.margin = margin(t=10,r=10,b=5,l=5,unit = "mm"),
    #       # axis.ticks.y = element_blank,
    #       axis.ticks.x = element_line(colour = "grey40",size = 0.5),
    #       axis.line = element_line(colour = "grey40",size = 0.5),
    #       panel.grid.major.y = element_line(colour = NA,size = 0.5),
    #       panel.grid.major.x = element_blank()) + 
    guides(size = guide_legend(title.position="top",
                               title.hjust = 0.5,
                               ncol = 3,
                               byrow = T,
                               override.aes = list(stroke = 0.4)),
           fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
  outdir = paste0('/data/rluo4/All/Output/res_fig/ANXA1_fig/')
  wid = length(unique(plotdata$Cell_Type))
  wid = ifelse(wid <=3, wid*2, ifelse(wid>4, wid*0.7, wid*1.2))
  ggsave(plot=p,height = 3.5, width=wid,
         filename=paste0(outdir, organ, "_Mesenchymal.png"), dpi = 500, device = "png")
  
}

# EpiSen_df <- NULL
# library(stringr)
# data_path = '/data/rluo4/All/Output'
# setwd(data_path)
# organ_all <- list.files('Cluster/')
# organ_all
# for (organ in organ_all) {
#   # organ = 'Cervix'
#   # save(newSTM, file = paste0('/data/rluo4/All/Output/Epi_Results/',organ, '_newSTM.Rdata'))
#   load(paste0('/data/rluo4/All/Output/Epi_Results/',organ, '_newSTM.Rdata'))
#   # newSTM$cell_class <- gsub('ADJ','Healthy', newSTM$cell_class)
#   Idents(newSTM) <- newSTM$cell_class
#   
# # for (mp in MPs) {
#   mp <-'Epithelial-senescence'
#   mp
#   markers <- MP_LR$MP_Gene[MP_LR$MPs==mp]
#   colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
#   colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
#   # dataTypeLevel <- c("CD4", "CD8", "MAIT", "Tgd", "NKT", "Proliferative")
#   gene <- intersect(markers, rownames(newSTM))
#   ggplotdata <- as.data.frame(t(test@assays$RNA@counts[gene, ]))
#   str(ggplotdata)
#   
#   p <- DotPlot(test, features = rev(gene))
#   STM_pdata <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
#   colnames(STM_pdata)  <- c('id','features.plot','Pct.exp','Avg.exp.scaled')
#   STM_pdata$organ <- organ
#   EpiSen_df <<- rbind(EpiSen_df, STM_pdata)
#   
# }
#   
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
# data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
## data$id <- factor(data$id, levels = rev(dataTypeLevel))
plotx <- ggplot(STM_pdata, aes(y = id, x = features.plot)) +        ## global aes
  geom_point(aes(fill = Avg.exp.scaled, size = Pct.exp),
             color = 'black',
             shape = 21,
             stroke = 0.005) +
  # scale_y_discrete(breaks=0:11, labels=paste0("Treg_c", 0:11), limits = rev) +
  scale_x_discrete(limits = rev) +
  xlab("") + ylab("") +
  scale_fill_gradientn(
    colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
  scale_size(range = c(0, 3.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
  ## scale_x_reverse() + scale_y_reverse() +
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
    axis.text.y=element_text(angle = -30, vjust = 1, hjust = 1),
    legend.position="bottom",
    legend.title=element_text(size=12)) +
  guides(size = guide_legend(title.position="top",
                             title.hjust = 0.5,
                             ncol = 3,
                             byrow = T,
                             override.aes = list(stroke = 0.4)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
plotx

outdir <- paste0('/data/rluo4/All/Output/Epi_Results/', organ, "_res_fig/MP_LR/")
mp = gsub('/','_',mp)
outfile = paste0(outdir, mp, "_subcluster_bubbleplot.pdf")
hei = length(unique(newSTM$cell_class))*15
# ggsave(file.path(outfile), plotx, width = 120, height = hei, units = "mm")
# }
library(ggplot2)
#features.plot <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")
colorForClass1 <- c("#C4DA5D", "#6DCCDD", "#F494BE", "#EDCAE0")
# dataTypeLevel <- c("CD4", "CD8", "MAIT", "Tgd", "NKT", "Proliferative")
gene <- intersect(markers, rownames(newSTM))
g <-  DotPlot(STM, features = rev(gene))#DotPlot(object = STM, features = features.plot, assay="RNA")
g$layers
# [[1]]
# mapping: size = ~pct.exp, colour = ~avg.exp.scaled 
# geom_point: na.rm = FALSE
# stat_identity: na.rm = FALSE
# position_identity 
g$data
#              avg.exp   pct.exp features.plot     id avg.exp.scaled
# LYZ     3.303911e+79 60.407407           LYZ pbmc3k            NaN
# CCL5    1.707925e+25 31.925926          CCL5 pbmc3k            NaN
# IL32    3.126373e+23 55.407407          IL32 pbmc3k            NaN
# PTPRCAP 1.608967e+12 64.037037       PTPRCAP pbmc3k            NaN
# FCGR3A  8.475464e+07 18.407407        FCGR3A pbmc3k            NaN
# PF4     3.125317e+23  1.592593           PF4 pbmc3k            NaN
g$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
g <- g + 
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  theme(axis.text.x = element_text(angle=90))
plot(g)
plots <- DotPlot(STM, features = rev(gene),assay=NULL,group.by ="Sample_Type",scale.by = "size") + RotatedAxis()
print(plots)
summary((colon_full@assays$RNA@data['ANXA1',colon_full$Sample_Type=='OSCC']))

# plot example line plots (Figure 3B)
colors <- fread("../CCLE_featurematrix_NK_PRISM/nk_crispr_colors.txt", data.table = F)
cols_example <- colors$color[colors$cancer_type %in% c("B-ALL", "AML")]
names(cols_example) <- c("THP1 (AML)", "697 (B-ALL)")


# data_path <- '/data/rluo4/All/Output/Cluster/'
# organ_all <- list.files(data_path)
# setwd(data_path)
# setwd('../')
# # for (organ in organ_all) {
# typeCCI <- list.files('CellChat')
# typeCCI <- typeCCI[grep('RData',typeCCI)]
# typeCCI <- gsub('_CCI.RData','',typeCCI)
# print(unique(typeCCI))
# 
# CCI_dir <- "CCI_figures"
# if (!dir.exists(paste0(CCI_dir))){
#   dir.create(paste0(CCI_dir))
# }
# CCI_dir = '/data/rluo4/All/Output/CCI_figures'
# setwd(CCI_dir)
# typeCCI =  vector("list", length(organ_all))
# names(typeCCI) <- organ_all
# ########################################################################
# for (organ in organ_all) {
#   # Define typeCCI
#   ########
#   print(organ)
#   if(organ == 'Breast'){#
#     typeCCI[[organ]] <- c('brca10_brca3_ctrl6_clusterC14', "brca1_brca2_brca3_clusterC9", 'brca10_brca2_brca3_DCIS2_M1_NCCBC14_NCCBC3_NCCBC5_P1_P2_clusterC4')#c("T1-C4", "T1-C9")
#   }
#   if(organ == 'Cervix'){#
#     typeCCI[[organ]] <- c('CA_HPV_1_CA_HPV_2_HSIL_HPV_1_HSIL_HPV_2_clusterC1', 'CA_HPV_1_CA_HPV_2_CA_HPV_3_HSIL_HPV_1_HSIL_HPV_2_clusterC5', 'CA_HPV_1_CA_HPV_2_CA_HPV_3_HSIL_HPV_1_clusterC8', 'CA_HPV_1_CA_HPV_2_CA_HPV_3_HSIL_HPV_1_HSIL_HPV_2_clusterC14',  'CA_HPV_1_clusterC6') #'CA_HPV_2_CA_HPV_3_HSIL_HPV_1_HSIL_HPV_2_clusterC4', 'HSIL_HPV_1_HSIL_HPV_2_clusterC10',
#   }
#   if(organ == 'Chen'){# 
#     typeCCI[[organ]] <- c('CRC1_8810_HTA11_10711_200000101113111_HTA11_99999974143_8462013111_clusterC1', 'HTA11_10711_200000101113111_HTA11_866_300476101113111_HTA11_99999971662_8245713111_HTA11_99999971662_8245713211_clusterC21', 'HTA11_1391_200000101113211_HTA11_3361_200000101113311_HTA11_347_200000101113111_HTA11_347_200000101113211_HTA11_696_200000101113111_HTA11_696_200000101113211_clusterC17','HTA11_1391_200000101113211_HTA11_3361_200000101113311_HTA11_696_200000101113111_HTA11_696_200000101113211_clusterC23')# 'A001-C-007_HTA11_99999974143_8462013111_clusterC26','HTA11_3410_200000101113111_HTA11_696_200000101113111_clusterC3',#c("T1-C1", "T1-C24", "T2-C26", "T2-C3", "T3-C17", "T3-C23")
#   }
#   if(organ == 'CRC'){# which is Becker
#     typeCCI[[organ]] <- c('A015-C-005_CRC3_11773_clusterC3', 'A002-C-010-R0_A002-C-114_A002-C-116_A002-C-201_A015-C-002_A015-C-006_A015-C-203_A015-C-204_CRC3_11773_clusterC19','A001-C-014_A001-C-207_A002-C-010-R0_A002-C-016_A002-C-114_A002-C-201_A002-C-205_A015-C-002_F034_clusterC17','A001-C-119_A002-C-114_A002-C-205_A015-C-002_A015-C-006_A015-C-104_A015-C-203_A015-C-204_A018-E-020_F034_clusterC16')
#   }
#   if(organ == 'Endometrium'){#
#     typeCCI[[organ]] <- c( 'AEH-subject5_EEC-subject4_clusterC24','AEH-subject5_clusterC19', 'AEH-subject5_clusterC4')#,'AEH-subject2_AEH-subject4_EEC-subject1_clusterC14','AEH-subject2_AEH-subject3_clusterC12')#
#   }
#   if(organ == 'Esophagus'){#
#     typeCCI[[organ]] <- c( 'P104T-E_P127T-E_P128T-E_P15T-E_P16T-E_P1T-E_P21T-E_P23T-E_P2T-E_P31T-E_P39T-E_P47T-E_P54T-E_P57T-E_P61T-E_P65T-E_P74T-E_P75T-E_P76T-E_P79T-E_P82T-E_P84T-E_P8T-E_clusterC10','LZE11D_LZE20T_LZE21D1_LZE21T_LZE22D1_LZE22D3_LZE22T_LZE24D1_LZE24T_LZE2D_LZE2T_LZE3D_LZE4T_LZE5T_LZE6T_LZE7T_LZE8T_P107T-E_P130T-E_P16T-E_P75T-E_P76T-E_clusterC12' ) 
#   }
#   if(organ == 'GC'){#
#     typeCCI[[organ]] <- c('EGC_Pt1_Superficial_SIM_1_SIM_2_SIM_4_clusterC3', 'Pt1_Superficial_SIM_1_SIM_2_clusterC7', 'Pat01-B_Pat02-B_Pat03-B_Pat16-B_Pat22-B_clusterC17', 'Pat02-B_Pat04-B_Pat05-B_Pat06-B_Pat09-B_Pat15-B_Pat18-B_Pat22-B_Pat24-B_Pat25-A_clusterC2')#c("T1-C13", "T1-C2", "T2-C3", "T2-C12")
#   }
#   if(organ == 'HNSCC'){#
#     typeCCI[[organ]] <- c('C08_C21_C43_C57_clusterC4', 'C09_C43_clusterC11', 'C43_clusterC19','C08_C21_C43_C46_C57_EOLP-2_clusterC3')#'EOLP-2_clusterC18', 'EOLP-1_EOLP-2_clusterC16')#   #c("T1-C11", "T1-C18", "T2-C19", "T2-C3")
#   }
#   if(organ=='Liver'){#
#     typeCCI[[organ]] <- c('cirrhotic1_cirrhotic2_cirrhotic3_HCC2_Meng_clusterC0','cirrhotic1_cirrhotic2_cirrhotic3_clusterC1','HCC1_Meng_Pt13.a_Pt13.b_Pt14.d_clusterC7','cirrhotic1_cirrhotic2_cirrhotic3_clusterC15')#'cirrhotic1_HCC1_Meng_HCC2_Meng_clusterC6',,'cirrhotic1_cirrhotic2_cirrhotic3_clusterC11')#  #c("T1-C0", "T1-C1", "T2-C7", "T2-C11")
#   }  
#   if(organ=='Lung'){#
#     typeCCI[[organ]] <- c('TD9_clusterC14','RNA-P10T2-P10T2-1_RNA-P10T2-P10T2-3_RNA-P25T1-P25T1-1_RNA-P25T1-P25T1-3_RNA-P25T1-P25T1-4_clusterC6','RNA-P10T1-P10T1-2_clusterC7','RNA-P7T1-P7T1-1_RNA-P7T1-P7T1-2_RNA-P7T1-P7T1-3_RNA-P7T1-P7T1-4_clusterC3')#'RNA-P23T2-P23T2-2_RNA-P23T2-P23T2-4_RNA-P24T2-P24T2-2_RNA-P6T1-P6T1-1_RNA-P6T1-P6T1-4_clusterC8',  'RNA-P5T2-P5T2-1_RNA-P5T2-P5T2-3_clusterC9',#'RNA-P17T-P17T-2_RNA-P17T-P17T-4_RNA-P17T-P17T-6_RNA-P17T-P17T-8_clusterC12', #c("T1-C12", "T1-C6", "T1-C4", "T1-C3")
#   }
#   if(organ == 'Pancreas'){#
#     typeCCI[[organ]] <- c('HTA12-9-1_HTA12-9-2_clusterC4','4347-EC_HTA12-9-1_clusterC9', 'HTA12-16-5_clusterC8', '4741-EC2_clusterC3')#c("T1-C8", "T1-C2", "T2-C4", "T2-C9")
#   }
#   if(organ == "Prostate"){#
#     typeCCI[[organ]] <- c('Dong_P1_Dong_P3_clusterC34', 'GSM5252126_BPH283PrGF_Via_GSM5252128_BPH327PrGF_Via_GSM5252130_BPH340PrGF_Via_GSM5252131_BPH340PrSF_Via_clusterC14', 'GSM5252130_BPH340PrGF_Via_GSM5252131_BPH340PrSF_Via_clusterC42','048752_1579-all-cells_052095_1628-all-cells_052097_1595-all-cells_GSM5252126_BPH283PrGF_Via_GSM5252130_BPH340PrGF_Via_GSM5252131_BPH340PrSF_Via_GSM5252132_BPH389PrGF_GSM5252134_BPH511PrG_Fcol_3GEX_GSM5252135_BPH511PrPUr_Fcol_3GEX_clusterC5')#c("T1-C6", "T1-C9", "T2-C34", "T2-C14")
#   }
#   if(organ == "Skin"){#
#     typeCCI[[organ]] <- c("P1_S1_AK_P2_S3_AK_clusterC11", 'P2_S4_SCCIS_clusterC1', 'P1_S1_AK_P2_S3_AK_P3_S6_AK_clusterC9')#,'P2_S4_SCCIS_clusterC14', 'P1_S1_AK_P2_S3_AK_P2_S4_SCCIS_P3_S6_AK_clusterC4', 'P1_S1_AK_P2_S3_AK_P2_S4_SCCIS_clusterC5', "P2_S3_AK_P2_S4_SCCIS_clusterC12",  'P2_S3_AK_P2_S4_SCCIS_clusterC17') # c("T1-C12", "T1-C9", "T2-C1", "T2-C11")
#   } 
#   if(organ=="THCA"){#
#     typeCCI[[organ]] <- c('Adj_PTCwithHT_6_PTCwithHT_6_PTCwithHT_8_clusterC5', 'Adj_PTCwithHT_6_Adj_PTCwithHT_8_PTCwithHT_6_clusterC2', 'Adj_PTCwithHT_8_PTCwithHT_6_PTCwithHT_8_clusterC1', 'Adj_PTCwithHT_6_Adj_PTCwithHT_8_PTCwithHT_6_clusterC15')
#   }
#   
# }
# 
# typeCCI_use <- unlist(typeCCI) 
# library(ArchR)
# for (TissueType in typeCCI_use) {
#   data_path = '/data/rluo4/All/Output/CCI_figures'
#   setwd(data_path)
#   cohort_directory <- paste0(data_path,'/', TissueType)
#   print(paste0(TissueType, ' from dir: ', cohort_directory))
#   # load(paste0(cohort_directory,'/CellChat/',TissueType, '_CCI.RData'))
#   load(paste0(data_path,'/../CellChat/',TissueType, '_CCI.RData'))
# 
#   if (!dir.exists(paste0(cohort_directory))){
#     dir.create(paste0(cohort_directory))
#   }
#   setwd(cohort_directory)  
#   if( file.exists('STM_source_net_bubble.png')){
#     print(paste0(TissueType, " is already over!"))
#     next;
#   }
#   #展示每个亚群作为source的信号传递
#   print(table(cellchat@idents))
#   
#   cell_color <- paletteDiscrete(values = unique(cellchat@meta$Cell_Type))
#   color_data <- data.frame(cell_color)
#   for (i in 1:length(cell_color)){
#     color1 <- c(cell_color[i], setdiff(cell_color, cell_color[i]))
#     color1 <- as.data.frame(color1)
#     color_data <- cbind(color_data, color1)
#   }
#   color_data <- color_data[, -1]
#   #展示每个亚群作为source的信号传递
#   mat <- cellchat@net$weight
#   print(dim(mat))
#   # par(mfrow = c(3,4), xpd=TRUE)
#   groupSize <- as.numeric(table(cellchat@idents))
#   # wid <- nrow(mat)*0.8
#   # hei <- nrow(mat)
#   # # par(mfrow = c(2,4), xpd=TRUE,mar=c(2,2,2,2))
#   # file = paste0(cohort_directory,'/CellChat/','CellChat_',TissueType,'.pdf')
#   # pdf(file)
#   for (i in 1:nrow(mat)) {
#     mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#     mat2[i, ] <- mat[i, ]
#     name <- rownames(mat2)[i]
#     mat2 <- as.data.frame(mat2)
#     data <- data.frame(mat2[i,])
#     data <- select(data,i,everything())
#     # remove celltype with weigth==0
#     logi <- table(data[1,]==0)
#     if( names(logi)[1]=='FALSE' | as.numeric(logi[1]) != ncol(data)){
#       data1 <- data.frame(mat2[-i, ])
#       data1 <- select(data1,i,everything())
#       
#       data2 <- rbind(data, data1)
#       data2 <- as.matrix(data2)
#       # name = gsub('A/X','A.X',name)
#       hei.wid = ifelse(nrow(mat)>25,ifelse(nrow(mat)>35,12,10),8)
#       file = paste0(name,'_net_circle.png')
#       png(file, width = hei.wid, height = hei.wid, res = 500,units = 'in')
#       netVisual_circle(data2,
#                        vertex.weight = groupSize,
#                        weight.scale = T,
#                        edge.weight.max = max(mat),
#                        title.name = rownames(data2)[1],
#                        color.use = color_data[, i])
#       dev.off()
#       # dir = paste0(cohort_directory,'/CellChat/','CellChat_',TissueType)
#       # ggsave(height=4,width=4, filename=file, dpi = 300, device = "png")
#     }
#   }
#   # table(cellchat@meta$Cell_Type)
#   print(table(cellchat@idents))
#   cell_use <- data.frame(table(cellchat@idents))
#   dim(cellchat@net$weight)
#   # rownames(cell_use) <- 1:nrow(cell_use)
#   n <- as.numeric(rownames(cell_use))
#   p <- NULL
#   plot <- lapply(setNames(n, n), function(nameindex) {
#     # for (nameindex in n) {
#     cell <- cell_use[rownames(cell_use) %in% nameindex,]$Var1
#     print(cell)
#     # cell <- gsub('A/X','A.X',cell)
#     # if(organ =='Esophagus' & TissueType=='Healthy' & cell =='PLA'| organ =='Pancreas' & TissueType=='HCC' & cell=='END'){
#     # Check if a file name exists in a directory
#     file_name <- paste0(cell,'_net_circle.png')
#     file_exists <- file_name %in% list.files()
#     if(! file_exists){
#       print(paste0("No interactions are detected for ",cell," in ",organ))
#       # stop(paste0("No interactions are detected in PLA of Esophagus"))
#     } else{
#       weight = apply(mat, 2, function(x){
#         y <- table(x==0)
#         logi <- (names(y)[1]=='FALSE' & length(y)==1) |as.numeric(y)[1] != length(x)
#         return(logi)
#       })
#       target_cell <- cell_use[weight,]
#       #netVisual_bubble(cellchat, remove.isolate = FALSE)
#       wid = ifelse(nrow(target_cell)>=25,nrow(target_cell)*0.25,6)
#       t <- as.numeric(rownames(target_cell))
#       
#       file = paste0(cell,'_source_net_bubble.png')
#       p_s <- netVisual_bubble(cellchat, sources.use = nameindex, targets.use = n, remove.isolate = TRUE)
#       # print(unique(p_s$data$interaction_name_2))
#       hei = length(unique(p_s$data$interaction_name_2))*0.15
#       hei = ifelse(hei>6,hei,6)
#       ggsave(file, width = wid, height = hei, dpi = 500, device = "png")
#       
#       if(nameindex %in% t){
#         file = paste0(cell,'_target_net_bubble.png')
#         p_t <- netVisual_bubble(cellchat, sources.use = n, targets.use = nameindex, remove.isolate = TRUE)
#         hei = length(unique(p_t$data$interaction_name_2))*0.15
#         hei = ifelse(hei>6,hei,6)
#         ggsave(file, width = wid, height = hei, dpi = 500, device = "png")
#       }
#       p <<- c(p,cell)#
#       return(p)
#     }
#   })
#   
# }
# 
# immune_cell_use <- c('STM', 'pDC', 'cDC', 'GDT', 'MDSC', 'M2MAC', 'PLA', 'MON', 'MAST', 'DC', 'INMON', 'CD8TEXP', 'NEUT', 'BN', 'TREG', 'CD8TEXINT', 'CD8TEREX', 'NK', 'KUP', 'M1MAC', 'MDSCs', 'LC', 'CD4TN', 'CD8TEX')
# 
# immune_cell_use <- c('STM','pDC', 'cDC',  'DC','MDSC', 'MDSCs',  'MON',  'INMON', 'MAST','NEUT', 'KUP', 'MAC', 'M1MAC','M2MAC', 'LC', 'ALVMAC'#,
#                     )# 'CD8TEXP', 'BN',  'PLA','TREG', 'CD8TEXINT', 'CD8TEREX', 'NK',  'CD4TN', 'CD8TEX')
# 
# stromal_cell_use <- unique(Lineage$Minor_type[Lineage$Major_type=='Mesenchymal'])#c("STM",'FIB','ICAF','INCAF','MYOFIB', 'PFIB','VFIB','CFIB', 'CAF', 'PERI','HSC', 'PSC', 'ECM')#'MVA', 'PVA','MSC.MVA','MSC.ADIPO',
# ANXA1 <- CCI_summary_All$DiseaseStage[CCI_summary_All$Pathway=='ANNEXIN']
# ANXA1
# 
# if(nrow(df.net.ANXA1) != 0){
#   targets <- as.numeric(rownames(cell_use[cell_use$Var1 %in% df.net.ANXA1$target,]))
# } else{
#   targets <- as.numeric(rownames(cell_use[cell_use$Var1 %in% df.net.ANXA1$target,]))
# }
# # cell <- gsub('A/X','A.X',cell)
# # if(organ =='Esophagus' & TissueType=='Healthy' & cell =='PLA'| organ =='Liver' & TissueType=='HCC' & cell=='END'){
# # Check if a file name exists in a directory
# file_name <- paste0(cell,'_net_circle.png')
# file_exists <- file_name %in% list.files()
# # if(! file_exists){
# #   print(paste0("No interactions are detected for ",cell," in ",organ))
# #   # stop(paste0("No interactions are detected in PLA of Esophagus"))
# # } else{
# weight = apply(mat, 2, function(x){
#   y <- table(x==0)
#   logi <- (names(y)[1]=='FALSE' & length(y)==1) |as.numeric(y)[1] != length(x)
#   return(logi)
# })
# target_cell <- cell_use[weight,]
# #netVisual_bubble(cellchat, remove.isolate = FALSE)
# # netVisual_aggregate(cellchat, signaling = 'ANNEXIN',  
# #                     vertex.receiver = vertex.receiver,layout = "hierarchy")
# wid = ifelse(nrow(target_cell)>=25,nrow(target_cell)*0.25,6)
# t <- as.numeric(rownames(target_cell))
# 
# file = paste0(cell,'_source_net_bubble.png')
# p_s <- netVisual_bubble(cellchat, sources.use = nameindex, targets.use = targets,  remove.isolate = TRUE)# signaling = c("ANNEXIN",'TGFb','WNT'),
# # print(unique(p_s$data$interaction_name_2))
# hei = length(unique(p_s$data$interaction_name_2))*0.15
# hei = ifelse(hei>6,hei,6)
# ggsave(file, width = wid, height = hei, dpi = 500, device = "png")
# 
# # if(nameindex %in% t){
# #   file = paste0(cell,'_target_net_bubble.png')
# #   p_t <- netVisual_bubble(cellchat, sources.use = n[!n %in% nameindex], targets.use = nameindex, signaling = c("ANNEXIN"), remove.isolate = TRUE)
# #   hei = length(unique(p_t$data$interaction_name_2))*0.15
# #   hei = ifelse(hei>6,hei,6)
# #   ggsave(file, width = wid, height = hei, dpi = 500, device = "png")
# # }
# # p <<- c(p,cell)#
# # return(p)
# # }
# # })
# 
# # }
# 


# cell_use <- cell_use[cell_use$Var1 %in% c("STM", 'MON', 'MAST','DC', 'M2MAC', 'INMON', 'CD8TEXP', 'NEUT', 'BN', 'PLA','TREG', 
#                                           'CD8TEXINT', 'CD8TEREX','GDT','NK'
#                                           #, 'MVA', 'PVA','MSC.MVA', 'PERI'
# ),]# Liver: C0
# 
# cell_use <- cell_use[cell_use$Var1 %in% c("STM", 'MON', 'CD8TEXP', 'KUP','BN','TREG', 'CD8TEREX'
#                                           #, 'MVA', 'PVA','MSC.MVA', 'PERI'
# ),]# Liver: C0
# 
# # cell_use <- cell_use[cell_use$Var1 %in% c("STM",'BN', 'PLA','DC', 'MON', 'KUP',
# #                                           'CD8TEXP','TREG','PVA','MSC.MVA'),]
# # 
# # cell_use <- cell_use[cell_use$Var1 %in% c("STM", 'DC', 'CD8TEREX', 'pDC',
# #                                           'PFIB','MVA', 'PERI','HSC'),]
# 
# cell_use <- cell_use[cell_use$Var1 %in% c("STM", 'INMON', 'CD8TEXP', 'M1MAC','M2MAC','TREG', 'CD8TEXP','CD8TEXINT', 'CD8TEREX','MDSCs'
#                                           #, 'MVA', 'PVA','MSC.MVA', 'PERI'
# ),]# Esophagus: C10
# 
# cell_use <- cell_use[cell_use$Var1 %in% c("STM", 'LC','cDC','pDC', 'INMON', 'CD8TEXP', 'M2MAC','BN','TREG', 'CD8TEXINT', 'CD8TEREX'
#                                           #, 'MVA', 'PVA','MSC.MVA', 'PERI'
# ),]# Lung: C12
# 
# cell_use <- cell_use[cell_use$Var1 %in% c("STM", 'CD4TN','NEUT', 'GDT', 'CD8TEX', 'M2MAC','MAST'
#                                           #, 'MVA', 'PVA','MSC.MVA', 'PERI'
# ),]# Pancreas: C4

# immune_cell_use <- c("STM", 'pDC','cDC', 'GDT', 'MDSC', 'M2MAC','PLA',"STM", 'MON', 'MAST','DC', 'M2MAC', 'INMON', 'CD8TEXP', 'NEUT', 'BN', 'PLA','TREG', 
#               'CD8TEXINT', 'CD8TEREX','GDT','NK',"STM", 'MON', 'CD8TEXP', 'KUP','BN','TREG', 'CD8TEREX', "STM", 'INMON', 'CD8TEXP', 'M1MAC','M2MAC','TREG', 'CD8TEXP','CD8TEXINT', 'CD8TEREX','MDSCs', "STM", 'LC','cDC','pDC', 'INMON', 'CD8TEXP', 'M2MAC','BN','TREG', 'CD8TEXINT', 'CD8TEREX', "STM", 'CD4TN','NEUT', 'GDT', 'CD8TEX', 'M2MAC','MAST',"STM", 'pDC','cDC', 'GDT', 'MDSC', 'M2MAC','PLA')
# paste0(unique(immune_cell_use), collapse = "', '")



##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) # load SCENIC object in ts860 not UTH
# data_path <- '/data/rluo4/database/'
# typeC <- c('Epi','Imm','Str')
# TF_summary <- NULL
# ntop1=5;ntop2=50
# for (organ in organ_all) {
#   print(organ)
#   inputDir=paste0(data_path,organ,'/Output/Epi/Results/res_fig')#'/data/rluo4/database/Liver/Output/Str/Results'
#   if(organ == 'Chen'){
#     inputDir = gsub('Chen','CRC',inputDir)
#   }
#   setwd(inputDir)
#   tfDir <- list.files()
#   tfDir <- tfDir[grep('TF_',tfDir)]
#   tfDir <- ifelse(organ %in% c('CRC','Chen'),ifelse(organ=='Chen', tfDir[1],tfDir[2]),
#                   tfDir)
#   print(tfDir)
#   inregulons_Type<-NULL
#   TissueType <- gsub('TF_','',tfDir)
#   TissueType
#   # TissueType <- 'AD'
#   # for (TissueType in Types) {
#   outputDir <- paste0(inputDir,'/TF_',TissueType)
#   if (dir.exists(paste0(outputDir))){
#     setwd(outputDir)
#   }
#   inregulons <- read.table('sd_regulon_RSS.list',header = T)
# 
#   ingrnPath <- ifelse(organ=='Chen',
#                       file.path(inputDir,'../SCENIC/Chen_ABN/grn.tsv'), file.path(inputDir,'../SCENIC/ABN/grn.tsv'))
# 
#   grn <- read.table(ingrnPath,sep='\t',header=T,stringsAsFactors=F)
#   inregulons1=gsub('[(+)]','', unique(inregulons$regulon))
#   inregulons$target <- inregulons$regulon
#   inregulons$regulon <- gsub('[(+)]','', (inregulons$regulon))
#   c1 <- which(grn$TF %in% inregulons1)
#   grn <- grn[c1,]
#   ntop <- ntop2
#   for (tf in unique(grn$TF)) {
#     tmp <- subset(grn,TF==tf)
#     if (dim(tmp)[1] > ntop) {
#       tmp <- tmp[order(tmp$importance,decreasing=T),]
#       tmp <- tmp[1:ntop,]
#     }
#     inregulons$target[inregulons$regulon==tf] <- paste(c(tmp$target[1:3],' etc.'),collapse = ',')
#     inregulons$target_t50[inregulons$regulon==tf] <- paste(c(tmp$target),collapse = ',')
#   }
#   inregulons$DiseaseStage <- TissueType
#   inregulons$Tissue <- organ
#   inregulons<- inregulons[,-(4:5)]
#   inregulons_Type <<- rbind(inregulons_Type,inregulons)
#   # }
#   TF_summary <<- rbind(TF_summary,inregulons_Type)
# }
# 
# TF_summary$Index <- 1:nrow(TF_summary)
# colnames(TF_summary)[1:5] <- c('TF','CellType','RSS','TargetGene','Target_Top50')
# save(TF_summary, file = paste0(data_path,'../TF_summary.rds'))
# write.table(TF_summary_All[,c(8,1:7)], '/data/rluo4/All/Output/Index_TF_All.txt',#paste0(TissueType,'_Index_CellChat.txt'),
#             sep='\t',quote=F,row.name=F,col.name=T)

# TF_summary_All <- rbind(TF_summary_Breast,TF_summary_Cervix,TF_summary_CRC, TF_summary_Endometrium,
#                         TF_summary_Esophagus, TF_summary_Liver,TF_summary_Lung,TF_summary_HNSCC,
#                         TF_summary_Prostate, TF_summary_Pancreas, TF_summary_Skin,TF_summary_GC,TF_summary_THCA
# )
# TF_summary_All$Index <- 1:nrow(TF_summary_All)
# colnames(TF_summary_All)[1:5] <- c('TF','CellType','RSS','TargetGene','Target_Top50')
# write.table(TF_summary_All[,c(8,1:7)], '/home/lorihan/lrh/All/Output/Index_TF_All.txt',#paste0(TissueType,'_Index_CellChat.txt'),
#             sep='\t',quote=F,row.name=F,col.name=T)

library(stringr)
setwd('/data/rluo4/All/Output')
TF_RSS <- read.csv('Index_TF_All.txt', sep = '\t')
disease_summary <- read.csv('disease_summary.txt',sep = '\t')
Path <- NULL
table(TF_RSS$DiseaseStage, TF_RSS$Tissue)
# TF_RSS$Tissue[TF_RSS$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
TF_RSS$DiseaseStage <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                                     gsub('Tumor','PRAD', TF_RSS$DiseaseStage)))
table(TF_RSS$DiseaseStage)
write.csv(TF_RSS, file = paste('/data/rluo4/All/Output/TF_RSS.csv'), quote = F, row.names = F) #Table S3
write.xlsx(TF_RSS[,-1], file = paste('/data/rluo4/All/Output/table_xlsx/TableS3.xlsx'),  rowNames = FALSE)


MajorCell <- apply(TF_RSS, 1, function(x){
  celltype <- x['CellType']
  search.1 <- grep(celltype,disease_summary$Epithelial_celltype)
  search.2 <- grep(celltype,disease_summary$Immune_celltype)
  search.3 <- grep(celltype,disease_summary$Stromal_celltype)
  dir <- ifelse(length(search.1)==0, ifelse(length(search.2)!=0 & length(search.3)==0, 
                                            'Imm','Str'),'Epi')
  Path <<- c(Path, dir)
})
sort(table(TF_RSS$TF[ TF_RSS$CellType %in% 'INCAF']))
precancer <- c('AAH','AD','AEH','AK',#'AIS',
               'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
               'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
               'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
               'SER','SIM','WIM')
TF_RSS_Precancer <- TF_RSS[TF_RSS$DiseaseStage %in% precancer & TF_RSS$CellType %in% 'TREG',]
Precancer_TF <- data.frame(sort(table(TF_RSS_Precancer$TF), decreasing = T))

TF_RSS_Cancer <- TF_RSS[ TF_RSS$DiseaseStage %ni% c(precancer,'ADJ','Healthy')  & TF_RSS$CellType %in% 'TREG',]
Cancer_TF <- data.frame(sort(table(TF_RSS_Cancer$TF), decreasing = T))

TF_RSS_ADJ <- TF_RSS[ TF_RSS$DiseaseStage %in% c('ADJ','Healthy')  & TF_RSS$CellType %in% 'TREG',]

TF_RSS_Precancer <- TF_RSS[TF_RSS$DiseaseStage %in% precancer & TF_RSS$CellType %in% c('MYOFIB','HSC','PSC'),]
Precancer_TF <- data.frame(sort(table(TF_RSS_Precancer$TF), decreasing = T))

TF_RSS_Cancer <- TF_RSS[ TF_RSS$DiseaseStage %ni% c(precancer,'ADJ','Healthy')  & TF_RSS$CellType %in% c('MYOFIB','HSC','PSC'),]
Cancer_TF <- data.frame(sort(table(TF_RSS_Cancer$TF), decreasing = T))

TF_RSS_ADJ <- TF_RSS[ TF_RSS$DiseaseStage %in% c('ADJ','Healthy')  & TF_RSS$CellType %in% c('MYOFIB','HSC','PSC'),]

ADJ_TF <- data.frame(sort(table(TF_RSS_ADJ$TF), decreasing = T))

intersect(Cancer_TF$Var1, Precancer_TF$Var1)

P <- setdiff(TF_RSS_Precancer$TF, c(TF_RSS_Cancer$TF, TF_RSS_ADJ$TF) )
C <- setdiff(TF_RSS_Cancer$TF, c(TF_RSS_Precancer$TF, TF_RSS_ADJ$TF) )

View(Precancer_TF[Precancer_TF$Var1 %in% P,])
View(Cancer_TF[Cancer_TF$Var1 %in% C,])
View(TF_RSS[TF_RSS$CellType %in% 'TREG',])



Index_TF_figures <- TF_RSS
Index_TF_figures$Path <- paste0('/data/rluo4/database/',
                                TF_RSS$Tissue,'/', Path, '/TF_',TF_RSS$DiseaseStage,'/regulons_activity_in_dotplot.png')
Index_TF_figures$Path <- gsub('Oral cavity','OralCavity',Index_TF_figures$Path)
# write.table(Index_TF_figures, '/data/rluo4/All/Output/Index_TF_All_figures.txt',#paste0(TissueType,'_Index_CellChat.txt'),
#             sep='\t',quote=F,row.name=F,col.name=T)

###################################################################
# remove duplicated genes (geneSymbol+Tissue)
###################################################################
gene_summary <- Patt
gene_summary$geneSymbol <- paste(gene_summary$Symbol,gene_summary$Organ,sep = ',')
# tmp=by(gene_summary,gene_summary$Symbol,function(x) rownames(x)[which.min(x$Pvalue_Adj)])
tmp=by(gene_summary,gene_summary$geneSymbol,function(x) rownames(x)[which.min(x$Pvalue_Adj)])
probes = as.character(tmp)
length(probes) # 59002
gene_summary=gene_summary[rownames(gene_summary) %in% probes,]
gene_summary$Index <- 1:nrow(gene_summary)
# table(gene_summary$Symbol %in% gene_summary_uniq$Symbol)
# gene_summary$geneID <- gene_summary_uniq$geneID[match(gene_summary$Symbol,gene_summary_uniq$Symbol)] 
# gene_summary$geneAlias <- gene_summary_uniq$geneAlias[match(gene_summary$Symbol,gene_summary_uniq$Symbol)] 
table(is.na(gene_summary$geneID))#3416 --> 1428 --> 1035
# table(Patt$Symbol %in% gene_summary$Symbol)
Na.geneid <- gene_summary[is.na(gene_summary$geneID),]
table(is.na(Patt$geneID))#71215
# Patt$geneID <- gene_summary_uniq$geneID[match(Patt$Symbol,gene_summary_uniq$Symbol)] 
table(is.na(Patt$geneID))#18501 

table(is.na(gene_summary$geneID))#1356-->567
miss.symbol <- gene_summary[is.na(gene_summary$geneID),]
################################
gene_summary <- gene_summary[!is.na(gene_summary$geneID),]
gene_summary$Index <- 1:nrow(gene_summary)
# Patt <- Patt[!is.na(Patt$geneID),]
# Patt$Index <- 1:nrow(Patt)
# library(AnnotationDbi, lib.loc = "/usr/local/lib/R/site-library")
library(org.Hs.eg.db)
eg2symbol=toTable(org.Hs.egSYMBOL)
eg2ensembl=toTable(org.Hs.egENSEMBL)
length(unique(eg2ensembl$ensembl_id))
eg2name=toTable(org.Hs.egGENENAME)
eg2alias=toTable(org.Hs.egALIAS2EG)
###############################################################################
eg2CHR <- toTable(org.Hs.egMAP)
eg2TYPE <- toTable(org.Hs.egGENETYPE)
eg2GO <- toTable(org.Hs.egGO2ALLEGS)
eg2UniProt <- toTable(org.Hs.egUNIPROT)
gene_summary$CHR <- eg2CHR$cytogenetic_location[match(gene_summary$geneID,eg2CHR$gene_id)]
gene_summary$geneAlias <- eg2alias$alias_symbol[match(gene_summary$geneID,eg2alias$gene_id)]
gene_summary$geneName <- eg2name$gene_name[match(gene_summary$geneID,eg2name$gene_id)]
gene_summary$geneType <- eg2TYPE$gene_type[match(gene_summary$geneID,eg2TYPE$gene_id)]
gene_summary$GO <- eg2GO$go_id[match(gene_summary$geneID,eg2GO$gene_id)]
gene_summary$UniProtAcc <- eg2UniProt$uniprot_id[match(gene_summary$geneID,eg2UniProt$gene_id)]
# ###############################################################################
write.table(gene_summary[,-c(3,12)],'gene-summary.txt',sep='\t',quote = F,row.names = F)

intersect(gene_summary$Symbol, Index_TF_figures$TF)
intersect(gene_summary$Symbol, Ligand_Receptor$Ligand_Receptor)
Index_TF_figures <- Index_TF_figures[Index_TF_figures$TF %in% gene_summary$Symbol,]
Index_TF_figures$Index <- 1:nrow(Index_TF_figures)
# write.table(Index_TF_figures, 'Index_TF_figures.txt', sep='\t',quote=F,row.name=F,col.name=T)
tar.all <- paste(unique(Index_TF_figures$Target_Top50),collapse = ',')
tar.all <-  unlist(strsplit(tar.all,','))

data_path = '/data/rluo4/All/Output/'
load(paste0(data_path, './TF_summary.rds'))
unique(intersect(TF_summary$TF, Trans_genes$MP_Gene))# 22 genes

View(Trans_genes[Trans_genes$MP_Gene %in% TF_summary$TF,])
View(TF_summary[TF_summary$TF %in% Trans_genes$MP_Gene,])
sort(table(Trans_genes[Trans_genes$MP_Gene %in% TF_summary$TF,]$MP_Gene))


###==================================================================================
### TF of QSC and PSP
dat = read.table("../TF/myeloid.auc_mtx.tsv",row.names=1,header=T,sep="\t",stringsAsFactors=F) 
dat2 = readRDS("/public/workspace/AML_Leukemia/AML_Leukemia.integrated.rds")

a1 = colnames(dat)[dat$cell_class %in% c("C2")]
a2 = colnames(dat)[dat$cell_class %in% c("C16.pre")]


group1 = dat[which(rownames(dat) %in% a1),]
group2 = dat[which(rownames(dat) %in% a2),]
cm = rbind(group1,group2)
p = apply(cm,2,function(x){wilcox.test(x[1:nrow(group1)],x[(nrow(group1)+1):nrow(cm)])[[3]]})
fc = apply(cm,2,function(x){mean(as.numeric(x[1:nrow(group1)]))/mean(as.numeric(x[(nrow(group1)+1):nrow(cm)]))})
mean.9 = apply(cm,2,function(x){mean(as.numeric(x[1:nrow(group1)]))})
mean.1_2 = apply(cm,2,function(x){mean(as.numeric(x[(nrow(group1)+1):nrow(cm)]))})

x = data.frame(label=names(fc),Mean.9=mean.9,Mean.1_2=mean.1_2, FC=as.numeric(fc), P.Value=as.numeric(p))
x$P.Value[x$P.Value==0]=1.381878e-300
rownames(x)=gsub("\\...","",rownames(x))

DefaultAssay(dat2) = "RNA"

tmp = dat2@assays$RNA@data[intersect(rownames(x),rownames(dat2)),a1]
tmp.n9 = apply(tmp,1,function(x){length(which(x>0))/ncol(tmp)})

tmp = dat2@assays$RNA@data[intersect(rownames(x),rownames(dat2)),a2]
tmp.n1_2 = apply(tmp,1,function(x){length(which(x>0))/ncol(tmp)})

select=intersect(names(which(tmp.n9<0.05)),names(which(tmp.n1_2<0.05)))
select=setdiff(names(tmp.n9),select)
x=x[select,]

logFCcut = 0.25
pvalCut = 1.30103 
logFCcut2 = 0.38
pvalCut2 = 12
logFCcut3=1
pvalCut3=20

n1 = length(x[, 1])
cols = rep("grey", n1)
names(cols)= rownames(x)
cols[-log10(x$P.Value) > pvalCut & log2(x$FC) >logFCcut]= "#9C9C9C"
cols[-log10(x$P.Value) > pvalCut2 & log2(x$FC) > logFCcut2]= "#ED4F4F"
cols[-log10(x$P.Value) > pvalCut & log2(x$FC) < -logFCcut]= "#B2DF8A"
cols[-log10(x$P.Value) > pvalCut2 & log2(x$FC) < -logFCcut2]= "#329E3F"
color_transparent = adjustcolor(cols, alpha.f = 0.5)
x$color_transparent = color_transparent

n1 = length(x[, 1])
size = rep(1, n1)
size[-log10(x$P.Value) > pvalCut & log2(x$FC) > logFCcut]= 2
size[-log10(x$P.Value) > pvalCut2 & log2(x$FC) > logFCcut2]= 4
size[-log10(x$P.Value) > pvalCut3 & log2(x$FC) > logFCcut3]= 6
size[-log10(x$P.Value) > pvalCut & log2(x$FC) < -logFCcut]= 2
size[-log10(x$P.Value) > pvalCut2 & log2(x$FC) < -logFCcut2]= 4
size[-log10(x$P.Value) > pvalCut3 & log2(x$FC) < -logFCcut3]= 6

p1 = ggplot(data=x, aes(log2(FC), -log10(P.Value), label = label)) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  scale_x_continuous(
    breaks = c( -1, -0.38, -logFCcut, 0, 0.38, logFCcut, 1),
    labels = c( -1, -0.38, -logFCcut, 0, 0.38, logFCcut, 1),
    limits = c(-1.5, 1.5) ) +
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey91", linetype="longdash", lwd = 0.5) + 
  geom_hline(yintercept = pvalCut, color="grey91", linetype="longdash", lwd = 0.5) +
  geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey91", linetype="longdash", lwd = 0.5) +
  geom_hline(yintercept = pvalCut2, color="grey91", linetype="longdash", lwd = 0.5)+
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank())


plot = p1 + geom_text_repel(aes(x = log2(FC), y = -log10(P.Value), label = ifelse(log2(FC) > log2(1.3) & -log10(P.Value) > 12 , rownames(x),"")),
                            colour="darkred", size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))+ 
  geom_text_repel(aes(x = log2(FC), y = -log10(P.Value), label = ifelse(log2(FC) < (-log2(1.3)) & -log10(P.Value) > 12 , rownames(x),"")),
                  colour="darkgreen", size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))

pdf("Figure3_B.pdf")
print(plot)
dev.off()


###################
library(stringr)
TF_RSS <- read.csv('/data/rluo4/All/Output/Index_TF_All.txt', sep = '\t')
disease_summary <- read.csv('/data/rluo4/All/Output/disease_summary.txt',sep = '\t')
Path <- NULL
MajorCell <- apply(TF_RSS, 1, function(x){
  celltype <- x['CellType']
  search.1 <- grep(celltype,disease_summary$Epithelial_celltype)
  search.2 <- grep(celltype,disease_summary$Immune_celltype)
  search.3 <- grep(celltype,disease_summary$Stromal_celltype)
  dir <- ifelse(length(search.1)==0, ifelse(length(search.2)!=0 & length(search.3)==0,
                                            'Imm','Str'),'Epi')
  Path <<- c(Path, dir)
})
Index_TF_figures <- TF_RSS
Index_TF_figures$Path <- paste0('/data/rluo4/database/',
                                TF_RSS$Tissue,'/', Path, '/TF_',TF_RSS$DiseaseStage,'/regulons_activity_in_dotplot.png')
Index_TF_figures$Path <- gsub('Oral cavity','OralCavity',Index_TF_figures$Path)
write.table(Index_TF_figures, '/data/rluo4/All/Output/Index_TF_All_figures.txt',#paste0(TissueType,'_Index_CellChat.txt'),
            sep='\t',quote=F,row.name=F,col.name=T)
# 
# intersect(gene_summary$Symbol, Index_TF_figures$TF)
# intersect(gene_summary$Symbol, Ligand_Receptor$Ligand_Receptor)
# Index_TF_figures <- Index_TF_figures[Index_TF_figures$TF %in% gene_summary$Symbol,]
# Index_TF_figures$Index <- 1:nrow(Index_TF_figures)
# write.table(Index_TF_figures, 'Index_TF_figures.txt', sep='\t',quote=F,row.name=F,col.name=T)
# tar.all <- paste(unique(Index_TF_figures$Target_Top50),collapse = ',')
# tar.all <-  unlist(strsplit(tar.all,','))

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) ##################---summarize the cell subtypes and the percentage by different stages---###############################
# results from database-sum.R: ts860
load('/data/rluo4/All/Output/pheno_all.rds')
load('/data/rluo4/All/Output/Imm_all.rds')
load('/data/rluo4/All/Output/Str_all.rds')

colnames(pheno_all$CRC_Epi)
colnames(Imm_all$CRC_Imm)[1]
colN <- c('barcode','orig.ident','batch', 'sample_name','Tissue','Cell_Type','SimplifiedSampleName')
colN

library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
data_path = '/data/rluo4/All/Output'
setwd(data_path)
organ_all <- list.files('Cluster/')
organ_all <- organ_all[organ_all!='Chen']
organ_all
cell_all <- c(pheno_all, Imm_all, Str_all)
pheno <- paste(organ_all, 'Epi', sep = '_')
All_ratio <- NULL
cell_summary <- NULL
for (obs in pheno) {
  Organ <- str_split(obs,'_', simplify = T)[,1]
  print(Organ)
  epi <- paste(Organ, 'Epi', sep = "_")
  str <- paste(Organ, 'Str', sep = "_")
  imm <- paste(Organ, 'Imm', sep = "_")
  
  Epi <- cell_all[[epi]]
  Epi$barcode = rownames(Epi)
  Epi$cell_type <- Epi$Cell_Type
  table(Epi$Organ,Epi$Tissue)
  Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
  Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                              gsub('Tumor','PRAD', Epi$Tissue)))
  Epi$Major_type <- 'Epithelial'
  print(table(Epi$Organ))
  Imm <- cell_all[[imm]]
  Imm$Major_type <- 'Immune'
  
  Str <- cell_all[[str]]
  Str$Major_type <- 'Stromal'
  
  common_columns <- c('barcode','orig.ident','batch','SimplifiedSampleName', 'sample_name','Tissue','Cell_Type','Major_type')#intersect(intersect(colnames(Epi),colnames(Imm)),colnames(Str))
  phe = dplyr::bind_rows(Epi, Imm, Str) %>%
    select(all_of(common_columns))
  # cohort_directory <- paste0(data_path,Organ,'/Output/Epi/')
  # setwd(cohort_directory)
  # setwd("/data/rluo4/All/Output/")
  # Cyto_dir <- paste0("/data/rluo4/All/Output/Cellratio/",Organ)
  # 
  # if (!dir.exists(paste0(Cyto_dir))){
  #   dir.create(paste0(Cyto_dir))
  # }
  # setwd(Cyto_dir)
  phe <- phe[! phe$sample_name %in% c('Cold'),]
  phe <- phe[! phe$orig.ident %in% c("p3", "p4", "p5"),]
  phe$Organ <- unique(Epi$Organ)
  table(phe$Tissue)
  phe <- phe[phe$Tissue!='UNC',]
  cell_summary <<- rbind(cell_summary, phe)
  # Cellratio <- prop.table(table(Cellratio$Cell_Type, Cellratio$Tissue), margin = 2)#计算各组样本不同细胞群比例
  # Cellratio <- as.data.frame(Cellratio)
  # colnames(Cellratio) <- c("CellType","Tissue","CellFraction") #对列名重命名
  Cellratio <- phe
  Cellratio <- prop.table(table(Cellratio$Cell_Type, Cellratio$orig.ident), margin = 2)#计算各组样本不同细胞群比例
  Cellratio <- as.data.frame(Cellratio)
  colnames(Cellratio) <- c("CellType","Sample","CellFraction") #对列名重命名
  Cellratio$Tissue <- unique(Epi$Organ)#phe$Organ[match(Cellratio$Sample,phe$orig.ident)]
  Cellratio$DiseaseStage <- phe$Tissue[match(Cellratio$Sample,phe$orig.ident)]
  table(Cellratio$DiseaseStage)
  # Cellratio$Major_type <- phe$Major_type[match(Cellratio$Sample,phe$orig.ident)]
  levels(factor(Cellratio$DiseaseStage))
  # Cellratio$DiseaseStage<- factor(Cellratio$DiseaseStage, ordered=T, levels = c('Healthy','ADJ','Precancer','GC'))  #调整画图的x轴坐标顺序
  All_ratio <<- rbind(All_ratio, Cellratio)
  # Cellratio$Cohort <- Cellratio$sample_name#[match(rownames(score_result),colnames(colon_full))]
  # head(Cellratio)
}
length(unique(All_ratio$Sample))
length(unique(cell_summary$orig.ident)) #1332 samples
datasets <- unique(cell_summary$sample_name)
setdiff(datasets,unique(tissue_summary$sample_name))
table(All_ratio$DiseaseStage)

index = cell_summary$Major_type=='Immune'
lym <- unique(cell_summary$Cell_Type[grepl('CD',cell_summary$Cell_Type) & index])
Tlym <- c(lym, 'TFH','TH17','TREG','TH1','TH2','MAIT',# For the CD8+ metacluster c16, nearly half of the cells harbored the semi-invariant TCR α chains of mucosal-associated invariant T cells (MAIT) 
          'ILC', 'NK', 'GDT', 'NKT')#Innate lymphoid cells(ILC): ILC1,2,3 and NK cells
#Innate-like lymphocytes ( ILLs): ILLs 主要包括 NKT 细胞, gamma-delta (γδ) T细胞, B1细胞, 其表面抗原识别受体( TCR 或 BCR) 由胚系基因直接编码产生, 为有限多样性抗原识别受体
Tlym
T_ratio <- All_ratio[ All_ratio$CellType %in% Tlym,]
unique(T_ratio$CellType)
table(cell_summary$Cell_Type[! cell_summary$Cell_Type %in% Tlym])

lym <- unique(cell_summary$Cell_Type[grepl('B',cell_summary$Cell_Type) & index])
lym
Blym <- c(lym, 'GC','PLA')
Blym
B_ratio <- All_ratio[ All_ratio$CellType %in% Blym,]
unique(B_ratio$CellType)

Mye <- unique(cell_summary$Cell_Type[index])
Mye <- setdiff(Mye, c(Tlym, Blym))
Mye
M_ratio <- All_ratio[ All_ratio$CellType %in% Mye,]

index = cell_summary$Major_type=='Stromal'
table(cell_summary$Cell_Type[index])
Mes <- unique(cell_summary$Cell_Type[grepl("CAF|FIB", cell_summary$Cell_Type) & index])#mesenchymal cells(间充质细胞)
Mes <- c(Mes, 'BSM', 'ECM','HSC', 'PSC', 'PERI','SMC')
Mes
End <- unique(cell_summary$Cell_Type[grepl("END|PVA|MVA|SEC", cell_summary$Cell_Type) & index])
End
Others <- unique(cell_summary$Cell_Type[index])
Others <- setdiff(Others, c(Mes, End))
Others

index = cell_summary$Major_type=='Epithelial'
table(cell_summary$Cell_Type[index])
Epi <- unique(cell_summary$Cell_Type[index])
table(cell_summary$Cell_Type[! cell_summary$Cell_Type %in% Epi])

Lineage <- data.frame(Minor_type = c(Epi, Tlym, Blym, Mye, Mes, End, Others),
                      Major_type = c(rep('Epithelia', length(Epi)),
                                     
                                     rep('Tcell', length(Tlym)),
                                     rep('Bcell', length(Blym)), 
                                     rep('Myeloid', length(Mye)),
                                     
                                     rep('Mesenchymal', length(Mes)), 
                                     rep('Endothelia', length(End)),
                                     rep('Others', length(Others))
                      )
)
# types <- unique(Lineage$Major_type)
# types
# for(i in types){
#   print(i)
#   df <-  cell_summary[cell_summary$Cell_Type %in% Lineage$Minor_type[Lineage$Major_type==i], ]
Lineage$Lineage <- Lineage$Major_type
Lineage$Lineage[Lineage$Major_type %in% c('Tcell','Bcell')] <- 'Lymphoid' 
write.csv(Lineage, file='organ13-celltypes.csv')

types <- unique(Lineage$Lineage)
types
library(ArchR)
for(i in types){
  print(i)
  df <-  cell_summary[cell_summary$Cell_Type %in% Lineage$Minor_type[Lineage$Lineage==i], ]
  precancer <- c('AAH','AD','AEH','AK',#'AIS',
                 'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
                 'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
                 'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
                 'SER','SIM','WIM')
  table(df$Tissue)
  df$Stage <- df$Tissue
  df$Stage[df$Tissue %in% precancer] <- 'Precancer'
  df$Stage[! df$Tissue %in% c(precancer,'Healthy','ADJ')] <- 'Cancer'
  print(table(df$Stage))
  # bar.df$Stage <- factor(bar.df$Stage, ordered=T, levels = c('Healthy','ADJ','Precancer','Cancer'))  #调整画图的x轴坐标顺序
  for(stage in unique(df$Stage)){
    bar.df <- df[df$Stage %in% stage,]
    table(bar.df$Organ)
    # bar.df <- mutate(bar.df,name=factor(bar.df$Organ))#
    color_cluster = paletteDiscrete(values = unique(df$Cell_Type))
    # names(color_cluster)=c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
    title <- paste0(stage)
    p <- ggplot(bar.df,aes(x=Organ))+
      geom_bar(aes(fill=Cell_Type),position = "fill",width = .7)+
      scale_x_discrete("")+
      scale_y_continuous(title,expand = c(0,0),labels = scales::label_percent(),position = "right")+
      scale_fill_manual("Cell Types",values = color_cluster)+
      theme_ArchR(baseSize = 10) + #ggtitle(TissueType) +
      theme(legend.position='right') +# theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
      theme(
        plot.title    = element_text(color = 'black', size   = 18, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 18,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 18,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 20, angle = 0),
        axis.text.y   = element_text(color = 'black', size = 20, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 30, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 20, angle = 90),
        legend.title  = element_text(color = 'black', size  = 18),
        legend.text   = element_text(color = 'black', size   = 18),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
        # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
      ) +   coord_flip() #让条形图横过来
    f = paste0(data_path,'/res_fig/Cellratio/', stage,'_',i,'_fraction','.png')
    f
    ggsave(plot=p,height=10,width=12, filename=f, dpi = 500, device = "png")
    
  }
}
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)
library(paletteer)
library(ggthemes)
data_path = '/data/rluo4/All/Output'
setwd(data_path)
LF <- All_ratio[All_ratio$CellType %in% Lineage$Minor_type[Lineage$Major_type=='Epithelia'],]
View(data.frame(table(LF$CellType)))
# organ = 'Colorectum'
celltype <- 'STM'
celltype <- 'GOB'
celltype <- 'BAS'
celltype <- 'EE'
celltype <- 'CILIA'
CellTypes <- c('STM','GOB','BAS','EE','CILIA','IMENT','ABS','ASC','PMC','GMC','HEP','CHO','CHIEF')
color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
color_cluster <- color_cluster[1:3]
names(color_cluster)= c('Healthy','Precancer','Cancer')#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
# for(organ in unique(All_ratio$Tissue)) {
# for(celltype in unique(All_ratio$CellType)) {
for(celltype in CellTypes) {
  # LF <- All_ratio[All_ratio$Tissue==organ & All_ratio$CellType %in% Lineage$Minor_type[Lineage$Major_type=='Epithelia'],]
  LF <- All_ratio#[All_ratio$CellType %in% Lineage$Minor_type[Lineage$Major_type=='Epithelia'],]
  LF$CellFraction <- LF$CellFraction*100
  # if(length(unique(LF$CellType %in% celltype))>0){
  LF <- LF[ LF$CellType==celltype,] 
  LF$Stage <- LF$DiseaseStage
  LF$DiseaseStage[LF$Stage %in% precancer] <- 'Precancer'
  LF$DiseaseStage[! LF$Stage %in% c(precancer,'Healthy','ADJ')] <- 'Cancer'
  table(LF$DiseaseStage)
  unique(LF$Sample)
  # LF_sum <- aggregate(CellFraction~Sample,data=LF ,FUN="sum")
  # LF$LF <- LF_sum$CellFraction[match(LF$Sample,LF_sum$Sample)]
  # LF <- LF[! duplicated(LF$Sample),]
  # LF$DiseaseStage <- paste(LF$organ, LF$DiseaseStage, sep = ':')
  # LF_sum <- aggregate(LF~DiseaseStage,data=LF ,FUN="median")
  # ylab <- "Stem-like cell fraction (%)"
  ylab <- paste0(celltype, "  fraction (%)")
  # ylab <- 'Goblet cell fraction (%)'
  # file = paste0(data_path,'/res_fig/Fraction_Compare/',  celltype,'_fraction_in_', organ, '.png')
  file = paste0(data_path,'/res_fig/Fraction_Compare/',  celltype,'_fraction.png')
  # ggsave(plot = p, filename = file, width = length(unique(LF$DiseaseStage)),height = 6, dpi = 500, device = "png")
  # LF <- LF[! LF$DiseaseStage %in%  'ADJ',]
  LF$DiseaseStage <- gsub('ADJ','Healthy',LF$DiseaseStage)
  LF$DiseaseStage<- factor(LF$DiseaseStage, ordered=T, levels = c('Healthy','Precancer','Cancer'))  #调整画图的x轴坐标顺序
  
  label.pos <- round(max(LF$CellFraction)[5],1)+0.05
  my_comparisons <- list(c("Healthy","Precancer"), c("Healthy", "Cancer"),
                         c("Precancer", "Cancer"))
  # ggerrorplot(LF, x = "DiseaseStage", y = "CellFraction",
  # desc_stat = "mean_sd", color = "black",fill = "DiseaseStage",
  # size = 0.5, 
  # add = "violin",#, add.params = list(color = "darkgray")
  # palette = color_cluster, 
  # position = position_dodge(0.3)  
  # )   + guides(fill = guide_legend(title = 'Disease Stage')) 
  
  p <- ggboxplot(LF, notch="TRUE",
                 "DiseaseStage", "CellFraction",fill = "DiseaseStage",
                 palette = color_cluster, 
                 size = 0.5, 
                 add = "")+xlab("")+ylab(ylab)+ 
    guides(fill = guide_legend(title = 'Disease Stage')) +
    theme(legend.position = "none",
          plot.title    = element_text(color = 'black', size   = 15, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 15,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
          axis.text.x   = element_blank(),#element_text(color = 'black', size = 15, angle = 270),#element_blank(),#
          axis.text.y   = element_text(color = 'black', size = 15, angle = 0),
          axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
          axis.title.y  = element_text(color = 'black', size = 18, angle = 90),
          #legend.title=element_blank(),
          legend.title  = element_text(color = 'black', size  = 14),
          legend.text   = element_text(color = 'black', size   = 14),
          axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
          axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
          panel.background = element_rect(fill='transparent')#, legend.position='right'
    )
  p1 <- p+ stat_compare_means(comparisons = my_comparisons)+ 
    stat_compare_means(label.y = label.pos, label.x = 1.5) # 添加全局p值#stat_compare_means(aes(group = DiseaseStage), label = "p.format")#,label.y = 0.6,label.x = 2)
  p1
  
  ggsave(plot = p1, filename = file, width = 3.5,height = 3.5, dpi = 500, device = "png")
  # }
}

compare_means(CellFraction~DiseaseStage, data=LF, method="wilcox.test",p.adjust.method="BH")
# ? INCAF, ICAF, SMC, PERI
CellTypes <- c('M1MAC','M2MAC','MDSC','INMON','CD4TN','TFH','TH1','TH2','TH17','GC',
               'TREG','CD8TEREX','PLA','NK','NEUT','ICAF','SMC','INCAF','PERI','MSC.MVA','LYMEND','HSC','PSC', 'MYOFIB')
for(celltype in CellTypes) {
  # celltype = 'GC'#'INCAF','PERI'#'NEUT' : use ADJ as healthy
  # LF <- All_ratio[All_ratio$Tissue==organ & All_ratio$CellType %in% Lineage$Minor_type[Lineage$Major_type=='Epithelia'],]
  LF <- All_ratio#[All_ratio$CellType %in% Lineage$Minor_type[Lineage$Major_type=='Epithelia'],]
  # if(length(unique(LF$CellType %in% celltype))>0){
  LF$CellFraction <- LF$CellFraction*100
  LF <- LF[ LF$CellType==celltype,] 
  LF$Stage <- LF$DiseaseStage
  LF$DiseaseStage[LF$Stage %in% precancer] <- 'Precancer'
  LF$DiseaseStage[! LF$Stage %in% c(precancer,'Healthy','ADJ')] <- 'Cancer'
  table(LF$DiseaseStage)
  unique(LF$Sample)
  # LF_sum <- aggregate(CellFraction~Sample,data=LF ,FUN="sum")
  # LF$LF <- LF_sum$CellFraction[match(LF$Sample,LF_sum$Sample)]
  # LF <- LF[! duplicated(LF$Sample),]
  # LF$DiseaseStage <- paste(LF$organ, LF$DiseaseStage, sep = ':')
  # LF_sum <- aggregate(LF~DiseaseStage,data=LF ,FUN="median")
  # ylab <- "Stem-like cell fraction (%)"
  ylab <- paste0(celltype, " fraction (%)")
  # ylab <- 'Goblet cell fraction (%)'
  file = paste0(data_path,'/res_fig/Fraction_Compare/',  celltype,'_fraction.png')
  # ggsave(plot = p, filename = file, width = length(unique(LF$DiseaseStage)),height = 6, dpi = 500, device = "png")
  if(celltype %in% c('INCAF','PERI','NEUT','MYOFIB')){
    LF <- LF[! LF$DiseaseStage %in%  'Healthy',]
    LF$DiseaseStage <- gsub('ADJ','Healthy',LF$DiseaseStage)
  } else{
    LF <- LF[! LF$DiseaseStage %in%  'ADJ',]
  }
  LF$DiseaseStage<- factor(LF$DiseaseStage, ordered=T, levels = c('Healthy','Precancer','Cancer'))  #调整画图的x轴坐标顺序
  # color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  # color_cluster <- color_cluster[1:length(unique(LF$DiseaseStage))]
  # names(color_cluster)=unique(LF$DiseaseStage)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
  # Healthy    Cancer Precancer 
  # "#f89e81" "#99a9cc" "#dd9bc5" 
  label.pos <-  ifelse(celltype %in% c('INMON','TREG','M2MAC','PSC'),round(summary(LF$CellFraction)[5],3)*1.5,
                       ifelse(celltype %in% c('TFH','SMC'),round(summary(LF$CellFraction)[5],3)*12, round(summary(LF$CellFraction)[5],3)*3.5))
  # label.pos <- mean(LF$CellFraction) + summary(LF$CellFraction)[5] # for PERY
  if(celltype %in% c('HSC')){
    label.pos <- mean(LF$CellFraction)*2 + summary(LF$CellFraction)[5] # for PERY
  }
  if(celltype %in% c('NEUT','GC')){
    label.pos <- mean(LF$CellFraction)*1.5 + summary(LF$CellFraction)[5] # for PERY
  }
  if(celltype %in% 'PERI'){
    label.pos <- mean(LF$CellFraction)*1.2 + summary(LF$CellFraction)[5] # for PERY
  }
  p = ggbarplot(LF, x = "DiseaseStage", y = "CellFraction",add = c("mean_se"),#, "jitter"),
                fill = "DiseaseStage",  palette = color_cluster, 
                size = 0.5,           #add = c("mean_se", "point"),
                alpha = 0.5,
                position = position_dodge(0.2))+
    #stat_compare_means() +    # Global p-value
    stat_compare_means(ref.group = 'Healthy', label = "p.format",label.y = label.pos)+
    # stat_compare_means(comparisons = my_comparisons)+ 
    # stat_compare_means(label.y = label.pos, label.x = 1.5) + # 添加全局p值#stat_compare_means(aes(group = DiseaseStage), label = "p.format")#,label.y = 0.6,label.x = 2)
    xlab("")+ylab(ylab)+                  # compare to ref.group
    guides(fill = guide_legend(title = 'Disease Stage')) +
    theme(legend.position = "none",
          plot.title    = element_text(color = 'black', size   = 15, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 15,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
          axis.text.x   = element_blank(),#element_text(color = 'black', size = 15, angle = 270),#element_blank(),#
          axis.text.y   = element_text(color = 'black', size = 15, angle = 0),
          axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
          axis.title.y  = element_text(color = 'black', size = 18, angle = 90),
          #legend.title=element_blank(),
          legend.title  = element_text(color = 'black', size  = 14),
          legend.text   = element_text(color = 'black', size   = 14),
          axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
          axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
          panel.background = element_rect(fill='transparent')#, legend.position='right'
    )
  p
  ggsave(plot = p, filename = file, width = 4,height = 4, dpi = 500, device = "png")
} 


celltype <- 'PSC'
celltype <- 'HSC'
celltype <- 'PERI'
# LF <- All_ratio[All_ratio$Tissue==organ & All_ratio$CellType %in% Lineage$Minor_type[Lineage$Major_type=='Epithelia'],]
LF <- All_ratio#[All_ratio$CellType %in% Lineage$Minor_type[Lineage$Major_type=='Epithelia'],]
LF <- All_ratio[All_ratio$Tissue=='Liver',]
LF <- LF[! LF$DiseaseStage %in%  c('Cyst','Healthy'),]
# if(length(unique(LF$CellType %in% celltype))>0){
LF <- LF[ LF$CellType==celltype,] 
LF$Stage <- LF$DiseaseStage
table(LF$DiseaseStage)
unique(LF$Sample)
# LF <- LF[LF$DiseaseStage!='Cyst',]
ylab <- paste0(celltype, " cell fraction (%)")
# ylab <- 'Goblet cell fraction (%)'
file = paste0(data_path,'/res_fig/Fraction_Compare/',  celltype,'_fraction.png')
if(celltype %in% c('INCAF','PERI','NEUT')){
  LF <- LF[! LF$DiseaseStage %in%  'Healthy',]
  LF$DiseaseStage <- gsub('ADJ','Healthy',LF$DiseaseStage)
} else{
  LF <- LF[! LF$DiseaseStage %in%  'ADJ',]
}
LF$DiseaseStage<- factor(LF$DiseaseStage, ordered=T, levels = c('Healthy','NAFLD','Cirrhotic','HCC'))  #调整画图的x轴坐标顺序
LF$DiseaseStage<- factor(LF$DiseaseStage, ordered=T, levels = c('Healthy','PanIN','PDAC'))  #调整画图的x轴坐标顺序
color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
color_cluster <- color_cluster[1:length(unique(LF$DiseaseStage))]
names(color_cluster)=unique(LF$DiseaseStage)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
# Healthy    Cancer Precancer 
# "#f89e81" "#99a9cc" "#dd9bc5" 
label.pos <-  ifelse(celltype %in% c('INMON','TREG','M2MAC'),round(summary(LF$CellFraction)[5],3)*1.5,
                     ifelse(celltype %in% c('GC','TFH','NEUT','SMC'),round(summary(LF$CellFraction)[5],3)*12, round(summary(LF$CellFraction)[5],3)*3))
compar <- split(t(combn(unique(LF$DiseaseStage), 2)), seq(nrow(t(combn(unique(LF$DiseaseStage), 2)))))

ggline(LF, x = "DiseaseStage", y = "CellFraction",
       combine = TRUE,
       ylab = "Expression", 
       color = "gray",                                     # Line color
       add = c("mean_sd", "violin"),                     
       add.params = list(color = "DiseaseStage"),
       palette = color_cluster#palette = "jco"
) + stat_compare_means(comparisons = compar)+
  stat_compare_means(label.y = label.pos, label.x = 1.5) + # 添加全局p值#stat_compare_means(aes(group = DiseaseStage), label = "p.format")#,label.y = 0.6,label.x = 2)
  xlab("")+ylab(ylab)+                  # compare to ref.group
  guides(fill = guide_legend(title = 'Disease Stage')) +
  theme(legend.position = "top",
        plot.title    = element_text(color = 'black', size   = 15, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 15,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
        axis.text.x   = element_blank(),#element_text(color = 'black', size = 15, angle = 270),#element_blank(),#
        axis.text.y   = element_text(color = 'black', size = 15, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 18, angle = 90),
        #legend.title=element_blank(),
        legend.title  = element_text(color = 'black', size  = 14),
        legend.text   = element_text(color = 'black', size   = 14),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
        panel.background = element_rect(fill='transparent')#, legend.position='right'
  )


Liver <- '/data/rluo4/All/Output/res_fig/Liver'
LF <- All_ratio[All_ratio$Tissue=='Liver',]
table(LF$DiseaseStage[LF$CellType %in% cell_summary$Cell_Type[cell_summary$Major_type=='Immune']])
# CellTypes <- c('M1MAC','M2MAC','MDSC','INMON','CD4TN','TFH','TH1','TH2','TH17','GC',
#                'TREG','CD8TEREX','PLA','NK','NEUT','ICAF','SMC','INCAF','PERI','MSC.MVA','LYMEND','HSC','PSC')
CellTypes <- unique(LF$CellType)
for(celltype in CellTypes) {
  celltype = 'PFIB'#'PERI'#'INCAF','PERI'#'NEUT' : use ADJ as healthy
  LF <- All_ratio[All_ratio$Tissue=='Liver',]
  LF <- LF[! LF$DiseaseStage %in%  c('Cyst'),]
  # if(length(unique(LF$CellType %in% celltype))>0){
  LF <- LF[ LF$CellType==celltype,] 
  LF$Stage <- LF$DiseaseStage
  # LF$DiseaseStage[LF$Stage %in% precancer] <- 'Precancer'
  # LF$DiseaseStage[! LF$Stage %in% c(precancer,'Healthy','ADJ')] <- 'Cancer'
  table(LF$DiseaseStage)
  unique(LF$Sample)
  ylab <- paste0(celltype, " cell fraction (%)")
  # ylab <- 'Goblet cell fraction (%)'
  file = paste0(Liver,'/Fraction_Compare/',  celltype,'_fraction.png')
  # ggsave(plot = p, filename = file, width = length(unique(LF$DiseaseStage)),height = 6, dpi = 500, device = "png")
  if(celltype %in% c('INCAF','PERI','NEUT','MYOFIB')){
    LF <- LF[! LF$DiseaseStage %in%  'Healthy',]
    LF$DiseaseStage <- gsub('ADJ','Healthy',LF$DiseaseStage)
  } else{
    LF <- LF[! LF$DiseaseStage %in%  'ADJ',]
  }
  LF$DiseaseStage<- factor(LF$DiseaseStage, ordered=T, levels = c('Healthy','NAFLD','Cirrhotic','HCC'))  #调整画图的x轴坐标顺序
  # LF$DiseaseStage<- factor(LF$DiseaseStage, ordered=T, levels = c('Cirrhotic','HCC'))  #调整画图的x轴坐标顺序
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  color_cluster <- color_cluster[1:length(unique(LF$DiseaseStage))]
  names(color_cluster)=unique(LF$DiseaseStage)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
  color_cluster
  label.pos <- mean(LF$CellFraction)*2 + summary(LF$CellFraction)[5] # for PERY
  maxgroup <- LF$Stage[which.max(LF$CellFraction)]
  index = LF$Stage==maxgroup
  print(summary(LF$CellFraction[index]))
  # label.pos <- max(LF$CellFraction[index])#
  label.pos <- mean(LF$CellFraction[index])*1.5 #+ se(LF$CellFraction[index])/2 # + 0.005#+ summary(LF$CellFraction[index])[5] # for PERY
  label.pos <- mean(LF$CellFraction[index])*2 # for PFIB
  p = ggbarplot(LF, x = "DiseaseStage", y = "CellFraction",add = c("mean_se"),#, "jitter"),
                fill = "DiseaseStage",  palette = color_cluster, 
                size = 0.5,           #add = c("mean_se", "point"),
                alpha = 0.5,
                position = position_dodge(0.2))+
    #stat_compare_means() +    # Global p-value
    stat_compare_means(ref.group = 'Healthy', label = "p.format",label.y = label.pos)+
    xlab("")+ylab(ylab)+                  # compare to ref.group
    guides(fill = guide_legend(title = 'Disease Stage')) +
    theme(legend.position = "top",
          plot.title    = element_text(color = 'black', size   = 15, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 15,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
          axis.text.x   = element_blank(),#element_text(color = 'black', size = 15, angle = 270),#element_blank(),#
          axis.text.y   = element_text(color = 'black', size = 15, angle = 0),
          axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
          axis.title.y  = element_text(color = 'black', size = 18, angle = 90),
          #legend.title=element_blank(),
          legend.title  = element_text(color = 'black', size  = 12),
          legend.text   = element_text(color = 'black', size   = 12),
          axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
          axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
          panel.background = element_rect(fill='transparent')#, legend.position='right'
    )
  p
  ggsave(plot = p, filename = file, width = 5.2, height = 4.5, dpi = 500, device = "png")
} 

# library(ggstatsplot, lib.loc = "/data/rluo4/lorihan/R/site-library")
# ggstatsplot::ggbetweenstats(
#   data = LF,
#   x = DiseaseStage,
#   y = CellFraction,
#   p.adjust.method = 'bonferroni',
#   notch = TRUE, # show notched box plot
#   mean.plotting = FALSE, # whether mean for each group is to be displayed
#   mean.ci = TRUE, # whether to display confidence interval for means
#   mean.label.size = 2.55, # size of the label for mean
#   type = "np", # which type of test is to be run
#   k = 2, # number of decimal places for statistical results
#   outlier.tagging = FALSE, # whether outliers need to be tagged
#   xlab = "Disease Stage", # label for the x-axis variable
#   ylab = "CIBERSORTx-inferred cell fractions (%)", # label for the y-axis variable
#   title = "C16.pre in TCGA LIHC cohort", # title text for the plot
#   package = "palettetown", #View(paletteer::palettes_d_names)
#   palette = "blastoise",# package= "ggsci",#提取调色板所需的包
#   messages = FALSE
# ) +  theme(legend.position='none')+
#   theme( #panel.border=element_rect(fill='transparent'),
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
# ggsave(height=6.5,width=8.5, filename=file, dpi = 500, device = "png")
# ggplot(LF, aes(DiseaseStage, CellFraction))+
#   #stat_summary(geom = "errorbar", fun.data = 'mean_sd', width = 0.3)+#误差棒
#   stat_summary(fun.y=mean, fun.args = list(mult=1),geom='bar',colour="black",fill=color_cluster,width=.7) +
#   stat_summary(fun.data = mean_sdl,fun.args = list(mult=1), geom='errorbar', color='black',width=.2) +
#   geom_jitter(aes(fill = DiseaseStage),position = position_jitter(0.2),shape=21, size = 0.35,alpha=0.9)+
#   scale_fill_manual(values=color_cluster)+
#   xlab('')+ylab(ylab)+ggtitle('')+
#   geom_signif(comparisons = my_comparisons, # 设置要对比的组
#               y_position = 1.6,#c(34,36,38), #设置3个显著性标记的高度
#               tip_length = c(0), #设置显著性那条横线两头向下的长度
#               map_signif_level = F,#T, #设置是否标记显著性的*号，还是直接标记数值
#               test =wilcox.test #设置显著性计算方式
#   ) +
#   theme_classic()+
#   theme(legend.position='none')+
#   theme(
#     plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
#     plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
#     plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#     axis.text.x   = element_text(color = 'black', size = 14, angle = 0),
#     axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
#     axis.title.x  = element_text(color = 'black', size = 14, angle = 0),
#     axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
#     legend.title  = element_text(color = 'black', size  = 16),
#     legend.text   = element_text(color = 'black', size   = 16),
#     axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#     axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#     #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
#   )

################################---Cell Dynamics---####################################################################################
data_path = '/data/rluo4/All/Output'
setwd(data_path)
organ_all <- list.files('Cluster/')
for (organ in organ_all) {
  print(organ)
  Organ <- gsub('THCA','Thyroid',gsub('HNSCC','Oral cavity',gsub('Chen','Colorectum',
                                                                 gsub('GC','Stomach', gsub('CRC','Colorectum',organ)))))
  
  LF <- cell_summary[cell_summary$Organ == Organ,]#[cell_summary$Major_type=='Epithelial',]#[! sdata.obs$Cell_Type %in% c('FIB','END'),]   
  #? A002-C-106 F007 ?
  # LF$orig.ident <- gsub("CRC1_8810",'CRC-1-8810',LF$orig.ident)
  # LF$orig.ident <- gsub("CRC3_11773",'CRC-3-11773',LF$orig.ident)
  # Count2Ratio<-function(x){
  LF<-(t(table(LF[c('SimplifiedSampleName', 'Cell_Type')])))
  LF<-as.data.frame(LF)
  LF<-dcast(LF, Cell_Type~SimplifiedSampleName)
  # write.table(LF, file='./count.xls', sep='\t', quote=F, row.names=F)
  for(col in colnames(LF)){
    if(col=='Cell_Type'){
      next
    }
    print(col)
    LF[col]<-prop.table(LF[col])
  }
  rownames(LF) <- LF$Cell_Type
  indir <- '/data/rluo4/All/Output/sdata_ABN/'
  pc_file = paste0(indir, organ,'_DEG.RData')
  if(organ !='Chen'){
    load(pc_file)
    if(organ=='CRC'){
      pc_df <-  pc_df_Becker
      Patt <- Patt_Becker
    }
  } else{
    pc_file <- gsub('Chen','CRC',pc_file)
    load(pc_file)
    pc_df <- pc_df_Chen
    Patt <- Patt_Chen
  }
  unique(pc_df$Sample)
  table(pc_df$Tissue)
  head(Patt)
  # table(pc_df$SimplifiedSampleName %in% colnames(TME))
  # STM_df <- pc_df[pc_df$SimplifiedSampleName %in% colnames(TME),]
  cell_LF <- data.frame(t(LF[,na.omit(match(pc_df$Sample, colnames(LF)))]*100))
  table(rownames(cell_LF) %in%  pc_df$Sample)
  cell_LF$Sample <- rownames(cell_LF)
  STM_df <- left_join(pc_df[,c("Sample",'nearest_spline_x_vals','Tissue','sample_name')], cell_LF,by='Sample')
  STM_df$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                                 gsub('Tumor','PRAD', STM_df$Tissue)))
  head(STM_df)
  colnames(STM_df)
  summary(STM_df$MAIT)
  STM_df[is.na(STM_df$BMEM),]
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  color_cluster <- color_cluster[1:length(unique(STM_df$Tissue))]
  names(color_cluster)=unique(STM_df$Tissue)#
  table(STM_df$Tissue)
  ######
  Sum_dir <- '/data/rluo4/All/Output/res_fig/Fraction_Dynamics/'
  if (!dir.exists(Sum_dir)){
    dir.create(paste0(Sum_dir))
  }
  setwd(Sum_dir)
  cell_folder <- paste0(Sum_dir,organ)
  if (!dir.exists(cell_folder)){
    dir.create(paste0(cell_folder))
  }
  setwd(cell_folder)
  getwd()
  # dim(test.diff)#4775 33
  n <- colnames(STM_df)[-(1:4)]#rownames(test.diff)
  p <- NULL
  plot <-  lapply(setNames(n, n), function(nameindex) {
    cell <- nameindex
    # cell_p <- pc_LF[match(colnames(test.diff),rownames(pc_LF)),]
    # cell_p$cell <- t(test.diff[cell,])
    ylab <- paste0('Cell fractions', ' (%)')
    p = ggplot(STM_df, aes(x=nearest_spline_x_vals, y= STM_df[,cell], color=Tissue)) +
      geom_point(size=4) + scale_color_manual(values=color_cluster) + ggtitle(cell) + 
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
    f <- paste0(cell,"_dynamics.png")
    ggsave(plot=p,height=4,width=6, filename=f, dpi = 300, device = "png")
    return(p)
  })
}

Index_cell_figures_all <- read.table('/data/rluo4/summary/Index_cell_figures_All.txt', sep = '\t', header = T)
Index_cell_figures <- read.table('/data/rluo4/database/Pancreas/Epi/Summary/Index_cell_figures.txt', sep = '\t', header = T)
Index_cell_figures$Tissue <- rep('Pancreas', nrow(Index_cell_figures))
Index_cell_figures_all <- rbind(Index_cell_figures_all, Index_cell_figures)
Index_cell_figures_all$Index <- 1:nrow(Index_cell_figures_all)
write.table(Index_cell_figures_all,'/data/rluo4/summary/Index_cell_figures.txt',quote = F,row.names = F,sep = "\t")
###############################################################################################################################
######------tissue_summary------######
###############################################################################################################################
# 7) ######summarize the sample information of 13 organs & save as  tissue_summary######
# ###############################################################################
identifier <- paste(cell_summary$orig.ident,cell_summary$sample_name,sep='_')
tissue_summary <- cell_summary[!duplicated(identifier),-1]
# tissue_summary <- cell_summary[!duplicated(cell_summary$orig.ident),]
index <- grepl('Adult',tissue_summary$sample_name) #Healthy controls
tissue_summary$PMID[index] <- '32214235'
tissue_summary$Source[index] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134355'
tissue_summary$Platform[index] <- 'Microwell-seq'
# ###############################################################################
# Breast
tissue_summary$PMID[tissue_summary$sample_name=='Nee'] <- '36914836'
tissue_summary$PMID[tissue_summary$sample_name=='Pal'] <- '33950524'
tissue_summary$PMID[tissue_summary$sample_name=='MJ'] <- '34031589'
tissue_summary$PMID[tissue_summary$sample_name=='Poornima'] <- '33763657'
tissue_summary$PMID[tissue_summary$sample_name=='Tokura'] <- '35852797'
tissue_summary$PMID[tissue_summary$sample_name=='Wei'] <- '35314812'
# Cervix
tissue_summary$PMID[tissue_summary$sample_name=='Guo'] <- '36967539'
tissue_summary$PMID[tissue_summary$sample_name=='Hua1'] <- '33996252'
tissue_summary$PMID[tissue_summary$sample_name=='Hua2'] <- '36357663'
tissue_summary$PMID[tissue_summary$sample_name=='Hua3'] <- '36569924'
tissue_summary$PMID[tissue_summary$sample_name=='Hua4'] <- 'NA'
tissue_summary$PMID[tissue_summary$sample_name=='SCP'] <- 'NA'
# CRC
tissue_summary$PMID[tissue_summary$sample_name=='Becker'] <- '35726067'
tissue_summary$PMID[tissue_summary$sample_name=='Chen'] <- '34910928'
# Esophagus
tissue_summary$PMID[tissue_summary$sample_name=='Lin'] <- '34489433'
tissue_summary$PMID[tissue_summary$sample_name=='Liu'] <- '35536873'
# GC
tissue_summary$PMID[tissue_summary$sample_name=='HY'] <- '34385296'
tissue_summary$PMID[tissue_summary$sample_name=='Kim'] <- '35087207'
tissue_summary$PMID[tissue_summary$sample_name=='Zhang'] <- '31067475'
# HNSCC
tissue_summary$PMID[tissue_summary$sample_name=='HG'] <- '36828832'
tissue_summary$PMID[tissue_summary$sample_name=='Xu'] <- 'NA'
tissue_summary$PMID[tissue_summary$sample_name=='Yu'] <- '35536873'
tissue_summary$PMID[tissue_summary$sample_name=='Williams'] <- '34129837'
# Endometrium
tissue_summary$PMID[tissue_summary$sample_name=='Berkley'] <- '35931863'
tissue_summary$PMID[tissue_summary$sample_name=='Matthew'] <- '34739872'
tissue_summary$PMID[tissue_summary$sample_name=='Ren'] <- '36273006'
tissue_summary$PMID[tissue_summary$sample_name=='Luz'] <- '36539619'
# Liver
tissue_summary$PMID[tissue_summary$sample_name=='Filliol'] <- '36198802'
tissue_summary$PMID[tissue_summary$sample_name=='Filliol.HCC'] <- '36198802'
tissue_summary$PMID[tissue_summary$sample_name=='Hassan'] <- '33332768'
tissue_summary$PMID[tissue_summary$sample_name=='Losic'] <- '31941899'
tissue_summary$PMID[tissue_summary$sample_name=='Meng'] <- '33619115'
tissue_summary$PMID[tissue_summary$sample_name=='Ramachandran'] <- '31597160'
tissue_summary$PMID[tissue_summary$sample_name=='Su'] <- '33531041'
tissue_summary$PMID[tissue_summary$sample_name=='Xinwei'] <- '36476645'
# Lung
tissue_summary$PMID[tissue_summary$sample_name=='Weimin'] <- '34764257'
tissue_summary$PMID[tissue_summary$sample_name=='Zhu'] <- '36434043'
tissue_summary$PMID[tissue_summary$sample_name=='Nayoung'] <- '32385277'
# Prostate
tissue_summary$PMID[tissue_summary$sample_name=='Daniel'] <- '35995947'
tissue_summary$PMID[tissue_summary$sample_name=='Marina'] <- '37021392'
tissue_summary$PMID[tissue_summary$sample_name=='Han'] <- '32988401'
tissue_summary$PMID[tissue_summary$sample_name=='JPeng'] <- '31273297'
# Prostate
tissue_summary$PMID[tissue_summary$sample_name=='Chen_PRAD'] <- '33420488'
tissue_summary$PMID[tissue_summary$sample_name=='Dong'] <- '33328604'
tissue_summary$PMID[tissue_summary$sample_name=='Joseph'] <- '34173975'
tissue_summary$PMID[tissue_summary$sample_name=='Song'] <- '35013146'
tissue_summary$PMID[tissue_summary$sample_name=='Vickman'] <- '35440548'
# Skin
tissue_summary$PMID[tissue_summary$sample_name=='He'] <- 'NA'
tissue_summary$PMID[tissue_summary$sample_name=='Ji'] <- '32579974'
tissue_summary$PMID[tissue_summary$sample_name=='Quan'] <- 'NA'
tissue_summary$PMID[tissue_summary$sample_name=='Boldo'] <- '32327715'
tissue_summary$PMID[tissue_summary$sample_name=='Xue'] <- '34042322'
# Stomach
tissue_summary$PMID[tissue_summary$sample_name=='HY'] <- '34385296'
tissue_summary$PMID[tissue_summary$sample_name=='Kim'] <- '35087207'
tissue_summary$PMID[tissue_summary$sample_name=='Zhang'] <- '31067475'
# THCA
tissue_summary$PMID[tissue_summary$sample_name=='Gao'] <- '33462507'
tissue_summary$PMID[tissue_summary$sample_name=='Lu'] <- '37053016'
tissue_summary$PMID[tissue_summary$sample_name=='Pan'] <- '34805166'
tissue_summary$PMID[tissue_summary$sample_name=='Peng'] <- '33588924'

##########################################################################################################
# Breast
tissue_summary$Source[tissue_summary$sample_name=='Nee'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174588'
tissue_summary$Source[tissue_summary$sample_name=='Pal'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529'
tissue_summary$Source[tissue_summary$sample_name=='MJ'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168660'
tissue_summary$Source[tissue_summary$sample_name=='Poornima'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164898'
tissue_summary$Source[tissue_summary$sample_name=='Tokura'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195861'
tissue_summary$Source[tissue_summary$sample_name=='Wei'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181254'
# Cervix
tissue_summary$Source[tissue_summary$sample_name=='Guo'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208653'
tissue_summary$Source[tissue_summary$sample_name=='Hua1'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168652'
tissue_summary$Source[tissue_summary$sample_name=='Hua2'] <- 'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11948'
tissue_summary$Source[tissue_summary$sample_name=='Hua3'] <- 'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12305'
tissue_summary$Source[tissue_summary$sample_name=='Hua4'] <- 'https://www.ebi.ac.uk/biostudies/studies/S-BSST1035?query=S-BSST1035'
tissue_summary$Source[tissue_summary$sample_name=='SCP'] <- 'https://singlecell.broadinstitute.org/single_cell/study/SCP1950/single-nucleus-rna-sequencing-and-deep-tissue-proteomics-reveal-distinct-tumour-microenvironment-in-stage-i-and-ii-cervical-cancer#study-download'
# CRC
tissue_summary$Source[tissue_summary$sample_name=='Becker'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201348'
tissue_summary$Source[tissue_summary$sample_name=='Chen'] <- 'https://data.humantumoratlas.org/explore?selectedFilters=%5B%7B%22value%22%3A%22HTAN+Vanderbilt%22%2C%22label%22%3A%22HTAN+Vanderbilt%22%2C%22group%22%3A%22AtlasName%22%2C%22count%22%3A4%2C%22isSelected%22%3Afalse%7D%5D&tab=file'
# Esophagus
tissue_summary$Source[tissue_summary$sample_name=='Lin'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160269'
tissue_summary$Source[tissue_summary$sample_name=='Liu'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199654'
# HNSCC
tissue_summary$Source[tissue_summary$sample_name=='HG'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919'
tissue_summary$Source[tissue_summary$sample_name=='Xu'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211630'
tissue_summary$Source[tissue_summary$sample_name=='Yu'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172577'
tissue_summary$Source[tissue_summary$sample_name=='Williams'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241'
# Endometrium
tissue_summary$Source[tissue_summary$sample_name=='Berkley'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203612'
tissue_summary$Source[tissue_summary$sample_name=='Matthew'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173682'
tissue_summary$Source[tissue_summary$sample_name=='Ren'] <- 'https://www.ncbi.nlm.nih.gov/sra/?term=SRP349751'
tissue_summary$Source[tissue_summary$sample_name=='Luz'] <- 'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10287'
# Liver
tissue_summary$Source[tissue_summary$sample_name=='Filliol'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174748'
tissue_summary$Source[tissue_summary$sample_name=='Filliol.HCC'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212046'
tissue_summary$Source[tissue_summary$sample_name=='Hassan'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146409'
tissue_summary$Source[tissue_summary$sample_name=='Losic'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112271'
tissue_summary$Source[tissue_summary$sample_name=='Meng'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166635'
tissue_summary$Source[tissue_summary$sample_name=='Ramachandran'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103'
tissue_summary$Source[tissue_summary$sample_name=='Su'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146115'
tissue_summary$Source[tissue_summary$sample_name=='Xinwei'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189903'
# Lung
tissue_summary$Source[tissue_summary$sample_name=='Weimin'] <- 'https://ngdc.cncb.ac.cn/gsa-human/browse/HRA001130'
tissue_summary$Source[tissue_summary$sample_name=='Zhu'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189357'
tissue_summary$Source[tissue_summary$sample_name=='Nayoung'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203360'
# Pancreas
tissue_summary$Source[tissue_summary$sample_name=='Daniel'] <- 'https://data.humantumoratlas.org/explore?selectedFilters=%5B%7B%22group%22%3A%22AtlasName%22%2C%22value%22%3A%22HTAN+WUSTL%22%7D%5D&tab=file'
tissue_summary$Source[tissue_summary$sample_name=='Han'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154778'
tissue_summary$Source[tissue_summary$sample_name=='JPeng'] <- 'https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA001063'
tissue_summary$Source[tissue_summary$sample_name=='Marina'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229413'

# Prostate
tissue_summary$Source[tissue_summary$sample_name=='Chen_PRAD'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141445'
tissue_summary$Source[tissue_summary$sample_name=='Dong'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137829'
tissue_summary$Source[tissue_summary$sample_name=='Joseph'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172357'
tissue_summary$Source[tissue_summary$sample_name=='Song'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176031'
tissue_summary$Source[tissue_summary$sample_name=='Vickman'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183676'
# Skin
tissue_summary$Source[tissue_summary$sample_name=='He'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193304'
tissue_summary$Source[tissue_summary$sample_name=='Ji'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144236'
tissue_summary$Source[tissue_summary$sample_name=='Quan'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221390'
tissue_summary$Source[tissue_summary$sample_name=='Boldo'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130973'
tissue_summary$Source[tissue_summary$sample_name=='Xue'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138669'
# Stomach
tissue_summary$Source[tissue_summary$sample_name=='HY'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167297'
tissue_summary$Source[tissue_summary$sample_name=='Kim'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150290'
tissue_summary$Source[tissue_summary$sample_name=='Zhang'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134520'
# THCA
tissue_summary$Source[tissue_summary$sample_name=='Gao'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148673'
tissue_summary$Source[tissue_summary$sample_name=='Lu'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193581'
tissue_summary$Source[tissue_summary$sample_name=='Pan'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163203'
tissue_summary$Source[tissue_summary$sample_name=='Peng'] <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158291'

##########################################################################################################
# Breast
tissue_summary$Platform[tissue_summary$sample_name=='Nee'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Pal'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='MJ'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Poornima'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Tokura'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Wei'] <- '10X Genomics'
# Cervix
tissue_summary$Platform[tissue_summary$sample_name=='Guo'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Hua1'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Hua2'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Hua3'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Hua4'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='SCP'] <- '10X Genomics'
# CRC
tissue_summary$Platform[tissue_summary$sample_name=='Becker'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Chen'] <- 'inDrop V2'
# Esophagus
tissue_summary$Platform[tissue_summary$sample_name=='Lin'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Liu'] <- '10X Genomics'
# GC
tissue_summary$Platform[tissue_summary$sample_name=='HY'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Kim'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Zhang'] <- '10X Genomics'
# HNSCC
tissue_summary$Platform[tissue_summary$sample_name=='HG'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Xu'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Yu'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Williams'] <- '10X Genomics'
# Endometrium
tissue_summary$Platform[tissue_summary$sample_name=='Berkley'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Matthew'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Ren'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Luz'] <- '10X Genomics'
# Liver
tissue_summary$Platform[tissue_summary$sample_name=='Filliol'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Filliol.HCC'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Hassan'] <- 'MARS-seq'
tissue_summary$Platform[tissue_summary$sample_name=='Losic'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Meng'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Ramachandran'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Su'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Xinwei'] <- '10X Genomics'
# Lung
tissue_summary$Platform[tissue_summary$sample_name=='Weimin'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Zhu'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Nayoung'] <- '10X Genomics'
# Pancreas
tissue_summary$Platform[tissue_summary$sample_name=='Daniel'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Han'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='JPeng'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Marina'] <- '10X Genomics'
# Prostate
tissue_summary$Platform[tissue_summary$sample_name=='Chen_PRAD'] <- 'BD Rhapsody'
tissue_summary$Platform[tissue_summary$sample_name=='Dong'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Joseph'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Song'] <- 'Seq-well'
tissue_summary$Platform[tissue_summary$sample_name=='Vickman'] <- '10X Genomics'
# Skin
tissue_summary$Platform[tissue_summary$sample_name=='He'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Ji'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Quan'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Boldo'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Xue'] <- '10X Genomics'
# Stomach
tissue_summary$Platform[tissue_summary$sample_name=='HY'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Kim'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Zhang'] <- '10X Genomics'
# THCA
tissue_summary$Platform[tissue_summary$sample_name=='Gao'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Lu'] <- '10X Genomics'
tissue_summary$Platform[tissue_summary$sample_name=='Pan'] <- 'Microwell-seq'
tissue_summary$Platform[tissue_summary$sample_name=='Peng'] <- 'BD Rhapsody'

##############################################################################################################################
# tissue_summary$Tissue[tissue_summary$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
# 
# tissue_summary$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
#                                                        gsub('Tumor','PRAD', tissue_summary$Tissue)))
unique(tissue_summary$orig.ident[tissue_summary$Tissue=='Healthy'])#133-->135
identifier <- paste(tissue_summary$orig.ident,tissue_summary$sample_name,sep='_')
length(unique(identifier))#1352
length(unique(identifier[tissue_summary$Tissue=='Healthy']))#135 Healthy amples
length(unique(identifier[tissue_summary$Tissue!='Healthy']))#1217 Diseased samples
identifier <- paste(tissue_summary$SimplifiedSampleName,tissue_summary$sample_name,sep='_')
length(unique(identifier[tissue_summary$Tissue!='Healthy']))#1192 patients
# length(unique(tissue_summary$SimplifiedSampleName[tissue_summary$Tissue!='Healthy'])) #1156-->1236
# length(unique(tissue_summary$orig.ident[tissue_summary$Tissue!='Healthy'])) #1261

tissue_summary$Index <- 1:nrow(tissue_summary)
colnames(tissue_summary)[1] <- 'Replicates'
tissue_summary$File <- paste0('/data/rluo4/database/',tissue_summary$Organ,'/Epi/sample_figures/uwot_projection_',tissue_summary$Replicates,'_3200_vargenes_all_32PCs.pdf')
index <- tissue_summary$Organ %in% c('Colorectum','Prostate') 
tissue_summary$File[index] <- gsub('3200_vargenes_all_32PCs','1600_vargenes_all_8PCs',tissue_summary$File[index])

unique(tissue_summary$sample_name)# 62 unique datasets

tissue_summary <- tissue_summary[,c('Index','Replicates','Organ','Tissue','File','PMID','Source','Platform')]
setwd('/data/rluo4/All/Output/')
colnames(tissue_summary)[3:4] <- c('Tissue','Disease Stage')
index <- tissue_summary$`Disease Stage` == "Healthy" 
tissue_summary$File[index] <- 'No Image'
tissue_summary$File <- gsub('Oral cavity','OralCavity',tissue_summary$File)
tissue_summary$Dataset <- tissue_summary$Source
tissue_summary$Dataset[grepl('ncbi',tissue_summary$Source)] <- str_split(tissue_summary$Dataset[grepl('ncbi',tissue_summary$Source)],'=',simplify = T)[,2]
tissue_summary$Dataset[grepl('HTAN',tissue_summary$Source)] <- 'HTA11'
tissue_summary$Dataset[grepl('SCP',tissue_summary$Source)] <- 'SCP1950'
tissue_summary$Dataset[grepl('ebi',tissue_summary$Source)] <- str_split(tissue_summary$Dataset[grepl('ebi',tissue_summary$Source)],'/studies/',simplify = T)[,2]
tissue_summary$Dataset[grepl('gsa',tissue_summary$Source)] <- str_split(tissue_summary$Dataset[grepl('gsa',tissue_summary$Source)],'/browse/',simplify = T)[,2]
tissue_summary$Index <- 1:nrow(tissue_summary)
length(unique(tissue_summary$Replicates))#1314 healthy + ADJ + diseased samples

table(tissue_summary$Tissue)#13 tissue types
# tissue_summary$Tissue <- gsub('THCA','Thyroid',gsub('HNSCC','OralCavity',
#                            gsub('GC','Stomach', gsub('CRC','Colorectum',tissue_summary$Tissue))))
table(tissue_summary$`Disease Stage`)
tissue_summary$Patients <- identifier
write.table(tissue_summary,'tissue_summary.txt',quote = F,row.names = F,sep = "\t")
table(cell_summary$Organ)
# cell_summary$Organ <- gsub('THCA','Thyroid',gsub('HNSCC','OralCavity',
#                            gsub('GC','Stomach', gsub('CRC','Colorectum',cell_summary$Organ))))
write.table(cell_summary,'cell_summary.txt',quote = F,row.names = F,sep = "\t")

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 8) ######------disease_summary------######
###############################################################################################################################
table(is.na(tissue_summary$Source))
disease <- paste(tissue_summary$`Disease Stage`,tissue_summary$Source)
disease_summary <- tissue_summary[!duplicated(disease),]
disease_summary$Dataset <- disease_summary$Source
disease_summary$Dataset[grepl('ncbi',disease_summary$Source)] <- str_split(disease_summary$Dataset[grepl('ncbi',disease_summary$Source)],'=',simplify = T)[,2]
disease_summary$Dataset[grepl('HTAN',disease_summary$Source)] <- 'HTA11'
disease_summary$Dataset[grepl('SCP',disease_summary$Source)] <- 'SCP1950'
disease_summary$Dataset[grepl('ebi',disease_summary$Source)] <- str_split(disease_summary$Dataset[grepl('ebi',disease_summary$Source)],'/studies/',simplify = T)[,2]
disease_summary$Dataset[grepl('gsa',disease_summary$Source)] <- str_split(disease_summary$Dataset[grepl('gsa',disease_summary$Source)],'/browse/',simplify = T)[,2]
disease_summary$Index <- 1:nrow(disease_summary)
# write.table(disease_summary[,c(1,3,4,9,7,6,8)],'disease_summary.txt',quote = F,row.names = F,sep = "\t")

# rm(list=ls())
data_path = '/data/rluo4/All/Output'
setwd(data_path)
cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
tissue_summary <- read.csv('tissue_summary.txt',sep = "\t",header = T)
unique(tissue_summary$Disease.Stage)
tissue_summary$Stage <- tissue_summary$Disease.Stage
precancer <- c('AAH','AD','AEH','AK',#'AIS',
               'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
               'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
               'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
               'SER','SIM','WIM')
unique(tissue_summary$Stage[! tissue_summary$Stage %in% c(precancer,'Healthy','ADJ')])
tissue_summary$Stage[tissue_summary$Stage %in% precancer] <- 'Precancer'
tissue_summary$Stage[! tissue_summary$Disease.Stage %in% c(precancer,'Healthy','ADJ')] <- 'Cancer'
table(tissue_summary$Stage)
tissue_summary$Stage <- factor(tissue_summary$Stage, ordered=T, levels = c('Healthy','ADJ','Precancer','Cancer'))  #调整画图的x轴坐标顺序
# tissue_summary$Tissue <- sort(tissue_summary$Tissue)
# tissue_summary$Tissue <- factor(tissue_summary$Tissue,ordered=F)
library(RColorBrewer)  
library(reshape2)
library(ggplot2)
RColorBrewer::brewer.pal(n = 12, name = 'Paired')
color_cluster=c("#1B9E77","#66A61E","skyblue","#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3")
color_cluster=c(RColorBrewer::brewer.pal(n = 12, name = 'Paired'),color_cluster)
color_cluster <- color_cluster[1:length(unique(tissue_summary$Tissue))]
names(color_cluster)=unique(tissue_summary$Tissue)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
dfm <- as.data.frame(table(tissue_summary$Tissue))
colnames(dfm) <- c('Tissue','Freq')
ggplot(dfm,aes(x = Freq, y = Tissue))+#, group = group)) +
  geom_col(aes(fill=Tissue)) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
  geom_text(aes(label = Freq), position = position_stack(vjust = .5), size = 6) + # labels inside the bar segments
  scale_fill_manual(values = color_cluster)+
  # scale_y_reverse("Datasets by Primary site",expand = c(0,0),position = 'left',limits=c(250,0))+#,position = "right",limits=c(5000,0))+  #yÖáÔÚÓÒ
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  # coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # ×ÖÌåÑùÊ½£¬¿ÉÑ¡ bold, plain, italic
    base_family = "serif", # ×ÖÌå¸ñÊ½£¬¿ÉÑ¡ serif, sans, mono, ArialµÈ
    base_size = 17,  # Í¼ÐÎµÄ×ÖÌå´óÐ¡
    base_line_size = 0.8, # ×ø±êÖáµÄ´ÖÏ¸
    axis_text_angle = 0) + theme(legend.position='none')
ggsave(file = paste0(data_path, '/res_fig/tissue_samples.png'),width = 9,height =7.5)

dfm <- as.data.frame(table(cell_summary$Organ))
colnames(dfm) <- c('Tissue','Freq')
ggplot(dfm,aes(x = Freq, y = Tissue))+#, group = group)) +
  geom_col(aes(fill=Tissue)) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
  geom_text(aes(label = Freq), position = position_stack(vjust = .5), size = 5) + # labels inside the bar segments
  scale_fill_manual(values = color_cluster)+
  # scale_y_reverse("Datasets by Primary site",expand = c(0,0),position = 'left',limits=c(250,0))+#,position = "right",limits=c(5000,0))+  #yÖáÔÚÓÒ
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  # coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # ×ÖÌåÑùÊ½£¬¿ÉÑ¡ bold, plain, italic
    base_family = "serif", # ×ÖÌå¸ñÊ½£¬¿ÉÑ¡ serif, sans, mono, ArialµÈ
    base_size = 17,  # Í¼ÐÎµÄ×ÖÌå´óÐ¡
    base_line_size = 0.8, # ×ø±êÖáµÄ´ÖÏ¸
    axis_text_angle = 0) + theme(legend.position='none')
ggsave(file = paste0(data_path, '/res_fig/tissue_cells.png'),width = 9,height =7.5)

library(ArchR)
bar.df <- tissue_summary# mutate(bar.df,name=factor(bar.df$Tissue))#
title <- paste0('Tissue type fraction')
color_cluster= color_cluster[c('Lung','Thyroid','Esophagus','Colorectum')]#RColorBrewer::brewer.pal(n = 5, name = 'Paired')
# color_cluster <- color_cluster[1:length(unique(bar.df$Stage))]
names(color_cluster)=levels(bar.df$Stage)
# color_cluster = paletteDiscrete(values = unique(bar.df$Stage))
p <- ggplot(bar.df,aes(x=Tissue))+
  geom_bar(aes(fill=Stage),position = "fill",width = .7)+
  scale_x_discrete("")+
  scale_y_continuous(title,expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("",values = color_cluster)+
  theme_ArchR(baseSize = 10) + #ggtitle('') +
  # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
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
    # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  )+   coord_flip() #让条形图横过来
# hei <- ceiling(length(unique(bar.df$Tissue))*1.5)
f = paste0(data_path,'/res_fig/tissue_stage_fraction.png')
p
ggsave(plot=p,width = 9,height =7.5, filename=f, dpi = 500, device = "png")


# from nbconvert import PythonExporter
# import nbformat
# notebook_path='/data/rluo4/lorihan/Pyscript/Ps_database/Liver-sum/'
# python_path='/data/rluo4/lorihan/Pyscript/script/'
# def notebook_to_python(notebook_path, python_path):
#   with open(notebook_path, 'r', encoding='utf-8') as notebook_file:
#     notebook_contents = nbformat.read(notebook_file, as_version=4)
# 
# python_exporter = PythonExporter()
# python_code, _ = python_exporter.from_notebook_node(notebook_contents)
# 
# with open(python_path, 'w', encoding='utf-8') as python_file:
#   python_file.write(python_code)
# # Usage
# notebook_to_python('/data/rluo4/lorihan/Pyscript/Ps_database/Liver-sum/scImm_pipeline.ipynb', '/data/rluo4/lorihan/Pyscript/script/scImm_pipeline.py')

celltype_summary <- read.csv('/data/rluo4/summary/celltype_summary.txt', sep = '\t') # Table S1; this file is modified by database work in June,2023
print(organ)
organ = 'Pancreas'

epi <- gsub('Chen','CRC',paste(organ, 'Epi', sep = "_"))
str <-  gsub('Chen','CRC',paste(organ, 'Str', sep = "_"))
imm <- gsub('Chen','CRC',paste(organ, 'Imm', sep = "_"))

Epi <- cell_all[[epi]]
Epi$barcode = rownames(Epi)
Epi$cell_type <- Epi$Cell_Type
table(Epi$Organ,Epi$Tissue)
Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                            gsub('Tumor','PRAD', Epi$Tissue)))
Epi$Major_type <- 'Epithelial'
print(table(Epi$Organ))
Imm <- cell_all[[imm]]
Imm$Major_type <- 'Immune'

Str <- cell_all[[str]]
Str$Major_type <- 'Stromal'
setwd('/data/rluo4/summary/Annotation')
epi_marker <- read.csv(paste0(epi, 'Marker.csv'), header = F, row.names = 1)
imm_marker <- read.csv(paste0(imm, 'Marker.csv'), header = F, row.names = 1)
str_marker <- read.csv(paste0(str, 'Marker.csv'), header = F, row.names = 1)
setdiff(unique(Str$Cell_Type), rownames(str_marker))
setdiff(unique(Imm$Cell_Type), rownames(imm_marker))

minor <- c(unique(Epi$Cell_Type), unique(Imm$Cell_Type), unique(Str$Cell_Type))
Pancreas_marker <- data.frame(Index = 1:length(minor), Tissue = rep(organ, length(minor)),
                              Major.Celltype = c(rep('Epi', length(unique(Epi$Cell_Type))), rep('Imm', length(unique(Imm$Cell_Type))), rep('Str', length(unique(Str$Cell_Type)))), 
                              Minor.Celltype = c(unique(Epi$Cell_Type), unique(Imm$Cell_Type), unique(Str$Cell_Type))
)

Pancreas_marker$FullName <- celltype_summary$FullName[match(Pancreas_marker$Minor.Celltype, celltype_summary$Minor.Celltype)]
Pancreas_marker$FullName[is.na(Pancreas_marker$FullName )] <- c('Duct-like1 cells', 'Acinar cells', 'Pancreatic islet cells',
                                                                'PanIN-like cells', 'Duct-like2 cells', 'Pancreatic stellate cells',
                                                                'Erythrocytes', 'Antigen-presenting cancer-associated fibroblasts')

epi.m <- apply(epi_marker,1,function(x){
  paste(x[x!=''], collapse=',') })
imm.m <- apply(imm_marker,1,function(x){
  paste(x[x!=''], collapse=',') })
str.m <- apply(str_marker,1,function(x){
  paste(x[x!=''], collapse=',') })
marker_genes <- c(epi.m, imm.m, str.m)
Pancreas_marker$Marker <- marker_genes[match(Pancreas_marker$Minor.Celltype, names(marker_genes))]
index = is.na(Pancreas_marker$Marker)
Pancreas_marker$Marker[index]  <- celltype_summary$Marker[match(Pancreas_marker$Minor.Celltype[index], celltype_summary$Minor.Celltype)]

celltype_summary <- rbind(celltype_summary, Pancreas_marker)
celltype_summary$Tissue <- gsub('THCA','Thyroid',gsub('HNSCC','Oral cavity',gsub('Chen','Colorectum',
                                                                                 gsub('GC','Stomach', gsub('CRC','Colorectum',celltype_summary$Tissue)))))
celltype_summary$Minor.Celltype <- gsub('MDSCs','MDSC',celltype_summary$Minor.Celltype)

celltype_summary <- celltype_summary[order(celltype_summary$Tissue),]
celltype_summary$Index <- 1:nrow(celltype_summary)
celltype_summary$Marker[celltype_summary$Tissue=='Liver' & celltype_summary$Minor.Celltype=='STM']<- 
  paste(c('EPCAM','CD24', 'ANPEP', 'SOX9', 'CD47'), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Liver' & celltype_summary$Minor.Celltype=='CHO']<- 
  paste(c("PROM1","HNF1B","ONECUT1","SCTR","CFTR","KRT7"), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Liver' & celltype_summary$Minor.Celltype=='HEP']<- 
  paste(c("TTR","ALB",'APOE','APOA1','CYP3A4','CPS1','G6PC'), collapse=',') 

celltype_summary$Marker[celltype_summary$Tissue=='Colorectum' & celltype_summary$Minor.Celltype=='STM']<- 
  paste(c('SMOC2','RGMB', 'SOX9','LGR5','OLFM4','ASCL2','CD44','CD166'), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Colorectum' & celltype_summary$Minor.Celltype=='EE']<- 
  paste(c('CHGA',"CHGB","TPH1"), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Colorectum' & celltype_summary$Minor.Celltype=='ABS']<- 
  paste(c('ALDOB',"CA2", "SI"), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Colorectum' & celltype_summary$Minor.Celltype=='SSC']<- 
  paste(c('GSDMB','GSDMD','IL18','RELB','MDK','RARA','RXRA','PDX1','S100P'), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Colorectum' & celltype_summary$Minor.Celltype=='ASC']<- 
  paste(c('AXIN2','RNF43','EPHB2','CDX2'), collapse=',') 

celltype_summary$Marker[celltype_summary$Tissue=='Prostate' & celltype_summary$Minor.Celltype=='STM']<- 
  paste(c('ITGB1', 'DLL1','LAMC2',"TP63",'VIM'), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Prostate' & celltype_summary$Minor.Celltype=='BAS']<- 
  paste(c("KRT14","ITGB4","KRT15","KRT5"), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Prostate' & celltype_summary$Minor.Celltype=='LUM']<- 
  paste(c("KLK3", 'KLK2',"ACPP", "MSMB", 'NKX3-1','AR','KRT8','KRT18','TACSTD2'), collapse=',') 
celltype_summary$Marker[celltype_summary$Tissue=='Prostate' & celltype_summary$Minor.Celltype=='NUER']<- 
  paste(c('SYP', 'KRT4', 'LY6D'), collapse=',') 

write.table(celltype_summary, file = paste0('/data/rluo4/All/Output/celltype_summary.csv'),quote = F, row.names = F)
write.table(celltype_summary, file = paste0('/data/rluo4/summary/Annotation/celltype_summary.txt'),quote = F, row.names = F, sep = '\t')

library(openxlsx)
hs <- createStyle(
  textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,
  fontName = "Arial Narrow", fgFill = "#4F80BD"
)
## Not run: 
write.xlsx(celltype_summary[,-1], file = paste0('/data/rluo4/All/Output/table_xlsx/TableS1.xlsx'), rowNames = FALSE)#,
          # colNames = TRUE, borders = "rows", headerStyle = hs)

setdiff(unique(cell_all$Cell_Type[cell_all$Major_type=='Stromal']), unique(celltype_summary$Minor.Celltype[celltype_summary$Major.Celltype=='Str']) )
# [1] "MSC.ADIPO" "MSC.MVA"   "MSC.SEC"   "MSC.PVA"  
bar.df <- disease_summary
#bar.df$AccessionID[duplicated(bar.df$AccessionID)] <- paste(bar.df$AccessionID[duplicated(bar.df$AccessionID)],'mus',sep = '_')
level <- unique(bar.df$Tissue)#unique(bar.df$AccessionID)
bar.df <- mutate(bar.df,name=factor(bar.df$Tissue, levels=level))
text.df <- as.data.frame(table(bar.df$name))
table(disease_summary$Platform)
color_cluster=c("#f89e81","#99a9cc","#acd485","#dd9bc5","#f6d573","#84c7b3")
names(color_cluster)=unique(disease_summary$Platform)#c('Human','Mus musculus')#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
dfm = data.frame(table(disease_summary$Tissue,disease_summary$Disease.Stage))
colnames(dfm)[1:2] <- c('Tissue','Platform')#bar.df <- bar.df[order(bar.df$Organ,decreasing = T),]
library(ggprism)
#colnames(bar.df)[3] <- 'State'
pdf(file = "Platform.pdf",height = 7,width = 7)
# bar.df <- bar.df[ordered(bar.df$Tissue),]
ggplot(bar.df,aes(x=Tissue))+
  geom_bar(aes(fill=Platform),position = "stack",width = .7)+
  scale_x_discrete("",position = "top")+  #x轴在上
  scale_fill_manual(values = color_cluster)+ ggtitle('Single Cell RNA-seq Platform')+
  scale_y_reverse("Num of Datasets by Tissue",expand = c(0,0),position = 'left',limits=c(15,0))+#,position = "right",limits=c(5000,0))+  #y轴在右
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # 字体样式，可选 bold, plain, italic
    base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
    base_size = 17,  # 图形的字体大小
    base_line_size = 0.8, # 坐标轴的粗细
    axis_text_angle = 0) + theme(legend.position='top')# 可选值有 0，45，90，270
# dev.off()
ggsave(paste0(data_path, 'Platform.png'), height = 6, width = 8, dpi = 500)

# 
# A <- data.frame(ratio, disease)#构建一个数据框
# library(ggplot2)
# library(ggforce)
# 
# ggplot()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_blank(),
#         legend.title=element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
#   xlab("")+ylab('')+#添加颜色
#   scale_fill_manual(values = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0',
#                                '#D6E7A3', '#57C3F3', '#476D87',
#                                '#E59CC4', '#AB3282', '#23452F', '#BD956A'))+
#   geom_arc_bar(data=,
#                stat = "pie",
#                aes(x0=0,y0=0,r0=0,r=2,
#                    amount=ratio,fill=disease)
#   )+#饼图
#   annotate("text",x=1.6,y=1.5,label="24.20%",angle=-50)+
#   annotate("text",x=1.6,y=-1.5,label="21.9%",angle=45)+
#   annotate("text",x=0,y=-2.2,label="7.6%",angle=0)+
#   annotate("text",x=-0.8,y=-2,label="5.2%",angle=-20)+
#   annotate("text",x=-1.3,y=-1.7,label="4.3%",angle=-40)+
#   annotate("text",x=-1.6,y=1.5,label="24.8%",angle=45)#手动注释，还是很麻烦

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 9) ######------cytotrace_summary------######
###############################################################################################################################
load('/data/rluo4/All/Output/pheno_all.rds')
table(pheno_all$Liver_Epi$tissue)
table(pheno_all$Cervix_Epi$Tissue)
library(ggpubr)
library(ggplot2)
library(ggprism)
library(ggthemes)
library(ggridges)
# cell_all <- c(pheno_all, Imm_all, Str_all)
pheno <- paste(organ_all, 'Epi', sep = '_')
# for (obs in pheno) {
#   # obs = 'HNSCC_Epi'#'CRC_Epi'#'Cervix_Epi'
#   print(obs)
#   organ <- str_split(obs,'_', simplify = T)[,1]
#   print(organ)
#   epi <- paste(organ, 'Epi', sep = "_")
#   
#   # score_result <- pheno_all$Liver_Epi
#   Epi <- pheno_all[[epi]]
#   Epi$barcode = rownames(Epi)
#   Epi$cell_type <- Epi$Cell_Type
#   # table(Epi$Organ,Epi$Tissue)
#   Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
#   Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
#                                               gsub('Tumor','PRAD', Epi$Tissue)))
#   Epi$Major_type <- 'Epithelial'
#   print(table(Epi$Organ))
#   common_columns <- c('barcode','orig.ident','batch','SimplifiedSampleName', 'sample_name','Tissue','Cell_Type','Major_type','CytoTRACE','Malignancy')#intersect(intersect(colnames(Epi),colnames(Imm)),colnames(Str))
#   phe = dplyr::bind_rows(Epi) %>%
#     select(all_of(common_columns))
#   
#   phe <- phe[! phe$sample_name %in% c('Cold'),]
#   phe <- phe[! phe$orig.ident %in% c("p3", "p4", "p5"),]
#   phe$Organ <- unique(Epi$Organ)
#   table(phe$Tissue)
#   phe <- phe[phe$Tissue!='UNC',]
#   score_result <- phe#[phe$Tissue!='Healthy',]
#   unique(score_result$SimplifiedSampleName)
#   # score_result$Tissue <- ordered(score_result$Tissue,levels=c("Cyst", "NAFLD","Cirrhotic","HCC"))
#   # score_result$Cell_Type <- ordered(score_result$Cell_Type,levels=c("STM", "CHO","HEP"))
#   Cyto_dir <- paste0("/data/rluo4/All/Output/res_fig/Cytotrace/",organ)
#   if (!dir.exists(paste0(Cyto_dir))){
#     dir.create(paste0(Cyto_dir))
#   }
#   setwd(Cyto_dir)
#  
#   ggplot(score_result) +   guides(fill = guide_legend(title = 'Disease State'))+
#     stat_compare_means(aes(x=Tissue,y=CytoTRACE,group = Tissue),
#                        label = "p.format",label.y =1,label.x = 1.25,# label.y = c(9.5,11,12),
#                        #method = "wilcox.test"
#     ) + 
#     geom_boxplot(data=score_result,aes(x=Tissue,y=CytoTRACE,fill=Tissue),notch="TRUE",
#                  width=0.35,position = position_dodge(0),size=0.5,outlier.size = 0) +
#     # geom_point(data=score_result,aes(x=Tissue,y=CytoTRACE,fill=Tissue),shape=23,colour="grey15",size=2)+
#     # geom_bezier(data=data.cells.1,aes(x=treatment,y=value,group=index,linetype = 'cubic'),
#     #             size=0.25,colour="grey20") + #需要调用ggforce包
#     scale_fill_manual(values=RColorBrewer::brewer.pal(n = 12, name = 'Paired')) +
#     facet_grid(.~Cell_Type) + ggtitle(paste0('Developmental potential in ',organ,' Tissues')) +
#     labs(x="",y="CytoTRACE Score")+theme_base(base_size = 18)+
#     theme(legend.position = "top",
#           plot.title    = element_text(color = 'black', size   = 15, hjust = 0.5),
#           plot.subtitle = element_text(color = 'black', size   = 15,hjust = 0.5),
#           plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#           axis.text.x   = element_blank(),#element_text(color = 'black', size = 15, angle = 310),
#           axis.text.y   = element_text(color = 'black', size = 15, angle = 0),
#           axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
#           axis.title.y  = element_text(color = 'black', size = 15, angle = 90),
#           #legend.title=element_blank(),
#           legend.title  = element_text(color = 'black', size  = 15),
#           legend.text   = element_text(color = 'black', size   = 15),
#           axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#           axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#           panel.background = element_rect(fill='transparent')#, legend.position='right'
#     )
#   file = paste0(Cyto_dir,'/Cytotrace_in_', organ,'.png')
#   ggsave(height=6.5,width=length(unique(score_result$Cell_Type))*1.5, filename=file, dpi = 500, device = "png")
#   
#   table(is.na(score_result$Malignancy))
#   library(ggpubr)
#   options(bitmapType = 'cairo')
#   file = paste0(Cyto_dir,'/Malignancy_in_', organ,'.png')
#   # pdf(file,width=6,height = 5.5)
#   ggscatter(score_result[!is.na(score_result$Malignancy),], y = "CytoTRACE", x = "Malignancy",
#             color = "#CAB2D6", size = 1,
#             add = "reg.line", add.params = list(color = "salmon", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
#             conf.int = TRUE, # 添加回归线的置信区间
#             cor.coef = TRUE, # 添加相关系数
#             cor.coeff.args = list(method = "pearson", label.x = 0,label.y=2, label.sep =", ") #"\n")#选择Pearson相关
#   )+ xlab('Malignancy continumn')+ylab("CytoTRACE Score")#+ stat_cor( label.x = 3)
#   ggsave(height=4.5,width=5, filename=file, dpi = 500, device = "png")
#   
#   # score_result <- score_result[! score_result$Tissue %in% c('ADJ'), ]#c("NAFLD","Cirrhotic",'HCC'),]
#   score_result$Tissue <- gsub('ADJ', 'Healthy',score_result$Tissue )
#   # score_result$Tissue <- gsub('CSG', 'Benign', gsub('Healthy', 'Benign',score_result$Tissue ))
#   # score_result$Tissue <- gsub('Cyst', 'Benign', gsub('N_HPV', 'Benign',score_result$Tissue ))
#   # score_result$Tissue <- gsub('Cyst', 'Benign', gsub('Goiters', 'Benign',score_result$Tissue ))
#   p0 <- ggplot(score_result, aes(x = CytoTRACE , y =Tissue  , fill = Tissue)) +
#     geom_density_ridges(alpha = 0.5) + xlab('') + ylab('') +#xlab('CytoTRACE') + ylab('Tissue Type') + (in Cytotrace_1)
#     theme_ridges() + 
#     theme_ArchR(baseSize = 10) + ggtitle(organ) +
#     theme(
#       plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
#       plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
#       plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#       axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
#       axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
#       axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
#       axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
#       legend.title  = element_text(color = 'black', size  = 16),
#       legend.text   = element_text(color = 'black', size   = 16),
#       axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#       axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#       panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
#     ) + theme(legend.position = "none")
#   p0
#   f = paste0(organ,'_TissueType_CytoTRACE.png')
#   ggsave(plot=p0,height=5,width=5, filename=f, dpi = 300, device = "png")
#   
#   print(getwd())
#   library(reshape2)
#   library(ggplot2)
#   library(ggridges)
#   table(score_result$Tissue)
#   table(score_result$Cell_Type)
#   #要用到geom_density_ridges()函数需要调用ggridges包
#   data_melt = melt(score_result, id.vars=c("Cell_Type","Tissue"),  measure.vars='CytoTRACE')
#   data_melt[1:5,]
#   table(data_melt$Tissue)
#   typeCyto <- data_melt$Cell_Type
#   print(unique(typeCyto))
#   for (CellType in unique(typeCyto)) {
#     print(paste0(CellType, ' in ', Cyto_dir))
#     setwd(Cyto_dir)
#     p1 <- ggplot(data_melt[data_melt$Cell_Type==CellType,], aes(x = value , y =Tissue  , fill = Tissue)) +
#       geom_density_ridges(alpha = 0.5) + xlab('') + ylab('') +#ylab('Cell Type') +
#       theme_ridges() + 
#       theme_ArchR(baseSize = 10) + ggtitle(CellType) +
#       # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#       # scale_color_manual(values=c("Reference"="#D5D5D5", 'ABS'='#D51F26','ASC'='#272E6A','CT'='#208A42','EE'='#89288F', 'GOB'='#F47D2B','IMENT'='#FEE500','IMGOB'='#8A9FD1','SSC'='#C06CAB','STM'='#D8A767','TAC'='#90D5E4','TUF'='#89C75F'))+
#       theme(
#         plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
#         plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
#         plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#         axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
#         axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
#         axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
#         axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
#         legend.title  = element_text(color = 'black', size  = 16),
#         legend.text   = element_text(color = 'black', size   = 16),
#         axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#         axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#         panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
#       ) + theme(legend.position = "none")
#     p1
#     f = paste0('developmental_potential_in_',CellType,'.png')
#     hei =   wid <- ifelse(length(unique(data_melt$Cell_Type))>10,length(unique(data_melt$Cell_Type))*0.6,
#                           ifelse(length(unique(data_melt$Cell_Type))>5,length(unique(data_melt$Cell_Type))*0.9,
#                                  length(unique(data_melt$Cell_Type))*1.2))
#     ggsave(plot=p1,height=hei,width=5, filename=f, dpi = 300, device = "png")
#   }
#   
# }
# 
# typeCyto <- data_melt$Tissue
# print(unique(typeCyto))
# typeCyto <- gsub('ADJ','Healthy',typeCyto)
# for (TissueType in unique(typeCyto)) {
#   print(paste0(TissueType, ' in ', cohort_directory))
#   setwd(Cyto_dir)
#   p1 <- ggplot(data_melt[data_melt$Tissue==TissueType,], aes(x = value , y =Cell_Type  , fill = Cell_Type)) +
#     geom_density_ridges(alpha = 0.5) + xlab('CytoTRACE') + ylab('Cell Type') +
#     theme_ridges() +
#     theme_ArchR(baseSize = 10) + ggtitle(TissueType) +
#     # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#     scale_color_manual(values=c("Reference"="#D5D5D5", 'ABS'='#D51F26','ASC'='#272E6A','CT'='#208A42','EE'='#89288F', 'GOB'='#F47D2B','IMENT'='#FEE500','IMGOB'='#8A9FD1','SSC'='#C06CAB','STM'='#D8A767','TAC'='#90D5E4','TUF'='#89C75F'))+
#     theme(
#       plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
#       plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
#       plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#       axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
#       axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
#       axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
#       axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
#       legend.title  = element_text(color = 'black', size  = 16),
#       legend.text   = element_text(color = 'black', size   = 16),
#       axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#       axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#       panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
#     ) + theme(legend.position = "none")
#   p1
#   f = paste0('developmental_potential_in_',TissueType,'.png')
#   hei =   wid <- ifelse(length(unique(data_melt$Cell_Type))>10,length(unique(data_melt$Cell_Type))*0.6,
#                         ifelse(length(unique(data_melt$Cell_Type))>5,length(unique(data_melt$Cell_Type))*0.9,
#                                length(unique(data_melt$Cell_Type))*1.2))
#   ggsave(plot=p1,height=hei,width=5, filename=f, dpi = 300, device = "png")
# }
# pheno <- c('HNSCC_Epi','Breast_Epi','Liver_Epi','THCA_Epi','CRC_Epi','Skin_Epi')
pheno <- paste(organ_all, 'Epi', sep = '_')
for (obs in pheno) {
  # obs = 'HNSCC_Epi'#'CRC_Epi'#'Cervix_Epi'
  print(obs)
  organ <- str_split(obs,'_', simplify = T)[,1]
  print(organ)
  epi <- paste(organ, 'Epi', sep = "_")
  
  # score_result <- pheno_all$Liver_Epi
  Epi <- pheno_all[[epi]]
  Epi$barcode = rownames(Epi)
  Epi$cell_type <- Epi$Cell_Type
  # table(Epi$Organ,Epi$Tissue)
  Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
  Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                              gsub('Tumor','PRAD', Epi$Tissue)))
  Epi$Major_type <- 'Epithelial'
  print(table(Epi$Organ))
  print(table(Epi$Tissue))
  common_columns <- c('barcode','orig.ident','batch','SimplifiedSampleName', 'sample_name','Tissue','Cell_Type','Major_type','CytoTRACE','Malignancy')#intersect(intersect(colnames(Epi),colnames(Imm)),colnames(Str))
  phe = dplyr::bind_rows(Epi) %>%
    select(all_of(common_columns))
  
  phe <- phe[! phe$sample_name %in% c('Cold'),]
  phe <- phe[! phe$orig.ident %in% c("p3", "p4", "p5"),]
  phe$Organ <- unique(Epi$Organ)
  table(phe$Tissue)
  phe <- phe[phe$Tissue!='UNC',]
  score_result <- phe#[phe$Tissue!='Healthy',]
  unique(score_result$SimplifiedSampleName)
  table(score_result$Tissue)
  
  Cyto_dir <- paste0("/data/rluo4/All/Output/res_fig/Cytotrace_1/",organ)
  if (!dir.exists(paste0(Cyto_dir))){
    dir.create(paste0(Cyto_dir))
  }
  setwd(Cyto_dir)
  table(score_result$Tissue)
  table(score_result$Cell_Type)
  score_result <- score_result[! score_result$Tissue %in% c('Healthy','ADJ','Cyst','CRC','FAP','Goiters','CSG'), ]#c("NAFLD","Cirrhotic",'HCC'),]
  print(getwd())
  if(obs == 'Liver_Epi'){
    score_result$Tissue <- ordered(score_result$Tissue,levels=c( "NAFLD","Cirrhotic","HCC"))
  }
  if(obs == 'Skin_Epi'){
    score_result$Tissue <- ordered(score_result$Tissue,levels=c( "AK","SCCIS","cSCC"))
  }
  if(obs == 'Cervix_Epi'){
    score_result$Tissue <- ordered(score_result$Tissue,levels=c( "N_HPV","HSIL_HPV","CC"))
  }
  if(obs == 'HNSCC_Epi'){
    score_result$Tissue <- ordered(score_result$Tissue,levels=c( "NEOLP","EOLP","LP","OSCC"))
  }
  # score_result$Cell_Type <- ordered(score_result$Cell_Type,levels=c("STM", "CHO","HEP"))
  library(reshape2)
  library(ggplot2)
  library(ggridges)
  table(score_result$Tissue)
  table(score_result$Cell_Type)
  #要用到geom_density_ridges()函数需要调用ggridges包
  data_melt = melt(score_result, id.vars=c("Cell_Type","Tissue"),  measure.vars='CytoTRACE')
  data_melt[1:5,]
  table(data_melt$Tissue)
  typeCyto <- data_melt$Cell_Type
  print(unique(typeCyto))
  for (CellType in unique(typeCyto)) {
    print(paste0(CellType, ' in ', Cyto_dir))
    setwd(Cyto_dir)
    p1 <- ggplot(data_melt[data_melt$Cell_Type==CellType,], aes(x = value , y =Tissue  , fill = Tissue)) +
      geom_density_ridges(alpha = 0.5) + xlab('') + ylab('') +#ylab('Cell Type') +
      theme_ridges() + 
      theme_ArchR(baseSize = 10) + ggtitle(CellType) +
      # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
      # scale_color_manual(values=c("Reference"="#D5D5D5", 'ABS'='#D51F26','ASC'='#272E6A','CT'='#208A42','EE'='#89288F', 'GOB'='#F47D2B','IMENT'='#FEE500','IMGOB'='#8A9FD1','SSC'='#C06CAB','STM'='#D8A767','TAC'='#90D5E4','TUF'='#89C75F'))+
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
      ) + theme(legend.position = "none")
    p1
    f = paste0('developmental_potential_in_',CellType,'.png')
    hei =   wid <- ifelse(length(unique(data_melt$Cell_Type))>10,length(unique(data_melt$Cell_Type))*0.6,
                          ifelse(length(unique(data_melt$Cell_Type))>5,length(unique(data_melt$Cell_Type))*0.9,
                                 length(unique(data_melt$Cell_Type))*1.2))
    ggsave(plot=p1,height=hei,width=5, filename=f, dpi = 300, device = "png")
  }
  score_result <- score_result[ score_result$Cell_Type %in% c('STM','HEP','KER','COL','ASC','ABS','SPI','COR','MLUM','TFC','BAS','LUM'), ]#c("NAFLD","Cirrhotic",'HCC'),]
  ggplot(score_result) +   guides(fill = guide_legend(title = ''))+
    stat_compare_means(aes(x=Tissue,y=CytoTRACE,group = Tissue),
                       label = "p.format",label.y =1,label.x = 1.25,# label.y = c(9.5,11,12),
                       #method = "wilcox.test"
    ) + 
    geom_boxplot(data=score_result,aes(x=Tissue,y=CytoTRACE,fill=Tissue),notch="TRUE",
                 width=0.35,position = position_dodge(0),size=0.5,outlier.size = 0) +
    # geom_point(data=score_result,aes(x=Tissue,y=CytoTRACE,fill=Tissue),shape=23,colour="grey20",size=2)+
    # geom_bezier(data=data.cells.1,aes(x=treatment,y=value,group=index,linetype = 'cubic'),
    #             size=0.25,colour="grey20") + #需要调用ggforce包
    scale_fill_manual(values=RColorBrewer::brewer.pal(n = 12, name = 'Paired')) +
    facet_grid(.~Cell_Type) + #ggtitle(paste0('Developmental potential in ',organ,' Tissues')) +
    labs(x="",y="")+theme_base(base_size = 22)+ #y="CytoTRACE Score"
    theme(legend.position = "top",
          plot.title    = element_text(color = 'black', size   = 20, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 20,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 20,face = 'italic', hjust = 1),
          axis.text.x   = element_blank(),#element_text(color = 'black', size = 20, angle = 310),
          axis.text.y   = element_text(color = 'black', size = 12, angle = 0),
          axis.title.x  = element_text(color = 'black', size = 20, angle = 0),
          axis.title.y  = element_text(color = 'black', size = 20, angle = 90),
          #legend.title=element_blank(),
          legend.title  = element_text(color = 'black', size  = 17),
          legend.text   = element_text(color = 'black', size   = 20),
          axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
          axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
          # panel.background = element_rect(fill='transparent')#, legend.position='right'
          panel.border = element_blank(),#element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
    ) + theme( panel.border = element_blank() )
  file = paste0(Cyto_dir,'/Cytotrace_in_', organ,'.png')
  wid = length(unique(score_result$Cell_Type)) + length(unique(score_result$Tissue))
  ggsave(height=4,width=wid, filename=file, dpi = 500, device = "png")
  
}
compare_means(CytoTRACE~Tissue, score_result[score_result$Cell_Type=='STM' & score_result$Tissue!='SCCIS',])
# score_result$cell_type <- paste(score_result$Cell_Type,score_result$Tissue,sep = ':')#sce_sub$cluster
# # score_result$cell_type <- paste(score_result$Tissue,score_result$Cell_Type,sep = ':')#sce_sub$cluster
# table(score_result$cell_type)
# # score_result$cell_type <- factor(score_result$cell_type, levels=c('NAFLD:CHO','NAFLD:HEP','NAFLD:STM',
# #                                                                   'Cirrhotic:HEP','Cirrhotic:STM',
# #                                                                   'HCC:HEP','HCC:STM'))
# score_result$cell_type <- factor(score_result$cell_type, levels=c('CHO:NAFLD','CHO:Cirrhotic','CHO:HCC',
#                                                                   'HEP:NAFLD','HEP:Cirrhotic','HEP:HCC', 'HCC:HEP','STM:NAFLD','STM:Cirrhotic','STM:HCC'))
ggstatsplot::ggwithinstats(#ggbetweenstats(
  data = score_result,
  x = cell_type,
  y = CytoTRACE,
  sort = "descending", # ordering groups along the x-axis based on
  sort.fun = median, # values of `y` variable
  pairwise.comparisons = FALSE,
  # pairwise.display = "s",
  pairwise.annotation = "p",
  # caption = "Data from: iris",
  # ggtheme = ggthemes::theme_fivethirtyeight(),
  ggstatsplot.layer = FALSE,
  # pairwise.display = 'all',
  # p.adjust.method = 'BH',#'bonferroni',
  notch = TRUE, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = TRUE, # whether to display confidence interval for means
  mean.label.size = 2.55,
  outlier.tagging = FALSE, # whether outliers need to be tagged
  xlab = "Cell:Sample Type", # label for the x-axis variable
  ylab = "CytoTRACE Score", # label for the y-axis variable
  # title = "Liver cohort", # title text for the plot
  # package = "palettetown", #View(paletteer::palettes_d_names)
  # palette = "elekid",# package= "ggsci",#提取调色板所需的包
  # palette = "uniform_startrek", # choosing a different color palette
  # messages = FALSE
) +   theme(
  plot.title    = element_text(color = 'black', size   = 13, hjust = 0.5),
  plot.subtitle = element_text(color = 'black', size   = 13,hjust = 0.5),
  plot.caption  = element_text(color = 'black', size   = 13,face = 'italic', hjust = 1),
  axis.text.x   = element_text(color = 'black', size =10, angle = 0),
  axis.text.y   = element_text(color = 'black', size = 13, angle = 0),
  axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
  axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
  legend.title  = element_text(color = 'black', size  = 13),
  legend.text   = element_text(color = 'black', size   = 13),
  axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
  axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
  panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
) + theme(legend.position='none')

# getwd()
# file = paste0(Cyto_dir,'/Cytotrace_in_', organ,'.png')
# ggsave(height=6.5,width=9, filename=file, dpi = 500, device = "png")
# 

# sce_subphenot <- paste(score_result$Tissue,score_result$Cell_Type,sep = ':')#sce_sub$cluster
# sce_subphenot <- as.character(sce_subphenot)
# names(sce_subphenot) <- rownames(score_result)
# sce_subemb <- read.csv(paste0(cohort_directory,'Epi/epithelial_results/cell_embeddings.csv'))#colon_full@reductions$umap@cell.embeddings
# rownames(sce_subemb) <- sce_subemb[,1]
# sce_subemb <- sce_subemb[,-1]
# dim(colon_full)
# # source('~/R/ts860_database/CytoTRACE_fixed.R')
# library(CytoTRACE)
# plotCytoTRACE(results, phenotype = phe, emb = sce_subemb, outputDir = Cyto_dir)
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 10) ######------scRNA_Mutation_summary------######
###############################################################################################################################
data_path <- '/data/rluo4/All/Output/Cluster/'
Mut_path <- paste0('/data/rluo4/All/Output/sdata_ABN/')
setwd(Mut_path)
organ_all <- list.files(data_path)
Mut_scRNA <- vector("list", length(organ_all))
names(Mut_scRNA) <- organ_all
Mut_Trans <- NULL
Mut_genes <- NULL
for (organ in organ_all) {
  # organ = 'THCA'
  file = paste0(Mut_path, organ, '_sdata_ABN.obs.csv' )
  phe <- read.csv(file)
  phe$Gene_count <- apply(phe,1,function(x){
    return(length(unique(strsplit(x["Mut_gene"],', ')[[1]])))
  })
  table(phe$Gene_count==0)
  index = phe$Gene_count!=0
  phe$Mut_count.gene[index] <- phe$Mut_count[index]/phe$Gene_count[index]
  #####change the phenotype-Tissue in the sdata_ABN.obs file#####
  if (organ=='CRC'|organ=='Chen') {
    # typeCyto <- typeCyto[typeCyto %ni% 'ADJ']
    print(table(phe$Tissue))
    phe$tissue <- phe$Tissue
    phe$Tissue <- gsub('Adenocarcinoma','CRC',phe$Tissue)
    phe$Tissue <- gsub('Unaffected','ADJ',phe$Tissue)
    phe$Tissue <- gsub('NL','ADJ',phe$Tissue)
    phe$Tissue <- gsub('Polyp','FAP',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    table(phe$Tissue)
    # phe <- phe[phe$sample_name=='Becker',]
    # table(phe$Polyp.Type,phe$Tissue)
    # table(phe$Tissue[phe$Polyp.Type==""])
    # phe <- phe[phe$sample_name=='Chen',]
    # phe <- phe[phe$Tissue!='UNC',]
  }
  if (organ=='Pancreas') {
    PanIN <- unique(pheno_all$Pancreas_Epi$orig.ident[pheno_all$Pancreas_Epi$lesion=="Lesion"])
    PanIN <- na.omit(PanIN)
    print(table(phe$Tissue))
    phe$tissue <- phe$Tissue
    phe$Tissue[phe$orig.ident %in% PanIN] <- 'PanIN'
    phe$cell_type <- phe$Cell_Type
    print(table(phe$Tissue))
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Skin') {
    # typeCyto <- typeCyto[typeCyto %ni% 'ADJ']
    print(table(phe$Tissue))
    phe$tissue <- phe$Tissue
    phe$Tissue <- str_split(phe$Tissue,' ',simplify = T)[,1]
    phe$Tissue <- gsub('normal','ADJ',phe$Tissue)
    phe$Tissue <- gsub('Old','Healthy',phe$Tissue)
    phe$Tissue <- gsub('Young','Healthy',phe$Tissue)
    phe$Tissue <- gsub('CONTROL','Healthy',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    table(phe$Tissue)
    
  }
  
  if (organ=='THCA') {
    phe$tissue <- phe$Tissue
    phe$Tissue[grepl("withHT",rownames(phe))] <- 'HT'
    phe$Tissue <- gsub('Tumor','ATC',phe$Tissue)
    phe$Tissue <- gsub('AdultThyroid','Healthy',phe$Tissue)
    phe$Tissue[grepl("PTC",(phe$Tissue))] <- 'PTC'
    phe$Tissue <- gsub('NOR','ADJ',phe$Tissue)
    # phe$Tissue <- gsub('goiters','Goiters',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
  }
  if (organ=='HNSCC') {
    phe$tissue <- phe$Tissue
    phe$Tissue <- str_split(phe$Tissue,' ', simplify = T)[,1]
    phe$Tissue <- gsub('Primary','OSCC',phe$Tissue)
    phe$Tissue <- gsub('Leukoplakia','LP',phe$Tissue)
    phe$Tissue <- gsub('Metastatic','OSCC',phe$Tissue)
    
    phe$Tissue <- gsub('Non-tumoral','ADJ',phe$Tissue)
    phe$Tissue <- gsub('normal','Healthy',phe$Tissue)
    phe$Tissue <- gsub('GM','Healthy',phe$Tissue)
    phe$Tissue <- gsub('BM','Healthy',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Esophagus') {
    # typeCyto <- typeCyto[typeCyto %ni% 'ADJ']
    phe$tissue <- phe$Tissue
    phe$Tissue <- gsub('High-grade intraepithelial neoplasia of esophagus','HGIN',phe$Tissue)
    phe$Tissue[grepl("carcinoma",(phe$Tissue))] <- 'ESCC'
    phe$Tissue <- gsub('Low-grade intraepithelial neoplasia of esophagus','LGIN',phe$Tissue)
    
    phe$Tissue <- gsub('Non-tumor esophgeal epithelium','ADJ',phe$Tissue)
    phe$Tissue <- gsub('adjacent normal','ADJ',phe$Tissue)
    phe$Tissue <- gsub('AdultEsophagus','Healthy',phe$Tissue)
    
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Lung') {
    # typeCyto <- typeCyto[typeCyto %ni% 'ADJ']
    phe$tissue <- phe$Tissue
    phe$Tissue <- gsub('IAC','IA',phe$Tissue)
    phe$Tissue <- gsub('IA','IAC',phe$Tissue)
    phe$Tissue <- gsub('tLung','IAC',phe$Tissue)#Early
    phe$Tissue <- gsub('tL/B','IAC',phe$Tissue)#Advanced
    phe$Tissue <- gsub('nLung','ADJ',phe$Tissue)#Early
    phe$Tissue <- gsub('Normal','ADJ',phe$Tissue)#Advanced
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Breast') {
    # typeCyto <- typeCyto[typeCyto %ni% 'ADJ']
    phe$tissue <- phe$Tissue
    phe$Tissue <- gsub('Prophylatic Mastectomy','Precancer',phe$Tissue)
    phe$Tissue <- gsub('Tumor','IDC',phe$Tissue)
    phe$Tissue <- gsub('Contralateral to DCIS','ADJ',phe$Tissue)
    phe$Tissue <- gsub('Contralateral to IDC','ADJ',phe$Tissue)
    phe$Tissue <- gsub('normal breast','ADJ',phe$Tissue)
    phe$Tissue <- gsub('Reduction Mammoplasty','Healthy',phe$Tissue)
    phe$Tissue <- gsub('Normal','Healthy',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Prostate') {
    # typeCyto <- typeCyto[typeCyto %ni% 'ADJ']
    phe$tissue <- phe$Tissue
    phe$Tissue[phe$sample_name=='Chen_PRAD'] <- 'Tumor'
    phe$Tissue[phe$sample_name=='Joseph'&phe$Tissue=='Tumor'] <- 'BPH'
    phe$Tissue <- gsub('Normal','ADJ',phe$Tissue)
    # phe$Tissue <- gsub('Tumor','PRAD',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Cervix') {
    phe$tissue <- phe$Tissue
    phe$Tissue <- gsub('cervical cancer tissue','CC',phe$Tissue)
    phe$Tissue <- gsub('CA_HPV','CC',phe$Tissue)
    phe$Tissue <- gsub('neoplasm','CC',phe$Tissue)
    phe$Tissue <- gsub('high','HSIL_HPV',phe$Tissue)
    phe$Tissue <- gsub('metastatic','CC',phe$Tissue)
    phe$Tissue <- gsub('normal adjacent tissue','ADJ',phe$Tissue)
    phe$Tissue <- gsub('normal','ADJ',phe$Tissue)
    phe$Tissue <- gsub('AdultCervix','Healthy',phe$Tissue)
    phe$Tissue <- gsub('NO_HPV','Healthy',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Endometrium') {
    phe$tissue <- phe$Tissue
    phe$Tissue <- gsub('Tumor of the endometrium','EEC',phe$Tissue)
    phe$Tissue <- gsub('tr','EEC',phe$Tissue)
    phe$Tissue <- gsub('nr','ADJ',phe$Tissue)
    phe$Tissue <- gsub('Normal','ADJ',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='GC') {
    phe$tissue <- phe$Tissue
    phe$Tissue <- gsub('NAG','CSG',phe$Tissue)
    phe$Tissue <- gsub('Superficial','AGC',phe$Tissue)
    phe$Tissue <- gsub('Deep','AGC',phe$Tissue)
    phe$Tissue <- gsub('NAG','CSG',phe$Tissue)
    phe$Tissue <- gsub('EGC','AGC',phe$Tissue)
    phe$Tissue <- gsub('AGC','GC',phe$Tissue)
    phe$Tissue <- gsub('Normal','ADJ',phe$Tissue)
    phe$Tissue <- gsub('AdultStomach','Healthy',phe$Tissue)
    phe$cell_type <- phe$Cell_Type
    # write.csv(phe,'/data/rluo4/database/CRC/Output/Epi/Results/sdata.obs.csv')
  }
  if (organ=='Liver') {
    phe$Tissue <- gsub('fibrotic','cirrhotic',phe$Tissue)
    # phe$Tissue[phe$sample_name=='Ramachandran'] <- 'cirrhotic'
    phe$Tissue <- gsub('Benign','Cyst',phe$Tissue)
    phe$Tissue <- gsub('cirrhotic','Cirrhotic',phe$Tissue)
    phe$Tissue <- gsub('tumor','HCC',phe$Tissue)
    phe$Tissue <- gsub('Tumor core','HCC',phe$Tissue)
    phe$Tissue <- gsub('Tumor border','HCC',phe$Tissue)
    phe$Tissue <- gsub('healthy','Healthy',phe$Tissue)
    phe$Tissue <- gsub('AdultLiver','Healthy',phe$Tissue)
    phe$tissue <- phe$Tissue
    phe$cell_type <- phe$Cell_Type
    #  names(results$CytoTRACE) <- substr(names(results$CytoTRACE),1,nchar(names(results$CytoTRACE))-2)
    # phe <- phe[rownames(phe) %in% names(results$CytoTRACE),]
  }
  phe <- phe[! phe$sample_name %in% c('Cold'),]
  phe <- phe[! phe$orig.ident %in% c("p3", "p4", "p5"),]
  table(phe$Tissue)
  phe$Tissue[phe$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
  phe$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                              gsub('Tumor','PRAD', phe$Tissue)))
  phe <- phe[! phe$Tissue %in% c('UNC','ADJ','Healthy'),]
  # Mutants <- list(phe[,c('barcode','orig.ident','SimplifiedSampleName','Sample_Type','cnv_score','Mut_count','Mut_gene','AAChange')])
  Mutants <- phe[,c('cell_id','orig.ident','SimplifiedSampleName','Tissue','sample_name','Sample_Type','Cell_Type','cnv_leiden','cnv_score','cnv_status','Mut_count','Mut_gene','AAChange','Gene_count')]
  Mutants$organ = organ
  Mut_scRNA[[organ]] = Mutants
  # Mutants <- list(Mutants)
  # names(Mutants) <- paste0(organ,'_SCOMATIC')
  # Mut_Trans <- c(Mut_Trans, Mutants)
  Mut_G <- na.omit(phe$Mut_gene)
  Mut_genes <- c(Mut_genes, paste0(Mut_G[Mut_G!='']))
  ##################################
}
save(Mut_scRNA, Mut_genes, file=paste0('/data/rluo4/All/Output/Mut_scRNA.rds'))
# Muttation and Trans_genes
##############################
load(file=paste0('/data/rluo4/All/Output/Mut_scRNA.rds'))
Mut_G <- paste0(unlist(strsplit(Mut_genes,', ')))
unique(Mut_G)
grep('TP53',Mut_G)
# 造血干细胞衰老，克隆造血标志物
table(Mut_G=='DNMT3A')#IDH1
table(Mut_G=='ASXL1')
table(Mut_G=='TET2')

table(Mut_G=='KRAS')
table(Mut_G=='NPM1')# high freq
table(Mut_G=='CTCF')
Mut_df <- as.data.frame(table(Mut_G))

getwd()
Mut_CNV <- read.csv('/data/rluo4/All/Mut_CNV_Tirosh.txt',sep = '\t')
# Mut_CNV <- Mut_CNV[!duplicated(Mut_CNV$Event),]
Mutants_TCGA <- unique(Mut_CNV$Event)
Mutants_TCGA <- unique(str_split(Mutants_TCGA, ' ',simplify = T)[,1])
Mutants_TCGA <- Mutants_TCGA[!grepl("^[0-9]+", Mutants_TCGA)]
Mutants_TCGA <- Mutants_TCGA[Mutants_TCGA %ni% c('Xp','Xq')]

Mut_genes <- NULL
Mut_all <- NULL
for (organ in organ_all) {
  # organ = 'Cervix'
  mut = Mut_scRNA[[organ]]
  library(ggpubr)
  table(mut$Tissue)
  Mut_G <- mut[,c('cell_id','Cell_Type','Mut_count','Mut_gene','AAChange','Gene_count','Tissue')]
  Mut_G$Organ <- organ
  Mut_G <- Mut_G[! is.na(Mut_G$Mut_gene),]
  Mut_G <- Mut_G[Mut_G$Mut_gene!='',]
  Mut_all <<- rbind(Mut_all, Mut_G)
  # Mut_G <- Mut_G[Mut_G$Cell_Type=='STM', ]
  # Mut_G <- Mut_G[Mut_G$Tissue %in% precancer, ]
  Mut_genes <- c(Mut_genes, paste0(Mut_G$Mut_gene))
}
Mut_all$Index = 1:nrow(Mut_all)
Mut_all <- Mut_all[, c('Index', colnames(Mut_all)[1:(ncol(Mut_all)-1)])]
n <- Mut_all$Mut_gene
m <- strsplit( Mut_all$Mut_gene,', ')
write.csv(Mut_all[, c('Index', 'Organ', 'Tissue','Cell_Type','Mut_gene'#,'AAChange'
                      )], file = paste('/data/rluo4/All/Output/SCOMATIC.csv'), quote = F, row.names = F) #Table S2
write.xlsx(Mut_all[, c( 'Organ', 'Tissue','Cell_Type','Mut_gene')], file = paste0('/data/rluo4/All/Output/table_xlsx/TableS2.xlsx'), rowNames = FALSE)#,

g <- lapply(m, function(x){
  index = unique(x %in% c('CDKN2A'))
  # print(x)
  return(index)
})
# View(Mut_all[g!='FALSE',])
table(Mut_all[g!='FALSE', 'Tissue'])
table(Mut_all[g!='FALSE', 'Cell_Type'])

g <- lapply(m, function(x){
  index = unique(x %in% c('CTNNB1'))
  # print(x)
  return(index)
})
# View(Mut_all[g!='FALSE',])

Mut_G <- paste0(unlist(strsplit(Mut_genes,', ')))
unique(Mut_G)
table(Mutants_TCGA %in% unique(Mut_G))
Mutants_TCGA[Mutants_TCGA %in%  unique(Mut_G)]#28 for STM; 30 for all celltypes
table(is.na(Mut_G))
grep('IDH1',Mut_G)
table(Mut_G=='TP53')#
table(Mut_G=='RB1')#
table(Mut_G=='TP53I3')#
table(Mut_G=='CTNNB1')##SMAD4
table(Mut_G=='RUNX1')##SMAD4

table(Mut_G=='DNMT3A')#IDH1
table(Mut_G=='ASXL1')
table(Mut_G=='KRAS')
table(Mut_G=='NPM1')# high freq
table(Mut_G=='CTCF')
Mut_df <- as.data.frame(table(Mut_G))
View(Mut_df[Mut_df$Mut_G %in% Mutants_TCGA,])
# NCOR1 5066
# NPM1  4100
# CDKN2A 2705
high_mut <- Mut_df[Mut_df$Mut_G %in% Mutants_TCGA,]
high_mut <- high_mut$Mut_G[order(high_mut$Freq,decreasing = T)]
high_mut <- head(high_mut, 10)
g <- lapply(m, function(x){
  index = unique(x %in% high_mut)
  # print(x)
  return(index)
})
# View(Mut_all[g!='FALSE',])
table(Mut_all[g!='FALSE', 'Tissue'])
table(Mut_all[g!='FALSE', 'Cell_Type'])
precancer <- c('AAH','AD','AEH','AK',#'AIS',
               'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
               'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
               'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
               'SER','SIM','WIM')
Mut_all$TissueType = ifelse(Mut_all$Tissue %in% precancer, 'Precancer', 'Cancer')

# high_mut <- as.data.frame(table(Mut_all[g!='FALSE', 'Cell_Type'], Mut_all[g!='FALSE', 'Tissue']))
top100 <- table(Mut_all[g!='FALSE', 'Cell_Type'])
top100 <- names(top100[top100>=100])
high_mut <- as.data.frame(table(Mut_all[g!='FALSE', 'Cell_Type'], Mut_all[g!='FALSE', 'TissueType']))
colnames(high_mut)[1:2] <- c('Cell_Type', 'TissueType')
high_mut$Cell_Type <- as.character(high_mut$Cell_Type)
high_mut <- high_mut[high_mut$Cell_Type %in% top100,]
high_mut <- high_mut[high_mut$Freq>=1,]
high_mut <- high_mut[order(high_mut$Freq,decreasing = F),]
str(high_mut)
table(high_mut$TissueType)
palette=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
# color_cluster=c(RColorBrewer::brewer.pal(n = 9, name = 'Paired'))#,color_cluster)
library(ggpubr)
ggdotchart(high_mut, x = 'Cell_Type', y = 'Freq', color = 'TissueType', palette = c("#f89e81","#99a9cc"),
           sorting = 'asc', sort.by.groups = TRUE, add = 'segments', add.params = list(color = 'lightgray', size = 2),
           group = 'TissueType', dit.size = 6, ggtheme = theme_pubclean()) + font('x.text', size = 8, vjust = 0.5)
ggdotchart(high_mut, x = 'Cell_Type', y = 'Freq',
           color = "TissueType",                                # Color by groups
           palette = palette,#c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           # sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 3), # Change segment color and size
           group = "TissueType",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(high_mut$Freq,1),                        # Add mpg values as dot labels
           font.label = list(color = "brown", size = 10,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
) + geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +  ylab('Frequency of cells with top10 mutants') + xlab('') +
  theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 11, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 12, angle = 90),
    # legend.title  = element_text(color = 'black', size  = 20),
    legend.text   = element_text(color = 'black', size   = 14),
    legend.title = element_blank()
    # axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    # axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  ) 
f = paste0('/data/rluo4/All/Output/res_fig/Top10_mutated_celltype.png')
ggsave(height=4,width=5, filename=f, dpi = 500, device = "png")


g <- lapply(m, function(x){
  index = unique(x %in% c('TP53','RB1','TP53I3','CDKN2A'))
  # print(x)
  return(index)
})
# which(g!='FALSE')
# View(Mut_all[g!='FALSE',])
table(Mut_all[g!='FALSE', 'Tissue'])
table(Mut_all[g!='FALSE', 'Cell_Type'])
high_mut <- as.data.frame(table(Mut_all[g!='FALSE', 'Cell_Type'], Mut_all[g!='FALSE', 'TissueType']))
colnames(high_mut)[1:2] <- c('Cell_Type', 'TissueType')
high_mut$Cell_Type <- as.character(high_mut$Cell_Type)
high_mut <- high_mut[high_mut$Freq>=1,]
high_mut <- high_mut[order(high_mut$Freq,decreasing = F),]
str(high_mut)
table(high_mut$TissueType)
palette=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
# color_cluster=c(RColorBrewer::brewer.pal(n = 9, name = 'Paired'))#,color_cluster)
library(ggpubr)
ggdotchart(high_mut, x = 'Cell_Type', y = 'Freq',
           color = "TissueType",                                # Color by groups
           palette = palette,#c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           # sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 3), # Change segment color and size
           group = "TissueType",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(high_mut$Freq,1),                        # Add mpg values as dot labels
           font.label = list(color = "brown", size = 10,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
) + geom_hline(yintercept = 0, linetype = 2, color = "lightgray") + xlab('') + ylab('Frequency of TP53,RB1,CDKN2A,TP53I3-mut') +
  theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 11, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 12, angle = 90),
    # legend.title  = element_text(color = 'black', size  = 20),
    legend.text   = element_text(color = 'black', size   = 14),
    legend.title = element_blank()
    # axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    # axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  ) + theme(legend.position='none')
f = paste0('/data/rluo4/All/Output/res_fig/TP53_mutated_celltype.png')
ggsave(height=3.5,width=5, filename=f, dpi = 500, device = "png")

##############################
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(ggpubr)
library(ggthemes)
themes <- theme(
  plot.title    = element_text(color = 'black', size   = 20, hjust = 0.5),
  plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
  plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
  axis.text.x   = element_text(color = 'black', size = 20, angle = 0),
  axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
  axis.title.x  = element_text(color = 'black', size = 20, angle = 0),
  axis.title.y  = element_text(color = 'black', size = 20, angle = 90),
  legend.title  = element_text(color = 'black', size  = 20),
  legend.text   = element_text(color = 'black', size   = 20),
  axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
  axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
  panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
) + theme(legend.position='none')

data_path <- '/data/rluo4/All/Output/Cluster/'
organ_all <- list.files(data_path)
velo_sample_name <- read.table('/data/rluo4/All/Output/velo_cohort.txt')
velo_sample_name <- unique(velo_sample_name$V1)
cohort_directory = '/data/rluo4/All/Output/SComatic/'
for (organ in organ_all) {
  # organ = 'Liver'
  mut = Mut_scRNA[[organ]]
  library(ggpubr)
  library(rstatix)
  print(table(mut$Tissue))
  print(compare_means(Mut_count~Tissue, data=mut, method="wilcox.test",p.adjust.method="BH"))
  # compar<-list(c(unique(mut$Tissue)))#(c("HCC","Cirrhotic", "NAFLD"))
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  color_cluster <- color_cluster[1:length(unique(mut$Tissue))]
  names(color_cluster)=unique(mut$Tissue)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
  ylab<- "Number of mutations"#expression(paste("FCGR3A  ", "log"["2"], "(SCNA)"))
  file = paste0(organ, ' Diseased Tissues')
  file <- gsub('THCA','Thyroid',gsub('HNSCC','OralCavity', gsub('Chen', 'Colorectum',
                                                                gsub('GC','Stomach', gsub('CRC','Colorectum',file)))))
  
  if(organ == 'Liver'){
    mut$Tissue <- ordered(mut$Tissue,levels=c( "NAFLD","Cirrhotic","HCC"))
  }
  if(organ == 'Chen'){
    mut$Tissue <- ordered(mut$Tissue,levels=c( "SER",'AD',"MSS","MSI-H"))
  }
  if(organ == 'Skin'){
    mut$Tissue <- ordered(mut$Tissue,levels=c( "AK","SCCIS","cSCC"))
  }
  if(organ == 'Cervix'){
    mut$Tissue <- ordered(mut$Tissue,levels=c( "N_HPV","HSIL_HPV","CC"))
  }
  if(organ == 'HNSCC'){
    mut$Tissue <- ordered(mut$Tissue,levels=c( "NEOLP","EOLP","LP","OSCC"))
  }
  if(organ == 'THCA'){
    mut <- mut[ ! mut$Tissue %in% c('Goiters'),]
    mut$Tissue <- ordered(mut$Tissue,levels=c( "HT","PTC","ATC"))
  }
  
  wid = ifelse(length(unique(mut$Tissue))>3, length(unique(mut$Tissue))*2, length(unique(mut$Tissue))*2.5)
  p = ggstatsplot::ggbetweenstats(
    data = mut,
    x = Tissue,
    y = cnv_score,
    pairwise.display = 'all',
    p.adjust.method = 'bonferroni',
    notch = TRUE, # show notched box plot
    mean.plotting = FALSE, # whether mean for each group is to be displayed
    mean.ci = TRUE, # whether to display confidence interval for means
    mean.label.size = 2.55, # size of the label for mean
    type = "np", # which type of test is to be run
    k = 2, # number of decimal places for statistical results
    outlier.tagging = FALSE, # whether outliers need to be tagged
    xlab = "", # label for the x-axis variable
    ylab = "",#"Infercnv score",#expression(paste("FCGR3A  ", "log"["2"], "(SCNA)"))
    title = "",#file, #"Liver cohort", # title text for the plot
    palette = "Dark2",#"Darjeeling1", # choosing a different color palette
    messages = FALSE
  )+   themes
  f = paste0(cohort_directory, organ, '_cnv.png')
  ggsave(plot=p,height=wid*0.75,width=wid*0.9, filename=f, dpi = 500, device = "png")
  
  mut <- mut[mut$sample_name %in% velo_sample_name,]
  mut <- mut[! mut$sample_name %in% c('Ji'),]
  print(table(mut$sample_name))
  p = ggstatsplot::ggbetweenstats(
    data = mut,
    x = Tissue,
    y = Mut_count,
    pairwise.display = 'all',
    p.adjust.method = 'bonferroni',
    notch = TRUE, # show notched box plot
    mean.plotting = FALSE, # whether mean for each group is to be displayed
    mean.ci = TRUE, # whether to display confidence interval for means
    mean.label.size = 2.55, # size of the label for mean
    type = "np", # which type of test is to be run
    k = 2, # number of decimal places for statistical results
    outlier.tagging = FALSE, # whether outliers need to be tagged
    xlab = "",#"Sample Type", # label for the x-axis variable
    ylab = "",# ylab = "No. of Mutants", # label for the y-axis variable
    title = '', #"Liver cohort", # title text for the plot
    package = "palettetown", #View(paletteer::palettes_d_names)
    palette = "phanpy",# package= "ggsci",#提取调色板所需的包
    messages = FALSE
  ) +   themes
  f = paste0(cohort_directory, organ, '_mut.png')
  ggsave(plot=p, height=wid*0.75, width=wid*0.9, filename=f, dpi = 500, device = "png")
}

mut = Mut_scRNA[[organ]]
library(ggpubr)
library(rstatix)
print(table(mut$Cell_Type))
mut$cell_type <- paste(mut$Tissue, mut$Cell_Type, sep = ': ')
print(table(mut$cell_type))
mut <- mut[! mut$Cell_Type %in% 'CHO',]
mut <- mut[ mut$Tissue %in% c('Cirrhotic','HCC'),]
print(compare_means(Mut_count~cell_type, data=mut, method="wilcox.test",p.adjust.method="BH"))
file = paste0(organ, ' Diseased Tissues')
file <- gsub('THCA','Thyroid',gsub('HNSCC','OralCavity', gsub('Chen', 'Colorectum',
                                                              gsub('GC','Stomach', gsub('CRC','Colorectum',file)))))

p = ggstatsplot::ggwithinstats(#ggbetweenstats(
  data = mut,
  x = cell_type,
  y = cnv_score,
  sort = "descending", # ordering groups along the x-axis based on
  sort.fun = median, # values of `y` variable
  pairwise.comparisons = FALSE,
  # pairwise.display = "s",
  pairwise.annotation = "p",
  # caption = "Data from: iris",
  # ggtheme = ggthemes::theme_fivethirtyeight(),
  ggstatsplot.layer = FALSE,
  # pairwise.display = 'all',
  # p.adjust.method = 'BH',#'bonferroni',
  notch = TRUE, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = TRUE, # whether to display confidence interval for means
  mean.label.size = 2.55, # size of the label for mean
  # type = "np", # which type of test is to be run
  # k = 2, # number of decimal places for statistical results
  # outlier.tagging = FALSE, # whether outliers need to be tagged
  xlab = "Sample:Cell Type", # label for the x-axis variable
  ylab = "Inferred CNV Score", # label for the y-axis variable
  # title = "Liver cohort", # title text for the plot
  package = "pals", #View(paletteer::palettes_d_names)
  palette = "polychrome",# package= "ggsci",#提取调色板所需的包
  # messages = FALSE
) + themes
f = paste0(cohort_directory, organ, '_cell_type_cnv.png')
ggsave(plot=p,height=5.5,width=9, filename=f, dpi = 500, device = "png")

mut <- mut[mut$sample_name %in% velo_sample_name,]
mut <- mut[! mut$sample_name %in% c('Ji'),]
print(table(mut$sample_name))
p = ggstatsplot::ggwithinstats(#ggbetweenstats(
  data = mut,
  x = cell_type,
  y = Mut_count,
  sort = "descending", # ordering groups along the x-axis based on
  sort.fun = median, # values of `y` variable
  pairwise.comparisons = FALSE,
  # pairwise.display = "s",
  pairwise.annotation = "p",
  # caption = "Data from: iris",
  # ggtheme = ggthemes::theme_fivethirtyeight(),
  ggstatsplot.layer = FALSE,
  # pairwise.display = 'all',
  # p.adjust.method = 'BH',#'bonferroni',
  notch = TRUE, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = TRUE, # whether to display confidence interval for means
  mean.label.size = 2.55, # size of the label for mean
  # type = "np", # which type of test is to be run
  # k = 2, # number of decimal places for statistical results
  # outlier.tagging = FALSE, # whether outliers need to be tagged
  xlab = "Sample:Cell Type", # label for the x-axis variable
  ylab = "No. of mutants", # label for the y-axis variable
  # title = "Liver cohort", # title text for the plot
  package = "pals", #View(paletteer::palettes_d_names)
  palette = "polychrome",# package= "ggsci",#提取调色板所需的包
  # messages = FALSE
) +   themes
getwd()
f = paste0(cohort_directory, organ, '_cell_type_mut.png')
ggsave(plot=p,height=5.5,width=9, filename=f, dpi = 500, device = "png")

# mkdir ~/lrh/All/Output/sdata_TME
# ls -d *|while read id;
# do echo $id;
# cp -p $id/Output/Imm/Results/sdata.obs.csv ../All/Output/sdata_TME/${id}_sdata_Imm.obs.csv;
# cp -p $id/Output/Str/Results/sdata.obs.csv ../All/Output/sdata_TME/${id}_sdata_Str.obs.csv;
# done
# ls -d */*/velocyto*/|while read id;
# do
# sample_name=`echo $id | awk -F '/' '{print $2}'`
# echo $sample_name
# done
# ls -d *|while read id;
# do echo $id;
# cp -p $id/Output/Epi/epithelial_results/*DEG.RData ../All/Output/sdata_ABN/${id}_DEG.RData;
# done

library(openxlsx)
CD8T_Signature <- read.xlsx('/data/rluo4/All/CD8T_Signature.xlsx')#,1)
colnames(CD8T_Signature) <- gsub('/','.', gsub(':','.',colnames(CD8T_Signature)))
colnames(CD8T_Signature) <- paste0('CD8-', colnames(CD8T_Signature))
CD4T_Signature <- read.xlsx('/data/rluo4/All/CD4T_Signature.xlsx')#,1)
colnames(CD4T_Signature) <- gsub('/','.', gsub(':','.',colnames(CD4T_Signature)))
colnames(CD4T_Signature) <- paste0('CD4-', colnames(CD4T_Signature))
Myeloid_Signature <- read.xlsx('/data/rluo4/All/Myeloid_Signature.xlsx')#,1)
NK_Signature <- read.xlsx('/data/rluo4/All/NK_Signature.xlsx')#,1)
# colnames(Myeloid_Signature) <- paste0('Myeloid-', colnames(Myeloid_Signature))
CAF_Signature <- read.xlsx('/data/rluo4/All/CAF_Signature.xlsx')#,1)
NeoTCR <- read.csv('/data/rluo4/All/NeoTCR.txt', sep = '\t')
# genesets <- subset(NeoTCR, select = c('NeoTCR8','NeoTCR4')) %>% as.data.frame()
NeoTCR4 <- data.frame(term = rep('CD4T-NeoTCR', length(NeoTCR$NeoTCR4[NeoTCR$NeoTCR4 !=''])) 
                      ,gene =NeoTCR$NeoTCR4[NeoTCR$NeoTCR4 !=''])
NeoTCR8 <- data.frame(term = rep('CD8T-NeoTCR', length(NeoTCR$NeoTCR8[NeoTCR$NeoTCR8 !=''])) 
                      ,gene =NeoTCR$NeoTCR8[NeoTCR$NeoTCR8 !=''])

sets <- c('CD8T', 'CD4T', 'Myeloid', 'NK', 'CAF')
Curated <- vector("list", length(sets))
names(Curated) <- sets
for(i in sets){
  Signature <- read.xlsx(paste0('/data/rluo4/All/', i , '_Signature.xlsx'))#, sheetIndex = 1)
  if(i %in% c('CD4T','CD8T')){
    colnames(Signature) <- gsub('/','.', gsub(':','.',colnames(Signature)))
    colnames(Signature) <- paste0(i,'-', colnames(Signature))
  }
  # Concatenate all values into a single vector
  all_values <- as.vector(t(t(Signature)))
  # Create a vector of corresponding column names
  all_column_names <- rep(colnames(Signature), each = nrow(Signature))
  # Combine into a new dataframe
  curated_signature <- data.frame(
    term = all_column_names,
    gene = all_values
  )
  curated_signature <- curated_signature[! is.na(curated_signature$gene),]
  
  if(i == 'CD4T'){
    curated_signature <- rbind(curated_signature, NeoTCR4)
  }
  if(i == 'CD8T'){
    curated_signature <- rbind(curated_signature, NeoTCR8)
  }
  curated_signature$term <- as.factor(curated_signature$term)
  
  Curated[[i]] <- curated_signature
}

table(curated_signature$term)
###############################---load results of GSVA-TME.R---###########################
load('/data/rluo4/All/Output/Imm_all.rds')
load('/data/rluo4/All/Output/Str_all.rds')
load('/data/rluo4/All/Output/pheno_all.rds')
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(stringr)
data_path = '/data/rluo4/All/Output'
setwd(data_path)
organ_all <- list.files('Cluster/')
organ_all <- organ_all[organ_all!='Chen']
organ_all
cell_all <- c(pheno_all, Imm_all, Str_all)
pheno <- paste(organ_all, 'Epi', sep = '_')
# All_ratio <- NULL
# cell_summary <- NULL
table(curated_signature$term)
Lineage <- read.csv('/data//rluo4/All/Output/organ13-celltypes.csv')
Tlym <- Lineage$Minor_type[Lineage$Major_type %in% 'Tcell']
CD8T <- Tlym[Tlym %in% c('CD8TRM', 'CD8TEXINT', 'CD8TCM', 'CD8TEREX', 'CD8TEFF', 'CD8TEXP', 'CD8TEX', 'CD8TN', 'MAIT')]
print(CD8T)
CD4T <- Tlym[! Tlym %in% c('CD8TRM', 'CD8TEXINT', 'CD8TCM', 'CD8TEREX', 'CD8TEFF', 'CD8TEXP', 'CD8TEX', 'CD8TN', 'MAIT', 'ILC', 'NK', 'GDT', 'NKT')]
print(CD4T)
Myeloid <- Lineage$Minor_type[Lineage$Major_type %in% 'Myeloid']
NK <- Lineage$Minor_type[Lineage$Major_type %in% 'NK']
CAF <- Lineage$Minor_type[Lineage$Major_type %in% 'Mesenchymal']
NK = 'NK'
cell_sets <- data.frame(Minor_type = c(CD8T, CD4T, Myeloid, NK, CAF),
                        Major_type = c(rep('CD8T', length(CD8T)),rep('CD4T', length(CD4T)),
                                       rep('Myeloid', length(Myeloid)),rep('NK', length(NK)),
                                       rep('CAF', length(CAF))
                        ))
pheno <- c(pheno, 'Chen_Epi')
# AUcell algorithm
# for (obs in pheno){
#   organ <- str_split(obs,'_', simplify = T)[,1]
#   print(organ)
#   # organ = 'Breast'
#   epi <- gsub('Chen','CRC',paste(organ, 'Epi', sep = "_"))
#   str <-  gsub('Chen','CRC',paste(organ, 'Str', sep = "_"))
#   imm <- gsub('Chen','CRC',paste(organ, 'Imm', sep = "_"))
#   
#   Epi <- cell_all[[epi]]
#   Epi$barcode = rownames(Epi)
#   Epi$cell_type <- Epi$Cell_Type
#   table(Epi$Organ,Epi$Tissue)
#   Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
#   Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
#                                               gsub('Tumor','PRAD', Epi$Tissue)))
#   Epi$Major_type <- 'Epithelial'
#   print(table(Epi$Organ))
#   Imm <- cell_all[[imm]]
#   Imm$Major_type <- 'Immune'
#   
#   Str <- cell_all[[str]]
#   Str$Major_type <- 'Stromal'
#   
#   data_path <- '/data/rluo4/All/Output/Enrichment/'
#   analysis_parent_folder <- paste0(data_path,organ,'_Enrichment')
#   analysis_parent_folder <- gsub('Chen','CRC',analysis_parent_folder)
#   if (!dir.exists(paste0(analysis_parent_folder))){
#     dir.create(paste0(analysis_parent_folder))
#   }
#   range = 4:12 # 迭代次数
#   gmin = 5 #最低关联数量
#   ncores = 5 #运行核数量
#   sets = unique(cell_sets$Major_type)
#   # sets <- c('CD8T','CD4T','Myeloid', 'NK', 'CAF')
#   for(s in sets){
#     # s = 'CAF'
#     setwd(analysis_parent_folder)
#     Func <- Curated[[s]]
#     # CD8_example <- data.frame(
#     #   features = c(colnames(CD8T_Signature)[1:3],'CD8-NeoTCR', colnames(CD8T_Signature)[4:19]),
#     #   functions = c(rep('Differentiation',3),rep('Function',12),rep('Metabolism',3),rep('Apoptosis',2))
#     # )
#     # CD8_example$features <- sub(".*?-","", CD8_example$features)#str_split(Func$features, '-', simplify = T)[,2]
#     # paste(CD8_example$features, collapse = "', '")
#     colnames(Func) <- c('features', 'functions')
#     Func$features <- as.character( Func$features)
#     Func$functions <- Func$features
#     Func <- dplyr::distinct(Func)
#     if(s %in% c('CD4T','CD8T')){
#       Func$features <- sub(".*?-","", Func$features)#str_split(Func$features, '-', simplify = T)[,2]
#     }
#     Differentiation <- c('Naïve', 'Activation.Effector.function', 'Exhaustion')
#     Function <-   c('NeoTCR', 'TCR.Signaling','TCR.signaling', 'Cytotoxicity', 'Cytokine.Cytokine.receptor', 'Chemokine.Chemokine.receptor', 'Senescence', 'Anergy', 'NFKB.Signaling', 'Stress.response', 'MAPK.Signaling', 'Adhesion', 'IFN.Response','IFN.response','Treg.signature','Costimulatory.molecules',
#                     'TaNK','Inflamatoty','Stress','HLA-dependent.inhibitory.receptor','HLA-independent.inhibitory.receptor','HLA-dependent.activating.receptor','HLA-independent.activating.receptor')
#     Metabolism <- c('Oxidative.phosphorylation', 'Glycolysis', 'Fatty.acid.metabolism','Lipid.metabolism','OXPHOS')
#     Apoptosis <- c('Pro-apoptosis', 'Anti-apoptosis')
#     
#     Func$functions[Func$features %in% Differentiation] <- 'Differentiation'
#     Func$functions[Func$features %in% Function] <- 'Function'
#     Func$functions[Func$features %in% Metabolism] <- 'Metabolism'
#     Func$functions[Func$features %in% Apoptosis] <- 'Apoptosis'
#     
#     Func$functions <- sub(".*?-","", Func$functions)#str_split(Func$features, '-', simplify = T)[,2]
#     Func$features <- sub("Inflammototy","Inflammotory", Func$features)#str_split(Func$features, '-', simplify = T)[,2]
#     
#     Func <- Func[ ! Func$features %in% c('NeoTCR','Senescence'), ]
#     print(s)
#     outfile_gsva <- paste0(organ, '_', s,'_ssgsea.rds')
#     outfile_gsva <- gsub('Chen', 'CRC', outfile_gsva)
#     
#     if (! file.exists(outfile_gsva)) {
#       print(paste0('There is no ', s, ' in ', organ,  " , so skip !"))
#       next;
#     }
#     
#     if (file.exists(outfile_gsva)) {
#       print(paste0('start gssea: ', s, " is already over!"))
#       # next;
#       load(outfile_gsva)
#     }
#     
#     outfile_auc <- paste0(organ, '_', s,'_aucell.rds')
#     outfile_auc <- gsub('Chen', 'CRC', outfile_auc)
#     
#     if (file.exists(outfile_auc)) {
#       print(paste0('start aucell: ', s, " is already over!"))
#       # next;
#       load(outfile_auc)
#     }
#     
#     if(s=='CAF'){
#       TMECell_obj <- Str[Str$cell_type %in% cell_sets$Minor_type[cell_sets$Major_type==s],]
#     } else{
#       TMECell_obj <- Imm[Imm$cell_type %in% cell_sets$Minor_type[cell_sets$Major_type==s],]
#     }   
#     print(table(colnames(AUC) == TMECell_obj$barcode))
#     AUC_colnames <- paste(TMECell_obj$SimplifiedSampleName,colnames(AUC), sep=':')
#     
#     TMECell_obj <- TMECell_obj[! TMECell_obj$sample_name %in% c('Cold'),]
#     TMECell_obj <- TMECell_obj[! TMECell_obj$orig.ident %in% c("p3", "p4", "p5"),]
#     # TMECell_obj$Organ <- unique(Epi$Organ)
#     table(TMECell_obj$Tissue)
#     TMECell_obj <- TMECell_obj[! TMECell_obj$Tissue %in% c('UNC','ADJ','Healthy'),]
#     print(table(TMECell_obj$Tissue))
#     precancer <- c('AAH','AD','AEH','AK',#'AIS',
#                    'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
#                    'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
#                    'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
#                    'SER','SIM','WIM')
#     samples <- unique(TMECell_obj$orig.ident)
#     print(samples)
#     setdiff(samples, unique(Epi$orig.ident))
#     setdiff(unique(Epi$orig.ident), samples)
#     table(TMECell_obj$orig.ident %in% Epi$orig.ident)
#     TMECell_obj <- TMECell_obj[TMECell_obj$orig.ident %in% Epi$orig.ident,]
#     sample.meta <- Epi[,c("SimplifiedSampleName","Tissue",'orig.ident')]
#     sample.meta <- distinct(sample.meta)
#     TMECell_obj$SimplifiedSampleName <- sample.meta$SimplifiedSampleName[match(TMECell_obj$orig.ident, sample.meta$orig.ident)]
#     table(TMECell_obj$SimplifiedSampleName == TMECell_obj$orig.ident)
#     TMECell_obj$Molecular.typing[TMECell_obj$Tissue %in% precancer] <- 'Precancer'
#     TMECell_obj$Molecular.typing[! TMECell_obj$Tissue %in% precancer] <- 'Cancer'
#     table(TMECell_obj$Molecular.typing)
#     table(TMECell_obj$SimplifiedSampleName)
#     print(table(colnames(AUC) %in% TMECell_obj$barcode))
#     # AUC <- AUC[,colnames(AUC) %in% TMECell_obj$barcode]
#     print(table(colnames(AUC) == TMECell_obj$barcode))
#     # AUC_colnames <- paste(TMECell_obj$SimplifiedSampleName,colnames(AUC), sep=':')
#     
#     TMEcells <- data.frame(table(TMECell_obj$SimplifiedSampleName))
#     TMEcells <- TMEcells$Var1[TMEcells$Freq>=10]
#     
#     TMECell_obj <- TMECell_obj[TMECell_obj$SimplifiedSampleName %in% TMEcells,]
#     # AUC <- ssgsea.res
#     # library(AUCell, lib.loc = "/data/rluo4/lorihan/R/site-library")
#     sample.meta <- unique(TMECell_obj$SimplifiedSampleName)
#     # AUC_colnames <- paste(TMECell_obj$SimplifiedSampleName,colnames(AUC), sep=':')
#     for(i in 1:length(sample.meta)){#   
#       sample.data <- AUC[,grep(sample.meta[i],AUC_colnames)]
#       # Print the dimensions of sample.data
#       # cat("Iteration:", i, "Dimensions of sample.data:", dim(sample.data), "\n")
#       medians <- apply(sample.data,1,median)
#       if(i == 1){
#         AUC.sample <- data.frame(medians)
#         colnames(AUC.sample)[i] <- sample.meta[i]
#       } else {
#         AUC.sample <- cbind(AUC.sample, data.frame(medians))
#         colnames(AUC.sample)[i] <- sample.meta[i]
#       }
#     }
#     sample.meta <- TMECell_obj[,c("SimplifiedSampleName","Tissue",'orig.ident')]
#     sample.meta <- distinct(sample.meta)
#     # obs <- paste(organ,'Epi',sep = '_')
#     # if(organ == 'Chen'){
#     #   obs <- gsub('Chen','CRC',obs)
#     # }
#     cytotrace.meta <- unique(Epi$SimplifiedSampleName)
#     cytotrace.meta <- cytotrace.meta[cytotrace.meta %in% sample.meta$SimplifiedSampleName]
#     Cyto_colnames <- paste(Epi$SimplifiedSampleName,rownames(Epi), sep=':')
#     # malignant.meta <- unique(Epi$SimplifiedSampleName)
#     # malignant.meta <- malignant.meta[malignant.meta %in% sample.meta$SimplifiedSampleName]
#     for(i in 1:length(cytotrace.meta)){#     sample.data <- Cyto[,grep(cytotrace.meta[i],colnames(Cyto))]
#       sample.data <- Epi[grep(cytotrace.meta[i],Cyto_colnames),]
#       sample.data <- sample.data['CytoTRACE',drop=FALSE]
#       medians <- apply(sample.data,2,median)
#       if(i == 1){
#         Cyto.sample <- data.frame(medians)
#         colnames(Cyto.sample)[i] <- cytotrace.meta[i]
#       } else {
#         Cyto.sample <- cbind(Cyto.sample, data.frame(medians))
#         colnames(Cyto.sample)[i] <- cytotrace.meta[i]
#       }
#     }
#     table(colnames(AUC.sample) %in% colnames(Cyto.sample))
#     AUC.sample <- AUC.sample[,colnames(AUC.sample) %in% colnames(Cyto.sample)]
#     dim(AUC.sample)#20 32
#     TME <- AUC.sample
#     # TME <- rbind(AUC.sample, Cyto.sample)
#     colnames(TME) 
#     indir <- '/data/rluo4/All/Output/sdata_ABN/'
#     pc_file = paste0(indir, organ,'_DEG.RData')
#     if(organ !='Chen'){
#       load(pc_file)
#       if(organ=='CRC'){
#         pc_df <-  pc_df_Becker
#         Patt <- Patt_Becker
#       }
#     } else{
#       pc_file <- gsub('Chen','CRC',pc_file)
#       load(pc_file)
#       pc_df <- pc_df_Chen
#       Patt <- Patt_Chen
#     }
#     
#     unique(pc_df$Sample)
#     table(pc_df$Tissue)
#     head(Patt)
#     table(pc_df$SimplifiedSampleName %in% colnames(TME))
#     STM_df <- pc_df[pc_df$SimplifiedSampleName %in% colnames(TME),]
#     STM_df$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
#                                                    gsub('Tumor','PRAD', STM_df$Tissue)))
#     
#     order_samples <- STM_df$SimplifiedSampleName
#     TME["Malignancy",] <- STM_df$nearest_spline_x_vals[match(colnames(TME),order_samples)]
#     TME <- TME[, ! is.na(TME["Malignancy",])]
#     table(colnames(TME) %in% order_samples)
#     
#     TME <- TME[, match(order_samples, colnames(TME))]
#     # save(Cor_MP_Malignancy, file=paste0(data_path,"../Cor_MP_Malignancy_new.rds"))
#     table(STM_df$Tissue)
#     color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
#     # color_cluster=c(RColorBrewer::brewer.pal(n = 9, name = 'Paired'))#,color_cluster)
#     color_cluster <- color_cluster[1:length(unique(STM_df$Tissue))]
#     names(color_cluster)=unique(STM_df$Tissue)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
#     col = list(Disease = color_cluster)#c("AK" = "#E41A1C","cSCC" = "#377EB8", "SCCIS" = "#984EA3"))
#     
#     setwd(analysis_parent_folder)
#     cell_folder <- 'res_fig'
#     if (!dir.exists(paste0(cell_folder))){
#       dir.create(paste0(cell_folder))
#     }
#     setwd(cell_folder)
#     getwd()
#     enrich_df <- as.data.frame(t(TME[-nrow(TME),]))
#     enrich_df$SimplifiedSampleName <- rownames(enrich_df)
#     STM_df <- left_join(STM_df[, c("SimplifiedSampleName",'nearest_spline_x_vals','Tissue','sample_name')], enrich_df,by='SimplifiedSampleName')
#     tcolors <- color_cluster#paletteDiscrete(values = unique(STM_df$Tissue))
#     n <- colnames(STM_df)[-(1:4)]#rownames(test.diff)
#     p <- NULL
#     # plot <-  lapply(setNames(n, n), function(nameindex) {
#     #   cell <- nameindex
#     #   # cell_p <- pc_df[match(colnames(test.diff),rownames(pc_df)),]
#     #   # cell_p$cell <- t(test.diff[cell,])
#     #   ylab <- paste0('AUCell Scores')
#     #   p = ggplot(STM_df, aes(x=nearest_spline_x_vals, y= STM_df[,cell], color=Tissue)) +
#     #     geom_point(size=2) + scale_color_manual(values=tcolors) + ggtitle(cell) +
#     #     xlab("Malignancy Continuum")+ylab(ylab)+theme(
#     #       plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
#     #       plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
#     #       plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
#     #       axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
#     #       axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
#     #       axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
#     #       axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
#     #       legend.title  = element_text(color = 'black', size  = 16),
#     #       legend.text   = element_text(color = 'black', size   = 16),
#     #       axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
#     #       axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
#     #       #panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
#     #     )
#     #   p
#     #   f <- paste0(cell,"_dynamics.png")
#     #   ggsave(plot=p,height=4,width=6.5, filename=f, dpi = 300, device = "png")
#     #   return(p)
#     # })
#     
#     library(pheatmap) # I like esoteric packages!
#     library(RColorBrewer)
#     library(iheatmapr)
#     library(datasets)
#     library(reshape2)
#     library(gplots)
#     library(ComplexHeatmap)
#     library(grid)
#     library(circlize)
#     # colnames(CD8T_Signature)
#     rownames(AUC)
#     
#     TME.matrix <- as.matrix(TME[-nrow(TME),])
#     if(s %in% c('CD4T','CD8T')){
#       rownames(TME.matrix) <-  sub(".*?-","", rownames(TME.matrix))
#     }
#     table(rownames(TME.matrix) %in% Func$features)
#     TME.matrix <- TME.matrix[Func$features, ]
#     ha1 <- HeatmapAnnotation(Disease=STM_df$Tissue,
#                              col = col,show_legend = TRUE)
#     # Create the heatmap annotation
#     options(bitmapType = 'cairo')
#     
#     set.seed(123)
#     par(mfrow = c(1, 1))
#     max <- round(max(TME.matrix),1)
#     max
#     min <- round(min(TME.matrix),1)
#     f = paste0(organ, '_',s,"_heatmap.png")
#     # png(height=hei,width=8, filename=f, units = "in", res = 400)
#     TME.matrix <- TME.matrix[!rownames(TME.matrix) %in% c('NeoTCR','Senescence'),]
#     p = Heatmap(TME.matrix,show_row_names = TRUE,show_column_names = FALSE,cluster_rows = FALSE,cluster_columns = FALSE,
#                 col=colorRamp2(c(min, (max(TME.matrix)-min(TME.matrix))/2, max),c("seagreen4", "lightyellow", "orange")), name = 'AUCell',
#                 row_split = Func[,ncol(Func)], #top_annotation = ha1, # rep(c("A","B","C"),each=6),
#                 #column_split = rep(c("A","B","C","D"),each=6),row_title = NULL,column_title = NULL
#                 clustering_method_columns = "complete",show_column_dend = TRUE,show_row_dend = F,row_names_gp = gpar(fontsize = 11))
#     # draw(p, heatmap_legend_side = "left", annotation_legend_side = "left", legend_grouping = "original")
#     # ggsave(height=8,width=8, filename=f, dpi = 400)#, device = "png")
#     # # Save as PNG
#     # f = paste0(organ, '_',s,"_heatmap.pdf")
#     # pdf(f, width = 10, height = 5)
#     hei <- ifelse(length(Func$features)<8, length(Func$features)*0.26,length(Func$features)*0.16)
#     wid <- length(unique(pc_df$SimplifiedSampleName))
#     wid <- ifelse(wid >30, wid*0.125,
#                   ifelse(wid<20, wid*0.3, wid*0.2))
# 
#     png(height=hei,width=wid, filename=f, units = "in", res = 500)
#     draw(p, heatmap_legend_side = "left", annotation_legend_side = "top", legend_grouping = "original")
#     dev.off()
#   }
# }
# 
# 

# Patt <- read.table('/data/rluo4/All/Output/organ13-deg.txt',sep = '\t',fill = TRUE,header = TRUE)
# GSVA algorithm
# for (obs in pheno[-5]){
for (obs in c('Liver_Epi','Cervix_Epi')){
  
  organ <- str_split(obs,'_', simplify = T)[,1]
  print(organ)
  # organ = 'Breast'
  epi <- gsub('Chen','CRC',paste(organ, 'Epi', sep = "_"))
  str <-  gsub('Chen','CRC',paste(organ, 'Str', sep = "_"))
  imm <- gsub('Chen','CRC',paste(organ, 'Imm', sep = "_"))
  
  Epi <- cell_all[[epi]]
  Epi$barcode = rownames(Epi)
  Epi$cell_type <- Epi$Cell_Type
  table(Epi$Organ,Epi$Tissue)
  Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
  Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                              gsub('Tumor','PRAD', Epi$Tissue)))
  Epi$Major_type <- 'Epithelial'
  print(table(Epi$Organ))
  Imm <- cell_all[[imm]]
  Imm$Major_type <- 'Immune'
  
  Str <- cell_all[[str]]
  Str$Major_type <- 'Stromal'
  
  data_path <- '/data/rluo4/All/Output/Enrichment/'
  analysis_parent_folder <- paste0(data_path,organ,'_Enrichment')
  analysis_parent_folder <- gsub('Chen','CRC',analysis_parent_folder)
  if (!dir.exists(paste0(analysis_parent_folder))){
    dir.create(paste0(analysis_parent_folder))
  }
  range = 4:12 # 迭代次数
  gmin = 5 #最低关联数量
  ncores = 5 #运行核数量
  sets = unique(cell_sets$Major_type)
  for(s in sets){
    # s = 'CAF'
    setwd(analysis_parent_folder)
    Func <- Curated[[s]]
    # CD8_example <- data.frame(
    #   features = c(colnames(CD8T_Signature)[1:3],'CD8-NeoTCR', colnames(CD8T_Signature)[4:19]),
    #   functions = c(rep('Differentiation',3),rep('Function',12),rep('Metabolism',3),rep('Apoptosis',2))
    # )
    # CD8_example$features <- sub(".*?-","", CD8_example$features)#str_split(Func$features, '-', simplify = T)[,2]
    # paste(CD8_example$features, collapse = "', '")
    colnames(Func) <- c('features', 'functions')
    Func$features <- as.character( Func$features)
    Func$functions <- Func$features
    Func <- dplyr::distinct(Func)
    if(s %in% c('CD4T','CD8T')){
      Func$features <- sub(".*?-","", Func$features)#str_split(Func$features, '-', simplify = T)[,2]
    }
    Differentiation <- c('Naïve', 'Activation.Effector.function', 'Exhaustion')
    Function <-   c('NeoTCR', 'TCR.Signaling','TCR.signaling', 'Cytotoxicity', 'Cytokine.Cytokine.receptor', 'Chemokine.Chemokine.receptor', 'Senescence', 'Anergy', 'NFKB.Signaling', 'Stress.response', 'MAPK.Signaling', 'Adhesion', 'IFN.Response','IFN.response','Treg.signature','Costimulatory.molecules',
                    'TaNK','Inflamatoty','Stress','HLA-dependent.inhibitory.receptor','HLA-independent.inhibitory.receptor','HLA-dependent.activating.receptor','HLA-independent.activating.receptor')
    Metabolism <- c('Oxidative.phosphorylation', 'Glycolysis', 'Fatty.acid.metabolism','Lipid.metabolism','OXPHOS')
    Apoptosis <- c('Pro-apoptosis', 'Anti-apoptosis')
    
    Func$functions[Func$features %in% Differentiation] <- 'Differentiation'
    Func$functions[Func$features %in% Function] <- 'Function'
    Func$functions[Func$features %in% Metabolism] <- 'Metabolism'
    Func$functions[Func$features %in% Apoptosis] <- 'Apoptosis'
    
    Func$functions <- sub(".*?-","", Func$functions)#str_split(Func$features, '-', simplify = T)[,2]
    Func$features <- sub("Inflammototy","Inflammotory", Func$features)#str_split(Func$features, '-', simplify = T)[,2]
    
    Func <- Func[ ! Func$features %in% c('NeoTCR','Senescence'), ]
    Func$features <- gsub('[.]','-',Func$features)
    print(s)
    outfile_gsva <- paste0(organ, '_', s,'_ssgsea.rds')
    outfile_gsva <- gsub('Chen', 'CRC', outfile_gsva)
    
    if (! file.exists(outfile_gsva)) {
      print(paste0('There is no ', s, ' in ', organ,  " , so skip !"))
      next;
    }
    
    if (file.exists(outfile_gsva)) {
      print(paste0('start gssea: ', s, " is already over!"))
      # next;
      load(outfile_gsva)
    }
    
    outfile_auc <- paste0(organ, '_', s,'_aucell.rds')
    outfile_auc <- gsub('Chen', 'CRC', outfile_auc)
    
    if (file.exists(outfile_auc)) {
      print(paste0('start aucell: ', s, " is already over!"))
      # next;
      load(outfile_auc)
    }
    ###### two options: ######
    for( method_use in c( 'AUCell','GSVA')){
      if(method_use=='GSVA'){
        AUC <- ssgsea.res # use GSVA algorithm
      }
      if(s=='CAF'){
        TMECell_obj <- Str[Str$cell_type %in% cell_sets$Minor_type[cell_sets$Major_type==s],]
      } else{
        TMECell_obj <- Imm[Imm$cell_type %in% cell_sets$Minor_type[cell_sets$Major_type==s],]
      }   
      print(table(colnames(AUC) == TMECell_obj$barcode))
      AUC_colnames <- paste(TMECell_obj$SimplifiedSampleName,colnames(AUC), sep=':')
      
      TMECell_obj <- TMECell_obj[! TMECell_obj$sample_name %in% c('Cold'),]
      TMECell_obj <- TMECell_obj[! TMECell_obj$orig.ident %in% c("p3", "p4", "p5"),]
      # TMECell_obj$Organ <- unique(Epi$Organ)
      table(TMECell_obj$Tissue)
      TMECell_obj <- TMECell_obj[! TMECell_obj$Tissue %in% c('UNC','ADJ','Healthy'),]
      print(table(TMECell_obj$Tissue))
      precancer <- c('AAH','AD','AEH','AK',#'AIS',
                     'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
                     'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
                     'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
                     'SER','SIM','WIM')
      samples <- unique(TMECell_obj$orig.ident)
      print(samples)
      setdiff(samples, unique(Epi$orig.ident))
      setdiff(unique(Epi$orig.ident), samples)
      table(TMECell_obj$orig.ident %in% Epi$orig.ident)
      TMECell_obj <- TMECell_obj[TMECell_obj$orig.ident %in% Epi$orig.ident,]
      sample.meta <- Epi[,c("SimplifiedSampleName","Tissue",'orig.ident')]
      sample.meta <- distinct(sample.meta)
      TMECell_obj$SimplifiedSampleName <- sample.meta$SimplifiedSampleName[match(TMECell_obj$orig.ident, sample.meta$orig.ident)]
      table(TMECell_obj$SimplifiedSampleName == TMECell_obj$orig.ident)
      TMECell_obj$Molecular.typing[TMECell_obj$Tissue %in% precancer] <- 'Precancer'
      TMECell_obj$Molecular.typing[! TMECell_obj$Tissue %in% precancer] <- 'Cancer'
      table(TMECell_obj$Molecular.typing)
      table(TMECell_obj$SimplifiedSampleName)
      print(table(colnames(AUC) %in% TMECell_obj$barcode))
      # AUC <- AUC[,colnames(AUC) %in% TMECell_obj$barcode]
      print(table(colnames(AUC) == TMECell_obj$barcode))
      # AUC_colnames <- paste(TMECell_obj$SimplifiedSampleName,colnames(AUC), sep=':')
      
      TMEcells <- data.frame(table(TMECell_obj$SimplifiedSampleName))
      TMEcells <- TMEcells$Var1[TMEcells$Freq>=10]
      
      TMECell_obj <- TMECell_obj[TMECell_obj$SimplifiedSampleName %in% TMEcells,]
      # library(AUCell, lib.loc = "/data/rluo4/lorihan/R/site-library")
      sample.meta <- unique(TMECell_obj$SimplifiedSampleName)
      # AUC_colnames <- paste(TMECell_obj$SimplifiedSampleName,colnames(AUC), sep=':')
      for(i in 1:length(sample.meta)){#   
        sample.data <- AUC[,grep(sample.meta[i],AUC_colnames)]
        # Print the dimensions of sample.data
        # cat("Iteration:", i, "Dimensions of sample.data:", dim(sample.data), "\n")
        medians <- apply(sample.data,1,median)
        if(i == 1){
          AUC.sample <- data.frame(medians)
          colnames(AUC.sample)[i] <- sample.meta[i]
        } else {
          AUC.sample <- cbind(AUC.sample, data.frame(medians))
          colnames(AUC.sample)[i] <- sample.meta[i]
        }
      }
      sample.meta <- TMECell_obj[,c("SimplifiedSampleName","Tissue",'orig.ident')]
      sample.meta <- distinct(sample.meta)
      # obs <- paste(organ,'Epi',sep = '_')
      # if(organ == 'Chen'){
      #   obs <- gsub('Chen','CRC',obs)
      # }
      cytotrace.meta <- unique(Epi$SimplifiedSampleName)
      cytotrace.meta <- cytotrace.meta[cytotrace.meta %in% sample.meta$SimplifiedSampleName]
      Cyto_colnames <- paste(Epi$SimplifiedSampleName,rownames(Epi), sep=':')
      # malignant.meta <- unique(Epi$SimplifiedSampleName)
      # malignant.meta <- malignant.meta[malignant.meta %in% sample.meta$SimplifiedSampleName]
      for(i in 1:length(cytotrace.meta)){#     sample.data <- Cyto[,grep(cytotrace.meta[i],colnames(Cyto))]
        sample.data <- Epi[grep(cytotrace.meta[i],Cyto_colnames),]
        sample.data <- sample.data['CytoTRACE',drop=FALSE]
        medians <- apply(sample.data,2,median)
        if(i == 1){
          Cyto.sample <- data.frame(medians)
          colnames(Cyto.sample)[i] <- cytotrace.meta[i]
        } else {
          Cyto.sample <- cbind(Cyto.sample, data.frame(medians))
          colnames(Cyto.sample)[i] <- cytotrace.meta[i]
        }
      }
      table(colnames(AUC.sample) %in% colnames(Cyto.sample))
      AUC.sample <- AUC.sample[,colnames(AUC.sample) %in% colnames(Cyto.sample)]
      dim(AUC.sample)#20 32
      TME <- AUC.sample
      # TME <- rbind(AUC.sample, Cyto.sample)
      colnames(TME) 
      # colEpi <- Epi[,c('orig.ident','SimplifiedSampleName')]
      # colEpi <- colEpi[! duplicated(colEpi$SimplifiedSampleName),]
      # TME <- TME[,colnames(TME) %in% colEpi$orig.ident]
      indir <- '/data/rluo4/All/Output/sdata_ABN/'
      pc_file = paste0(indir, organ,'_DEG.RData')
      if(organ !='Chen'){
        load(pc_file)
        if(organ=='CRC'){
          pc_df <-  pc_df_Becker
          Patt <- Patt_Becker
        }
      } else{
        pc_file <- gsub('Chen','CRC',pc_file)
        load(pc_file)
        pc_df <- pc_df_Chen
        Patt <- Patt_Chen
      }
      unique(pc_df$Sample)
      table(pc_df$Tissue)
      head(Patt)
      table(pc_df$SimplifiedSampleName %in% colnames(TME))
      STM_df <- pc_df[pc_df$SimplifiedSampleName %in% colnames(TME),]
      STM_df$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                                     gsub('Tumor','PRAD', STM_df$Tissue)))
      
      order_samples <- STM_df$SimplifiedSampleName
      TME["Malignancy",] <- STM_df$nearest_spline_x_vals[match(colnames(TME),order_samples)]
      TME <- TME[, ! is.na(TME["Malignancy",])]
      table(colnames(TME) %in% order_samples)
      
      TME <- TME[, match(order_samples, colnames(TME))]
      # save(Cor_MP_Malignancy, file=paste0(data_path,"../Cor_MP_Malignancy_new.rds"))
      table(STM_df$Tissue)
      color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
      # color_cluster=c(RColorBrewer::brewer.pal(n = 9, name = 'Paired'))#,color_cluster)
      color_cluster <- color_cluster[1:length(unique(STM_df$Tissue))]
      names(color_cluster)=unique(STM_df$Tissue)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
      col = list(Disease = color_cluster)#c("AK" = "#E41A1C","cSCC" = "#377EB8", "SCCIS" = "#984EA3"))
      
      setwd(analysis_parent_folder)
      # cell_folder <- ifelse(organ != 'Chen','res_fig_gsva', 'res_fig_gsva_Chen')
      cell_folder <- ifelse(method_use == 'GSVA', 'res_fig_gsva', 'res_fig')
      if (!dir.exists(paste0(cell_folder))){
        dir.create(paste0(cell_folder))
      }
      setwd(cell_folder)
      getwd()
      enrich_df <- as.data.frame(t(TME[-nrow(TME),]))
      enrich_df$SimplifiedSampleName <- rownames(enrich_df)
      STM_df <- left_join(STM_df[, c("SimplifiedSampleName",'nearest_spline_x_vals','Tissue','sample_name')], enrich_df,by='SimplifiedSampleName')
      tcolors <- color_cluster#paletteDiscrete(values = unique(STM_df$Tissue))
      n <- colnames(STM_df)[-(1:4)]#rownames(test.diff)
      
      library(pheatmap) # I like esoteric packages!
      library(RColorBrewer)
      library(iheatmapr)
      library(datasets)
      library(reshape2)
      library(gplots)
      library(ComplexHeatmap)
      library(grid)
      library(circlize)
      # colnames(CD8T_Signature)
      rownames(AUC)
      
      TME.matrix <- as.matrix(TME[-nrow(TME),])
      if(s %in% c('CD4T','CD8T')){
        rownames(TME.matrix) <-  sub(".*?-","", rownames(TME.matrix))
      }
      table(rownames(TME.matrix) %in% Func$features)
      setdiff(rownames(TME.matrix) , Func$features)
      rownames(TME.matrix) <-  gsub("[.]","-", rownames(TME.matrix))
      TME.matrix <- TME.matrix[match(Func$features, rownames(TME.matrix)), ]
      setdiff(rownames(TME.matrix) , Func$features)
      ha1 <- HeatmapAnnotation(Disease=STM_df$Tissue,
                               col = col,show_legend = TRUE)
      # Create the heatmap annotation
      options(bitmapType = 'cairo')
      
      set.seed(123)
      par(mfrow = c(1, 1))
      max <- round(max(TME.matrix),1)
      max
      min <- round(min(TME.matrix),1)
      f = paste0(organ, '_',s,"_heatmap.png")
      # hei <- ifelse(length(Func$features)>8, length(Func$features)*0.5,length(Func$features)*0.6)
      # png(height=hei,width=8, filename=f, units = "in", res = 400)
      # method_use <- 'AUCell'
      if( method_use == 'GSVA'){
        min <- ifelse(min<0, 0, min)
        range_col <- c(min, (max(TME.matrix)+min(TME.matrix))*1/4,  (max(TME.matrix)+min(TME.matrix))*3/4 , max)
      }else{
        range_col <- c(0, max*1/4, max*3/4, max) 
      } 
      # p = Heatmap(TME.matrix,show_row_names = TRUE,show_column_names = FALSE,cluster_rows = FALSE,cluster_columns = FALSE,
      #             col=colorRamp2(range_col,c("seagreen4", "lightyellow", "#CAB2D6")), name = 'GSVA',
      #             row_split = Func[,ncol(Func)], #top_annotation = ha1, # rep(c("A","B","C"),each=6),
      #             #column_split = rep(c("A","B","C","D"),each=6),row_title = NULL,column_title = NULL
      #             clustering_method_columns = "complete",show_column_dend = TRUE,show_row_dend = F,row_names_gp = gpar(fontsize = 12))
      #"#f89e81","#99a9cc","#dd9bc5"
      TME.matrix <- TME.matrix[!rownames(TME.matrix) %in% c('NeoTCR','Senescence'),]
      rownames(TME.matrix)  <- gsub('-Chemokine','', rownames(TME.matrix) )
      rownames(TME.matrix)  <- gsub('-Cytokine','', rownames(TME.matrix) )
      
      hei <- ifelse(length(Func$features)<8, length(Func$features)*0.28,length(Func$features)*0.18)
      wid <- length(unique(pc_df$SimplifiedSampleName))
      wid <- ifelse( wid >30, wid*0.123, ifelse(wid<20, wid*0.3, wid*0.2) )
      if(s %in% c('CAF','Myeloid')){
        wid <- wid*0.67
      }
      
      p = Heatmap(TME.matrix,show_row_names = TRUE,show_column_names = FALSE,cluster_rows = FALSE,cluster_columns = FALSE,
                  col=colorRamp2(range_col, c('white',  "lightyellow", "#CAB2D6", "pink")), name = method_use,#'AUCell',
                  row_split = Func[,ncol(Func)], #top_annotation = ha1, # rep(c("A","B","C"),each=6),
                  #column_split = rep(c("A","B","C","D"),each=6),row_title = NULL,column_title = NULL
                  clustering_method_columns = "complete",show_column_dend = TRUE,show_row_dend = F,row_names_gp = gpar(fontsize = 11.5))
      
      if(organ == 'Endometrium'){
        wid <- wid*1.3
      }
      png(height=hei,width=wid, filename=f, units = "in", res = 500)
      draw(p, heatmap_legend_side = "left", annotation_legend_side = "top", legend_grouping = "original")
      dev.off()
      
      
      if(s %in% c('NK')){
        wid <- wid*0.7
        # hei = hei*2
        if(organ == 'Endometrium'){
          wid <- wid*1.4
        }
        png(height=hei,width=wid, filename=f, units = "in", res = 500)
        p = Heatmap(TME.matrix,show_row_names = TRUE,show_column_names = FALSE,cluster_rows = FALSE,cluster_columns = FALSE,
                    col=colorRamp2(range_col, c('white',  "lightyellow", "#CAB2D6", "pink")), name = method_use,#'AUCell',
                    row_split = Func[,ncol(Func)], #top_annotation = ha1, # rep(c("A","B","C"),each=6),
                    #column_split = rep(c("A","B","C","D"),each=6),row_title = NULL,column_title = NULL
                    clustering_method_columns = "complete",show_column_dend = TRUE,show_row_dend = F,row_names_gp = gpar(fontsize = 7.5))
        
        draw(p, heatmap_legend_side = 'left', annotation_legend_side = "top", legend_grouping = "original")
        dev.off()
      }    
      
    }
  }
}


indir = '/data/rluo4/All/Output/sdata_TME/'
for (organ in organ_all[-5]) {
  pandas.Imm <- read.csv(paste0(indir, organ, '_imm_res_fig/', organ, '_Imm_sdata.obs.csv') )
  table(pandas.Imm$Sample_Type)
  if(organ == 'Breast'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "BRCA1-mut","DCIS","IDC"))
  }
  if(organ == 'Liver'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "NAFLD","Cirrhotic","HCC"))
  }
  if(organ == 'Chen'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "SER",'AD',"MSS","MSI-H"))
  }
  if(organ == 'Skin'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "AK","SCCIS","cSCC"))
  }
  if(organ == 'Cervix'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "N_HPV","HSIL_HPV","CC"))
  }
  if(organ == 'HNSCC'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "NEOLP","EOLP","LP","OSCC"))
  }
  if(organ == 'Endometrium'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "AEH","EEC"))
  }
  if(organ == 'Lung'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "AAH","AIS","MIAC",'IAC'))
  }
  if(organ == 'Prostate'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "BPH",'PRAD'))
  }
  if(organ == 'Pancreas'){
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "PanIN",'PDAC'))
  }
  if(organ == 'THCA'){
    pandas.Imm <- pandas.Imm[ ! pandas.Imm$Sample_Type %in% c('Goiters'),]
    pandas.Imm$Sample_Type <- ordered(pandas.Imm$Sample_Type,levels=c("Healthy", "HT","PTC","ATC"))
  }
  
  library(ggpubr)
  library(rstatix)
  colnames(pandas.Imm)
  stat.test <- pandas.Imm %>%
    # group_by(CellType) %>%
    t_test(CD8..NeoTCR.Signature ~ Sample_Type) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
  stat.test 
  # Make facet and add p-values
  stat.test <- stat.test %>% add_xy_position(x = "Sample_Type")
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  color_cluster <- rev(color_cluster)
  color_cluster <- color_cluster[1:length(unique(pandas.Imm$Sample_Type))]
  names(color_cluster)=unique(pandas.Imm$Sample_Type)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
  p = ggboxplot(pandas.Imm,"Sample_Type", "CD4..NeoTCR.Signature", fill = "Sample_Type",   notch="TRUE",
               palette = color_cluster, add = c("mean_se"),#, "jitter"),
               size = 0.3)+xlab("")+ylab("")+
    
    stat_pvalue_manual(
      stat.test, bracket.nudge.y = 0.64, hide.ns = TRUE,
      label = "{p.adj.signif}"#"{p.adj}{p.adj.signif}"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    theme(legend.position = 'none',
          legend.key.width = unit(0.3, "cm"), legend.box = "vertical", legend.direction = "horizontal",
          # legend.key.width = unit(0.3, "cm"), legend.box = "vertical", legend.direction = "vertical",
          plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 13, angle = 15),#element_blank(),#
          axis.text.y   = element_text(color = 'black', size = 14, angle = 0),
          axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
          axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
          legend.title =  element_text(color = 'black', size  = 13), #element_blank(),#
          legend.text   = element_text(color = 'black', size   = 14),
          axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
          axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    )# + facet_wrap(CellType ~ ., scales = "free")
  wid = length(unique(pandas.Imm$Sample_Type))
  wid = ifelse(wid < 4, wid*2, 
               ifelse(wid >6,  wid*0.6, wid*2.2))
  hei = length(unique(pandas.Imm$Sample_Type))
  hei = ifelse(hei < 3, hei*1.6, 
               ifelse(hei >6,  hei*0.4, hei))
  ggsave(plot = p, dpi = 500, file = paste0(indir,'../res_fig/NeoTCR_fig/', organ, '-NeoTCR4.png'), width = hei , height = hei)
  
  color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
  # color_cluster <- rev(color_cluster)
  color_cluster <- color_cluster[1:length(unique(pandas.Imm$Sample_Type))]
  names(color_cluster)=unique(pandas.Imm$Sample_Type)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
  
  p<-ggboxplot(pandas.Imm,"Sample_Type", "CD8..NeoTCR.Signature", fill = "Sample_Type",   notch="TRUE",
               palette = color_cluster, add = c("mean_se"),#, "jitter"),
               size = 0.3)+xlab("")+ylab("")+
    
    stat_pvalue_manual(
      stat.test, bracket.nudge.y = 0.02, hide.ns = TRUE,
      label = "{p.adj.signif}"#"{p.adj}{p.adj.signif}"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    theme(legend.position = 'none',
          legend.key.width = unit(0.3, "cm"), legend.box = "vertical", legend.direction = "horizontal",
          # legend.key.width = unit(0.3, "cm"), legend.box = "vertical", legend.direction = "vertical",
          plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 13, angle = 15),#element_blank(),#
          axis.text.y   = element_text(color = 'black', size = 14, angle = 0),
          axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
          axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
          legend.title =  element_text(color = 'black', size  = 13), #element_blank(),#
          legend.text   = element_text(color = 'black', size   = 14),
          axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
          axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    )# + facet_wrap(CellType ~ ., scales = "free") 
  
  wid = length(unique(pandas.Imm$Sample_Type))
  wid = ifelse(wid < 4, wid*2, 
               ifelse(wid >6,  wid*0.6, wid*2.2))
  hei = length(unique(pandas.Imm$Sample_Type))
  hei = ifelse(hei < 3, hei*1.6, 
               ifelse(hei >6,  hei*0.4, hei))
  ggsave(plot = p, dpi = 500, file = paste0(indir,'../res_fig/NeoTCR_fig/', organ, '-NeoTCR8.png'), width = hei , height = hei)
  
  
}




print(organ)
# organ = 'Breast'
epi <- gsub('Chen','CRC',paste(organ, 'Epi', sep = "_"))
str <-  gsub('Chen','CRC',paste(organ, 'Str', sep = "_"))
imm <- gsub('Chen','CRC',paste(organ, 'Imm', sep = "_"))

Epi <- cell_all[[epi]]
Epi$barcode = rownames(Epi)
Epi$cell_type <- Epi$Cell_Type
table(Epi$Organ,Epi$Tissue)
Epi$Tissue[Epi$orig.ident %in% c('PTCwithHT_1','PTCwithHT_6', 'PTCwithHT_8')] <- 'PTC'
Epi$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                            gsub('Tumor','PRAD', Epi$Tissue)))
Epi$Major_type <- 'Epithelial'
print(table(Epi$Organ))
Imm <- cell_all[[imm]]
Imm$Major_type <- 'Immune'

Str <- cell_all[[str]]
Str$Major_type <- 'Stromal'

data_path <- '/data/rluo4/All/Output/Enrichment/'
analysis_parent_folder <- paste0(data_path,organ,'_Enrichment')
analysis_parent_folder <- gsub('Chen','CRC',analysis_parent_folder)
if (!dir.exists(paste0(analysis_parent_folder))){
  dir.create(paste0(analysis_parent_folder))
}
range = 4:12 # 迭代次数
gmin = 5 #最低关联数量
ncores = 5 #运行核数量
sets = unique(cell_sets$Major_type)
# sets <- c('CD8T','CD4T','Myeloid', 'NK', 'CAF')
# for(s in sets){
# s = 'CAF'
setwd(analysis_parent_folder)
Func <- Curated[[s]]
colnames(Func) <- c('features', 'functions')
Func$features <- as.character( Func$features)
Func$functions <- Func$features
Func <- dplyr::distinct(Func)
if(s %in% c('CD4T','CD8T')){
  Func$features <- sub(".*?-","", Func$features)#str_split(Func$features, '-', simplify = T)[,2]
}
Differentiation <- c('Naïve', 'Activation.Effector.function', 'Exhaustion')
Function <-   c('NeoTCR', 'TCR.Signaling','TCR.signaling', 'Cytotoxicity', 'Cytokine.Cytokine.receptor', 'Chemokine.Chemokine.receptor', 'Senescence', 'Anergy', 'NFKB.Signaling', 'Stress.response', 'MAPK.Signaling', 'Adhesion', 'IFN.Response','IFN.response','Treg.signature','Costimulatory.molecules',
                'TaNK','Inflamatoty','Stress','HLA-dependent.inhibitory.receptor','HLA-independent.inhibitory.receptor','HLA-dependent.activating.receptor','HLA-independent.activating.receptor')
Metabolism <- c('Oxidative.phosphorylation', 'Glycolysis', 'Fatty.acid.metabolism','Lipid.metabolism','OXPHOS')
Apoptosis <- c('Pro-apoptosis', 'Anti-apoptosis')

Func$functions[Func$features %in% Differentiation] <- 'Differentiation'
Func$functions[Func$features %in% Function] <- 'Function'
Func$functions[Func$features %in% Metabolism] <- 'Metabolism'
Func$functions[Func$features %in% Apoptosis] <- 'Apoptosis'

Func$functions <- sub(".*?-","", Func$functions)#str_split(Func$features, '-', simplify = T)[,2]
Func$features <- sub("Inflammototy","Inflammotory", Func$features)#str_split(Func$features, '-', simplify = T)[,2]
print(s)
outfile_gsva <- paste0(organ, '_', s,'_ssgsea.rds')
outfile_gsva <- gsub('Chen', 'CRC', outfile_gsva)

if (! file.exists(outfile_gsva)) {
  print(paste0('There is no ', s, ' in ', organ,  " , so skip !"))
  next;
}

if (file.exists(outfile_gsva)) {
  print(paste0('start gssea: ', s, " is already over!"))
  # next;
  load(outfile_gsva)
}

outfile_auc <- paste0(organ, '_', s,'_aucell.rds')
outfile_auc <- gsub('Chen', 'CRC', outfile_auc)

if (file.exists(outfile_auc)) {
  print(paste0('start aucell: ', s, " is already over!"))
  # next;
  load(outfile_auc)
}

if(s=='CAF'){
  TMECell_obj <- Str[Str$cell_type %in% cell_sets$Minor_type[cell_sets$Major_type==s],]
} else{
  TMECell_obj <- Imm[Imm$cell_type %in% cell_sets$Minor_type[cell_sets$Major_type==s],]
}   
print(table(colnames(AUC) == TMECell_obj$barcode))
AUC_colnames <- paste(TMECell_obj$SimplifiedSampleName,colnames(AUC), sep=':')

TMECell_obj <- TMECell_obj[! TMECell_obj$sample_name %in% c('Cold'),]
TMECell_obj <- TMECell_obj[! TMECell_obj$orig.ident %in% c("p3", "p4", "p5"),]
# TMECell_obj$Organ <- unique(Epi$Organ)
table(TMECell_obj$Tissue)
TMECell_obj <- TMECell_obj[! TMECell_obj$Tissue %in% c('UNC','ADJ','Healthy'),]
print(table(TMECell_obj$Tissue))
precancer <- c('AAH','AD','AEH','AK',#'AIS',
               'BPH','CAG','CAG with IM','Cirrhotic','CSG','Cyst',#'DCIS',
               'EOLP','FAP','Goiters', 'HGIN','HSIL_HPV','HT','LGIN','LP',#'MIAC',
               'N_HPV','NAFLD','NEOLP','PanIN','BRCA1-mut',#'SCCIS',
               'SER','SIM','WIM')
samples <- unique(TMECell_obj$orig.ident)
print(samples)
setdiff(samples, unique(Epi$orig.ident))
setdiff(unique(Epi$orig.ident), samples)
table(TMECell_obj$orig.ident %in% Epi$orig.ident)
TMECell_obj <- TMECell_obj[TMECell_obj$orig.ident %in% Epi$orig.ident,]
sample.meta <- Epi[,c("SimplifiedSampleName","Tissue",'orig.ident')]
sample.meta <- distinct(sample.meta)
TMECell_obj$SimplifiedSampleName <- sample.meta$SimplifiedSampleName[match(TMECell_obj$orig.ident, sample.meta$orig.ident)]
table(TMECell_obj$SimplifiedSampleName == TMECell_obj$orig.ident)
TMECell_obj$Molecular.typing[TMECell_obj$Tissue %in% precancer] <- 'Precancer'
TMECell_obj$Molecular.typing[! TMECell_obj$Tissue %in% precancer] <- 'Cancer'
table(TMECell_obj$Molecular.typing)
table(TMECell_obj$SimplifiedSampleName)
print(table(colnames(AUC) %in% TMECell_obj$barcode))
AUC <- AUC[,colnames(AUC) %in% TMECell_obj$barcode]
print(table(colnames(AUC) == TMECell_obj$barcode))

TME <- AUC[grep('TCR', rownames(AUC)), ]
colnames(TME) 
indir <- '/data/rluo4/All/Output/sdata_ABN/'
table(TMECell_obj$barcode %in% colnames(TME))
STM_df <- TMECell_obj[TMECell_obj$barcode %in% colnames(TME),]
STM_df$Tissue <- gsub('goiters','Goiters',gsub('Precancer','BRCA1-mut',
                                               gsub('Tumor','PRAD', STM_df$Tissue)))

table(STM_df$Tissue)
setwd(analysis_parent_folder)
cell_folder <- 'res_fig'
if (!dir.exists(paste0(cell_folder))){
  dir.create(paste0(cell_folder))
}
setwd(cell_folder)
getwd()
enrich_df <- as.data.frame(t(TME))
enrich_df$barcode <- rownames(enrich_df)
# STM_df <- left_join(STM_df[, c("SimplifiedSampleName",'nearest_spline_x_vals','Tissue','sample_name')], enrich_df,by='SimplifiedSampleName')
enrich_df$Tissue <- STM_df$Tissue[match(enrich_df$barcode,STM_df$barcode)]
enrich_df$DiseaseStage <- gsub('ADJ','Healthy',enrich_df$Tissue)
table(enrich_df$Tissue)
# enrich_df$DiseaseStage<- factor(enrich_df$DiseaseStage, ordered=T, levels = c('Healthy','Precancer','Cancer'))  #调整画图的x轴坐标顺序

# label.pos <- round(max(LF$CellFraction)[5],1)+0.05
my_comparisons <- list(c("CC","HSIL_HPV"), c("CC", "N_HPV"),
                       c("HSIL_HPV", "N_HPV"))
colnames(enrich_df)[1:2] <- c('NeoTCR8','CD8TCR')
color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
# color_cluster=c(RColorBrewer::brewer.pal(n = 9, name = 'Paired'))#,color_cluster)
color_cluster <- color_cluster[1:length(unique(enrich_df$DiseaseStage))]
names(color_cluster)=unique(enrich_df$DiseaseStage)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
col = list(Disease = color_cluster)#c("AK" = "#E41A1C","cSCC" = "#377EB8", "SCCIS" = "#984EA3"))

p <- ggboxplot(enrich_df, notch="TRUE",
               "DiseaseStage", "NeoTCR8",fill = "DiseaseStage",
               palette = color_cluster, 
               size = 0.5, 
               add = "")+xlab("")+ylab('NeoTCR8 signature score')+ 
  guides(fill = guide_legend(title = 'Disease Stage')) +
  theme(legend.position = "top",
        plot.title    = element_text(color = 'black', size   = 15, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 15,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
        axis.text.x   = element_blank(),#element_text(color = 'black', size = 15, angle = 270),#element_blank(),#
        axis.text.y   = element_text(color = 'black', size = 15, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 18, angle = 90),
        #legend.title=element_blank(),
        legend.title  = element_text(color = 'black', size  = 14),
        legend.text   = element_text(color = 'black', size   = 14),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
        panel.background = element_rect(fill='transparent')#, legend.position='right'
  )
p
p1 <- p+ stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 0.9, label.x = 1.5) # 添加全局p值#stat_compare_means(aes(group = DiseaseStage), label = "p.format")#,label.y = 0.6,label.x = 2)
p1
summary(enrich_df$NeoTCR8[enrich_df$DiseaseStage=='N_HPV'])

sdata_path = '/data/rluo4/All/Output/sdata_Imm/'
file = paste0(sdata_path, organ, '_Imm_sdata.obs.csv')
Imm.obs <- read.csv(file)
# Imm.obs <- read.csv('/data/rluo4/All/Output/sdata_TME/Cervix_sdata_Imm.obs.csv')#/data/rluo4/database/Cervix/Output/Imm/Results/idata.obs.csv
colnames(Imm.obs)
table(Imm.obs$cell_type, Imm.obs$Sample_Type)
CD8T
CD4T
table(Imm.obs$Sample_Type)
enrich_df <- Imm.obs[Imm.obs$cell_type %in% CD8T,]
enrich_df <- enrich_df[ !enrich_df$Sample_Type %in% c('Healthy','ADJ'),]
enrich_df$DiseaseStage <- enrich_df$Sample_Type
color_cluster=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
color_cluster <- color_cluster[1:length(unique(enrich_df$DiseaseStage))]
names(color_cluster)=unique(enrich_df$DiseaseStage)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
col = list(Disease = color_cluster)#c("AK" = "#E41A1C","cSCC" = "#377EB8", "SCCIS" = "#984EA3"))
# For Wilcoxon test
wilcox_test_results <- pairwise.wilcox.test(enrich_df$CD8..NeoTCR.Signature, enrich_df$Sample_Type, p.adjust.method = "bonferroni")
wilcox_test_results
# ggplot(plot_data, aes(x = Comparison, y = Value, fill = Comparison)) +
#   geom_boxplot() +
#   labs(title = "Pairwise Comparisons", x = "Comparison", y = "Value") +
#   theme_minimal()
p <- ggboxplot(enrich_df, notch="TRUE",
               "DiseaseStage", "CD8..NeoTCR.Signature",fill = "DiseaseStage",
               palette = color_cluster, 
               size = 0.5, 
               add = "")+xlab("")+ylab('NeoTCR8 signature score')+ 
  guides(fill = guide_legend(title = 'Disease Stage')) +
  theme(legend.position = "top",
        plot.title    = element_text(color = 'black', size   = 15, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 15,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
        axis.text.x   = element_blank(),#element_text(color = 'black', size = 15, angle = 270),#element_blank(),#
        axis.text.y   = element_text(color = 'black', size = 15, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 15, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 18, angle = 90),
        #legend.title=element_blank(),
        legend.title  = element_text(color = 'black', size  = 14),
        legend.text   = element_text(color = 'black', size   = 14),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
        panel.background = element_rect(fill='transparent')#, legend.position='right'
  )
p
p1 <- p+ stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 0.9, label.x = 1.5) # 添加全局p值#stat_compare_means(aes(group = DiseaseStage), label = "p.format")#,label.y = 0.6,label.x = 2)
p1

summary(enrich_df$CD8..NeoTCR.Signature[enrich_df$DiseaseStage=='N_HPV'])
summary(enrich_df$CD8..NeoTCR.Signature[enrich_df$DiseaseStage=='HSIL_HPV'])
summary(enrich_df$CD8..NeoTCR.Signature[enrich_df$DiseaseStage=='CC'])
