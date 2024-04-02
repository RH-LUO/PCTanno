##############################################################################################################################
# Guo cervix data
##############################################################################################################################
library(data.table)
library(stringr)
library(tidyverse) 
library(GEOquery)
library(reshape2)
library(Matrix)
library(Seurat)
setwd("/data/rluo4/All/Output/GEO/Guo/")
gset<-getGEO('GSE208654',destdir=".",
             AnnotGPL=F,
             getGPL=F)
pdata <- pData(gset[[1]])
pdata.ST <- pdata
setwd("/data/rluo4/lorihan/datasets/")
ascp_ena_link <- read.table("filereport_read_run_PRJNA860567_tsv.txt",sep="\t",header = T, fill = T)

table(ascp_ena_link$sample_alias %in% pdata.ST$geo_accession)
setdiff(pdata$geo_accession, ascp_ena_link$sample_alias)#SS2 samples
pdata.ST$SRA <- ascp_ena_link$run_accession[match(pdata.ST$geo_accession,
                                                   ascp_ena_link$sample_alias) ]
save(gset,pdata.ST,ascp_ena_link,file = "/data/rluo4/All/Output/GEO/Guo/Guo_ST_anno.RData")
write.table(pdata.ST[,c(1,2,47,50)],"/data/rluo4/All/Output/GEO/Guo/Guo_ST_pdata.txt",
                        quote = F,sep = ";",row.names = F,col.names = F)
asp <- ascp_ena_link[ascp_ena_link$run_accession %in% pdata.ST$SRA,c("run_accession","fastq_aspera")]
colnames(asp)[1] <- "Run"
asp_sc=asp$fastq_aspera#fastq_aspera
asp_sc=unlist(str_split(asp_sc,";"))
write.table(asp_sc,file="asp_ST_Guo.link",sep = "\n",quote = F, col.names = F, row.names = F)

# path = '/data/rluo4/All/Output/GEO/Guo/GSE208654_RAW'
# save(pdata, file='/data/rluo4/All/Output/GEO/Wan/Wan_anno.RData')
# for (i in unique(pdata$geo_accession)) {
#   path_new = paste(path,i,sep = '/')
#   dir.create(path_new)
#   setwd(paste(path_new,sep = ''))
#   cmd = paste('mv ../',i,'*barcodes.tsv.gz ','./barcodes.tsv.gz',sep = '')
#   system(cmd)
#   cmd = paste('mv ../',i,'*features.tsv.gz ','./features.tsv.gz',sep = '')
#   system(cmd)
#   cmd = paste('mv ../',i,'*mtx.gz ','./matrix.mtx.gz',sep = '')
#   system(cmd)
#   cmd = paste('mv ../',i,'*.jpg.gz ','./',sep = '')
#   system(cmd)
# }

# for (i in unique(pdata$geo_accession)) {
#   path_new = paste(path,i,sep = '/')
#   dir.create(path_new)
#   setwd(paste(path_new,sep = ''))
#   path_new = paste0(path_new,'/filtered_feature_bc_matrix/')
#   dir.create(path_new)
#   setwd(paste(path_new,sep = ''))
#   cmd = paste('mv ../','*barcodes.tsv.gz ','./barcodes.tsv.gz',sep = '')
#   system(cmd)
#   cmd = paste('mv ../','*features.tsv.gz ','./features.tsv.gz',sep = '')
#   system(cmd)
#   cmd = paste('mv ../','*mtx.gz ','./matrix.mtx.gz',sep = '')
#   system(cmd)
# }
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 0) #######------ prepare for scRNA_reference in ts860 ----#######
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
library(dplyr)
library(RColorBrewer)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)
library(future)
plan(multisession, workers=40)
availableCores() #12 #查看几个核可用
nbrOfWorkers() #4 当前可用的核有多少个
pid <- Sys.getpid()
pid #2109512
library(tidyverse)
library(rhdf5)
library(data.table)
library(SeuratDisk, lib.loc='/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(data.table)
library(stringr)
library(dplyr)
load('/home/lorihan/lrh/All/Output/pheno_all.rds')
outdir <- '/home/lorihan/lrh/All/Output/Cluster/'
setwd(outdir)
organ = 'Cervix'
data_path = '/data/rluo4/All/Output'
setwd(data_path)
cohort_directory <- paste0(data_path, '/sdata_ABN/')
setwd(cohort_directory)
#############################---import rdata_STM.obs from scEpi_pipleline---###########################
refer_file <- paste0(cohort_directory, organ, '_rdata_STM.obs.csv')
rdata.obs <- read.csv(refer_file)
# # for (obs in pheno){
#   # obs <- 'Esophagus_Epi'
#   # organ <- str_split(obs,'_', simplify = T)[,1]
print(organ)
if(organ == 'Chen'){
  # colon <- readRDS(paste0(cohort_directory, "CRC_diet_clustered_normal_colon_proj_seurat.rds"))
  colon_full <- readRDS(paste0(cohort_directory, "CRC_diet_clustered_all_samples_epithelial_colon_proj_seurat.rds"))
} else{
  # colon <- readRDS(paste0(cohort_directory, organ, "_diet_clustered_normal_colon_proj_seurat.rds"))
  colon_full <- readRDS(paste0(cohort_directory, organ, "_diet_clustered_all_samples_epithelial_colon_proj_seurat.rds"))
}
colon_full$tissue <- colon_full$orig.ident
table(colon_full$Tissue)
table(colon_full$Cell_Type)
if(organ =='Chen'){
  colon_full$SimplifiedSampleName[colon_full$sample_name=='Chen'] <- colon_full$HTAN.Specimen.ID[colon_full$sample_name=='Chen']
} else{
  colon_full$SimplifiedSampleName <- str_split(colon_full$orig.ident, '_cd45',simplify = T)[,1]#colon_full$orig.ident
  if(organ=='CRC'){
    colon_full$SimplifiedSampleName <- gsub('CRC','CRC-',colon_full$SimplifiedSampleName)
    colon_full$SimplifiedSampleName[grep("CRC",colon_full$SimplifiedSampleName)] <- 
      gsub('_','-',colon_full$SimplifiedSampleName[grep("CRC",colon_full$SimplifiedSampleName)])
  }
}
table(colon_full$Tissue)
table(colon_full$Tissue, colon_full$SimplifiedSampleName)
setwd(analysis_parent_folder)
table(colon_full$Cell_Type)
colon_full
extract_first_element <- function(x) {
  split_string <- strsplit(x, "-")
  first_element <- sapply(split_string, function(y) head(y, n = 1))
  return(first_element)
}
extract_last_element <- function(x) {
  split_string <- strsplit(x, "-")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
colon_full$cell_id <- extract_last_element(colnames(colon_full))
colon_full$cell_id <- paste(colon_full$orig.ident, colon_full$cell_id,sep = '_')
table(rdata.obs$barcode %in% colnames(colon_full))
table(rdata.obs$barcode %in% colon_full$cell_id)
setdiff(rdata.obs$barcode, colnames(colon_full))
# colon_full$cell_id[! colon_full$cell_id %in% rdata.obs$barcode] <- colnames(colon_full)[! colon_full$cell_id %in% rdata.obs$barcode]#extract_last_element(colnames(colon_full))
colon_full$cell_id <- colnames(colon_full)
table(rdata.obs$barcode %in% colon_full$cell_id)
identical(rdata.obs$barcode,colon_full$cell_id)
# subcluster <- data.frame(table(rdata.obs$cell_subtype))
# subcluster <- subcluster[subcluster$Freq>10,]
# rdata_subset <- rdata.obs#[rdata.obs$cell_subtype %in% subcluster$Var1,]
rdata <- subset(colon_full, subset = cell_id %in% rdata.obs$barcode)
rdata_Epi <- pheno_all$Cervix_Epi[colnames(rdata),! colnames(pheno_all$Cervix_Epi) %in% colnames(rdata.obs)]
intersect(colnames(rdata_Epi), colnames(rdata.obs))
table(colnames(rdata_Epi) %in% colnames(rdata.obs))
identical(rownames(rdata_Epi), colnames(rdata))
rdata_Epi$cell_id <- rdata$cell_id
table(rownames(rdata_Epi) == colnames(rdata))
rdata_Epi <- left_join(rdata.obs, rdata_Epi,by="cell_id")
rownames(rdata_Epi) <- colnames(rdata)[match(rdata_Epi$barcode, rdata$cell_id)]
rdata@meta.data <- rdata_Epi
table(rdata$SimplifiedSampleName)
keep_gene <- rowSums(rdata@assays$RNA@counts != 0) > ncol(rdata@assays$RNA@counts)/1000 # 0.618
table(keep_gene)
table(rdata$sample_name)
# rdata <- rdata[keep_gene,]
table(rdata$cell_type)
rdata$Tissue <- as.character(rdata$Sample_Type)
rdata$cell_subtype <- paste0("C",rdata$leiden_res2)#paste(group_ABN$Sample_Type,group_ABN$Tissue,sep = ':')
# index = rdata$leiden_res1.5==16
table(rdata$cell_subtype)

# next step is GSVA-Cervix.R
path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
file = paste0(path_PCT, 'PCT-',organ,'.rds')
load(file)
unique(rdata_filter$cell_class)
rdata_filter$cell_subtype <- paste0("C",rdata_filter$leiden_res2)#paste(group_ABN$Sample_Type,group_ABN$Tissue,sep = ':')

###################################subset the cell_types#######################################
subtype <- table(rdata$SimplifiedSampleName, rdata$cell_subtype)
subtype <- data.frame(subset = subtype)
subtype <- subtype[subtype$subset.Freq !=0,]
subtype <- subtype[subtype$subset.Freq >=10,]

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)
library(future)

# subset the expression data of signaling genes for saving computation cos
################################################################################
# subcluster = 'C14'
# patient = unique(subtype$subset.Var1[subtype$subset.Var2==subcluster])
# # TissueType
# future::plan("multicore", workers = 80)
# run_cellchat(object=scRNA_sub,Cluster = subcluster,TissueType = patient)
load('~/lrh/All/Output/Transition_MP.rds')
t <- unique(DEG_list_all$Cervix_all$cluster)
c <- str_split(t, '-', simplify = T)[,2]
patient = unique(subtype$subset.Var1[subtype$subset.Var2 %in% c])
patient
# scRNA_sub$orig.ident <- str_split(scRNA_sub$orig.ident, '_cd45',simplify = T)[,1]#colon_full$orig.ident
# count = subset(scRNA_sub, subset = orig.ident %in% patient)
# print(table(count$orig.ident))
# print(table(count$Tissue))
# # barcode = colnames(rdata)[rdata$leiden_res1.5==Cluster]
# # print(table(colnames(count) %in% barcode))#1737
# print(table(count$Cell_Type))
patient = paste(patient,collapse = '.')
# save(count, file = paste0('~/lrh/All/Output/sdata_TME/', patient, '_count.rds'))# saved in ts860 
count = readRDS( file = paste0('/data/rluo4/All/Output/sdata_TME/', patient, '_count.rds'))
# index <- ! is.na(count$Sample_Type) & count$cell_subtype %in% c
# print(table(count$Cell_Type[index]))
# print(table(count$Cell_Type[! index]))
# table(count$cell_subtype)
# count = count[,! index]
print(table(count$Tissue))
Idents(count) <- count$Cell_Type
print(table(count$Cell_Type))
table(count$cell_subtype)

table(rdata_filter$cell_subtype)
exclude_cluster = na.omit(setdiff(unique(count$cell_subtype), unique(rdata_filter$cell_subtype)))
table(count$cell_subtype[count$cell_subtype %in% exclude_cluster])

count = count[, ! count$cell_subtype %in% exclude_cluster]#subset(count, ! cell_subtype %in% exclude_cluster)
count
table(count$cell_subtype)
# count$cell_subtype[count$cell_subtype %in% unique(rdata_filter$cell_subtype)]
count$cell_subtype[is.na(count$cell_subtype)] <- count$Cell_Type[is.na(count$cell_subtype)]
table(count$cell_subtype) # Cervix - Epi + TME
subtype <- data.frame(subset =  table(count$cell_subtype))
subtype <- subtype[subtype$subset.Freq !=0,]
subtype <- subtype[subtype$subset.Freq >=10,]
count <- count[, count$cell_subtype %in% subtype$subset.Var1]
table(count$cell_subtype)

# Seurat - 空间spot注释 
# 1, FindTransferAnchors
# 首先使用SCTransform 处理下数据，如果前面是SCTransform处理的则可以省略该步骤。
scRNA_reference <- SCTransform(count, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# the annotation is stored in the 'celltype' column of object metadata
DimPlot(scRNA_reference, group.by = "cell_subtype", label = TRUE)
save(scRNA_reference, file = '/data/rluo4/database/Cervix/Guo/scRNA_reference.RData')
# 空转是SCT处理的，因此此处无需再次处理 
# FindTransferAnchors 结合
scRNA_reference$cell_subtype[scRNA_reference$cell_subtype=='C14'] = 'C1'
table(scRNA_reference$cell_subtype)

################################################################################
########------ prepare Healthy count data for scRNA_reference in UTH ----#######
# Cervix.H_SpaCET.R
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) #######------scRNA & ST standard pipeline ----#######
for (sample in unique(pdata.ST$geo_accession)) {
  indir = '/data/rluo4/database/Cervix/Guo/spaceranger_count/'
  # sample = 'GSM6360690'
  st <- Load10X_Spatial(paste0(indir, sample, '/outs'))
  ##### Normaliztion #####
  VlnPlot(st, features = 'nCount_Spatial', pt.size = 0.1) + NoLegend()
  SpatialFeaturePlot(st, features = 'nCount_Spatial') + theme(legend.position = 'right')
  SpatialFeaturePlot(st, features = 'nFeature_Spatial') + theme(legend.position = 'right')
  st <- SCTransform(st, assay = 'Spatial', verbose = FALSE)
  ##### Gene Expression Visualiztion #####
  SpatialFeaturePlot(st, features = c('EPCAM', 'CD47', 'ANXA1', 'SAA1'), ncol = 2)
  SpatialFeaturePlot(st, features = c('ANXA1'), pt.size.factor = 1) # 表示spot原有大小
  SpatialFeaturePlot(st, features = c('ANXA1'), alpha = c(0.1, 1)) # 表达越低越透明
  ##### Dimential reduction and clustering ######
  st <- RunPCA(st, assay = 'SCT', verbose = FALSE)
  st <- FindNeighbors(st, reduction = 'pca', dims = 1:30)
  st <- FindClusters(st, verbose = FALSE)
  st <- RunUMAP(st, reduction = 'pca', dims = 1:30)
  
  DimPlot(st, reduction = 'umap', pt.size = 1.5, label = TRUE, repel = T, label.size = 5) + guides(colour = guide_legend(override.aes = list(size=7)))
  SpatialDimPlot(st, label = TRUE, repel = T, label.size = 8)
  SpatialDimPlot(st, cells.highlight = CellsByIdentities(object = st, idents = c(0:8)), facet.highlight = TRUE, ncol = 5)
  
  # load('/data/rluo4/database/Cervix/Guo/scRNA_reference.RData')
  # 空转是SCT处理的，因此此处无需再次处理 
  # FindTransferAnchors 结合
  anchors <- FindTransferAnchors(reference = scRNA_reference, 
                                 query = st, 
                                 normalization.method = "SCT")
  # 2, TransferData
  # 使用TransferData函数进行注释，其中refdata 即为单细胞转录组中的参考列
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = scRNA_reference$cell_subtype, 
                                    prediction.assay = TRUE,
                                    weight.reduction = st[["pca"]], 
                                    dims = 1:30)
  #查看一下结果  
  class(predictions.assay)
  #[1] "Assay"
  #attr(,"package")
  #[1] "SeuratObject"
  predictions.assay@data[,1:4]
  
  # 可以看到结果的 每行 是单细胞celltype列中的细胞类型，列为barcode的ID ，会给出每个spot 的celltype 占比，介绍2种保存方式
  # （1）可以像之前一样添加至metadata.
  # （2）因为predictions.assay 本身就是SeuratObject ，因此可以单独作为一个新的slot 。
  #保存结果 
  #（1） metadata 
  predictions.res <- predictions.assay@data
  rownames(predictions.res) <- make.names(rownames(predictions.res))
  cell_types <- rownames(predictions.res)
  for(cell_type in cell_types){
    st <- AddMetaData(object= st,metadata = predictions.res[cell_type,],col.name = cell_type)
  }
  #（2）slot 
  st[["predictions"]] <- predictions.assay
  # 以上就简单的完成了Seurat中结合单细胞和空转的分析。
  # 三 空转区域差异 
  # 得到了每个spot的每种celltype的预测结果，可以查看一下关心的细胞类型分布
  DefaultAssay(st) <- "predictions"
  SpatialFeaturePlot(st, features = c("NEUT", "C1"),
                     pt.size.factor = 1.6, ncol = 2, crop = TRUE)
  # Image
  # 基于这些预测分数，还可以预测其位置受空间限制的细胞类型。使用每个spot 细胞类型的预测分数替代基因表达来作为“marks” 得到spatially variable features 
  st <- FindSpatiallyVariableFeatures(st,  assay = "predictions", 
                                      selection.method = "moransi",
                                      features = rownames(st), 
                                      r.metric = 5, 
                                      slot = "data")
  top.clusters <- head(SpatiallyVariableFeatures(st, selection.method = "moransi"), 4)
  SpatialPlot(object = st, features = top.clusters, ncol = 2)
  
  saveRDS(st, file = paste0(indir, sample, "_vis.seu.rds"))
}
# mkdir /data/rluo4/database/Cervix/Guo/vis_seu
# mv /data/rluo4/database/Cervix/Guo/spaceranger_count/*_vis.seu.rds ../vis_seu
indir = '/data/rluo4/database/Cervix/Guo/vis_seu/'
st <- readRDS(paste0(indir, sample, "_vis.seu.rds"))
SpatialDimPlot(st, label = TRUE, label.size = 3)
st@meta.data$region = as.character(st@meta.data$seurat_clusters)
st@meta.data$region = paste0('reg', st@meta.data$region)
st@meta.data$region = factor(st@meta.data$region, levels = sort(unique(st@meta.data$region)))

assignLabels <- function(object, prediction = "predictions") {
  pred <- object[[prediction]]@data
  pred <- pred[1:(nrow(pred)-1), ]
  # label each spot based on the maximum prediction probability
  labels = rownames(pred)[apply(pred, 2, which.max)]
  names(labels) <- colnames(pred)
  object$labels <- factor(labels)
  Idents(object) <- "labels"
  return(object)
}
st <- assignLabels(st, prediction = "predictions")
table(st$labels)

##############################################################################################################################
anchors <- FindTransferAnchors(reference = scRNA_reference, 
                               query = st, 
                               normalization.method = "SCT")
# 2, TransferData
# 使用TransferData函数进行注释，其中refdata 即为单细胞转录组中的参考列
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = scRNA_reference$cell_subtype, 
                                  prediction.assay = TRUE,
                                  weight.reduction = st[["pca"]], 
                                  dims = 1:30)
#查看一下结果  
class(predictions.assay)
#[1] "Assay"
#attr(,"package")
#[1] "SeuratObject"
predictions.assay@data[,1:4]

# 可以看到结果的 每行 是单细胞celltype列中的细胞类型，列为barcode的ID ，会给出每个spot 的celltype 占比，介绍2种保存方式
# （1）可以像之前一样添加至metadata.
# （2）因为predictions.assay 本身就是SeuratObject ，因此可以单独作为一个新的slot 。
#保存结果 
#（1） metadata 
predictions.res <- predictions.assay@data
rownames(predictions.res) <- make.names(rownames(predictions.res))
cell_types <- rownames(predictions.res)
for(cell_type in cell_types){
  st <- AddMetaData(object= st,metadata = predictions.res[cell_type,],col.name = cell_type)
}
#（2）slot 
st[["predictions"]] <- predictions.assay
# 以上就简单的完成了Seurat中结合单细胞和空转的分析。
# 三 空转区域差异 
# 得到了每个spot的每种celltype的预测结果，可以查看一下关心的细胞类型分布
DefaultAssay(st) <- "predictions"
SpatialFeaturePlot(st, features = c("NEUT", "C1"),
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE)
# Image
# 基于这些预测分数，还可以预测其位置受空间限制的细胞类型。使用每个spot 细胞类型的预测分数替代基因表达来作为“marks” 得到spatially variable features 
st <- FindSpatiallyVariableFeatures(st,  assay = "predictions", 
                                    selection.method = "moransi",
                                    features = rownames(st), 
                                    r.metric = 5, 
                                    slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(st, selection.method = "moransi"), 4)
SpatialPlot(object = st, features = top.clusters, ncol = 2)
# 
# # Image
# # 以上分析就完成了，最后给一个彩蛋，使用grid.arrange 函数自定义空转图形的分布，参考2022年Immunity文献。
# # 四 彩蛋- 空转主图可视化
# # 可以学习以下几个点
# # （1）外层设定统一的theme ， 如hide_axis
# # （2）多张图根据重要性分布 ，HE染色图可以大一些 ，CA9是肿瘤的marker，也可以重点展示
# # （3）lapply 得到图形的list结果，grid.arrange 来拼图，可以随意调整位置和比例，是直接FeaturePlot得不到的。
# # （4）图像相关参数随意改，适合文章其他图的风格即可
# #hide spatial axes for seurat plots 
# hide_axis <- theme(axis.title.x=element_blank(),
#                    axis.text.x=element_blank(),
#                    axis.ticks.x=element_blank(),
#                    axis.title.y=element_blank(),
#                    axis.text.y=element_blank(),
#                    axis.ticks.y=element_blank())
# Idents(st)  <- as.numeric(Idents(st)) 
# # H&E 
# he_st <- SpatialPlot(st,repel = F,label = F,
#                      image.alpha=1,
#                      alpha = c(0,0), 
#                      pt.size.factor =0.000001) + 
#   geom_point(alpha=0)+
#   NoLegend() +
#   DarkTheme() +
#   hide_axis + 
#   ggtitle("st")+  
#   theme(text=element_text(size=14))+ 
#   theme(text=element_text(face = "bold"))
# # celltype score 
# cell_types <- c("C5" , "ICAF","CD8TEREX" , "INMON" ,"NK" ,"TREG", "PLA","M2MAC")
# 
# plot_list_mcp <- lapply(cell_types[c(5,1,2,4,6,3,7,8)],function(x){
#   plot_list_mcp <- SpatialPlot(st,
#                                features = x,
#                                image.alpha=0,
#                                pt.size.factor = 1.8) +
#     DarkTheme() +
#     hide_axis +    
#     theme(text=element_text(size=14))+ 
#     theme(text=element_text(face = "bold"))+
#     theme(legend.text=element_text(size=7))
#   
# })
# #CA9 重点marker
# ca9 <- SpatialPlot(st,
#                    features = "C4",
#                    image.alpha=0,
#                    pt.size.factor = 1.8) +
#   DarkTheme() +
#   hide_axis +    
#   theme(text=element_text(size=14))+ 
#   theme(text=element_text(face = "bold"))+
#   theme(legend.text=element_text(size=7))
# 
# #########拼图############
# plot_list_mcp[[9]] <- ca9
# plot_list_mcp[[10]] <- he_st
# #图形摆放
# lay <- rbind(c(10,9,1,2,3,4),
#              c(10,9,5,6,7,8))
# 
# grid.arrange(grobs = plot_list_mcp, layout_matrix = lay)
assignLabels <- function(object, prediction = "predictions") {
  pred <- object[[prediction]]@data
  pred <- pred[1:(nrow(pred)-1), ]
  # label each spot based on the maximum prediction probability
  labels = rownames(pred)[apply(pred, 2, which.max)]
  names(labels) <- colnames(pred)
  object$labels <- factor(labels)
  Idents(object) <- "labels"
  return(object)
}
st <- assignLabels(st, prediction = "predictions")
table(st$labels)
color.use <- scPalette(nlevels(st)); names(color.use) <- levels(st)
Seurat::SpatialDimPlot(st, label = F, label.size = 3, cols = color.use)

p1 <- Seurat::SpatialDimPlot(seu1, label = F, label.size = 3, cols = color.use)
color.use <- scPalette(nlevels(seu2)); names(color.use) <- levels(seu2)
p2 <- Seurat::SpatialDimPlot(st, label = F, label.size = 3, cols = color.use) + NoLegend()
p2

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) #######------SPOTlight pipeline ----#######
library(data.table)
library(stringr)
library(tidyverse) 
library(GEOquery)
library(reshape2)
library(Matrix)
library(Seurat)
# library(devtools)
# install_github("https://github.com/MarcElosua/SPOTlight")
library(SPOTlight)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(scran)
load(file = "/data/rluo4/All/Output/GEO/Guo/Guo_ST_anno.RData")
print(table(pdata.ST$geo_accession))
load('/data/rluo4/database/Cervix/Guo/scRNA_reference.RData')
sce <- as.SingleCellExperiment(scRNA_reference)
### Feature selection
sce <- logNormCounts(sce)
### Variance modelling
# 去掉核糖体和线粒体基因
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(sce))
dec <- modelGeneVar(sce , subset.row = genes)
# 计算高变基因
hvg <- getTopHVGs(dec, n = 3000)
# 加上细胞注释信息
colLabels(sce) <- colData(sce)$cell_subtype
##################################################################################
indir = '/data/rluo4/database/Cervix/Guo/spaceranger_count/'
outdir = '/data/rluo4/database/Cervix/Guo/vis_seu/'
SPOTlight_obj <- NULL
assignLabels <- function(object, prediction = "predictions") {
  pred <- object[[prediction]]@data
  # pred <- pred[1:(nrow(pred)-1), ] This line is for {st <- readRDS(paste0(indir, sample, "_vis.seu.rds"))} predicted by TransferData() in Seurat
  # label each spot based on the maximum prediction probability
  labels = rownames(pred)[apply(pred, 2, which.max)]
  names(labels) <- colnames(pred)
  object$labels <- factor(labels)
  Idents(object) <- "labels"
  return(object)
}
for (sample in pdata.ST$geo_accession) {
  # sample  = 'GSM6360689'
  # save(st, res , mgs_df, file = paste0(indir, sample, "_spotlight.RData"))
  load(file = paste0(indir, sample, "_spotlight.RData"))
  st <- assignLabels(st, prediction = "SPOTlight")
  print(table(st$labels))
  ################################ 三 SPOTlight 结果可视化 ########################
  # pdf( paste0(outdir, sample, '_SPOTlight.pdf'), height = 10, width = 10 )
  # # 1，提取SPOTlight结果
  # head(mat <- res$mat)[, seq_len(length(unique(sce$cell_subtype)))]
  # mod <- res$NMF
  # res.data <- (mat <- res$mat)[, seq_len(length(unique(sce$cell_subtype)))]
  # # 2，topic对细胞类型的拟合情况
  # # 使用plotTopicProfiles 函数绘制拟合情况图
  # ## 检查topic对细胞类型的拟合情况
  # plotTopicProfiles(
  #   x = mod,
  #   y = sce$cell_subtype,  #
  #   facet = FALSE,
  #   min_prop = 0.01,
  #   ncol = 1) +
  #   theme(aspect.ratio = 1)
  # 
  # # 3，Correlation Matrix 和 Co-localization
  # ## Spatial Correlation Matrix
  # p1 <- plotCorrelationMatrix(mat)
  # ## Co-localization
  # p2 <- plotInteractions(mat, "heatmap")
  # plotInteractions(mat, "network")
  # # 4，Scatterpie 饼图
  # #
  # # 这时SPOTlight 注释spot后的核心图，将每个spot中的各celltype比例绘制为饼图，可以绘制到切片tiff背景上（左图），也可以同样的形状绘制在白板上（右图）。
  # ct <- colnames(mat)
  # mat[mat < 0.1] <- 0
  # #自定义颜色
  # paletteMartin <- c(
  #   "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
  #   "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
  #   "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
  # pal <- colorRampPalette(paletteMartin)(length(ct))
  # names(pal) <- ct
  # p3 <- plotSpatialScatterpie(
  #   x = st,
  #   y = mat,
  #   cell_types = colnames(mat),
  #   img = T, #以tiff为背景
  #   scatterpie_alpha = 1,
  #   pie_scale = 0.4,
  #   # Rotate the image 90 degrees counterclockwise
  #   degrees = -90,
  #   # Pivot the image on its x axis
  #   axis = "h") +
  #   scale_fill_manual(
  #     values = pal,
  #     breaks = names(pal))
  # 
  # p4 <- plotSpatialScatterpie(
  #   x = st,
  #   y = mat,
  #   cell_types = colnames(mat),
  #   img = FALSE,
  #   scatterpie_alpha = 1,
  #   pie_scale = 0.4) +
  #   scale_fill_manual(
  #     values = pal,
  #     breaks = names(pal))
  # p3 + p4
  # # 注意：（1）可以通过img 是否添加背景 ；
  # # （2）pal 自定义颜色，注意长度 与celltype个数一致 ；
  # # （3）degrees 调整翻转角度 和 tiff图片一致
  # # 5，批量绘制解析后的结果
  # # 前面Seurat的 四 彩蛋- 空转主图可视化 部分了介绍了lapply 得到list然后自定义拼图的方式，这里介绍一下SpatialFeaturePlot进行绘制的方式。
  # #
  # # 注意要用 & 而不是 + ，否则后续的theme等设置只会对最后一张图有效果。
  # celltypes = rownames(st)
  # SpatialFeaturePlot(st, features = celltypes,
  #                    pt.size.factor = 1.6,
  #                    ncol = 6,
  #                    crop = TRUE) &
  #   DarkTheme() &
  #   theme(text=element_text(size=14)) &
  #   theme(text=element_text(face = "bold")) &
  #   theme(legend.text=element_text(size=7))
  # dev.off()
  ################################################################
  st_obj <-  vector("list", 3)
  names(st_obj) <- c('obj', 'res', 'mgs_df')
  st_obj$obj = st
  st_obj$res <- res
  st_obj$mgs_df <- mgs_df
  names(st_obj) <- paste0(sample,'_', names(st_obj))
  SPOTlight_obj <- c(SPOTlight_obj, st_obj)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) #######------SpaCET pipeline ----#######
indir = '/data/rluo4/database/Cervix/Guo/spaceranger_count/'
outdir = '/data/rluo4/database/Cervix/Guo/vis_seu/'
organ = 'Cervix'

library(SpaCET)
load(file = "/data/rluo4/All/Output/GEO/Guo/Guo_ST_anno.RData")
print(table(pdata.ST$geo_accession))
##############################################################################################################################
load(file = '/data/rluo4/database/Cervix/Guo/scRNA_reference.RData')
# show count matrix
scRNA_reference <- scRNA_reference[, !scRNA_reference$cell_subtype %in% c('C10','C4')]
sc_counts <- as.data.frame(scRNA_reference@assays$RNA$counts)#read.csv(file.path(refdir,"counts.csv")) 
sc_counts[1:6,1:5]
# show cell annotation matrix
sc_annotation <- scRNA_reference@meta.data[,c('barcode', 'cell_subtype')]
sc_annotation[1:6,]
colnames(sc_annotation) <- c('cellID','bio_celltype')
table(sc_annotation$bio_celltype)
# show cell type lineage tree
sc_lineageTree <-  vector("list", 8)
names(sc_lineageTree) <- c('Epithelia', 'B', 'CD8T','CD4T','Myeloid','Lymphoid','Mesenchymal','Endothelia')
# sc_lineageTree$Epithelia = c('C1','C10','C14','C16','C4','C5','C6','C8')
sc_lineageTree$Epithelia = c('C1','C14','C16','C5','C6','C8')
sc_lineageTree$B <- c('BN','PLA')
sc_lineageTree$CD8T <- c( 'CD8TCM',   'CD8TEFF',  'CD8TEREX', 'CD8TEXINT',   'CD8TEXP', 'CD8TRM' )
sc_lineageTree$CD4T <- c('CD4TN','TFH','TREG')
sc_lineageTree$Myeloid <- c('DC','INMON','M2MAC','MAST','MON','NEUT')
sc_lineageTree$Lymphoid <- c('GDT','NK')
sc_lineageTree$Mesenchymal <- c('ICAF','INCAF')
sc_lineageTree$Endothelia <- 'END'
head(sc_lineageTree)
# i) 依据匹配的单细胞数据集分析空间数据
# saveRDS(SpaCET_obj, file.path( paste0(outdir, sample,'_SpaCET_obj.rds')))
# PDAC_Path <- system.file("extdata", 'oldST_PDAC', package = 'SpaCET')
# load(paste0(PDAC_Path,"/st_PDAC.rda"))
# st <- readRDS(paste0(outdir, sample, "_vis.seu.rds"))
# counts <- as.matrix(st@assays$Spatial$counts)#
# coords <- read.csv(file.path(paste0(indir, sample, '/outs/spatial/tissue_positions.csv')))#
# rownames(coords) <- coords[,1]; coords <- coords[,3:4]#coords[,1] <- NULL
# colnames(coords) <- c('x','y')
# spotCoordinates <- coords[match(colnames(counts), rownames(coords)),]
# # show count matrix
# counts[1:6,1:5]
# # show coordinate matrix
# spotCoordinates[1:5,]

# 
# # load ST data to create an SpaCET object.
# SpaCET_obj <- create.SpaCET.object(
#   counts=counts,
#   spotCoordinates=spotCoordinates,
#   imagePath=NA,
#   platform = "oldST"
# )

# show this object.
# str(SpaCET_obj)
# 解卷积
# load sc data
# PDAC_Path <- system.file("extdata", 'oldST_PDAC', package = 'SpaCET')
# load(paste0(PDAC_Path,"/sc_PDAC.rda"))
# show count matrix
# sc_counts <- read.csv(file.path(refdir,"counts.csv")) 
# sc_counts[1:6,1:5]
# show cell annotation matrix
# table(sc_annotation$bio_celltype)
# show cell type lineage tree
# head(sc_lineageTree)
# SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
#   SpaCET_obj, 
#   sc_counts=sc_counts, 
#   sc_annotation=sc_annotation, 
#   sc_lineageTree=sc_lineageTree, 
#   coreNo=8
# )
assignLabels <- function(object, prediction = "prediction") {
  pred <- prediction@results$deconvolution$propMat#object[[prediction]]@data
  pred <- pred[! rownames(pred) %in% c('Epithelia', 'B', 'CD8T','CD4T','Myeloid','Lymphoid','Mesenchymal','Endothelia'),]
  # label each spot based on the maximum prediction probability
  labels = rownames(pred)[apply(pred, 2, which.max)]
  names(labels) <- prediction@input$spotCoordinates$barcode
  object$labels <- factor(labels)
  Idents(object) <- "labels"
  return(object)
}
assignLabels_normal <- function(object, prediction = "prediction") {
  pred <- prediction@results$deconvolution$propMat
  pred <- pred[ rownames(pred) %in% c('Epithelia','Mesenchymal','B', 'CD8T','CD4T','Myeloid','Lymphoid','Endothelia'),]
  # label each spot based on the maximum prediction probability
  labels = rownames(pred)[apply(pred, 2, which.max)]
  names(labels) <- prediction@input$spotCoordinates$barcode
  object$labels <- factor(labels)
  Idents(object) <- "labels"
  return(object)
}
SPOTlight_obj <- NULL
for (sample in pdata.ST$geo_accession) {
  # sample  = 'GSM6360689'
  load(file = paste0(indir, sample, "_spotlight.RData"))
  SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj.rds'))
  if(sample != 'GSM6360689'){
    st <- assignLabels(st, prediction = SpaCET_obj)
  } else{
    # SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj_Hea.rds'))
    st <- assignLabels_normal(st, prediction = SpaCET_obj)
  }
  print(table(st$labels))
  
  st_obj <-  vector("list", 3)
  names(st_obj) <- c('obj', 'res', 'mgs_df')
  st_obj$obj = st
  st_obj$res <- res
  st_obj$mgs_df <- mgs_df
  names(st_obj) <- paste0(sample,'_', names(st_obj))
  SPOTlight_obj <- c(SPOTlight_obj, st_obj)
}
SPOTlight_obj_Cervix <- SPOTlight_obj


# plan(multisession, workers=4)
library(ArchR)
organ = 'Cervix'
CellChat <- NULL
for (sample in pdata.ST$geo_accession) {
  data_path = '/data/rluo4/database/Cervix/Guo/vis_seu/CCI_figures'
  if(! dir.exists(data_path)){
    dir.create(data_path)
  }
  setwd(data_path)
  file = paste0(outdir,sample, "_SpaCET_cellchat", ".rds")
  if(file.exists(file)){
    cellchat <- readRDS(file)
  } else{
    print(paste0(organ, ' ', Sample, " doesn't exist!"))
  }
  # if (!dir.exists(paste0(sample))){
  #   dir.create(paste0(sample))
  # }
  print(table(cellchat@idents))
  signaling.name <- cellchat@netP$pathways
  # print(table(pathways.show.all %in% pathways.all))
  # print(pathways.all[!pathways.all %in% pathways.show.all])
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  print(signaling.name)
  data.cci=NULL
  data.cci <- subsetCommunication(cellchat)
  data.cci <- data.cci[, c(1:7, 9)]#data.cci[,c(3,4,1,2)]
  # colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
  data.cci$Tissue <- organ
  data.cci$DiseaseStage <- sample
  CellChat <<- rbind(CellChat, data.cci)
}
CellChat_Cervix <- CellChat

#######################################################################################################
# plan(multisession, workers=4)
# library(CellChat)
# organ = 'Cervix'
# CellChat <- NULL
# for (sample in pdata.ST$geo_accession) {
#   file = paste0(outdir,sample, "_SpaCET_cellchat", ".rds")
#   if(file.exists(file)){
#     cellchat <- readRDS(file)
#   } else{
#     print(paste0(organ, ' ', sample, " doesn't exist!"))
#   }
#   
#   print(table(cellchat@idents))
#   signaling.name <- cellchat@netP$pathways
#   # print(table(pathways.show.all %in% pathways.all))
#   # print(pathways.all[!pathways.all %in% pathways.show.all])
#   if (is.null(signaling.name)) {
#     signaling.name <- signaling
#   }
#   # CCI <- vector("list",length(signaling.name))
#   # names(CCI) <- signaling.name
#   data.cci=NULL
#   thresh <- 0.05
#   net <- cellchat@net
#   pairLR.use.name <- dimnames(net$prob)[[3]]
#   for (i in signaling.name) {
#     pairLR <- searchPair(signaling = i, pairLR.use = cellchat@LR$LRsig,
#                          key = "pathway_name", matching.exact = T, pair.only = T)
#     pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
#     pairLR <- pairLR[pairLR.name, ]
#     prob <- net$prob
#     pval <- net$pval
#     prob[pval > thresh] <- 0
#     if (length(pairLR.name) > 1) {
#       pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 3, sum) != 0]
#     }else{
#       pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 0]
#     }
#     print(setdiff(pairLR.name, pairLR.name.use))
#     if (length(pairLR.name.use) == 0) {
#       stop(paste0("There is no significant communication of ",
#                   signaling.name))
#     }else {
#       pairLR <- pairLR[pairLR.name.use, ]
#     }
#     
#     data.cci <<- rbind(data.cci,pairLR)
#   }
#   data.cci <- data.cci[,c(3,4,1,2)]
#   colnames(data.cci) <- c('Ligand','Receptor','LRpair','Pathway')
#   data.cci$Tissue <- organ
#   data.cci$DiseaseStage <- sample
#   CellChat <<- rbind(CellChat, data.cci)
# }
# 
# library(ArchR)
# for (sample in pdata.ST$geo_accession[-1]) {
#   data_path = '/data/rluo4/database/Cervix/Guo/vis_seu/CCI_figures'
#   if(! dir.exists(data_path)){
#     dir.create(data_path)
#   }
#   setwd(data_path)
#   file = paste0(outdir,sample, "_SpaCET_cellchat", ".rds")
#   if(file.exists(file)){
#     cellchat <- readRDS(file)
#   } else{
#     print(paste0(organ, ' ', Sample, " doesn't exist!"))
#   }
#   if (!dir.exists(paste0(sample))){
#     dir.create(paste0(sample))
#   }
#   setwd(sample)  
#   if( file.exists('C1_source_net_bubble.png')){
#     print(paste0(sample, " is already over!"))
#     next;
#   }
#   #展示每个亚群作为source的信号传递
#   print(table(cellchat@idents))
#   
#   cell_color <- paletteDiscrete(values = unique(cellchat@meta$labels))
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
#       ggsave(file, width = wid, height = hei, dpi = 500, device = "png", limitsize =
#                FALSE)
#       
#       if(nameindex %in% t){
#         file = paste0(cell,'_target_net_bubble.png')
#         p_t <- netVisual_bubble(cellchat, sources.use = n, targets.use = nameindex, remove.isolate = TRUE)
#         hei = length(unique(p_t$data$interaction_name_2))*0.15
#         hei = ifelse(hei>6,hei,6)
#         ggsave(file, width = wid, height = hei, dpi = 500, device = "png", limitsize =
#                  FALSE)
#       }
#       p <<- c(p,cell)#
#       return(p)
#     }
#   })
#   
# }

##############################################################################################################################
# ii)  normal cervix
# pdf( paste0(outdir, 'Guo_SpaCET.pdf'), height = 12, width = 12, onefile = TRUE )

for (sample in pdata.ST$geo_accession) {
# if(sample == 'GSM6360689')  {
#   # sample = 'GSM6360689'
#   # load(file = '/data/rluo4/database/Cervix/Guo/scRNA_reference_Hea.RData')
#   # table(scRNA_reference$Cell_Type)
#   # SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj_HeaAdj.rds'))
#   SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj_Adj.rds'))
#   # SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj_Hea.rds'))
# } else{
  
  # sample = 'GSM6360691'
  SpaCET_obj <- readRDS( paste0(outdir, sample,'_SpaCET_obj.rds'))
# }
# str(SpaCET_obj)
# visiumPath = file.path('/data/rluo4/database/Cervix/Guo/spaceranger_count/',sample, '/outs/')
# SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
# # show this object.
# str(SpaCET_obj)
# # show this object.
# SpaCET_obj@input$counts[1:8,1:6]
# # calculate the QC metrics
# SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
# # plot the QC metrics
# SpaCET.visualize.spatialFeature(
#   SpaCET_obj,
#   spatialType = "QualityControl",
#   spatialFeatures=c("UMI","Gene")
# )
# 
# # show this object.
# str(SpaCET_obj)
# # 解卷积
# # load sc data
# SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
#   SpaCET_obj, 
#   sc_counts=sc_counts, 
#   sc_annotation=sc_annotation, 
#   sc_lineageTree=sc_lineageTree, 
#   coreNo=8
# )
  # par(mar=c(3,3,1,0), mfrow=c(1,2))
pdf( paste0(outdir, sample, '_SpaCET.pdf'), height = 10, width = 10, onefile = TRUE )
if(sample!='GSM6360689'){
 a <-  SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction",
  spatialFeatures = c("C5",'C1',"C14",'C8','C6','M2MAC','INMON',"MON",'NEUT','ICAF','INCAF','TREG'),
  nrow=4
)
}
if(sample=='GSM6360689'){
 a <- SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "CellFraction",
  spatialFeatures = c('Epithelia','M1MAC','M2MAC','NEUT','INMON','TREG','Mesenchymal'),#,'INMON','INCAF','TREG','END'),
  nrow=3
)
}
print(a)
# 基因空间表达可视化
# Markers for cancer clone A and B, acinar cell, and centroacinar like ductal cell
b <- SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneExpression",
  spatialFeatures =  c('CD47','THBS2',"ANXA1","FPR1","FPR2",'FPR3',"SAA1",'TGFB1','TGFBR1', 
                       'LGALS9', 'COL6A3','MMP1','ACTA2','CTGF', 'MYH11','FAP'),#c('CD47','THBS2',"ANXA1","FPR1","FPR2",'FPR3',"SAA1",'TGFB1','TGFBR1','TGFBR2','VEGFA','NECTIN2','TIGIT',"CTLA4",'MRC1'),#c("TM4SF1","S100A4","PRSS1","CRISP3"),
  nrow=4
)
print(b)
# ii) CCI and Co-locolization
# # set the path to the in-house breast cancer ST data. User can set the paths to their own data.
# visiumPath <- file.path(system.file(package = "SpaCET"), "extdata/Visium_BC")
# # load ST data to create an SpaCET object.
# SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
# show this object.
# str(SpaCET_obj)
# show this object.
# SpaCET_obj@input$counts[1:8,1:6]
# calculate the QC metrics
# SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
# plot the QC metrics
# SpaCET.visualize.spatialFeature(
#   SpaCET_obj,
#   spatialType = "QualityControl",
#   spatialFeatures=c("UMI","Gene")
# )

# deconvolve ST data
# SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType="BRCA", coreNo=8)
# show the ST deconvolution results
SpaCET_obj@results$deconvolution$propMat[1:13,1:6]
pred <- SpaCET_obj@results$deconvolution$propMat#object[[prediction]]@data
rownames(pred)
if(sample!='GSM6360689'){
pred <- pred[! rownames(pred) %in% c('Epithelia','B', 'CD8T','CD4T','Myeloid','Lymphoid','Mesenchymal','Endothelia'),]
} else{
  pred <- pred[ rownames(pred) %in% c('Epithelia','Mesenchymal','B', 'CD8T','CD4T','Myeloid','Lymphoid','Endothelia'),]
}
# label each spot based on the maximum prediction probability
labels = rownames(pred)[apply(pred, 2, which.max)]
names(labels) <- colnames(pred)
print(table(labels))

# print(table(SPOTlight_obj_Cervix$GSM6360691_obj$labels))
# print(table(SPOTlight_obj_Cervix$GSM6360689_obj$labels))
# show the spatial distribution of malignant cells and macrophages.
# SpaCET.visualize.spatialFeature(
#   SpaCET_obj,
#   spatialType = "CellFraction",
#   spatialFeatures=c("C14","NEUT",'ICAF')
# )

# show the spatial distribution of all cell types.
nrow <- ifelse(sample=='GSM6360689', 5, 6)
c <- SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "CellFraction",
  spatialFeatures="All",
  pointSize = 0.1,
  nrow=nrow
)
print(c)
# calculate the cell-cell colocalization.
# SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)
# 
# # visualize the cell-cell colocalization.
d <- SpaCET.visualize.colocalization(SpaCET_obj)
print(d)
# 
# # calculate the L-R network score across ST spots.
# SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj,coreNo=8)

# visualize the L-R network score.
e <- SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "LRNetworkScore",
  spatialFeatures=c("Network_Score","Network_Score_pv")
)
table(labels)
print(e)
dev.off()

}
# dev.off()
# Ligand-Receptor analysis for a co-localized cell-type pair
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C14","ICAF"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C14","NEUT"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C14","INMON"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C14","INCAF"))

SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C1","ICAF"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C1","NEUT"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C1","INMON"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C1","INCAF"))

SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C5","ICAF"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C5","NEUT"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C5","INCAF"))

SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C8","ICAF"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C8","NEUT"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C8","INCAF"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C8","INMON"))

SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C6","ICAF"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C6","NEUT"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("C6","INCAF"))

SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("INCAF","NEUT"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("INCAF","INMON"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("ICAF","NEUT"))
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("ICAF","INMON"))


## [1] "CAF and Macrophage M2 have potential intercellular interaction in the current tissue."
# Visualize the interaction analysis of a co-localized cell-type pair.
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("C14","INCAF"))
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("C14","ICAF"))
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("INCAF","NEUT"))
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("ICAF","NEUT"))
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("ICAF","INMON"))
if(sample =='GSM6360689'){
  SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("Epithelia","Mesenchymal"))
  SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("Epithelia","Mesenchymal"))
  pairs <- c("Epithelia","Mesenchymal")
  source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
  SI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
  ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.STM_Stromal.Pair.png'), width=12, height=6,bg="white")
  
}

pairs <- c("C5","NEUT")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
SI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C5_NEUT.Pair.png'), width=11, height=6,bg="white")
# pairs <- c("C1","DC")
# source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
# SI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
# ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C1_DC.Pair.png'), width=11, height=6,bg="white")
pairs <- c("C1","NEUT")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
SI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C1_NEUT.Pair.png'), width=11, height=6,bg="white")
pairs <- c("C1","INMON")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
SI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C1_INMON.Pair.png'), width=11, height=6,bg="white")

pairs <- c("C14","NEUT")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
SI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C14_NEUT.Pair.png'), width=11, height=6,bg="white")
pairs <- c("C14","INMON")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
SI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C14_INMON.Pair.png'), width=11, height=6,bg="white")

pairs <- c("C8","NEUT")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
MI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C8_NEUT.Pair.png'), width=11, height=6,bg="white")
pairs <- c("C8","INMON")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
MI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C8_INMON.Pair.png'), width=11, height=6,bg="white")

pairs <- c("C6","NEUT")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
MI.visualize.cellTypePair(SpaCET_obj, cellTypePair= pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.C6_NEUT.Pair.png'), width=11, height=6,bg="white")

pairs <- c("ICAF","NEUT")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
FI.visualize.cellTypePair(SpaCET_obj, cellTypePair=pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.ICAF_NEUT.Pair.png'), width=11, height=6,bg="white")
pairs <- c("ICAF","INMON")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
FI.visualize.cellTypePair(SpaCET_obj, cellTypePair=pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.ICAF_INMON.Pair.png'), width=11, height=6,bg="white")

pairs <- c("INCAF","NEUT")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
FI.visualize.cellTypePair(SpaCET_obj, cellTypePair=pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.INCAF_NEUT.Pair.png'), width=11, height=6,bg="white")

pairs <- c("INCAF","INMON")
source('~/Rcode/UTH/ST.visualize.cellTypePair.r')
FI.visualize.cellTypePair(SpaCET_obj, cellTypePair=pairs)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.visualize.INCAF_INMON.Pair.png'), width=11, height=6,bg="white")

# Enrich cell-cell interactions at the tumor-immune interface
# # Interestingly, besides the interaction significance, we found an enrichment of CAF-M2 interactions close to boundaries between tumor and immune/stromal regions. The distance between CAF-M2 and the tumor-immune border was calculated by averaging the distances between each CAF-M2 interaction spot and its nearest tumor border spot. We randomly selected the same number of spots as CAF-M2 spots from the non-malignant regions and calculated their distances to the border as the null distribution. The result showed that CAF-M2 interaction spots are significantly closer to the tumor-immune boundaries.
# Identify the Tumor-Stroma Interface
propMat <- SpaCET_obj@results$deconvolution$propMat 
Malignant <- apply(propMat[c('C1','C14','C8','C6'),], 2, function(x){
  mal <- sum(x)
  return(mal)
})
Malignant <- as.matrix(t(Malignant))
rownames(Malignant) <- 'Malignant'
SpaCET_obj@results$deconvolution$propMat <- rbind(propMat, Malignant)

SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj, Malignant = 'Malignant')
# Visualize the Interface
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface")
# Combine the interface and interaction spots
SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("C14","NEUT"))
SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("C14","ICAF"))
SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("INCAF","NEUT"))
SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("C1","NEUT"))
SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("C1","ICAF"))
# Visualize the Interface. The spatialFeature should be formated as "Interface&celltype1_celltype2". Celltype 1 and 2 is in the order of alphabet.
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface&INCAF_NEUT")
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface&C14_NEUT")
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface&C14_ICAF")
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface&C1_NEUT")
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface&C1_ICAF")

# Compute the distance of CAF-M2 to tumor border
# SpaCET.distance.to.interface(SpaCET_obj, cellTypePair=c("PSC", "INMON"))
# saveRDS(SpaCET_obj, file.path( paste0(outdir, sample,'_SpaCET_obj.rds')))


# # further deconvolve malignant cell states
# SpaCET_obj <- SpaCET.deconvolution.malignant(SpaCET_obj, malignantCutoff = 0.7, coreNo = 8)
# 
# # show cancer cell state fraction of the first five spots
# SpaCET_obj@results$deconvolution$propMat[c("Malignant cell state A","Malignant cell state B"),1:6]
# 
# ##                           50x102        59x19        14x94     47x13        73x43     61x97
# ## Malignant cell state A 0.2295498 9.999900e-01 6.845962e-02 0.2038680 9.608802e-01 0.6517794
# ## Malignant cell state B 0.0565137 1.239573e-11 3.921715e-08 0.1861075 2.340661e-09 0.2675332
# 
# SpaCET.visualize.spatialFeature(
#   SpaCET_obj,
#   spatialType = "CellFraction",
#   spatialFeatures=c("Malignant","Malignant cell state A","Malignant cell state B"),
#   nrow=1
# )

#######################################################################################################
samples_Cervix = pdata.ST$geo_accession[-1] # Cervix_ABN_ST_CCI.R
outdir = '/data/rluo4/database/Cervix/Guo/vis_seu/'
cellchat <- readRDS(file = paste0(outdir, paste0(samples_Cervix, collapse = '_'), "_ST_cellchat", ".rds"))
table(cellchat@idents)
table(cellchat@meta$labels)
subtype <- data.frame(subset =  table(cellchat@meta$labels))
subtype <- subtype[subtype$subset.Freq !=0,]
subtype <- subtype[subtype$subset.Freq >=5,]
cellchat <-  subsetCellChat(cellchat, idents.use = subtype$subset.Var1)
table(cellchat@meta$labels)
cellchat@meta$labels <- as.character(cellchat@meta$labels)
data.cci.ST_Cervix <- subsetCommunication(cellchat)
# execution.time = Sys.time() - ptm
# print(as.numeric(execution.time, units = "secs"))
# # We can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot or heatmap plot.
# ptm = Sys.time()
pathways.show <- c("ANNEXIN") 
# Circle plot
par(mfrow=c(1,1))
pdf(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/Cervix_SpaCET.ANNEXIN.pdf'), height = 5, width = 5)                    
netVisual_aggregate(cellchat, signaling = pathways.show, layout = 'circle', alpha.image = 0.15, point.size = 2,
                    edge.width.max = 4,  vertex.size.max = 4, vertex.label.cex = 1, arrow.size = 1)
# ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', 'Cervix_SpaCET.ANNEXIN.png'), width=5, height=4,  dpi = 500, limitsize = FALSE)
dev.off()
png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/Cervix_SpaCET.ANNEXIN.png'), height = 5, width = 5, units = 'in',  res = 500)
par(mfrow = c(1,1))
# netVisual_individual(cellchat, signaling = pathways.show,   arrow.size = 0.5, #vertex.weight = 15, edge.weight.max = 10,
#                      layout = "circle" , vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = 'circle', alpha.image = 0.15, point.size = 2,
                    edge.width.max = 4,  vertex.size.max = 4, vertex.label.cex = 1, arrow.size = 1)
dev.off()

pathways.show <- c("RA") 
pdf(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/Cervix_SpaCET.RA.pdf'), height = 5, width = 5)                    
netVisual_aggregate(cellchat, signaling = pathways.show, layout = 'circle', alpha.image = 0.15, point.size = 2,
                    edge.width.max = 4,  vertex.size.max = 4, vertex.label.cex = 1, arrow.size = 1)
dev.off()
png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/Cervix_SpaCET.RA.png'), height = 5, width = 5, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = 'circle', alpha.image = 0.15, point.size = 2,
                    edge.width.max = 4,  vertex.size.max = 4, vertex.label.cex = 1, arrow.size = 1)
dev.off()

# pdf(paste0(outdir, samples_Cervix, '_cellchat.pdf'), height = 10, width = 10)
pdf(paste0(outdir, paste(samples_Cervix, collapse = '.'), '_CellChat.pdf'))
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 120000 * 1024^2)
library(future)
availableCores() #12 #查看几个核可用
# plan(multisession, workers=40)
nbrOfWorkers() #4 当前可用的核有多少个
pid <- Sys.getpid()
pid #2109512
levels(cellchat@idents)
levels <- c(extract_last_element(PCT_Clusters[[organ]]), 
            unique(marker_Cervix$Cell_Type[marker_Cervix$Lineage=='Myeloid']),
            unique(marker_Cervix$Cell_Type[marker_Cervix$Lineage=='Mesenchymal']))

sources.use = na.omit( match(levels, levels(cellchat@idents))  )
vertex.receiver = na.omit( match( levels, levels(cellchat@idents))  )

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
#> Do heatmap based on a single object
#netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
# Part III: Visualization of cell-cell communication network
# Upon infering the cell-cell communication network, CellChat provides various functionality for further data exploration, analysis, and visualization. Here we only showcase the circle plot and the new spatial plot.
# Visualization of cell-cell communication at different levels: One can visualize the inferred communication network of signaling pathways using netVisual_aggregate, and visualize the inferred communication networks of individual L-R pairs associated with that signaling pathway using netVisual_individual.
# Here we take input of one signaling pathway as an example. All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways.

# Compute and visualize the network centrality scores:
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# USER can show this information on the spatial transcriptomics when visualizing a signaling network, e.g., bigger circle indicates larger incoming signaling
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, slice.use = samples_Cervix[1], layout = "spatial", edge.width.max = 3, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 2, vertex.label.cex = 4)
netVisual_aggregate(cellchat, signaling = pathways.show, slice.use = samples_Cervix[2], layout = "spatial", edge.width.max = 3, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 2, vertex.label.cex = 4)
netVisual_aggregate(cellchat, signaling = pathways.show, slice.use = samples_Cervix[3], layout = "spatial", edge.width.max = 3, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 2, vertex.label.cex = 4)

png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', samples_Cervix[1], '_SpaCET.ANNEXIN.png'), height = 6, width = 6, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  title.space = 0.1, show.legend = FALSE,sources.use = sources.use, targets.use = vertex.receiver,
                    slice.use = samples_Cervix[1], layout = "spatial", edge.width.max = 1, arrow.size = 0.25, point.size = 2.5, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 1, vertex.label.cex = 4)
dev.off()

png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', samples_Cervix[2], '_SpaCET.ANNEXIN.png'), height = 6, width = 6, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = sources.use, targets.use = vertex.receiver, title.space = 0.1, show.legend = FALSE,
                    slice.use = samples_Cervix[2], layout = "spatial", edge.width.max = 1, arrow.size = 0.25, point.size = 2.5, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 1, vertex.label.cex = 4)
dev.off()
png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', samples_Cervix[3], '_SpaCET.ANNEXIN.png'), height = 6, width = 6, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = sources.use, targets.use = vertex.receiver, title.space = 0.1, show.legend = FALSE, 
                    slice.use = samples_Cervix[3], layout = "spatial", edge.width.max = 1, arrow.size = 0.25, point.size = 2.5, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 1, vertex.label.cex = 4)
dev.off()

pathways.show <- c("RA") 
png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', samples_Cervix[1], '_SpaCET.RA.png'), height = 6, width = 6, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  title.space = 0.1, show.legend = FALSE,sources.use = sources.use, targets.use = vertex.receiver,
                    slice.use = samples_Cervix[1], layout = "spatial", edge.width.max = 1, arrow.size = 0.25, point.size = 2.5, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 1, vertex.label.cex = 4)
dev.off()

png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', samples_Cervix[2], '_SpaCET.RA.png'), height = 6, width = 6, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = sources.use, targets.use = vertex.receiver, title.space = 0.1, show.legend = FALSE,
                    slice.use = samples_Cervix[2], layout = "spatial", edge.width.max = 1, arrow.size = 0.25, point.size = 2.5, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 1, vertex.label.cex = 4)
dev.off()
png( paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', samples_Cervix[3], '_SpaCET.RA.png'), height = 6, width = 6, units = 'in',  res = 500)
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = sources.use, targets.use = vertex.receiver, title.space = 0.1, show.legend = FALSE, 
                    slice.use = samples_Cervix[3], layout = "spatial", edge.width.max = 1, arrow.size = 0.25, point.size = 2.5, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 1, vertex.label.cex = 4)
dev.off()

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
# When visualizing gene expression distribution on tissue using spatialFeaturePlot, users also need to provide the slice.use as an input.

# Take an input of a few genes
spatialFeaturePlot(cellchat, features = c("ANXA1","FPR3"), slice.use = samples_Cervix[1], point.size = 0.8, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("ANXA1","FPR3"), slice.use = samples_Cervix[2], point.size = 0.8, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("ANXA1","FPR3"), slice.use = samples_Cervix[3], point.size = 0.8, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("TGFB1","TGFBR1_TGFBR2"), slice.use = samples_Cervix[1],point.size = 0.8, color.heatmap = "Reds", direction = 1)
# Take an input of a ligand-receptor pair
spatialFeaturePlot(cellchat, pairLR.use = "TGFB1_TGFBR1_TGFBR2", slice.use = samples_Cervix[1], point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "TGFB1_TGFBR1_TGFBR2", slice.use = samples_Cervix[2], point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "TGFB1_TGFBR1_TGFBR2", slice.use = samples_Cervix[3], point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
#> Applying a cutoff of  0.05 to the values...
# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat, pairLR.use = "ANXA1_FPR3", slice.use = samples_Cervix[1], point.size = 1.5, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "ANXA1_FPR3", slice.use = samples_Cervix[2], point.size = 1.5, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "ANXA1_FPR3", slice.use = samples_Cervix[3], point.size = 1.5, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
dev.off()

saveRDS(cellchat, file = paste0(outdir, paste0(samples_Cervix, collapse = '_'), "_ST_cellchat", ".rds"))

# NB: Upon infering the intercellular communication network from spatial transcriptomics data, CellChat’s various functionality can be used for further data exploration, analysis, and visualization. Please check other functionality in the basic tutorial named CellChat-vignette.html
# Part V: Save the CellChat object
# saveRDS(cellchat, file = paste0(outdir,sample, "_ST_cellchat", ".rds"))
# 1，celltype之间通讯结果
sample = 'GSM6360689'
# 1）根据使用netVisual_circle显示任意两个celltype之间的通讯次数（左）或总通讯强度(右)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight),
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# 2）根据使用netVisual_heatmap显示任意两个celltype之间的通讯次数（左）或总通讯强度(右)
p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

p1 + p2
# 2，单个信号通路可视化
# 首先根据cellchat@netP$pathways展示当前有哪些通路结果，选择感兴趣的进行展示，此处示例展示SPP1通路。
levels(cellchat@idents) # 查看当前的celltype顺序，然后可以通过vertex.receiver指定target 的细胞类型。
# 1）层级图
# 绘制层级图的话 ，需要指定layout为hierarchy ，当前版本默认下出来的是circle图。
sort(cellchat@netP$pathways)
spatialFeaturePlot(cellchat, features = c("ANXA1","FPR3"), slice.use = sample, point.size = 0.8, color.heatmap = "Reds", direction = 1)
ggsave(file = paste0('/data/rluo4/database/Cervix/Guo/vis_seu/', sample, '_SpaCET.CellChat.Pair.png'), width=7, height=6,bg="white")

pathways.show <- c("ANNEXIN")
levels(cellchat@idents)
#[1] "C1"   "C14"  "C8"   "MAST" "MON"  "NEUT"
vertex.receiver = c(2:6)  #选择的是levels(cellchat@idents) 中的
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
# 2）和弦图
# 可以额外绘制空间转录组版本的和弦图，添加layout = "spatial" 。
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1,
                    alpha.image = 0.2, vertex.label.cex = 3.5)
# 3） network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 8, height = 2.5, font.size = 10)

# 3，绘制配体受体气泡图
# 1）指定受体-配体细胞类型
# 绘制指定受体-配体细胞类型中的全部配体受体结果的气泡图，通过sources.use 和 targets.use指定。

#指定受体-配体细胞类型
netVisual_bubble(cellchat, sources.use = c(3,5),
                 targets.use = c(1,2,4,6), remove.isolate = FALSE)

# 2）指定受体-配体细胞类型 且 指定通路
# 同时通过signaling指定展示通路

netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6),
                 signaling = c("TGFb","SPP1",'ANNEXIN'), remove.isolate = FALSE)

# 3）ligand-receptor pair 表达
# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat, pairLR.use = "ANXA1_FPR3", point.size = 1.5,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)



##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) #######------RCTD pipeline ----#######
library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)
indir = '/data/rluo4/database/Cervix/Guo/spaceranger_count/'
outdir = '/data/rluo4/database/Cervix/Guo/vis_seu/'

organ = 'Cervix'
path_PCT <- paste0('/data/rluo4/All/Output/PCT/')#('/data/rluo4/database/',organ,'/Output/')
file = paste0(path_PCT, 'PCT-',organ,'.rds')
load(file)
unique(rdata_filter$cell_class)
rdata_filter$cell_subtype <- paste0("C",rdata_filter$leiden_res2)#paste(group_ABN$Sample_Type,group_ABN$Tissue,sep = ':')

# saveRDS(myRCTD, file.path( paste0(outdir, sample,'_RCTD_visium_full.rds')))
myRCTD <- readRDS( paste0(outdir, sample,'_RCTD_visium_full.rds'))
# 看看myRCTD的结构：
# 
str(myRCTD)
# ## result
# myRCTD <- readRDS(file.path(savedir,'myRCTD_visium_full.rds'))
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)

# 看一下norm_weights：每一行代表一个spots，每一列代表一个细胞类型，值表示这个spot中细胞类型所占的比例。每一行的和为1。
# observe weight values
cell_types <- c('C5','C14','C8','C6','NEUT','INMON','ICAF')#c('Denate', 'Neurogenesis','Cajal_Retzius')
print(head(norm_weights[,cell_types]))

head(norm_weights)
assignLabels <- function(object, prediction = "predictions") {
  weights <- object@results$weights
  pred <- normalize_weights(weights)
  pred <- t(pred)#pred[1:(nrow(pred)-1), ]
  # label each spot based on the maximum prediction probability
  labels = rownames(pred)[apply(pred, 2, which.max)]
  names(labels) <- colnames(pred)
  table(labels)
  
  object$labels <- factor(labels)
  Idents(object) <- "labels"
  return(object)
}
st <- assignLabels(st, prediction = "predictions")
table(st$labels)
# 有了这个结果，再加上空间对象的rds文件，可以借助其他反卷积的绘图函数绘制出其他反卷积结果中一样的那种饼图了，比如SPOTlight软件可以绘制的饼图：
#
# 教程中可以展示某一个细胞类型在空间数据中反卷积出来的比例可视化：
#
# 总共有17中细胞类型，这里展示Denate

# plot Dentate weights
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'ICAF'], ylimit = c(0,0.5),
                          title ='plot of Dentate weights', size=4.5, alpha=0.8)
p
# ggsave(paste0(savedir, "/Spaital_weights.png"), width=8, height=6, plot=p,bg="white")
coords <- read.csv(file.path(paste0(indir, sample, '/outs/spatial/tissue_positions.csv')))#
rownames(coords) <- coords[,1]; coords <- coords[,3:4]#coords[,1] <- NULL
colnames(coords) <- c('x','y')
# 结果图：
# wget -c https://www.bioconductor.org/packages/release/bioc/src/contrib/STdeconvolve_1.6.0.tar.gz
# conda activate py_env
# install.packages(c(‘topicmodels’, ‘liger’))
# R CMD INSTALL --configure-vars='ICUDT_DIR=/home/rluo4/R/R-source'  STdeconvolve_1.6.0.tar.gz --library=/home/rluo4/R/x86_64-pc-linux-gnu-library/4.1

# 借用STdeconvolve绘制一下饼图看看：
library(STdeconvolve)
library(ggplot2)
library(ggsci)
packageVersion("STdeconvolve")

m <- as.matrix(norm_weights)
p <- coords

plt <- vizAllTopics(theta = m,
                    pos = p,
                    topicOrder=seq(ncol(m)),
                    topicCols=rainbow(ncol(m)),
                    groups = NA,
                    group_cols = NA,
                    r = 56, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0.3,
                    showLegend = TRUE,
                    plotTitle = "scatterpies")

## function returns a `ggplot2` object, so other aesthetics can be added on:
plt <- plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))
plt
# ggsave(paste0(savedir, "/Spaital_scatterpies.png"), width=12, height=6, plot=plt, bg="white")

# 
sessionInfo()
# # 
# # options(stringsAsFactors = FALSE)
# # devtools::install_github("immunogenomics/presto")
# # devtools::install_github("jinworks/CellChat")
# library(CellChat)
# library(Seurat)
# library(tidyverse)
# library(viridis)
# library(RColorBrewer)
# #载入数据
# #查看数据情况
# # st@meta.data$celltype <- paste("C",st$seurat_clusters,sep = "")
# Idents(st) <- st["region"]
# head(st)
# #可定义颜色
# color.use <- scPalette(nlevels(st))
# names(color.use) <- levels(st)
# SpatialDimPlot(st, label = TRUE, label.size = 3, cols = color.use)
# # 和单细胞区别之一：空间转录组除矩阵和meta外，还需要额外输入空间图像信息 以及 Scale factor信息 。
# #矩阵信息
# data.input = Seurat::GetAssayData(st, slot = "data", assay = "SCT") 
# #meta信息
# meta = data.frame(labels = Idents(st), #名字自定义
#                   row.names = names(Idents(st))) # manually create a dataframe consisting of the cell labels
# unique(meta$labels)
# 
# # 空间图像信息
# spatial.locs = Seurat::GetTissueCoordinates(st, scale = NULL,  cols = c("imagerow", "imagecol")) 
# # spatial.locs <- rbind(spatial.locs1, spatial.locs2)
# # rownames(spatial.locs) <- colnames(data.input)
# # same as head(dt@images$slice1@coordinates, 2)
# # Scale factors and spot diameters 信息 
# # scale.factors = jsonlite::fromJSON(txt = file.path(paste0(indir, sample, '/outs/spatial/scalefactors_json.json')))
# # 
# # scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
# #                      fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
# # )
# scalefactors1 = jsonlite::fromJSON(txt = file.path(paste0(indir, sample, '/outs/spatial/scalefactors_json.json')))#jsonlite::fromJSON(txt = file.path("/Users/suoqinjin/Library/CloudStorage/OneDrive-Personal/works/CellChat/tutorial/spatial_imaging_data-intestinalA1", 'scalefactors_json.json'))
# spot.size = 65 # the theoretical spot size (um) in 10X Visium
# conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
# spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)
# 
# # scalefactors2 = jsonlite::fromJSON(txt = file.path("/Users/suoqinjin/Library/CloudStorage/OneDrive-Personal/works/CellChat/tutorial/spatial_imaging_data-intestinalA2", 'scalefactors_json.json'))
# # conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
# # spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)
# 
# # spatial.factors <- rbind(spatial.factors1, spatial.factors2)
# # rownames(spatial.factors) <- c("A1", "A2")
# 
# spatial.factors <- spatial.factors1
# rownames(spatial.factors) <- sample#c("A1", "A2")
# 
# #### create a CellChat object ####
# cellchat <- createCellChat(object = data.input, 
#                            meta = meta, 
#                            group.by = "labels", #前面的meta ，定义的名字是labels
#                            datatype = "spatial", ###
#                            coordinates = spatial.locs, 
#                            spatial.factors = spatial.factors)
# # 注意datatype 要选择 "spatial" 。
# # 因为空转数据是人的，这里直接选择CellChatDB.human（鼠的话选择 CellChatDB.mouse） 。
# 
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# # use a subset of CellChatDB for cell-cell communication analysis
# #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# # use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# 
# # set the used database in the object
# cellchat@DB <- CellChatDB.use
# # 注1：使用全部的用于cellchat分析，也可以不进行subsetDB，直接指定cellchat@DB <- CellChatDB 即可。
# # 
# # 注2：如果你有关心的配受体对 且 不在该数据库中，也可以自行添加上。大概步骤就是下载对应的csv（数据库），在对应的列上添加上你的配受体对信息，保存后重新读取新的csv即可，详细见https://htmlpreview./?https://github.com/sqjin/CellChat/blob/master/tutorial/Update-CellChatDB.html。
# # 
# # 5，CellChat预处理
# # 可以使用subsetData选择进行cellchat的子集，注意使用全集的话也要subsetData一下。
# 
# # subset the expression data of signaling genes for saving computation cost
# cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# future::plan("multisession", workers = 1) #笔记本可以选1
# ##识别过表达基因
# cellchat <- identifyOverExpressedGenes(cellchat)
# #识别过表达配体受体对
# cellchat <- identifyOverExpressedInteractions(cellchat)
# 
# 
# # 1，推断细胞通讯网络
# # 使用表达值推测细胞互作的概率
# cellchat <- computeCommunProb(cellchat, 
#                               type = "truncatedMean", trim = 0.1, 
#                               distance.use = TRUE, 
#                               scale.distance = 0.01)
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# # 注1：type 默认为triMean，producing fewer but stronger interactions; 当设置'type = "truncatedMean"', 需要跟 'trim'参数 , producing more interactions.
# # 
# # 注2：distance.use = FALSE 会过滤掉较远空间距离之间的交互
# # 
# # 注3：空转默认的是truncatedMean 和 trim组合，根据经验进行适当的调整。
# # 
# # 2，计算cell-cell communication
# # 使用computeCommunProbPathway计算每个信号通路的所有配体-受体相互作用的通信结果，结存存放在net 和 netP中 。
# # 
# # 然后使用aggregateNet计算细胞类型间整合的细胞通讯结果。
# 
# #计算每个信号通路相关的所有配体-受体相互作用的通信结果
# cellchat <- computeCommunProbPathway(cellchat)
# #计算整合的细胞类型之间通信结果
# cellchat <- aggregateNet(cellchat)
# 
# 
# # 1，celltype之间通讯结果
# # 1）根据使用netVisual_circle显示任意两个celltype之间的通讯次数（左）或总通讯强度(右)
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
#                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
#                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# 
# # 2）根据使用netVisual_heatmap显示任意两个celltype之间的通讯次数（左）或总通讯强度(右)
# p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
# 
# p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
# 
# p1 + p2
# # 2，单个信号通路可视化
# # 首先根据cellchat@netP$pathways展示当前有哪些通路结果，选择感兴趣的进行展示，此处示例展示SPP1通路。
# levels(cellchat@idents) # 查看当前的celltype顺序，然后可以通过vertex.receiver指定target 的细胞类型。
# # 1）层级图
# # 绘制层级图的话 ，需要指定layout为hierarchy ，当前版本默认下出来的是circle图。
# cellchat@netP$pathways
# pathways.show <- c("SPP1")
# levels(cellchat@idents)
# #[1] "C9" "C3" "C5" "C4" "C0" "C7" "C6" "C8" "C1" "C2"
# vertex.receiver = c(1,2,4,5)  #选择的是levels(cellchat@idents) 中的
# netVisual_aggregate(cellchat, signaling = pathways.show,
#                     vertex.receiver = vertex.receiver,layout = "hierarchy")
# # 2）和弦图
# # 可以额外绘制空间转录组版本的和弦图，添加layout = "spatial" 。
# # Circle plot
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# # Spatial plot
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
#                     edge.width.max = 2, vertex.size.max = 1,
#                     alpha.image = 0.2, vertex.label.cex = 3.5)
# # 3） network centrality scores
# # Compute the network centrality scores
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# par(mfrow=c(1,1))
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
#                                   width = 8, height = 2.5, font.size = 10)
# 
# # 3，绘制配体受体气泡图
# # 1）指定受体-配体细胞类型
# # 绘制指定受体-配体细胞类型中的全部配体受体结果的气泡图，通过sources.use 和 targets.use指定。
# 
# #指定受体-配体细胞类型
# netVisual_bubble(cellchat, sources.use = c(3,5),
#                  targets.use = c(1,2,4,6), remove.isolate = FALSE)
# 
# # 2）指定受体-配体细胞类型 且 指定通路
# # 同时通过signaling指定展示通路
# 
# netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6),
#                  signaling = c("TGFb","SPP1",'ANNEXIN'), remove.isolate = FALSE)
# 
# # 3）ligand-receptor pair 表达
# # Take an input of a ligand-receptor pair and show expression in binary
# spatialFeaturePlot(cellchat, pairLR.use = "ANXA1_FPR3", point.size = 1.5,
#                    do.binary = TRUE, cutoff = 0.05, enriched.only = F,
#                    color.heatmap = "Reds", direction = 1)
# 
# 


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) #######------CARD pipeline ----#######
### 安装CARD包
# conda activate py_env
#####w conda install -c conda-forge udunits2 # locate libudunits2.so.0 # /data/rluo4/lorihan/R/site-library/units/libs
config <- c(units="--with-udunits2-lib=/data/rluo4/bin/miniconda3/envs/py_env/lib --with-udunits2-include=/data/rluo4/bin/miniconda3/envs/py_env/include")
# config <- c(units="--with-udunits2-lib=/data/rluo4/bin/anaconda3/envs/r_env/lib --with-udunits2-include=/data/rluo4/bin/anaconda3/envs/r_env/include")
install.packages("sf", configure.args = config )
# install.packages("RcppGSL", configure.args = config )
devtools::install_github('YingMa0107/CARD',  force = TRUE)

library(dplyr,lib.loc='/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(CellChat,lib.loc='/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData,lib.loc='/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(dplyr)
library(RColorBrewer)
library(ArchR, lib.loc='/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)
library(future)
plan(multisession, workers=40)
availableCores() #12 #查看几个核可用
nbrOfWorkers() #4 当前可用的核有多少个
pid <- Sys.getpid()
pid #2109512
library(tidyverse)
library(rhdf5)
library(data.table)
library(SeuratDisk, lib.loc='/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')
library(data.table)
library(stringr)
library(dplyr)
load(file = "/home/lorihan/lrh/database/Cervix/Guo/Guo_ST_anno.RData")
print(table(pdata.ST$geo_accession))
indir = '/home/lorihan/lrh/database/Cervix/Guo/vis_seu/'
load('/home/lorihan/lrh/database/Cervix/Guo/scRNA_reference.RData')
library(CARD)
library(Seurat)
library(patchwork)
library(tidyverse)
for (sample in unique(pdata.ST$geo_accession)) {
  setwd(indir)
  sample = 'GSM6360691'
  load(paste0(indir, sample, '_CARD.RData'))
  assignLabels <- function(object, prediction = "predictions") {
    pred <- t(CARD_obj@Proportion_CARD)#object[[prediction]]@data
    # pred <- pred[1:(nrow(pred)-1), ]
    # label each spot based on the maximum prediction probability
    labels = rownames(pred)[apply(pred, 2, which.max)]
    names(labels) <- colnames(pred)
    print(table(labels))
    CARD_obj$labels <- factor(labels)
    Idents(CARD_obj) <- "labels"
    return(CARD_obj)
  }
  
  CARD_obj <- assignLabels(CARD_obj, prediction = "Proportion_CARD")
  table(CARD_obj$labels)
  subtype <-as.data.frame(table(seu1$labels))
  subtype <- subtype[subtype$Freq >=10,]
  seu1 <- seu1[, seu1$labels %in% subtype$Var1]
  
  
  
pdf( paste0(sample, '_CARD.pdf'), height = 10, width = 14 )
# 3.1，spot比例可视化
# 使用CARD.visualize.pie 函数绘制spot的细胞类型分布饼图
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,
                         spatial_location = CARD_obj@spatial_location, 
                         colors = colors)
p1
# 3.2，展示感兴趣的细胞类型
# 
# 重点关注各细胞类型在空间中的分布，这也是空转特色

## select the cell type that we are interested
ct.visualize = c('C5','C1','C14','C8','C6','C10','C4',"ICAF","INCAF","TREG",'NEUT','INMON','M2MAC')
## visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4)                                 ### number of columns in the figure panel
print(p2)


# 3.3，可视化两种细胞类型
# 
# 可以同时展示2种细胞类型，更清晰的对比空间位置。

## visualize the spatial distribution of two cell types on the same plot
p3 = CARD.visualize.prop.2CT(
  proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
  spatial_location = CARD_obj@spatial_location,                  ### two cell types you want to visualize
  ct2.visualize = c("ICAF","C5"),
  colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales                             
p3
# 3.4，细胞类型比例相关图
p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
p4

# 4.1，imputation 函数推断精度
#1. Imputation on the newly grided spatial locations #
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)

# 这里的NumGrids 自定义，值越大推断的精度越高
# 
# 提升前后的空间位置网格轮廓

p5 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))

p50 <- ggplot(CARD_obj@spatial_location, 
              aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
p50 + p5

# 4.2，分辨率增强-细胞类型比例可视化
p6 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("lightblue","lightyellow","red"),    
  NumCols = 4)                                  
p6

# 4.3，分辨率增强-marker gene 可视化
#增强后空间位置 
interest_gene = c('CD47','THBS2',"ANXA1","FPR1","FPR2",'FPR3',"SAA1",'TGFB1','TGFBR1','TGFBR2','VEGFA','NECTIN2','TIGIT',"CTLA4",'MRC1')
genes.common <- intersect(rownames(CARD_obj@refined_expression), interest_gene)  
p7 <- CARD.visualize.gene(spatial_expression = CARD_obj@refined_expression,
                          spatial_location = location_imputation,
                          gene.visualize = genes.common,
                          colors = NULL,
                          NumCols = 6)

# 原始空间位置 
p77 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = genes.common,
  colors = NULL,
  NumCols = 6)
p7 / p77
dev.off()
}