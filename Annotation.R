
# Load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(SingleR)
library(celldex)
library(gridExtra)
library(Cairo)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(celldex)
library(gridExtra)
library(tidyr)
library(data.table)
library(Matrix)
library(Cairo)
library(assorthead)
library(SingleR)
Dataset  <- readRDS("/data_storage/datasets/all_combo_with_UMAP_PCs.rds")
dim(Dataset@meta.data)
head(Dataset@meta.data)


p1 = DimPlot(Dataset,split.by="Source", repel=T,raster=FALSE)
png("/data_storage/basic_dimplot_by_source_Dataset.png",res=300,units="in",width=8,height=8)
p1
dev.off()

###############################################################################
table(Dataset@meta.data$ann_blueprint)

p1 = DimPlot(Dataset,group.by="ann_blueprint", repel=T,raster=FALSE)
png("/data_storage/basic_dimplot_by_ann_blueprint_Dataset.png",res=300,units="in",width=8,height=8)
p1
dev.off()

###############################################################################
# read in celltypist preds
lowres_preds_list = read_csv("data_storage/Cell_typist/celltypist_predictions_low.csv")
highres_preds_list = read_csv("data_storage/Cell_typist/celltypist_predictions_high.csv")

# add to main Seurat object metadata
Dataset[['ann_celltypist_highres']] = highres_preds_list$predicted_labels
Dataset[['ann_celltypist_lowres']] = lowres_preds_list$predicted_labels
Dataset[['ann_celltypist_highres_conf_score']] = highres_preds_list$conf_score
Dataset[['ann_celltypist_lowres_conf_score']] = lowres_preds_list$conf_score

######################################

# Celltypist highres
totals = Dataset@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_celltypist_highres) %>% summarise(total = sum(n))
ct_calls = Dataset@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_celltypist_highres) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_celltypist_highres))%>% slice_max(n=1,order_by=prop)
write_csv(ct_calls,"/data_featherstone/Mohammad/SCRNA/Batch_11/QC/basic_qc/celltypist_highres_cluster_calls_Dataset.csv")
p=ggplot(ct_calls,aes(seurat_clusters,prop,fill=ann_celltypist_highres,label=ann_celltypist_highres))+geom_col()+geom_text(angle=90)+theme_bw()
png("data_storage/celltypist_highres_calls.png",res=300,units="in",width=15,height=8)
p
dev.off()

###################################################
# Plot celltypist annotation
p1 = DimPlot(Dataset,group.by="ann_celltypist_highres", repel=T,raster=FALSE) 
png("/data_storage/basic_dimplot_celltypist.png",res=300,units="in",width=25,height=8) #CG
p1
dev.off()

p1 = DimPlot(Dataset,group.by="ann_celltypist_lowres", repel=T,raster=FALSE)
png("data_storage/basic_dimplot_celltypist_low.png",res=300,units="in",width=25,height=8)
p1
dev.off()

####################################
# SingleR annotation
####################################

blueprint = readRDS("/data_storage/SingleR_Reference/blueprint.RDS")
#monaco = readRDS("/data_storage/SingleR_Reference/monaco.RDS")
#hpca = readRDS("/data_storage/SingleR_Reference/hpca.RDS")
#dice = readRDS("/data_storage/SingleR_Reference/dice.RDS")

blueprint_annotations = SingleR(test = all_combo@assays$RNA@data, ref = blueprint, labels = blueprint$label.fine)
all_combo[["ann_blueprint"]] = blueprint_annotations$pruned.labels
#monaco_annotations = SingleR(test = all_combo@assays$RNA@data, ref = monaco, labels = monaco$label.fine)
#all_combo[["ann_monaco"]] = monaco_annotations$pruned.labels
#hpca_annotations = SingleR(test = all_combo@assays$RNA@data, ref = hpca, labels = hpca$label.fine)
#all_combo[["ann_hpca"]] = hpca_annotations$pruned.labels
#dice_annotations = SingleR(test = all_combo@assays$RNA@data, ref = dice, labels = dice$label.fine)
#all_combo[["ann_dice"]] = dice_annotations$pruned.labels


####################################
# Plot different annotations
####################################

# Celltypist highres
totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_celltypist_highres) %>% summarise(total = sum(n))
ct_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_celltypist_highres) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_celltypist_highres))%>% slice_max(n=1,order_by=prop)
write_csv(ct_calls,"celltypist_highres_cluster_calls.csv")
p=ggplot(ct_calls,aes(seurat_clusters,prop,fill=ann_celltypist_highres,label=ann_celltypist_highres))+geom_col()+geom_text(angle=90)+theme_bw()
png("../cluster_plots/celltypist_highres_calls.png",res=300,units="in",width=15,height=8)
p
dev.off()

# Celltypist lowres
totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_celltypist_lowres) %>% summarise(total = sum(n))
ct_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_celltypist_lowres) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_celltypist_lowres))%>% slice_max(n=1,order_by=prop)
write_csv(ct_calls,"celltypist_lowres_cluster_calls.csv")
p=ggplot(ct_calls,aes(seurat_clusters,prop,fill=ann_celltypist_lowres,label=ann_celltypist_lowres))+geom_col()+geom_text(angle=90)+theme_bw()
png("../cluster_plots/celltypist_lowres_calls.png",res=300,units="in",width=15,height=8)
p
dev.off()

# Blueprint
totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_blueprint) %>% summarise(total = sum(n))
singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_blueprint) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_blueprint))%>% slice_max(n=1,order_by=prop)
write_csv(singler_calls,"blueprint_cluster_calls.csv")
p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_blueprint,label=ann_blueprint))+geom_col()+geom_text(angle=90)+theme_bw()
png("../cluster_plots/blueprint_calls.png",res=300,units="in",width=15,height=8)
p
dev.off()

#Monaco
#totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_monaco) %>% summarise(total = sum(n))
#singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_monaco) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_monaco))%>% slice_max(n=1,order_by=prop)
#write_csv(singler_calls,"monaco_cluster_calls.csv")
#p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_monaco,label=ann_monaco))+geom_col()+geom_text(angle=90,hjust=1)+theme_bw()
#png("../cluster_plots/monaco_calls.png",res=300,units="in",width=20,height=12)
#p
#dev.off()

#HPCA
#totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_hpca) %>% summarise(total = sum(n))
#singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_hpca) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_hpca))%>% slice_max(n=1,order_by=prop)
#write_csv(singler_calls,"hpca_cluster_calls.csv")
#p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_hpca,label=ann_hpca))+geom_col()+geom_text(angle=90)+theme_bw()
#png("../cluster_plots/hpca_calls.png",res=300,units="in",width=20,height=8)
#p
#dev.off()

#Dice
#totals = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_dice) %>% summarise(total = sum(n))
#singler_calls = all_combo@meta.data %>% group_by(seurat_clusters) %>% dplyr::count(ann_dice) %>% left_join(totals,by="seurat_clusters") %>% mutate(prop = n/total*100) %>% filter(prop > 0.05) %>% filter(!is.na(ann_dice))%>% slice_max(n=1,order_by=prop)
#write_csv(singler_calls,"dice_cluster_calls.csv")
#p=ggplot(singler_calls,aes(seurat_clusters,prop,fill=ann_dice,label=ann_dice))+geom_col()+geom_text(angle=90)+theme_bw()
#png("../cluster_plots/dice_calls.png",res=300,units="in",width=20,height=8)
#p
#dev.off()


# Plot celltypist annotation
p1 = DimPlot(all_combo,group.by="ann_celltypist_highres", repel=T) 
png("../cluster_plots/basic_dimplot_celltypist.png",res=300,units="in",width=25,height=8) 
p1 
dev.off() 
p1 = DimPlot(all_combo,group.by="ann_celltypist_lowres", repel=T) 
png("../cluster_plots/basic_dimplot_celltypist_low.png",res=300,units="in",width=25,height=8) 
p1 
dev.off() 

####################################
# cluster biomarkers
####################################

# Find markers
all_combo_markers = FindAllMarkers(all_combo, min.pct=0.25, logfc.threshold = 0.25,only.pos=TRUE,recorrect_umi=FALSE)
write_csv(all_combo_markers,"cluster_biomarkers.csv")
#all_combo_markers = read_csv("cluster_biomarkers.csv")

# Prepare for manual labelling of cluster IDs
cluster_info = all_combo@meta.data[,c("seurat_clusters","ann_celltypist_highres","ann_celltypist_lowres")]
cluster_info1 = table(cluster_info[,1:2])
cluster_info2 = table(cluster_info[,c(1,3)])
write.table(cluster_info1,"cluster_info1.tsv",sep="\t",quote=F)
write.table(cluster_info2,"cluster_info2.tsv",sep="\t",quote=F)

