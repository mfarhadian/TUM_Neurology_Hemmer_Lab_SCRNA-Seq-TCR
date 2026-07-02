
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


# set WD
setwd(paste0(wd,"/datasets/")) 

# read in outputs from deconvolution & integration step
all_combo = readRDS("individual_post_qc_datasets_new_v5_integrated.rds")

# harmony vs pcs
p1=DimPlot(all_combo,reduction="pca",group.by="Source",raster=FALSE) 
p2=DimPlot(all_combo,reduction="harmony",group.by="Source",raster=FALSE) 
png(res=300,units="in",height=8,width=16,file="harmony_vs_pcs_Source.png")
grid.arrange(p1,p2,nrow=2)
dev.off()

p3=DimPlot(all_combo,reduction="pca",group.by="Batch_id",pt.size = 0.5,raster=FALSE) + NoLegend()
p4=DimPlot(all_combo,reduction="harmony",group.by="Batch_id",pt.size = 0.5,raster=FALSE) + NoLegend()
png(res=300,units="in",height=8,width=16,file="harmony_vs_pcs_BatchID.png")
grid.arrange(p3,p4,nrow=1)
dev.off()

# Elbow plot
CairoPNG(res=300,units="in",height=8,width=20,file="elbow_harmony.png")
ElbowPlot(all_combo,reduction="harmony")
dev.off()

# get var explained
eigenvalues = all_combo@reductions$harmony@stdev^2

centiles = c(50,60,70,80,90,95,99)
sapply(centiles,function(x){
  no_pcs = which(cumsum(eigenvalues) / sum(eigenvalues) * 100 > x)[1]
  return(no_pcs)
})

# see how PCs / Harmony PCs are influenced by batch and source
harmony_embeddings = data.frame(all_combo@reductions$harmony@cell.embeddings)
harmony_embeddings = harmony_embeddings %>% mutate(cell_name = rownames(harmony_embeddings))
pc_embeddings = data.frame(all_combo@reductions$pc@cell.embeddings)
pc_embeddings = pc_embeddings %>% mutate(cell_name = rownames(pc_embeddings))
meta_data = data.frame(all_combo@meta.data)
meta_data = meta_data %>% mutate(cell_name = rownames(meta_data))
joint_df = meta_data %>% left_join(harmony_embeddings,by="cell_name") %>% left_join(pc_embeddings,by="cell_name")

p1=ggplot(joint_df,aes(Batch_id,harmony_1,fill=Source))+geom_violin()+theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_brewer(palette="Set2")
p2=ggplot(joint_df,aes(Batch_id,PC_1,fill=Source))+geom_violin()+theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_brewer(palette="Set2")
png(res=300,units="in",height=8,width=16,file="harmony_vs_pcs_embeddings.png")
grid.arrange(p1,p2,ncol=1)
CairoPNG(res=300,units="in",height=8,width=8,file="harmony_vs_pcs_embeddings.png")
dev.off()

#################################
# remove RBCs
combo <- all_combo

DefaultAssay(combo) = "RNA"
# step1: Calculate % of RBC gene expression per cell
rbc_genes=c("HBB","HBA1","HBA2")
combo$percentage.rbc= PercentageFeatureSet(combo, features = rbc_genes)
summary(combo$percentage.rbc)
table(combo$percentage.rbc > 0.0000005)

# Step3: Filter out RBC-Contamination cells

combo = subset(combo, subset=percentage.rbc < 0.0000005)
dim(combo@meta.data)
#############################

# UMAP over a range of PCs and resolutions
# repeat with different parameters to test robustness of clusters
umap_change_pcs = function(x,y){
  message(paste0("Start time:",date()))
  message(paste0("Running UMAP & clustering with ",x," PCs and a resolution of ",y))
  all_combo_test = all_combo %>% RunUMAP(reduction="harmony", dims=1:x) %>%
    FindNeighbors(reduction="harmony",dims=1:x) %>%
    FindClusters(resolution=y,group.singletons = FALSE)
  message(paste0("Finished UMAP & clustering with ",x,"PCs and a resolution of",y,". Moving on."))
  message(paste0("End time:",date()))
  CairoPNG(res=300,units="in",height=8,width=12,file=paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,".png"))
  print(DimPlot(all_combo_test,label=T,raster=FALSE))
  dev.off()
  saveRDS(all_combo_test,paste0("all_combo_with_UMAP_PCs_",x,"resolution",y,".rds"))
}

mat = matrix(c(15,15,15,15,15,10,10,10,10,10,
               2.5,2.0,1.5,1.0,0.5,2.5,2.0,1.5,1.0,0.5),ncol=2)
mapply(umap_change_pcs,x=mat[,1],y=mat[,2]) # takes several hours - not routinely run

