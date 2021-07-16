#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(data.table)
pony_colors<-fread("~/papers/TWAS/pony_palette")

xCell_results<-args[1]
tissue_sampleIDs<-args[2]

###############################################3
prep_columns<-function(dat) {
  CD4<-colMeans(dat[c("CD4+ memory T-cells","CD4+ naive T-cells","CD4+ T-cells","CD4+ Tcm","CD4+ Tem"),])
  CD8<-colMeans(dat[c("CD8+ naive T-cells","CD8+ T-cells","CD8+ Tcm","CD8+ Tem"),])
  DC<-colMeans(dat[c("aDC","cDC","DC","iDC","pDC"),])
  Endothelial<-colMeans(dat[c("Endothelial cells","ly Endothelial cells","mv Endothelial cells"),])
  Macrophage<-colMeans(dat[c("Macrophages","Macrophages M1","Macrophages M2"),])
  Bcell<-colMeans(dat[c("B-cells","Memory B-cells","naive B-cells","Class-switched memory B-cells","Plasma cells","pro B-cells"),])
  
  remove<-c("CD4+ memory T-cells","CD4+ naive T-cells",
            "CD4+ T-cells","CD4+ Tcm","CD4+ Tem",
            "CD8+ naive T-cells","CD8+ T-cells","CD8+ Tcm","CD8+ Tem",
            "aDC","cDC","DC","iDC","pDC",
            "Endothelial cells","ly Endothelial cells","mv Endothelial cells",
            "Macrophages","Macrophages M1","Macrophages M2",
            "B-cells","Memory B-cells","naive B-cells","Class-switched memory B-cells","pro B-cells",
            "Plasma cells","StromaScore","MicroenvironmentScore",
            "HSC","CLP","CMP","GMP","MEP","MPP",
            "Megakaryocytes","Erythrocytes","Platelets",
            "Smooth muscle","ImmuneScore")
  
  dat<-dat[which(!(rownames(dat)%in%remove)),]
  
  dat["CD4+ cells",]<-CD4
  dat["CD8+ cells",]<-CD8
  dat["Dendritic cells",]<-DC
  dat["Endothelial cells",]<-Endothelial
  dat["Macrophages",]<-Macrophage
  dat["B cells",]<-Bcell
  
  #Scale proportions so they sum to 1
  scaled_dat<-dat
  for (i in 1:ncol(dat)) {
    scaled_dat[,i]<-dat[,i]/sum(dat[,i])
  }
  return(scaled_dat)
}
###############################################

celltypes<-fread("~/midway/xCell/CMC_GTEx_xCell_tpm_results.txt",data.table=FALSE)
prep_columns(celltypes_tissue)

# Get sex information

if (tissue=="cortex") {
  
} else {
  sex<-fread("~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt")
  sex<-sex[,c('SUBJID',"SEX")]
  subjs<-paste("GTEX",sapply(strsplit(rownames(celltypes_tissue),'-'),'[[',2),sep='-')
  sex<-sex$SEX[match(subjs,sex$SUBJID)]
}

#Get tissue information
tissue_samples<-c("V1",scan("~/midway/expression/WHR/adipose_subcutaneous/adipose_subcutaneous_sample_IDs",what='character',sep='\n'))
celltypes_tissue<-celltypes[,colnames(celltypes)%in%tissue_samples]
rownames(celltypes_tissue)<-celltypes_tissue$V1
celltypes_tissue$V1<-NULL
celltypes_tissue<-prep_columns(celltypes_tissue)

props_F<-as.data.frame(t(rbind(rownames(celltypes_tissue),
                                  rowMeans(celltypes_tissue[,sex==2]),
                                  apply(celltypes_tissue[,sex==2],1,sd))),stringsAsFactors = FALSE)
colnames(props_F)<-c("celltype","mean_prop","sd")
props_F$mean_prop<-as.numeric(props_F$mean_prop)
props_F$sd<-as.numeric(props_F$sd)
props_F$celltype<-factor(props_F$celltype,levels=unique(props_F$celltype))
props_F$sex="F"

props_M<-as.data.frame(t(rbind(rownames(celltypes_tissue),
                               rowMeans(celltypes_tissue[,sex==1]),
                               apply(celltypes_tissue[,sex==1],1,sd))),stringsAsFactors = FALSE)
colnames(props_M)<-c("celltype","mean_prop","sd")
props_M$mean_prop<-as.numeric(props_M$mean_prop)
props_M$sd<-as.numeric(props_M$sd)
props_M$celltype<-factor(props_M$celltype,levels=unique(props_M$celltype))
props_M$sex="M"

props<-rbind(props_F,props_M)

C<-ggplot(props,aes(x=celltype,y=mean_prop))+
  geom_bar(aes(fill=sex),position="dodge",stat="identity")+
  geom_errorbar(aes(ymin=mean_prop-sd, ymax=mean_prop+sd), width=.2,
                position=position_dodge(.9))+
  theme_minimal()+
  scale_y_continuous("Estimated cell type proportion",limits=c(-0.05,0.75))+
  scale_fill_manual(values=c(rgb(pony_colors[8,1:3]),rgb(pony_colors[10,1:3])))+
  theme(
    plot.title = element_text(size=24),
    axis.title=element_text(size=15),
    legend.position = "none",
    axis.title.x=element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  ggtitle("Frontal cortex RNA-seq")
C
