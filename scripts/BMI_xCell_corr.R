library(tidyverse)
library(data.table)
library(gridExtra)

#Correlation between BMI and high adipocyte cell proportion in xCell (GTEx data)

#Get BMI from GTEx subject data, combine it with sample names
phenos<-fread("~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt")
bmicol<-which(colnames(phenos)=="BMI") #for QC later
phenos<-phenos[,c("SUBJID","BMI")]

subjs_samples<-fread("~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt",stringsAsFactors = FALSE)
subjs_samples<-as.data.frame(subjs_samples$SAMPID,stringsAsFactors=FALSE)
subjs_samples$subj<-paste("GTEX-",sapply(strsplit(subjs_samples$`subjs_samples$SAMPID`,'-'),'[[',2),sep='')
colnames(subjs_samples)<-c("SAMPLE","SUBJID")

bmis<-merge(subjs_samples,phenos,by="SUBJID")

#Load xCell data
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

groups<-c("Adipose","Neural","Immune","Stromal","Immune"
          ,rep("Stromal",4),"Immune"
          ,rep("Stromal",2),"Immune","Stromal","Stromal",
          "Neural",rep("Immune",3),rep("Stromal",2),
          "Adipose","Stromal","Stromal",rep("Immune",7),
          "Stromal","Immune","Immune")

## Adipose
celltypes<-fread("~/midway/xCell/CMC_GTEx_xCell_tpm_results.txt",data.table=FALSE)
adipose_samples<-c("V1",scan("~/midway/expression/WHR/adipose_subcutaneous/adipose_subcutaneous_sample_IDs",what='character',sep='\n'))
celltypes_adipose<-celltypes[,colnames(celltypes)%in%adipose_samples]
rownames(celltypes_adipose)<-celltypes_adipose$V1
celltypes_adipose$V1<-NULL
celltypes_adipose<-prep_columns(celltypes_adipose)
celltypes_adipose<-as.data.frame(t(celltypes_adipose),stringsAsFactors=FALSE)

#Merge xCell and bmi data
dat<-merge(bmis,celltypes_adipose,by="SAMPLE")

#Plot correlation for cell types from Science direct paper: 
# https://www.sciencedirect.com/science/article/pii/S0002929719301211
# Cell types from their paper: adipocytes, macrophages, CD4+ T cells, and micro-vascular endothelial cells
# Cell types that I included as covariates: Adipocytes	Chondrocytes	Epithelial_cells	Fibroblasts	Preadipocytes	Dendritic_cells
Adipo_cor=cor(dat$BMI,dat$Adipocytes,method="pearson")
Chondo_cor=cor(dat$BMI,dat$Chondrocytes,method="pearson")
Epi_cor=cor(dat$BMI,dat$`Epithelial cells`,method="pearson")
Fibro_cor=cor(dat$BMI,dat$Fibroblasts,method="pearson")
Preadipo_cor=cor(dat$BMI,dat$Predipocytes,method="pearson")
Dendro_cor=cor(dat$BMI,dat$`Dendritic cells`,method="pearson")

##Adipocytes
r=cor(dat$BMI,dat$Adipocytes,method="pearson")
A<-ggplot(data=dat,aes(BMI,Adipocytes))+
  geom_point(color="gray40")+
  ggtitle("Correlation between BMI and adipocyte cell proportion in \n GTEx subcutaneous adipose samples")+
  scale_y_continuous(name="Cell proportion")+
  geom_smooth(method="lm",se=FALSE,color="steelblue4",width=3)+
  annotate("text",label=paste("Pearson's r^2=",round(r^2,digits=6),sep=''),x=20,y=0.22,size=6)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        title=element_text(size=14))

#Macrophages
r=cor(dat$BMI,dat$Macrophages,method="pearson")
M<-ggplot(data=dat,aes(BMI,Macrophages))+
  geom_point(color="gray40")+
  ggtitle("Correlation between BMI and macrophage cell proportion in \n GTEx subcutaneous adipose samples")+
  scale_y_continuous(name="Cell proportion")+
  geom_smooth(method="lm",se=FALSE,color="steelblue4",width=3)+
  annotate("text",label=paste("Pearson's r^2=",round(r^2,digits=4),sep=''),x=25,y=0.1,size=6)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        title=element_text(size=14))

#CD4 T-cells
r=cor(dat$BMI,dat$`CD4+ cells`,method="pearson")
CD4T<-ggplot(data=dat,aes(BMI,`CD4+ cells`))+
  geom_point(color="gray40")+
  ggtitle("Correlation between BMI and CD4+ T-cell proportion in \n GTEx subcutaneous adipose samples")+
  scale_y_continuous(name="Cell proportion")+
  geom_smooth(method="lm",se=FALSE,color="steelblue4",width=3)+
  annotate("text",label=paste("Pearson's r^2=",round(r^2,digits=4),sep=''),x=25,y=0.018,size=6)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        title=element_text(size=14))

#Endothelial cells
r=cor(dat$BMI,dat$`Endothelial cells`,method="pearson")
E<-ggplot(data=dat,aes(BMI,`Endothelial cells`))+
  geom_point(color="gray40")+
  ggtitle("Correlation between BMI and endothelial cell proportion in \n GTEx subcutaneous adipose samples")+
  scale_y_continuous(name="Cell proportion")+
  geom_smooth(method="lm",se=FALSE,color="steelblue4",width=3)+
  annotate("text",label=paste("Pearson's r^2=",round(r^2,digits=4),sep=''),x=25,y=0.25,size=6)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        title=element_text(size=14))

jpeg("~/midway/xCell/xCell_BMI_corrs.jpg",width=800,height=800)
grid.arrange(A,M,CD4T,E,ncol=2,nrow=2)
dev.off()


#QC:
  #Check BMI:subject:sample
  t1<-strsplit(system(paste("awk '$1==\"'",dat$SUBJID[1],"'\" {print $0}' ~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt",sep=''),intern=TRUE),'\t')[[1]]
  t1[bmicol]==dat[1,c("BMI")]
  
  t1<-strsplit(system(paste("awk '$1==\"'",dat$SUBJID[5000],"'\" {print $0}' ~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt",sep=''),intern=TRUE),'\t')[[1]]
  t1[bmicol]==dat[5000,c("BMI")]
  
  t1<-strsplit(system(paste("awk '$1==\"'",dat$SUBJID[17000],"'\" {print $0}' ~/midway/GTEx_Analysis_2017-06-05_v8_Annotations%2FGTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt",sep=''),intern=TRUE),'\t')[[1]]
  t1[bmicol]==dat[17000,c("BMI")]
  
  #Check cell props: sample
  t1<-strsplit(system(paste("awk '$1==\"'",dat$SAMPLE[1],"'\" {print $0}' ~/medusa/papers/TWAS/supplement/supplementary_tables/Table2_xCellresults.txt",sep=''),intern=TRUE),'\t')[[1]]
  t1[which(colnames(xcell)=="Adipocytes")]==dat[1,c("Adipocytes")]
  
  t2<-strsplit(system(paste("awk '$1==\"'",dat$SAMPLE[5000],"'\" {print $0}' ~/medusa/papers/TWAS/supplement/supplementary_tables/Table2_xCellresults.txt",sep=''),intern=TRUE),'\t')[[1]]
  t2[which(colnames(xcell)=="Adipocytes")]==dat[5000,c("Adipocytes")]
  
  t3<-strsplit(system(paste("awk '$1==\"'",dat$SAMPLE[17000],"'\" {print $0}' ~/medusa/papers/TWAS/supplement/supplementary_tables/Table2_xCellresults.txt",sep=''),intern=TRUE),'\t')[[1]]
  t3[which(colnames(xcell)=="Adipocytes")]==dat[17000,c("Adipocytes")]
  