#!/usr/bin/Rscript

######### Author: Grace Hansen #########
#This script runs xCell to estimate cell type proportions on GTEx and CMC data.

args=commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(xCell)
library(data.table)

# Load expression dtaa
dat<-fread(args[1],data.table=FALSE)
rownames(dat)<-make.names(dat[,1],unique=TRUE)
dat<-dat[,-1]

#Prepare naming conventions
working_dir<-dirname(args[1])
prefix<-strsplit(basename(args[1]),'[.]')[[1]][1]

print(working_dir)
print(prefix)
if (len(args)>1) {
  celltypes<-strsplit(args[2],',')[[1]]
} else {
  celltypes<-c("Adipocytes","Epithelial cells","Hepatocytes","Keratinocytes","Myocytes","Neurons","Neutrophils")
}

#Run xCell
out<-xCellAnalysis(dat)
out<-out[which(rownames(out)%in%celltypes),]
write.table(out,paste(working_dir,"/",prefix,"_results_7celltypes.txt",sep=''),quote=FALSE,sep='\t')
