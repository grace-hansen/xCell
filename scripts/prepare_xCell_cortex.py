import os, gzip
import pandas as pd


############## Author: Grace Hansen #################
#This script combines GTex data with obesity cortex data so that xCell can be run on both together.
###################################################################
CMC=pd.read_csv("~/midway/expression/obesity/cortex/cortex_CMC_tpm.txt.gz",compression="gzip",sep='\t')
CMC=CMC.drop(['coord'], axis=1)
CMC=CMC.drop(['chr'], axis=1)

GTEx=pd.read_csv("/project2/yangili1/GTEx_v8/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",sep='\t',compression="gzip",skiprows=2)
#Convert gene ID from ENSG to gene symbol
GTEx=GTEx.drop(['Name'], axis=1)
GTEx=GTEx.rename(columns={"Description":"gene"})

#Merge GTEx and cCMC data
out=pd.merge(CMC,GTEx,on="gene")
out.to_csv("~/midway/expression/obesity/cortex/CMC_GTEx_xCell_tpm.txt",sep='\t',index=False)