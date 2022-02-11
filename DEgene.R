## ���ڶ� series_marix ���ӱ�ǩ������ݣ����� T ����
## ȷ�� T ����� p ֵ������ p/pdf <0.05 ��gene������


## clear
rm(list = ls())

## package
library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       # ��>�� �ܵ������ĵ��ã�����
library(tidyr)
library(tidyverse)   # tibble �ĵ���
library(fdrtool)     # fdrУ��


## load function
source("C:\\Users\\LiLingyu\\Desktop\\LogReg\\R\\Ttest.R")


## load data
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15")
x1 = read.table("Data/GSE59491_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)
GSE59491 <- t(x1)
dim(x1)    # 24479   326


############################  GSE59491 ���������t ���飩 #############################
# p < 0.05 ��   2,809 
test1 <- as.matrix(my_test_1(GSE59491))
test0 <- as.matrix(my_test_0(GSE59491))
p59491 <- my_p(GSE59491)                    
# write.csv(p59491,"p59491.csv")

p59491_BH_fdr <- my_BH_fdr(p59491)
# write.csv(p59491_BH_fdr,"p59491_BH_fdr.csv")

# p59491_fdr <- my_fdr(p59491,p59491)
# View(p59491_fdr)
# write.csv(p59491_fdr,"p59491_fdr.csv")

GSE59491_T <- my_p_data(p59491_BH_fdr, GSE59491)
# GSE59491_T <- my_p_data(p59491, GSE59491)   # GSE73685
dim(GSE59491_T)    # 326 360
# write.table(GSE59491_T,"GSE59491_DE.txt",quote=F,sep="\t") 

## ��׼��,scale������ʱdata.fram
data59491 <- my_scale(data.frame(t(GSE59491_T)))
dim(data59491)    # 326 360
View(data59491[,1:10])
# write.table(data59491,"GSE59491_scale_DE.txt",quote=F,sep="\t") 




############################  GSE46510 ���������t ���飩 #############################

## ��������
x1 = read.table("GSE46510_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)

## ��׼��,ȥNA
x2 <- my_scale(x1)
# dim(x2) # 18419   154

x2[!complete.cases(x2),] 
x3 <- na.omit(x2) 
# dim(x3) # 17171   154

# ��ȡgeneԭ����ֵ ------------------------------------------------------------------
x4 <- as.matrix(rownames(x3)[-1])
colnames(x4) <- c('gene')

x5 <- cbind(rownames(x2)[-1], x1[-1,])
colnames(x5) <- c('gene', colnames(x1))
x6 <- merge(x4, x5, by = "gene") %>% tibble::column_to_rownames(colnames(.)[1])
# dim(x6)
x7 <- rbind(t(x3[1,]),x6)
dim(x7)
rownames(x7) <- c('Lable', rownames(x6))
# write.table(x7, "GSE46510_NA.txt", quote=F, sep="\t") 

# Test --------------------------------------------------------------------

x1 = read.table("GSE46510_NA.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)
# dim(x1)
GSE59491 <- t(x1)

## p < 0.05 ��   2,809 
test1 <- as.matrix(my_test_1(GSE59491))
test0 <- as.matrix(my_test_0(GSE59491))
p59491 <- my_p(GSE59491)                    
# dim(p59491)
# write.csv(p59491,"p59491.csv")

p59491_BH_fdr <- my_BH_fdr(p59491)
# View(p59491_BH_fdr)
# write.csv(p59491_BH_fdr,"p59491_BH_fdr.csv")

# p59491_fdr <- my_fdr(p59491,p59491)
# View(p59491_fdr)
# write.csv(p59491_fdr,"p59491_fdr.csv")

GSE59491_T <- my_p_data(p59491_BH_fdr, GSE59491)
# GSE59491_T <- my_p_data(p59491, GSE59491)
dim(GSE59491_T)    # 154  62
# write.table(GSE59491_T,"GSE59491_DE.txt",quote=F,sep="\t") 

## ��׼��,scale������ʱdata.fram
data59491 <- my_scale(data.frame(t(GSE59491_T)))
dim(data59491)    # 62 154
# write.table(data59491,"GSE46510_scale_pDE.txt",quote=F,sep="\t") 