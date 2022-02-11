# 2020.6.29


## clear
rm(list = ls())


## package
library(dplyr)  
library(tidyr)
library(tidyverse)  


## load data
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data")
Data1 = read.table("GSE73685_scale_blood_zf17_zhafter.txt", header = T, check.names = FALSE)
dim(Data1)    # 20910    17

gene = read.csv("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37\\result\\gene_overlap_scale20.csv", header=TRUE, sep = ',') [,-1]
gene <- as.matrix(gene)
colnames(gene) <- c('gene')


##  extract data from independent dataset
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))
genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
# write.table(genedata2,"GSE73685_17.txt",quote=F,sep="\t")  # 2020.6.29


##  extract data from origina discovery dataset
data1 = read.table("GSE59491_scale_DE.txt", header = T, check.names = FALSE)
dim(data1)    # 360 326
gene <- as.matrix(genedata[,1])
colnames(gene) <- c('gene')
data2 <- cbind(rownames(data1), data1)
colnames(data2) <- c('gene', colnames(data1))

genedata <- merge(gene, data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
# write.table(genedata2,"GSE59491_17.txt",quote=F,sep="\t") 