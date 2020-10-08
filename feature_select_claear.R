library(dplyr)        
library(tidyr)
library(tidyverse)   

setwd('D:\\E\\博士\\R_程序\\SVM_RFE\\Data\\limma')

Data1 = read.table("GSE73685_scale_blood_zf17_new.txt", header = T, check.names = FALSE)
gene = read.csv("gene_feature.csv", header=TRUE, sep = ',')
dim(Data1)

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
# write.table(genedata2,"GSE73685_46.txt",quote=F,sep="\t") 

View(genedata[,1])


# 在原数据集提取特征 ---------------------------------------------------------------

data1 = read.table("GSE59491_outcome_scale.txt", header = T, check.names = FALSE)
dim(data1)
gene <- as.matrix(genedata[,1])

colnames(gene) <- c('gene')
data2 <- cbind(rownames(data1), data1)
colnames(data2) <- c('gene', colnames(data1))

genedata <- merge(gene, data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
# write.table(genedata2,"GSE59491_46.txt",quote=F,sep="\t") 