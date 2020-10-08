library(dplyr)     
library(tidyr)
library(tidyverse)  

setwd('D:\\E\\博士\\R_程序\\SVM_RFE\\Data\\limma')



# 提取 feature gene 表达值 -----------------------------------------------------------------

Data1 = read.table("GSE59491_outcome_scale.txt", header = T, check.names = FALSE)
gene = read.csv("gene_feature.csv", header=TRUE, sep = ',')
dim(Data1)

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
dim(genedata2)    #  55 326
# write.table(genedata2,"GSE59491_feature.txt",quote=F,sep="\t") 

