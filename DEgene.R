library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用

setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37\\DE")
# load(".RData")

Data1 = read.table("D:\\E\\博士\\R_程序\\GSE59491_15\\Data\\GSE73685_outcome_blood17.txt", header = T, check.names = FALSE)
# Data1 = read.table("D:\\E\\博士\\R_程序\\GSE59491_15\\Data\\GSE59491_scale_DE.txt", header = T, check.names = FALSE) #scale以后的 gene*sample
gene = read.csv("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37\\result\\gene_scale.csv", header=TRUE, sep = ',')
# View(gene)

colnames(gene) <- c('no','gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
# write.table(genedata2,"genedata.txt",quote=F,sep="\t") 

# 找DEgene -----------------------------------------------------------------

# 函数，区分0-1，标签在第1列 ---------------------------------------------------------

## 全为 1 的
my_test_1 <- function(x){
  arry <- x[,1]
  test1 <- x[arry == 1, ]
  return(test1)
}
## 全为 0 的
my_test_0 <- function(x){
  arry <- x[,1]
  test0 <- x[arry == 0, ]
  return(test0)
}

# my_p T检验，函数，标签在第1列 -----------------------------------------------------------

my_p <- function(x){
  p <- matrix(data=NA, nrow = dim(x)[2]-1, ncol = 1, byrow = FALSE, dimnames=list(c(colnames(x[,-1])),c("pvalue")))
  for(i in 2:ncol(x)){
    p[i-1] <- t.test(as.numeric(test0[,i]), as.numeric(test1[, i]))$p.value  # 未scale的数据
  }
  return(p)
}

# BH校正提取 ------------------------------------------------------------------

my_data <- function(x){
  x_hat <- t(x)[-1,]  
  x_hat_Right <- x_hat[Right, ]
  dim(x_hat_Right)    
  y <- x[,1]
  y <- as.matrix(y)
  colnames(y) <- 'Lable'
  x_lab <- cbind(y, t(x_hat_Right))
  return(x_lab)
}  

my_BH_fdr <- function(x){
  BH_fdr <- as.matrix(p.adjust(x, "BH"))
  p_BH_fdr <- matrix(data=BH_fdr, nrow = dim(x)[1], ncol = 1, byrow = FALSE, dimnames=list(c(as.matrix(rownames(x))),c("p_BH_fdr")))
  return(p_BH_fdr)
}

# 检验 -----------------------------------------------------------------------

# gene  <- read.table("D:\\E\\博士\\R_程序\\GSE59491_14\\Data\\genedata.txt", header = T, check.names = FALSE)
genedata3 <- t(genedata2)   
test0 <- my_test_0(genedata3)
test1 <- my_test_1(genedata3)
View(test1)
dim(genedata3)

p <- my_p(genedata3)
p_fdr <- my_BH_fdr(p)
Right <- which(p_fdr <= 1)

p59491 <- my_data(genedata3)
dim(p59491)  
write.table(t(p59491),"p73685_17_orig.txt",quote=F,sep="\t") 
# write.table(t(p59491),"p59491_20_orig.txt",quote=F,sep="\t") 

