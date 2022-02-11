
## clear
rm(list = ls())

## package
library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       
library(tidyr)
library(tidyverse)    
library(fdrtool)      


## load data
setwd("D:\\E\\博士\\R_程序\\GSE59491_15")
expSet <- read.table("Data/GSE59491_series_matrix.txt", header=T, sep='\t', fill=TRUE, strip.white = T)
expSet$ID_REF <- as.character(expSet$ID_REF)  # 将ID_REF列全部转换成符号,为了同anno2合并


## load annotation
anno <- read.table("Data/GPL18964.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  # 探针和基因ID的对应文件
anno2 <- anno[,c('ID','ENTREZ_GENE_ID')]     # 提取这两列标签
colnames(anno2) <- c('ID_REF','EntrezID')    # 将这两列的标签替换, 'GeneSymbol''Gene.Symbol'
anno2$ID_REF <- as.character(anno2$ID_REF)   # 将ID_REF列全部转换成符号,为了同expSet合并


## 将基因表达数据与芯片注释文件的探针名进行对应
expset2 <- expSet %>%                      # ％>％来自dplyr包的管道函数，
  inner_join(anno2,by='ID_REF') %>%        # 作用是将前一步的结果直接传参给下一步的函数，省略中间的赋值步骤
  select(ID_REF,EntrezID, everything())    # %>%# 重新排列
# View(expset2[,1:3])                      # expset2 与 expset 相比，第2列多了基因号


## 整理芯片注释文件，把其中一个探针对应多个基因的拆分开
expset3 <- expset2
a <- tibble(expset3[,1:2])                 # 把第 1 和第 2 列提取出来，放在 a 中  
test1 <- apply(a,1, function(x){
  str_split(x[2], '///', simplify=T)       # 把 a 的第2列提取出，给test1
} )


test2 <- apply(a, 1, function(x){          # 将探针号和基因号，进行---的链接
  paste(x[1],str_split(x[2], '///', simplify=T), sep = "---")
})


unlist(test2)                              # 将 list 数据变成字符串向量或者数字向量的形式
x <- tibble(unlist(test2))                 # tibble，取代传统data.frame，读取并自动添加列名：unlist(test2)
colnames(x) <- "lala"                      # 改变 x 的列名：将 unlist(test2) 定义为 lala
x2 <- separate(x,lala,c("id","entrezID"),sep = '---')     # 识别 lala 中的 ---，将数据分离，单独成列并附新标签
x3 <- merge(x2,expset3,by.x = "id", by.y="ID_REF", all=FALSE)  #  将两个文件按顺序合并为一个，
x4<-x3[,-c(1,3)]                           # 将 第1 和第3 两列删除, 剩下的数据还是字符型的，带着“ "
zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))
XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))


## 用基因id对整理好的芯片注释文件进行基因名的更新
homo<-read.table("Data/homo.txt",header=T,sep='\t')
x5 <- merge(homo, XX1, by.x="GeneID", by.y = "entrezID", all=FALSE) 
# 合并， x5 从25088，降为24478，基因号与基因名进行匹配


## 探针名匹配基因名，取出多个探针对应一个基因的数据计算IQR，保留IQR最大的探针数据
expset4 <- x5 %>% 
  dplyr::select(-GeneID) %>%              # 去掉多余信息
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>% # 计算每行的IQR
  arrange(desc(rowIQR)) %>%               # 把表达量的平均值按从大到小排序
  distinct(Symbol,.keep_all = T) %>%      # symbol留下第一个
  dplyr::select(-rowIQR)   %>%            # 反向选择去除rowIQR这一列
  tibble::column_to_rownames(colnames(.)[1]) # 把第一列变成行名并删除
View(expset4[1:10,1:10])
dim(expset4)    # 24478   326
  

# 标签 ----------------------------------------------------------------------
lable2 = read.csv("Data/GSE59491_all.csv", header = T, sep=',')
dim(lable2)    # 326   2
data = rbind(as.matrix(t(lable2[,2])), as.matrix(expset4))
rownames(data) <- c('Lable', rownames(expset4))
View(data[1:10,1:10])
dim(data)    #  24479   326
# write.table(data,"Data/GSE59491_outcome.txt",quote=F,sep="\t")


