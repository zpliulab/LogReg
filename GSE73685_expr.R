library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用
library(fdrtool)     # fdr校正
library(data.table)  # 用 fread 读入.soft 文件
library(fdrtool)     # 处理 p 值
library(data.table)  # 使用 fread 读入数据
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data")

## 读入基因表达数据
exprSet <- read.table("GSE73685_series_matrix.txt",header=T,sep='\t',fill=TRUE,strip.white = T)
exprSet$ID_REF <- as.character(exprSet$ID_REF)

## 读入芯片探针注释文件
anno <- read.table("GPL6244.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  # 探针和基因ID的对应文件

library(dplyr)
library(tidyr)
anno2 <- anno %>%                   # GSE73685_family.soft 文件的内容 
  select(ID,GeneID) %>%             # 提取文件中 ID,Gene_ID 两列 
  filter(GeneID != '') 
colnames(anno2) <- c('ID_REF','EntrezID')    # 将这两列的标签替换, 'GeneSymbol''Gene.Symbol'
anno2$ID_REF <- as.character(anno2$ID_REF)   # 将ID_REF列全部转换成符号,为了同expSet合并

## 将基因表达数据与芯片注释文件的探针名进行对应
exprset2 <- exprSet %>%                      # ％>％来自dplyr包的管道函数，
  inner_join(anno2,by='ID_REF') %>%          # 作用是将前一步的结果直接传参给下一步的函数，省略中间的赋值步骤
  select(ID_REF,EntrezID, everything())      # %>%# 重新排列
# View(exprset2[,1:2])                       # expset2 与 expset 相比，第2列多了基因号


## 整理芯片注释文件，把其中一个探针对应多个基因的拆分开
exprset3 <- exprset2
a <- tibble(exprset3[,1:2])

test1 <- apply(a,1, function(x){
  str_split(x[2],'///',simplify=T)       # 把 a 的第2列提取出，给test1
} )

test2 <- apply(a, 1, function(x){          # 将探针号和基因号，进行---的链接
  paste(x[1],str_split(x[2],'///', simplify=T), sep = "---")
})

unlist(test2)                              # 将 list 数据变成字符串向量或者数字向量的形式
x <- tibble(unlist(test2))                 # tibble，取代传统data.frame，读取并自动添加列名：unlist(test2)
colnames(x) <- "lala"                      # 改变 x 的列名：将 unlist(test2) 定义为 lala

x2 <- separate(x,lala,c("id","entrezID"),sep = '---')     # 识别 lala 中的 ---，将数据分离，单独成列并附新标签
x2[1:10,1]
expset3[1:10,1]
x3 <- merge(x2,exprset3,by.x = "id",by.y="ID_REF",all=FALSE)  #  将两个文件按顺序合并为一个，

x4<-x3[,-c(1,3)]                           # 将 第1 和第3 两列删除, 剩下的数据还是字符型的，带着“ "
View(x4[,1:3])
x4[1:6,1]

zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))
View(zz)

XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))

## 用基因id对整理好的芯片注释文件进行基因名的更新
homo<-read.table("homo.txt",header=T,sep='\t')
homo[1:6,1]
x5 <- merge(homo,XX1,by.x="GeneID",by.y = "entrezID",all=FALSE) 
# 合并， x5 从25088，降为24478，基因号与基因名进行匹配
# View(x5[,1:10])
# dim(x5)    # 28468   185

## 探针名匹配基因名，取出多个探针对应一个基因的数据计算IQR，保留IQR最大的探针数据
expset4 <- x5 %>%
  dplyr::select(-GeneID) %>%              # 去掉多余信息
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>% # 计算每行的IQR
  arrange(desc(rowIQR)) %>%               # 把表达量的平均值按从大到小排序
  distinct(Symbol,.keep_all = T) %>%      # symbol留下第一个
  dplyr::select(-rowIQR) %>%                 # 反向选择去除rowIQR这一列
  tibble::column_to_rownames(colnames(.)[1]) # 把第一列变成行名并删除
View(expset4[1:10.1:10])
dim(expset4)   # 20909   183
# write.table(expset4,"GSE73685_expr.txt",quote=F,sep="\t")  

# 添加 0/1 标签 ---------------------------------------------------------------
# 此处需要将 tibble 注销

lable = read.csv("GSE73685_all.csv", header = T, sep=',')
# View(lable2)
expset5 <- expset4[,which(lable$Source.name == "Maternal Blood")]
dim(expset5)          # 20909    24
# View(expset5[1:10,])


lable2 = read.csv("GSE73685_blood.csv", header = T, sep=',')
dim(lable2)    # 154   2
data = rbind(as.matrix(t(lable2[,2])), as.matrix(expset5))
rownames(data) <- c('Lable', rownames(expset5))
View(data[1:10,1:10])
dim(data)    # 20910    24
# write.table(data,"GSE46510_outcome_blood.txt",quote=F,sep="\t")   

## 1个gene
lable = read.csv("GSE73685_all.csv", header = T, sep=',')
lable_1 <- lable[which(lable$Source.name == "Maternal Blood"),]
lable_2 <- lable_1[-which(lable_1$Outcome == "PNL preterm no labor"),]
lable_3 <- lable_2[-which(lable_2$Outcome == "pPROM with labor"),]
lable_3[,1]
# View(lable_3)
expset5 <- expset4[,as.character(lable_3$Accession)]
dim(expset5)          # 20909    24
# View(expset5[1:10,])

lable2 = read.csv("GSE73685_blood17.csv", header = T, sep=',')
dim(lable2)    # 154   2
data = rbind(as.matrix(t(lable2[,2])), as.matrix(expset5))
rownames(data) <- c('Lable', rownames(expset5))
View(data[1:10,1:10])
dim(data)    # 20910    24
# write.table(data,"GSE73685_outcome_blood17.txt",quote=F,sep="\t")   


# 样本变化，scale一下 ------------------------------------------------------------

setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data")

## 读入数据
data = read.table("GSE73685_outcome_blood17.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)

# my_scale -------------------------------------------------------------------

my_scale <- function(x){
  x1 <- cbind(t(x[1,]), scale(t(x[-1,])))
  x2 <- t(x1)
  return(x2)
}

# 重组 ----------------------------------------------------------------------

data1 <- my_scale(data)
dim(data1)
write.table(data1,"GSE73685_scale_blood17.txt",quote=F,sep="\t")   
