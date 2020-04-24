library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       # ��>�� �ܵ������ĵ��ã�����
library(tidyr)
library(tidyverse)   # tibble �ĵ���
library(fdrtool)     # fdrУ��
library(data.table)  # �� fread ����.soft �ļ�
library(fdrtool)     # ���� p ֵ
library(data.table)  # ʹ�� fread ��������
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data")

## ��������������
exprSet <- read.table("GSE73685_series_matrix.txt",header=T,sep='\t',fill=TRUE,strip.white = T)
exprSet$ID_REF <- as.character(exprSet$ID_REF)

## ����оƬ̽��ע���ļ�
anno <- read.table("GPL6244.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  # ̽��ͻ���ID�Ķ�Ӧ�ļ�

library(dplyr)
library(tidyr)
anno2 <- anno %>%                   # GSE73685_family.soft �ļ������� 
  select(ID,GeneID) %>%             # ��ȡ�ļ��� ID,Gene_ID ���� 
  filter(GeneID != '') 
colnames(anno2) <- c('ID_REF','EntrezID')    # �������еı�ǩ�滻, 'GeneSymbol''Gene.Symbol'
anno2$ID_REF <- as.character(anno2$ID_REF)   # ��ID_REF��ȫ��ת���ɷ���,Ϊ��ͬexpSet�ϲ�

## ���������������оƬע���ļ���̽�������ж�Ӧ
exprset2 <- exprSet %>%                      # ��>������dplyr���Ĺܵ�������
  inner_join(anno2,by='ID_REF') %>%          # �����ǽ�ǰһ���Ľ��ֱ�Ӵ��θ���һ���ĺ�����ʡ���м�ĸ�ֵ����
  select(ID_REF,EntrezID, everything())      # %>%# ��������
# View(exprset2[,1:2])                       # expset2 �� expset ��ȣ���2�ж��˻����


## ����оƬע���ļ���������һ��̽���Ӧ�������Ĳ�ֿ�
exprset3 <- exprset2
a <- tibble(exprset3[,1:2])

test1 <- apply(a,1, function(x){
  str_split(x[2],'///',simplify=T)       # �� a �ĵ�2����ȡ������test1
} )

test2 <- apply(a, 1, function(x){          # ��̽��źͻ���ţ�����---������
  paste(x[1],str_split(x[2],'///', simplify=T), sep = "---")
})

unlist(test2)                              # �� list ���ݱ���ַ�����������������������ʽ
x <- tibble(unlist(test2))                 # tibble��ȡ����ͳdata.frame����ȡ���Զ�����������unlist(test2)
colnames(x) <- "lala"                      # �ı� x ���������� unlist(test2) ����Ϊ lala

x2 <- separate(x,lala,c("id","entrezID"),sep = '---')     # ʶ�� lala �е� ---�������ݷ��룬�������в����±�ǩ
x2[1:10,1]
expset3[1:10,1]
x3 <- merge(x2,exprset3,by.x = "id",by.y="ID_REF",all=FALSE)  #  �������ļ���˳��ϲ�Ϊһ����

x4<-x3[,-c(1,3)]                           # �� ��1 �͵�3 ����ɾ��, ʣ�µ����ݻ����ַ��͵ģ����š� "
View(x4[,1:3])
x4[1:6,1]

zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))
View(zz)

XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))

## �û���id�������õ�оƬע���ļ����л������ĸ���
homo<-read.table("homo.txt",header=T,sep='\t')
homo[1:6,1]
x5 <- merge(homo,XX1,by.x="GeneID",by.y = "entrezID",all=FALSE) 
# �ϲ��� x5 ��25088����Ϊ24478������������������ƥ��
# View(x5[,1:10])
# dim(x5)    # 28468   185

## ̽����ƥ���������ȡ�����̽���Ӧһ����������ݼ���IQR������IQR����̽������
expset4 <- x5 %>%
  dplyr::select(-GeneID) %>%              # ȥ��������Ϣ
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>% # ����ÿ�е�IQR
  arrange(desc(rowIQR)) %>%               # �ѱ�������ƽ��ֵ���Ӵ�С����
  distinct(Symbol,.keep_all = T) %>%      # symbol���µ�һ��
  dplyr::select(-rowIQR) %>%                 # ����ѡ��ȥ��rowIQR��һ��
  tibble::column_to_rownames(colnames(.)[1]) # �ѵ�һ�б��������ɾ��
View(expset4[1:10.1:10])
dim(expset4)   # 20909   183
# write.table(expset4,"GSE73685_expr.txt",quote=F,sep="\t")  

# ���� 0/1 ��ǩ ---------------------------------------------------------------
# �˴���Ҫ�� tibble ע��

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

## 1��gene
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


# �����仯��scaleһ�� ------------------------------------------------------------

setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data")

## ��������
data = read.table("GSE73685_outcome_blood17.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)

# my_scale -------------------------------------------------------------------

my_scale <- function(x){
  x1 <- cbind(t(x[1,]), scale(t(x[-1,])))
  x2 <- t(x1)
  return(x2)
}

# ���� ----------------------------------------------------------------------

data1 <- my_scale(data)
dim(data1)
write.table(data1,"GSE73685_scale_blood17.txt",quote=F,sep="\t")   