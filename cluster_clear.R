source("http://bioconductor.org/biocLite.R")
biocLite('clusterProfiler')
biocLite('biocLite')

library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(carData)
library(org.Hs.eg.db)

# 数据读入 --------------------------------------------------------------------
setwd('D:\\E\\博士\\R_程序\\SVM_RFE\\Data\\limma')

coef_Bridge <- read.csv("gene_feature.csv", header=TRUE, sep = ',')
x1 <- as.matrix(coef_Bridge)
x1 <- as.character(x1)

View(x1)
eg1 <- bitr(x1, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db"); 
head(eg1)

genelist <- as.matrix(eg1[,2])
colnames(genelist) <- c("ENTREZID")

genelist[duplicated(genelist)]

go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP', pAdjustMethod = 'BH', pvalueCutoff = 0.2, 
                 qvalueCutoff = 0.2, keyType = 'ENTREZID')
# write.csv(go, file = "goBP6.csv")
head(go)

dim(go)   # 419  10
dim(go[go$ONTOLOGY=='BP',])    #298  10
dim(go[go$ONTOLOGY=='CC',])    #59  10
dim(go[go$ONTOLOGY=='MF',])    #62  10

# 进行简单的可视化
# pdf(file = "go6.pdf",width = 8,height = 6)
barplot(go)
dotplot(go)
# dev.off()

barplot(go,showCategory=20,drop=T)
dotplot(go,showCategory=50)
# dev.off()

go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
# pdf(file = "gonet.pdf")
plotGOgraph(go.BP)
# dev.off()

adjustedp <- read.csv("adjustedp.csv", header=TRUE, sep = ',')

gene1 <- as.matrix(coef_Bridge[,2])
colnames(gene1) <- c("gene1")

colnames(adjustedp) <- c("gene2","p")

x3 <- merge(gene1, adjustedp, by.x = "gene1", by.y="gene2",all=FALSE) 

# write.csv(x3, file = "coef_p.csv")
