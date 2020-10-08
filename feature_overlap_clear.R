setwd('D:\\E\\博士\\R_程序\\SVM_RFE\\Data\\limma')


# 读入数据 --------------------------------------------------------------------

svm <- read.csv("gene_ref_svm.csv", header=TRUE, sep = ',')
knn <- read.csv("gene_ref_knn.csv", header=TRUE, sep = ',')
rf <- read.csv("gene_ref_rf.csv", header=TRUE, sep = ',')
ab <- read.csv("gene_ref_ab.csv", header=TRUE, sep = ',')
nn <- read.csv("gene_ref_nnet.csv", header=TRUE, sep = ',')

svm <- as.matrix(svm)
knn <- as.matrix(knn)
rf <- as.matrix(rf)
ab <- as.matrix(ab)
nn <- as.matrix(nn)
  
gene <- intersect(svm, intersect(knn, intersect(rf, ab)))


svm_1 <- intersect(svm, knn)
svm_2 <- intersect(svm, rf)
svm_3 <- intersect(svm, ab)
svm_4 <- intersect(svm, nn)

knn_1 <- intersect(knn, rf)
knn_2 <- intersect(knn, ab)
knn_3 <- intersect(knn, nn)

rf_1 <- intersect(rf, ab)
rf_2 <- intersect(rf, nn)

ab_1 <- intersect(ab, nn)

gene <- union(svm_1, union(svm_2, union(svm_3, union(svm_4, union(knn_1, union(knn_2, union(knn_3, union(rf_1, union(rf_2, ab_1)))))))))


# write.csv(gene, file = "gene_feature.csv", row.names = F)
