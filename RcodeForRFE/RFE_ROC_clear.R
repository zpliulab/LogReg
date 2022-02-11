setwd('D:\\E\\博士\\R_程序\\SVM_RFE\\Data\\limma')
library(e1071)
library(limma)
library(ggplot2)
library(reshape2)
library(gmodels)


# 读入数据 --------------------------------------------------------------------

eset <- read.table("matrix.txt",header = TRUE,sep = "\t") 
eset <- read.table("GSE59491_feature.txt",header = TRUE,sep = "\t") 
eset <- as.matrix(eset)

# 读入排序文件 ------------------------------------------------------------------

rfe_svm <- read.table("ranklist_rf.txt",stringsAsFactors = FALSE)
colnames(rfe_svm) <- c("rank","chipID")
rfe_svm_dif<- head(rfe_svm,n = 40L)
eset_svm <- eset[rfe_svm_dif$chipID,]

# 差异表达gene ----------------------------------------------------------------

library(limma)

data1 <- eset
disease = read.table("phe_zh.txt",header = TRUE, sep = "\t")
disease <- factor(disease[,"class"])

design <- model.matrix(~-1+disease)
contrast.matrix <- makeContrasts(contrasts = "diseaseSPTB - diseaseTerm", levels = design)
fit_feature <- lmFit(data1,design)
fit1_feature <- contrasts.fit(fit_feature,contrast.matrix)
fit2_feature <- eBayes(fit1_feature)
dif <- topTable(fit2_feature,coef = "diseaseSPTB - diseaseTerm",n = nrow(fit2_feature))
dif_feature <- dif[dif[,"adj.P.Val"]<0.05,]  
# write.csv(dif_feature, file = "feature_svm.csv", row.names = F)

matrix_DE <- eset[rownames(dif_feature),]
# write.table(matrix_DE, file = "matrix_DE.txt", quote=F, sep="\t")

# 训练---测试 -----------------------------------------------------------------

x <- t(eset[rownames(dif_feature),])
gene_ref <- colnames(x)
# write.csv(gene_ref, file = "gene_ref_svm.csv", row.names = F)

y1 <- as.matrix(c(rep(1,98), rep(0,228)))
# colnames(y1) <- c("Label")
# x_hat <- t(x)
# x_hat1 <- rbind(t(y1),x_hat)
# write.table(x_hat1, file = "x_hat32.txt", quote=F, sep="\t")

disease <- as.factor(c(rep("SPTB",98), rep("Term",228)))
y <- disease


# Acc ---------------------------------------------------------------------

svm_acc <- list()
cost <- c(rep(10,4),rep(100,4))
gamma <- rep(c(0.1,0.01,0.001,0.0001),2)
parameters <- data.frame(cost,gamma)
for(k in 1:dim(parameters)[1]){
  svm_acc_tmp <- c()
  costtmp <- parameters[k,"cost"]
  gammatmp <- parameters[k,"gamma"]
  all_test_label <- c()
  all_pred_label <- c()
  
  set.seed(666) # for reproducing results
  rowIndices <- 1 : nrow(x) # prepare row indices
  sampleSize <- 0.70 * length(rowIndices) # training sample size
  trainingRows <- sample (rowIndices, sampleSize) # random sampling
  trainingData <- x[trainingRows, ] # training data
  testData <- x[-trainingRows, ] # test data
  trainingLabel <- y[trainingRows]
  testLabel <- y[-trainingRows]
  # trainingLabel <- y1[trainingRows]
  # testLabel <- y1[-trainingRows]
  # for(i in 1:dim(eset)[2]){
  #   ind <- rep(1,dim(eset)[2])
  #   ind[i] <- 2
    # trainingData <-x[ind==1,]
    # testData <-t(as.matrix(x[ind==2,]))
    # trainingLabel <- disease[ind==1]
    # testLabel <- disease[ind==2]
    svmfit <- svm (trainingData,trainingLabel, kernel = "radial", cost = costtmp, gamma=gammatmp, scale = FALSE) # radial svm, scaling turned OFF
    eset_pred<- predict(svmfit, testData)
    all_test_label <- c(all_test_label,as.vector(testLabel))
    all_pred_label <- c(all_pred_label,as.vector(eset_pred))
  # }
  svm_acc <- c(svm_acc,mean(all_test_label == all_pred_label))
}

library(gmodels) 
CrossTable(x=all_test_label,y=all_pred_label, prop.chisq=FALSE)

parametertypes <- c()
for(k in 1:dim(parameters)[1]){
  costtmp <- parameters[k,"cost"]
  costtmp <- paste("cost:",costtmp,sep = "")
  gammatmp <- parameters[k,"gamma"]
  gammatmp <- paste("gamma:",gammatmp,sep = "")
  parametertmp <- paste(costtmp,gammatmp)
  parametertypes <- c(parametertypes,parametertmp)
}
# View(parametertypes)
names(svm_acc) <- parametertypes
svm_acc<- data.frame(svm_acc)
# View(svm_acc)
library(ggplot2)
library(reshape2)
svm_melt<- melt(svm_acc)
colnames(svm_melt) <- c("Parameter","Accuracy")
svm_melt$Accuracy <- round(svm_melt$Accuracy,3)
# pdf(file = "KNN-RFE+SVM.pdf",width = 7,height = 5)
ggplot(data = svm_melt,aes(x = Parameter,y = Accuracy,fill = Parameter))+
  geom_bar(stat = 'identity', width = 0.6)+
  geom_text(aes(label = Accuracy),vjust=-0.5)+
  labs(title = "Accuracy of SVM in different parameter") +
  theme(axis.text.x = element_text(angle=30,size=10))
# dev.off()


# 画图 ----------------------------------------------------------------------

trainingLabel <- y1[trainingRows]
testLabel <- y1[-trainingRows]

svmfit <- svm (trainingData,trainingLabel, kernel = "radial", cost = 100, gamma=1e-3, scale = FALSE) # radial svm, scaling turned OFF
print(svmfit)
pred <- predict(svmfit, testData, decision.values = TRUE)
compareTable <- table(testLabel, pred) 
compareTable

library(caret)
confusionMatrix(compareTable)


library(pROC)
# pdf(file = "ROC_SVM_knn.pdf",width = 5,height = 5)
plot.roc(testLabel, pred, print.auc=T, main="pAUC")
# Roc_svm <- smooth(roc, method="binormal")
# legend("bottomright", legend=c("Accuracy=0.942", "Precision=0.750 ", "Sensitivity=0.360", "Specificity=0.990", "F-measure=0.486"))
# dev.off()



# 性能指标 --------------------------------------------------------------------

predict = ifelse(pred > 0.5, 1, 0)
predict_value = predict 
true_value = testLabel
error = predict_value-true_value
table(true_value, predict_value) 

data <- x
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure= 2*precision*recall/(precision+recall)    
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 


# 保存结果 --------------------------------------------------------------------
A_roc <- cbind(testLabel, pred)
result_svm <- c(accuracy, precision, recall, specificity, F_measure)
# write.csv(A_roc, file = "A_roc_rf.csv", row.names = F)
# write.csv(result_svm, file = "result_rf.csv", row.names = F)
# write.csv(dif_feature, file = "feature_rf.csv", row.names = F)
# write.csv(gene_ref, file = "gene_ref_rf.csv", row.names = F)


# library(pROC)
# # pdf(file = "ROC_SVM_knn.pdf",width = 5,height = 5)
# plot.roc(testLabel, pred, print.auc=T, main="pAUC")
# legend("bottomright", legend=c("Accuracy=0.942", "Precision=0.750 ", "Sensitivity=0.360", "Specificity=0.990", "F-measure=0.486"))
# # dev.off()




