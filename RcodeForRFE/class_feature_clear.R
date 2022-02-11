setwd('D:\\E\\博士\\R_程序\\SVM_RFE\\Data\\limma')


# GSE73685 7个gene (17个样本)-----------------------------------------------------------------
Data  = read.table("GSE73685_46.txt", header = T, check.names = FALSE)
Data2 = read.table("GSE59491_46.txt", header = T, check.names = FALSE)


# 独立数据集 -------------------------------------------------------------------
x.train <- data.frame(t(Data2)[,-1])
y.train <- t(Data2)[,1]
x.test <- data.frame(t(Data)[,-1])
y.test <- t(Data)[,1]


# SVM ---------------------------------------------------------------------

library(e1071)
set.seed(666) 
tuned <- tune.svm(x.train,y.train, gamma = 10^(-6:-1), cost = 10^(1:2)) # tune
summary (tuned) # to select best gamma and cost

model <- svm(x.train, y.train, kernel = "radial", cost = tuned$best.parameters$cost, gamma=tuned$best.parameters$gamma,  scale = FALSE)
summary(model)

##################################### 利用svm() 函数建立的模型进行预测
p_test <- predict(model, x.test, type = "response")
pre <- cbind(p_test,y.test)
colnames(pre) <- c('y.test', 'Lable')
# write.table(pre,"pre73685_glm37.txt",quote=F,sep="\t")


# test --------------------------------------------------------------------
library(pROC)
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, y.test)
names(A_test)<- c("p", "outcome")
# pdf(file = "ROC_independent46.pdf",width = 5,height = 5)
plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
# show_25 <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
# s1 <- smooth(show_25,method="binormal")
# plot(s1)

# plot.roc(A_test$outcome, A_test$p)
legend("bottomright", legend=c("Accuracy = 0.882", "Precision = 1.000 ", "Sensitivity = 0.600", "Specificity = 1.000", "F-measure = 0.750", "AUC = 0.817"))
# dev.off()


## 性能
predict = ifelse(pred_glm[,1] > 0.40, 1, 0)
predict_value = predict
true_value = pred_glm[,2]
error = predict_value-true_value

data <- t(Data)
accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure = 2*precision*recall/(precision+recall)   
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))

accuracy
precision
recall
F_measure
specificity

table(true_value, predict_value) 
