
Data = read.table("E:\\GSE59491_15\\Data37\\DE\\p73685.txt", header = T, check.names = FALSE)
Data2 = read.table("E:\\GSE59491_15\\Data37\\DE\\p59491_20.txt", header = T, check.names = FALSE)



# �������ݼ� -------------------------------------------------------------------
x.train <- data.frame(t(Data2)[,-1])
y.train <- t(Data2)[,1]
x.test <- data.frame(t(Data)[,-1])
y.test <- t(Data)[,1]

glm.fit <- glm(y.train~., data = x.train, family = binomial)
# glm.fit <- glm(y.train~., data = x.train, family = binomial, control = list(maxit = 100))
summary(glm.fit)


setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\DE")

p_test <- predict(glm.fit, x.train, type = "response")

pred_glm <- cbind(p_test, y.train)
colnames(pred_glm) <- c('y.test', 'Lable')
p <- pred_glm[,1]
p_glm <- cbind(log(p/(1-p)), pred_glm[,2])
colnames(p_glm) <- c('y.test', 'Lable')

library(pROC)
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, y.test)
names(A_test)<- c("p", "outcome")
plot.roc(A_test$outcome, A_test$p)
legend("bottomright", legend=c("Acc=0.875", "Pre=1.000 ", "Sn=0.75", "F-measure=0.857", "Sp=1.000", "AUC=0.875"))

## ����
predict = ifelse(pred_glm[,1] > 0.3, 1, 0)
predict_value = predict
true_value = pred_glm[,2]
error = predict_value-true_value

data <- t(Data)
# ����ģ��׼ȷ�ԣ�accuracy������ȷ�ȣ�Precision�����ٻ��ʣ�Recall��-- �����ԣ�sensitivity��-- �������ʣ�TPR���� F��ȣ�F-measure���ͻ�������
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

## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) 