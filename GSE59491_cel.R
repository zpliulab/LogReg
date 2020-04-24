############################# ���� R �� #############################
library(tidyverse)
library(caret)     
library(glmnet)
library("l0ara")
library(ggplot2)
library(caTools)
library(magrittr)
library(ROCR)
library(pROC)
library(glmnetUtils)
library(ncpen)
library(stargazer) #ת��latex 
library(broom) #�ع�������
library(ncvreg)
library(plyr)
library(pROC)


setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37")
# load(".RData")

############################# ��������  #############################

x = read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data\\GSE59491_scale_pDE.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))

## �洢Ԥ����
coef_ridge <- matrix()       # �洢ϵ�����
pred_ridge <- matrix()       # �洢Ԥ����

## ���� ���� ѵ����+���Լ�
set.seed(1234)
training.samples <- data$Lable %>% createDataPartition(p = 0.8, list = FALSE)
# View(training.samples)
train.data  <- data[training.samples, ]
test.data <- data[-training.samples, ]


x <- model.matrix(Lable ~., train.data)[,-1]   # ɾ����Lable~.
y <- train.data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable

############################ glmnet ģ�� Ridge �ͷ� ##########################

set.seed(1234)
## ������֤���õ�ģ�Ͳ���
cv.ridge = cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 10)
lambda.min <- cv.ridge$lambda.min
lambda.1se <- cv.ridge$lambda.1se
print("***lambda.min��lambda.1se***")
print(lambda.min)
print(lambda.1se)

## ���
ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.1se)
coef_ridge <- as.matrix(coef(ridge.model))
## Ԥ��
p <- predict(ridge.model, newx = x.test, type = "response")
pred_ridge <- cbind(y.test, p)

## ���ô洢·��
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\ridge")
write.csv(coef_ridge, file = "coef_ridge.csv")
write.csv(pred_ridge, file = "pred_ridge.csv")

## ROC curve
jpeg(file = "pAUC_ridge.jpg")
plot.roc(pred_ridge[,1], pred_ridge[,2], print.auc=T, main="pAUC")
dev.off()

## ����ָ��
predict = ifelse(pred_ridge[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_ridge[,1]
error = predict_value-true_value

# ����ģ��׼ȷ�ԣ�accuracy������ȷ�ȣ�Precision�����ٻ��ʣ�Recall��-- �����ԣ�sensitivity��-- �������ʣ�TPR���� F��ȣ�F-measure���ͻ�������
# Precision����������б���������item��TP+FP����,"Ӧ�ñ���������item��TP����ռ�ı�����
# Recall����������м�������item��TP��ռ����"Ӧ�ñ���������item��TP+FN��"�ı���
# һ����˵��Precision���Ǽ�����������Ŀ�����磺�ĵ�����ҳ�ȣ��ж�����׼ȷ�ģ�Recall��������׼ȷ����Ŀ�ж��ٱ�����������
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #��ʵֵԤ��ֵȫΪ1 / Ԥ��ֵȫΪ1 --- ��ȡ������ȷ��Ϣ����/��ȡ������Ϣ����
recall = sum(predict_value & true_value)/sum(true_value)        #��ʵֵԤ��ֵȫΪ1 / ��ʵֵȫΪ1 --- ��ȡ������ȷ��Ϣ���� /�����е���Ϣ����
# P��Rָ����ʱ�����ֵ�ì�ܵ��������������Ҫ�ۺϿ������ǣ�����ķ�������F-Measure���ֳ�ΪF-Score��
F_measure= 2*precision*recall/(precision+recall)    #F-Measure��Precision��Recall��Ȩ����ƽ������һ���ۺ�����ָ��
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_ridge[,2], pred_ridge[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "table.csv")     

setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")
## ������
result <- matrix(0, 7, 6)
colnames(result) <- c("accuracy", "precision", "recall", "F_measure", "specificity", "AUC")
rownames(result) <- c("ridge", "lasso", "elatic net", "L1/2", "L0", "SCAD", "MCP")

i <- 1
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result.csv") 

######################### glmnet ģ�� Lasso �ͷ� #############################

## ��������
set.seed(12345)

## ������֤���õ�ģ�Ͳ���
cv.lasso = cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 10 )
lambda.min <- cv.lasso$lambda.min
lambda.1se <- cv.lasso$lambda.1se
print("***lambda.min��lambda.1se***")
print(lambda.min)
print(lambda.1se)

## ���
lasso.model <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.1se)
coef_lasso <- as.matrix(coef(lasso.model))        # ���Ǿ��󲻿���

## Ԥ��
p <- predict(lasso.model, newx = x.test, type = "response")
pred_lasso <- cbind(y.test, p)   #temp���к�pred�ϲ�

## ���ô洢·��
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\lasso")
write.csv(coef_lasso, file = "coef_lasso.csv")
write.csv(pred_lasso, file = "pred_lasso.csv")

## ROC curve
jpeg(file = "pAUC_lasso.jpg")
plot.roc(pred_lasso[,1], pred_lasso[,2], print.auc=T, main="pAUC")
dev.off()

## ����ָ��
predict = ifelse(pred_lasso[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_lasso[,1]
error = predict_value-true_value

# ����ģ��׼ȷ�ԣ�accuracy������ȷ�ȣ�Precision�����ٻ��ʣ�Recall��-- �����ԣ�sensitivity��-- �������ʣ�TPR���� F��ȣ�F-measure���ͻ�������
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)   
recall = sum(predict_value & true_value)/sum(true_value)       
F_measure= 2*precision*recall/(precision+recall)     
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_lasso[,2], pred_lasso[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "table.csv")     

## ������
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")
i <- 2
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result.csv") 

######################### glmnet ģ�� Elastic Net �ͷ� ##########################
#################################################################################

set.seed(12345)

## ������֤���õ�ģ�Ͳ���
cv.elastic = cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = 10)
lambda.min <- cv.elastic$lambda.min
lambda.1se <- cv.elastic$lambda.1se
print("***lambda.min��lambda.1se***")
print(lambda.min)
print(lambda.1se)

## ���
elastic.model <- glmnet(x, y, alpha = 0.5, family = "binomial", lambda = cv.elastic$lambda.1se)
coef_elastic <- as.matrix(coef(elastic.model))        # ���Ǿ��󲻿���

## Ԥ��
p <- predict(elastic.model, newx = x.test, type = "response")
pred_elastic <- cbind(y.test, p)   #temp���к�pred�ϲ�

## ���ô洢·��
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\elastic_net")
write.csv(coef_elastic, file = "coef_elastic.csv")
write.csv(pred_elastic, file = "pred_elastic.csv")


## ROC curve
jpeg(file = "pAUC_elastic.jpg")
plot.roc(pred_elastic[,1], pred_elastic[,2], print.auc=T, main="pAUC")
dev.off()

## ����ָ��
predict = ifelse(pred_elastic[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_elastic[,1]
error = predict_value-true_value


accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)         
F_measure= 2*precision*recall/(precision+recall)     
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_elastic[,2], pred_elastic[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "table.csv")     

## ������
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")
i <- 3
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result.csv") 

######################### ncpen ģ�� Bridge 1/2 �ͷ�#############################
#################################################################################

## ��������
set.seed(123456)
 
## ��ϣ��õ�ģ�Ͳ���
Bridge.model <- ncpen(y.vec = y, x.mat = x, family="binomial", penalty="mbridge")
opt.lambda <- gic.ncpen(Bridge.model, pch="*", type="b")$opt.lambda

print("*** optional lambda ***")
print(opt.lambda)

## ��ȡϵ��
coef <- coef(Bridge.model)
coef_Bridge <- as.matrix(coef[, dim(coef)[2]])

## Ԥ��
p <- predict(Bridge.model, "prob", new.x.mat = x.test)
p <- p[,dim(p)[2]]
pred_Bridge <- cbind(y.test, p)   #temp���к�pred�ϲ�

## ���ô洢·��
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\L0.5")
write.csv(coef_Bridge, file = "coef_Bridge.csv")
write.csv(pred_Bridge, file = "pred_Bridge.csv")

## ROC curve
jpeg(file = "pAUC_Bridge.jpg")
plot.roc(pred_Bridge[,1], pred_Bridge[,2], print.auc=T, main="pAUC")
dev.off()

## ����ָ��
predict = ifelse(pred_Bridge[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_Bridge[,1]
error = predict_value-true_value

# ����ģ��׼ȷ�ԣ�accuracy������ȷ�ȣ�Precision�����ٻ��ʣ�Recall��-- �����ԣ�sensitivity��-- �������ʣ�TPR���� F��ȣ�F-measure���ͻ�������
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)   
recall = sum(predict_value & true_value)/sum(true_value)         
F_measure= 2*precision*recall/(precision+recall)     
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_Bridge[,2], pred_Bridge[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "table.csv")     

## ������
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")
i <- 4
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result.csv") 


######################### ncvreg ģ�� SCAD �ͷ� #################################
#################################################################################

X <- model.matrix(Lable ~., train.data)[,-1]   # ɾ����Lable~.
y <- train.data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable

library(ncvreg)

## ��������
set.seed(12345)

## ������֤���õ�ģ�Ͳ���
cv.SCAD <- cv.ncvreg(X, y, family ="binomial", penalty="SCAD") 
lambda.min <- cv.SCAD$lambda.min 

print("*** lambda.min ***")
print(lambda.min)

## ���
SCAD.model <- ncvreg(X, y, family ="binomial", lambda = cv.SCAD$lambda.min, penalty="SCAD")
coef_SCAD <- coef(SCAD.model)

## Ԥ��
p <- predict(SCAD.model, x.test, type = "response")
pred_SCAD <- cbind(y.test, p)   #temp���к�pred�ϲ�


## ���ô洢·��
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\SCAD")

## �洢 ������ y.test��p �� ���� 
write.csv(coef_SCAD, file = "coef_SCAD.csv")
write.csv(pred_SCAD, file = "pred_SCAD.csv")

## ROC curve
jpeg(file = "pAUC_SCAD.jpg")
plot.roc(pred_SCAD[,1], pred_SCAD[,2], print.auc=T, main="pAUC")
dev.off()
                         
predict = ifelse(pred_SCAD[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_SCAD[,1]
error = predict_value-true_value

 
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)   
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure= 2*precision*recall/(precision+recall)    
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_SCAD[,2], pred_SCAD[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "table.csv")     

## ������
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")
i <- 6
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result.csv")   


######################### ncvreg ģ�� MCP �ͷ� #############################
############################################################################

library(ncvreg)

## ��������
set.seed(12345)

## ������֤���õ�ģ�Ͳ���
cv.MCP <- cv.ncvreg(X, y, family ="binomial", penalty="MCP") 
lambda.min <- cv.MCP$lambda.min 

print("*** lambda.min ***")
print(lambda.min)

## ���
MCP.model <- ncvreg(X, y, family ="binomial", lambda = cv.MCP$lambda.min, penalty="MCP")
coef_MCP <- coef(MCP.model)

## Ԥ��
p <- predict(MCP.model, x.test, type = "response")
pred_MCP <- cbind(y.test, p)   #temp���к�pred�ϲ�


## ���ô洢·��
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\MCP")
write.csv(coef_MCP, file = "coef_MCP.csv")
write.csv(pred_MCP, file = "pred_MCP.csv")

## ROC curve
jpeg(file = "pAUC_MCP.jpg")
plot.roc(pred_MCP[,1], pred_MCP[,2], print.auc=T, main="pAUC")
dev.off()

## ����ָ��
predict = ifelse(pred_MCP[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_MCP[,1]
error = predict_value-true_value


accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)   
recall = sum(predict_value & true_value)/sum(true_value)         
F_measure= 2*precision*recall/(precision+recall)  
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_MCP[,2], pred_MCP[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "table.csv")     

## ������
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")
i <- 7
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result.csv")


###################################################################################################################
###################################################################################################################

stargazer(result)
# 
###################### l0ara ģ�� Logistic �ع� + L0 �ͷ� #########################################################
###################################################################################################################

library("l0ara")

## �洢Ԥ����
coef_L0 <- matrix()       # �洢ϵ�����
pred_L0 <- data.frame()   # �洢Ԥ����

for(i in 1:10){
  # i = 1
  test.data <- data[folds[[i]],]     # ȡfolds[[i]]��Ϊ���Լ�
  train.data <- data[-folds[[i]],]   # ʣ�µ�������Ϊѵ����

  print("***���***")
  print(i)

  x <- model.matrix(Lable ~., train.data)[,-1]
  y <- train.data$Lable
  x.test <- model.matrix(Lable ~., test.data)[,-1]
  y.test <- test.data$Lable

  ## ��������
  set.seed(12345)

  ## ������֤���õ�ģ�Ͳ���
  # lam <- seq(1,0.05,-0.05)
  # lam <- seq(1,0.5,-0.5)
  # lam <- seq(1,0.005,-0.005)
  # lam <- seq(0.1, 2, 0.1)
  lam <- seq(1e-2, 1, 1e-2)
  # lam <- seq(0.1, 1, 0.1)
  cv.l0 <- cv.l0ara(x, y, family="logit", lam, measure = "mse")
  lambda.min <- cv.l0$lam.min
  print("*** lambda.min ***")
  print(lambda.min)

  ## ���
  l0.model <- l0ara(x, y, family = "logit", lambda.min)

  ## ��ȡϵ��
  coef_l0 <- coef(l0.model)
  coef_l0 = as.matrix(coef_l0)

  new_coef_l0 <- matrix()
  for (j in 2:nrow(coef_l0)){
    new_coef_l0[j+1] <- coef_l0[j]
  }
  for (k in 1:2){
    new_coef_l0[k] <- coef_l0[k]
  }

  coef <- as.matrix(new_coef_l0)

  tcross <- rep(i, length(coef))                 # i�ǵڼ���ѭ�����棬��K��
  lable <- as.matrix(rownames(coef_lasso))
  step_L0 <- data.frame(cbind(coef, tcross))
  coef_L0 <- cbind(coef_L0, lable, step_L0)      #temp���к�pred�ϲ�

  ## Ԥ��
  p <- predict(l0.model, newx = x.test, type = "response")

  ## AUC
  p = as.matrix(p)
  pred <- prediction(p, y.test)
  auc <- performance(pred,'auc')
  auc <- unlist(slot(auc,'y.values'))
  print("***auc***")
  print(auc)

  ## �ۼ�
  kcross <- rep(i, length(p))
  temp_L0 <- data.frame(cbind(y.test, p, kcross)) # ��ʵֵ��Ԥ��ֵ�����ɭ��������Ԥ������������һ������µ����ݿ�tenp
  pred_L0 <- rbind(pred_L0, temp_L0)              # temp���к�pred�ϲ�


}

View(pred_L0)
View(coef_L0)


## ���ô洢·��
setwd("C:\\Users\\LiLingyu\\Desktop\\GSE59491_5\\L0")

## �洢 ������ y.test��p �� ����
write.csv(coef_L0, file = "coef_L0.csv")
write.csv(pred_L0, file = "pred_L0.csv")


folds[[1]]
## ROC curve
jpeg(file = "pAUC_L01.jpg")
plot.roc(pred_L0[,1], pred_L0[,2], print.auc=T, main="pAUC")
dev.off()

## ����ָ��
predict = ifelse(pred_L0[,2] > 0.5, 1, 0)
predict_value = predict
true_value = pred_L0[,1]
error = predict_value-true_value

accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
precision = sum(true_value & predict_value)/sum(predict_value)  
recall = sum(predict_value & true_value)/sum(true_value)        
F_measure = 2*precision*recall/(precision+recall)    
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))
# AUC
pred <- prediction(pred_L0[,2], pred_L0[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## ����������ʾ�������ΪTP��FN��FP��TN
sink("C:\\Users\\LiLingyu\\Desktop\\GSE59491_5\\L0\\table.txt")
table(true_value, predict_value)
sink()

## ������
setwd("C:\\Users\\LiLingyu\\Desktop\\GSE59491_5\\result")
i <- 5
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result.csv")


########################################  ����ϲ� ###########################################################
##############################################################################################################

setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")

## ���� gene ����

coef_ridge <- read.csv("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\ridge\\coef_ridge.csv", header=TRUE, sep = ',')
coef_lasso <- read.csv("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\lasso\\coef_lasso.csv", header=TRUE, sep = ',')
coef_elastic <- read.csv("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\\\elastic_net\\coef_elastic.csv", header=TRUE, sep = ',')
coef_Bridge <- read.csv("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\l0.5\\coef_Bridge.csv", header=TRUE, sep = ',')
coef_l0 <- read.csv("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\l0\\coef_L0.csv", header=TRUE, sep = ',')
coef_SCAD <- read.csv("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\SCAD\\coef_SCAD.csv", header=TRUE, sep = ',')
coef_MCP <- read.csv("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\MCP\\coef_MCP.csv", header=TRUE, sep = ',')


sum(coef_ridge[,2] != 0)
sum(coef_lasso[,2] != 0)
sum(coef_elastic[,2] != 0)
sum(coef_l0[,2] != 0)
sum(coef_Bridge[,2] != 0)
sum(coef_SCAD[,2] != 0)
sum(coef_MCP[,2] != 0)


ridge <- coef_ridge[which(coef_ridge[,2] != 0),1]
lasso <- coef_ridge[which(coef_lasso[,2] != 0),1]
elastic <- coef_ridge[which(coef_elastic[,2] != 0),1]
ridge <- coef_ridge[which(coef_l0[,2] != 0),1]
Bridge <- coef_ridge[which(coef_Bridge[,2] != 0),1]
SCAD <- coef_ridge[which(coef_SCAD[,2] != 0),1]
MCP <- coef_ridge[which(coef_MCP[,2] != 0),1]


setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")

# gene <- intersect(intersect(intersect(intersect(coef_ridge, coef_lasso), coef_Elastic), coef_SCAD), coef_MCP)
gene <- intersect(intersect(intersect(intersect(intersect(ridge, lasso), elastic), Bridge), SCAD), MCP)
# gene <-intersect(intersect(intersect(intersect(intersect(intersect(coef_ridge, coef_lasso), coef_Elastic), coef_Bridge), coef_SCAD), coef_MCP), coef_l0)
# write.csv(gene, "gene6.csv")

View(gene)

stargazer(gene) #ʹ��stargazer����LaTeX�����µ�������ͳ�Ʊ���

###################################################################################################################

## ��������

setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\result")


##################################    ���ƹ����� ROC ����    ######################################################
###################################################################################################################
## ���� A ����
setwd("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\result")
jpeg(file = "pAUCall.jpg")

coef_ridge <- read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\ridge\\pred_ridge.csv", header=TRUE, sep = ',')
coef_lasso <- read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\lasso\\pred_lasso.csv", header=TRUE, sep = ',')
coef_Elastic <- read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\elastic_net\\pred_elastic.csv", header=TRUE, sep = ',')
coef_Bridge <- read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\l0.5\\pred_Bridge.csv", header=TRUE, sep = ',')
coef_l0 <- read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\l0\\pred_L0.csv", header=TRUE, sep = ',')
coef_SCAD <- read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\SCAD\\pred_SCAD.csv", header=TRUE, sep = ',')
coef_MCP <- read.table("D:\\E\\��ʿ\\R_����\\GSE59491_15\\Data37\\MCP\\pred_MCP.csv", header=TRUE, sep = ',')


library(pROC)

roc_ridge <- plot.roc( coef_ridge[,2], coef_ridge[,3], col="1" )
roc_lasso <- lines.roc( coef_lasso[,2], coef_lasso[,3], col="2" )
roc_Elastic <- lines.roc( coef_Elastic[,2], coef_Elastic[,3], col="3" )
roc_Bridge <- lines.roc( coef_Bridge[,2], coef_Bridge[,3], col="4" )
roc_l0 <- lines.roc( coef_l0[,2], coef_l0[,3], col="5" )
roc_SCAD <- lines.roc( coef_SCAD[,2], coef_SCAD[,3], col="6" )
roc_MCP <- lines.roc( coef_MCP[,2], coef_MCP[,3], col="7" )


legend("bottomright", legend=c("ridge", "lasso", "elastic net", "Bridge", "SCAD", "MCP"), col=c("1", "2", "3", "4",  "6", "7"), lwd=6)

dev.off()

