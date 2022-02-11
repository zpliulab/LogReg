##  clear
rm(list = ls())


############################# 下载 R 包 #############################
library(tidyverse)
library(caret)   # 十折交叉验证     
library(glmnet)
library("l0ara")
library(ggplot2)
library(caTools)
library(magrittr)
library(ROCR)
library(pROC)
library(glmnetUtils)
library(ncpen)
library(stargazer) #转换latex 
library(broom) #回归结果保存
library(ncvreg)
library("l0ara")
# install.packages("caret")
library(plyr)
library(pROC)


############################# 读入数据  #############################

x = read.table("D:\\E\\博士\\R_程序\\GSE59491_15\\Data\\GSE59491_scale_DE.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))


## Data split —— training data + testing data
set.seed(1234)
training.samples <- data$Lable %>% createDataPartition(p = 0.7, list = FALSE)
train.data  <- data[training.samples, ]
test.data <- data[-training.samples, ]


x <- model.matrix(Lable ~., train.data)[,-1]   # 删除了Lable~.
y <- train.data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable
# sum(y == 0)


# pathways ----------------------------------------------------------------
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37")
set.seed(123)

ridge.fit <- glmnet(x, y, family="binomial", alpha = 0, lambda = NULL) 
jpeg(file = "ridge_fit.jpg")
# postscript("ridge_fit.eps")
plot(ridge.fit, xvar = "lambda")
dev.off()
cv.ridge = cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 10, type.measure = "class")
# jpeg(file = "cv_ridge.jpg")
postscript("cv_ridge.eps")
plot(cv.ridge)
dev.off()

lasso.fit <- glmnet(x, y, family="binomial", alpha = 1, lambda = NULL) 
jpeg(file = "lasso_fit.jpg")
# postscript("lasso_fit.eps")
plot(lasso.fit, xvar = "lambda")
dev.off()
cv.lasso = cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 10, type.measure = "class")
jpeg(file = "cv_lasso.jpg")
# postscript("cv_lasso.eps")
plot(cv.lasso)
dev.off()

elastic.fit <- glmnet(x, y, family="binomial", alpha = 0.5, lambda = NULL) 
jpeg(file = "elastic_fit.jpg")
# postscript("elastic_fit.eps")
plot(elastic.fit, xvar = "lambda")
dev.off()
cv.elastic = cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = 10, type.measure = "class")
jpeg(file = "cv_elastic.jpg")
# postscript("cv_elastic.eps")
plot(cv.elastic)
dev.off()

cv.Bridge <- cv.ncpen(y.vec=y, x.mat=x, family="binomial", penalty="mbridge")
jpeg(file = "cv_Bridge.jpg")
# postscript("cv_Bridge.eps")
plot(cv.Bridge, type = "rmse", log.scale = T)
dev.off()
Bridge.model <- ncpen(y.vec = y, x.mat = x, family="binomial", penalty="mbridge")
jpeg(file = "Bridge_fit.jpg")
# postscript("Bridge_fit.eps")
plot(Bridge.model)
dev.off()

lam <- seq(1,0.05,-0.05)
cv.l0 <- cv.l0ara(x, y, family="logit", lam, measure = "class")
lambda.min <- cv.l0$lam.min
jpeg(file = "cv_l0.jpg")
# postscript("cv_l0.eps")
plot(cv.l0, col = 2)
dev.off()

l0.fit <- l0ara(x, y, family = "logit", 0.4)
jpeg(file = "l0_fit.jpg")
# postscript("l0_fit.eps")
plot(l0.fit, auc = F, split = F, col = 4) # 描绘 解 路径的非局部凸起的区域
dev.off()



## X, y 
X <- model.matrix(Lable ~., train.data)[,-1]   # 删除了Lable~.
y <- train.data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable

cv.SCAD <- cv.ncvreg(X, y, family ="binomial", penalty="SCAD") 
jpeg(file = "cv_SCAD.jpg")
# postscript("cv_SCAD.eps")
plot(cv.SCAD)
dev.off()
SCAD.fit <- ncvreg(X, y, family ="binomial", penalty="SCAD")
jpeg(file = "SCAD_fit.jpg")
# postscript("SCAD_fit.eps")
plot(SCAD.fit)
dev.off()

cv.MCP <- cv.ncvreg(X, y, family ="binomial", penalty="MCP")
jpeg(file = "cv_MCP.jpg")
# postscript("cv_MCP.eps")
plot(cv.MCP)
dev.off()
MCP.fit <- ncvreg(X, y, family ="binomial", penalty="MCP")
jpeg(file = "MCP_fit.jpg")
# postscript("MCP_fit.eps")
plot(MCP.fit)
dev.off()

############################ glmnet 模拟 Ridge 惩罚 ##########################
## 存储预测结果
coef_ridge <- matrix()       # 存储系数结果
pred_ridge <- matrix()   # 存储预测结果


for(i in 1:30){
  
  ## 设置种子
  set.seed(i)
  
  ## 交叉验证，得到模型参数
  cv.ridge = cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 10, type.measure = "class")
  lambda.min <- cv.ridge$lambda.min
  lambda.1se <- cv.ridge$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  
  ## 拟合
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.1se)
  coef <- as.matrix(coef(ridge.model))
  tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
  step_ridge <- data.frame(cbind(coef, tcross))
  coef_ridge <- cbind(coef_ridge, step_ridge)   #temp按行和pred合并
  
  ## 预测
  p <- predict(ridge.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_ridge <- data.frame(cbind(y.test, p, kcross))
  pred_ridge <- cbind(pred_ridge,temp_ridge)   #temp按行和pred合并
  
  print(paste("第：",i)) 
}

## 设置存储路径
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37")
# write.csv(coef_ridge, file = "ridge\\coef_ridge.csv")
# write.csv(pred_ridge, file = "ridge\\pred_ridge.csv")

# my_pred 函数 输入p_ridge
my_pred <- function(x){
  p_ridge1 <- x[,-2]
  p_ridge2 <- matrix(data=0, nrow = dim(p_ridge1)[1], ncol = 1, byrow = FALSE, dimnames=list(c(as.character(x[,1])),c("prob")))
  for (j in 1:10){
    p_ridge2 <- p_ridge2[] + p_ridge1[,3*j]
  }
  p_ridge3 <- p_ridge2/10
  return(p_ridge3)
}


# compute average
p_ridge <- read.table("ridge\\pred_ridge.csv", header=TRUE, sep = ',')
pred_ridge <- cbind(y.test,my_pred(p_ridge))
# write.csv(pred_ridge, file = "ridge\\pred_ridge0.csv")

## ROC curve
# jpeg(file = "ridge\\pAUC_ridge.jpg")
plot.roc(pred_ridge[,1], pred_ridge[,2], print.auc=T, main="pAUC")
# dev.off()

## 性能指标
predict = ifelse(pred_ridge[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_ridge[,1]
error = predict_value-true_value

# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
# Precision计算的是所有被检索到的item（TP+FP）中,"应该被检索到的item（TP）”占的比例；
# Recall计算的是所有检索到的item（TP）占所有"应该被检索到的item（TP+FN）"的比例
# 一般来说，Precision就是检索出来的条目（比如：文档、网页等）有多少是准确的，Recall就是所有准确的条目有多少被检索出来了
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
# P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
F_measure= 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_ridge[,2], pred_ridge[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "ridge\\table.csv")     


## 输出结果
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
# result %>% write.csv(file = "result\\result.csv") 


######################### glmnet 模拟 Lasso 惩罚 #############################
## 存储预测结果
coef_lasso <- matrix()       # 存储系数结果
pred_lasso <- matrix()   # 存储预测结果

for(i in 1:30){
  
  ## 设置种子
  set.seed(i)
  
  ## 交叉验证，得到模型参数
  cv.lasso = cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 10, type.measure = "class" )
  lambda.min <- cv.lasso$lambda.min
  lambda.1se <- cv.lasso$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  
  ## 拟合
  lasso.model <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.1se)
  coef <- as.matrix(coef(lasso.model))
  tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
  step_lasso <- data.frame(cbind(coef, tcross))
  coef_lasso <- cbind(coef_lasso, step_lasso)   #temp按行和pred合并
  
  ## 预测
  p <- predict(lasso.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_lasso <- data.frame(cbind(y.test, p, kcross)) # 真实值、预测值、随机森林树数、预测组编号捆绑在一起组成新的数据框tenp
  pred_lasso <- cbind(pred_lasso,temp_lasso)   #temp按行和pred合并

  print(paste("第：",i)) 
}
  
## 设置存储路径
# write.csv(coef_lasso, file = "lasso\\coef_lasso.csv")
# write.csv(pred_lasso, file = "lasso\\pred_lasso.csv")

# 求平均
p_lasso <- read.table("lasso\\pred_lasso.csv", header=TRUE, sep = ',')
pred_lasso <- cbind(y.test,my_pred(p_lasso))
# write.csv(pred_lasso, file = "lasso\\pred_lasso0.csv")

## ROC curve
# jpeg(file = "lasso\\pAUC_lasso.jpg")
plot.roc(pred_lasso[,1], pred_lasso[,2], print.auc=T, main="pAUC")
# dev.off()

## 性能指标
predict = ifelse(pred_lasso[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_lasso[,1]
error = predict_value-true_value

# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
# P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
F_measure= 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_lasso[,2], pred_lasso[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "lasso\\table.csv")     

## 输出结果
i <- 2
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv") 

######################### glmnet 模拟 Elastic Net 惩罚 #############################

## 存储预测结果
coef_elastic <- matrix()       # 存储系数结果
pred_elastic <- matrix()   # 存储预测结果

for(i in 1:30){
  
  ## 设置种子
  set.seed(i)
  
  ## 交叉验证，得到模型参数
  cv.elastic = cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = 10, type.measure = "class")
  lambda.min <- cv.elastic$lambda.min
  lambda.1se <- cv.elastic$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  
  ## 拟合
  elastic.model <- glmnet(x, y, alpha = 0.5, family = "binomial", lambda = cv.elastic$lambda.1se)
  coef <- as.matrix(coef(elastic.model))
  tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
  step_elastic <- data.frame(cbind(coef, tcross))
  coef_elastic <- cbind(coef_elastic, step_elastic)   #temp按行和pred合并
  
  ## 预测
  p <- predict(elastic.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_elastic <- data.frame(cbind(y.test, p, kcross))
  
  pred_elastic <- cbind(pred_elastic,temp_elastic)   #temp按行和pred合并
  
  print(paste("第：",i)) 
}
  
## 设置存储路径
# write.csv(coef_elastic, file = "elastic_net\\coef_elastic.csv")
# write.csv(pred_elastic, file = "elastic_net\\pred_elastic.csv")

# 求平均 
p_elastic <- read.table("elastic_net\\pred_elastic.csv", header=TRUE, sep = ',')
pred_elastic <- cbind(y.test,my_pred(p_elastic))
# write.csv(pred_elastic, file = "elastic_net\\pred_elastic0.csv")

## ROC curve
# jpeg(file = "elastic_net\\pAUC_elastic.jpg")
plot.roc(pred_elastic[,1], pred_elastic[,2], print.auc=T, main="pAUC")
# dev.off()

## 性能指标
predict = ifelse(pred_elastic[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_elastic[,1]
error = predict_value-true_value

# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
# Precision计算的是所有被检索到的item（TP+FP）中,"应该被检索到的item（TP）”占的比例；
# Recall计算的是所有检索到的item（TP）占所有"应该被检索到的item（TP+FN）"的比例
# 一般来说，Precision就是检索出来的条目（比如：文档、网页等）有多少是准确的，Recall就是所有准确的条目有多少被检索出来了
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
# P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
F_measure= 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_elastic[,2], pred_elastic[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "elastic_net\\table.csv")     

## 输出结果
i <- 3
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv") 

######################### ncpen 模拟 Bridge 1/2 惩罚#############################

## 存储预测结果
coef_Bridge <- matrix()       # 存储系数结果
pred_Bridge <- matrix()   # 存储预测结果

for(i in 1:30){

  ## 设置种子
  set.seed(i)
  
  ## 交叉验证，得到模型参数
  # cv.Bridge <- cv.ncpen(y.vec=y, x.mat=x, family="binomial", penalty="mbridge")
  # plot(cv.Bridge, type = "rmse", log.scale = T)
  # coef(cv.Bridge)$lambda
  
  ## 拟合，得到模型参数
  Bridge.model <- ncpen(y.vec = y, x.mat = x, family="binomial", penalty="mbridge")
  opt.lambda <- gic.ncpen(Bridge.model, pch="*", type="b")$opt.lambda
  
  print("*** optional lambda ***")
  print(opt.lambda)
  
  ## 提取系数
  coef <- coef(Bridge.model)
  # coef_Bridge <- as.matrix(coef[, dim(coef)[2]])
  coef <- as.matrix(coef[, dim(coef)[2]])
  tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
  step_Bridge <- data.frame(cbind(coef, tcross))
  coef_Bridge <- cbind(coef_Bridge, step_Bridge)   #temp按行和pred合并
  ## 预测
  p <- predict(Bridge.model, "prob", new.x.mat = x.test)
  p <- p[,dim(p)[2]]
  kcross <- rep(i, length(p)) 
  temp_Bridge <- data.frame(cbind(y.test, p, kcross))
  pred_Bridge <- cbind(pred_Bridge,temp_Bridge)   #temp按行和pred合并
  
  print(paste("第：",i)) 
}


## 设置存储路径
# write.csv(coef_Bridge, file = "L0.5\\coef_Bridge.csv")
# write.csv(pred_Bridge, file = "L0.5\\pred_Bridge.csv")
# 求平均
p_Bridge <- read.table("L0.5\\pred_Bridge.csv", header=TRUE, sep = ',')
pred_Bridge <- cbind(y.test,my_pred(p_Bridge))
# write.csv(pred_Bridge, file = "L0.5\\pred_Bridge0.csv")

## ROC curve
# jpeg(file = "L0.5\\pAUC_Bridge.jpg")
plot.roc(pred_Bridge[,1], pred_Bridge[,2], print.auc=T, main="pAUC")
# dev.off()

## 性能指标
predict = ifelse(pred_Bridge[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_Bridge[,1]
error = predict_value-true_value

# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
F_measure= 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_Bridge[,2], pred_Bridge[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "L0.5\\table.csv")     

## 输出结果
i <- 4
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv") 

######################### ncvreg 模拟 SCAD 惩罚 #############################

X <- model.matrix(Lable ~., train.data)[,-1]   # 删除了Lable~.
y <- train.data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)
x.test <- model.matrix(Lable ~., test.data)[,-1] 
y.test <- test.data$Lable


library(ncvreg)

## 存储预测结果
coef_SCAD <- matrix()       # 存储系数结果
pred_SCAD <- matrix()   # 存储预测结果

for(i in 1:30){
  
  set.seed(i)
  ## 交叉验证，得到模型参数
  cv.SCAD <- cv.ncvreg(X, y, family ="binomial", penalty="SCAD") 
  lambda.min <- cv.SCAD$lambda.min 
  
  print("*** lambda.min ***")
  print(lambda.min)
  
  ## 拟合
  SCAD.model <- ncvreg(X, y, family ="binomial", lambda = cv.SCAD$lambda.min, penalty="SCAD")
  coef <- as.matrix(coef(SCAD.model))
  tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
  step_SCAD <- data.frame(cbind(coef, tcross))
  coef_SCAD <- cbind(coef_SCAD, step_SCAD)   #temp按行和pred合并
  
  ## 预测
  p <- predict(SCAD.model, x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_SCAD <- data.frame(cbind(y.test, p, kcross))
  pred_SCAD <- cbind(pred_SCAD,temp_SCAD)   #temp按行和pred合并
  
  print(paste("第：",i)) 
  
}


## 存储 整个的 y.test、p 和 折数 
# write.csv(coef_SCAD, file = "SCAD\\coef_SCAD.csv")
# write.csv(pred_SCAD, file = "SCAD\\pred_SCAD.csv")
# 求平均
p_SCAD <- read.table("SCAD\\pred_SCAD.csv", header=TRUE, sep = ',')
pred_SCAD <- cbind(y.test,my_pred(p_SCAD))
# write.csv(pred_SCAD, file = "SCAD\\pred_SCAD0.csv")

## ROC curve
# jpeg(file = "SCAD\\pAUC_SCAD.jpg")
plot.roc(pred_SCAD[,1], pred_SCAD[,2], print.auc=T, main="pAUC")
# dev.off()

predict = ifelse(pred_SCAD[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_SCAD[,1]
error = predict_value-true_value

# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
# Precision计算的是所有被检索到的item（TP+FP）中,"应该被检索到的item（TP）”占的比例；
# Recall计算的是所有检索到的item（TP）占所有"应该被检索到的item（TP+FN）"的比例
# 一般来说，Precision就是检索出来的条目（比如：文档、网页等）有多少是准确的，Recall就是所有准确的条目有多少被检索出来了
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
# P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
F_measure= 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_SCAD[,2], pred_SCAD[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "SCAD\\table.csv")     

## 输出结果
i <- 6
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
result %>% write.csv(file = "result\\result.csv")

######################### ncvreg 模拟 MCP 惩罚 #############################

library(ncvreg)
## 存储预测结果
coef_MCP <- matrix()       # 存储系数结果
pred_MCP <- matrix()   # 存储预测结果

for(i in 1:30){
  
  ## 设置种子
  set.seed(i)
  ## 交叉验证，得到模型参数
  cv.MCP <- cv.ncvreg(X, y, family ="binomial", penalty="MCP") 
  lambda.min <- cv.MCP$lambda.min 
  
  print("*** lambda.min ***")
  print(lambda.min)
  
  ## 拟合
  MCP.model <- ncvreg(X, y, family ="binomial", lambda = cv.MCP$lambda.min, penalty="MCP")
  coef <- as.matrix(coef(MCP.model))
  tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
  step_MCP <- data.frame(cbind(coef, tcross))
  coef_MCP <- cbind(coef_MCP, step_MCP)   #temp按行和pred合并
  
  ## 预测
  p <- predict(MCP.model, x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_MCP <- data.frame(cbind(y.test, p, kcross)) # 真实值、预测值、随机森林树数、预测组编号捆绑在一起组成新的数据框tenp
  pred_MCP <- cbind(pred_MCP,temp_MCP)   #temp按行和pred合并
  
  print(paste("第：",i)) 
}

## 设置存储路径
# write.csv(coef_MCP, file = "MCP\\coef_MCP.csv")
# write.csv(pred_MCP, file = "MCP\\pred_MCP.csv")
# 求平均 
p_MCP <- read.table("MCP\\pred_MCP.csv", header=TRUE, sep = ',')
pred_MCP <- cbind(y.test,my_pred(p_MCP))
# write.csv(pred_MCP, file = "MCP\\pred_MCP0.csv")

## ROC curve
# jpeg(file = "MCP\\pAUC_MCP.jpg")
plot.roc(pred_MCP[,1], pred_MCP[,2], print.auc=T, main="pAUC")
# dev.off()

## 性能指标
predict = ifelse(pred_MCP[,2] > 0.5, 1, 0)
predict_value = predict 
true_value = pred_MCP[,1]
error = predict_value-true_value

# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
# Precision计算的是所有被检索到的item（TP+FP）中,"应该被检索到的item（TP）”占的比例；
# Recall计算的是所有检索到的item（TP）占所有"应该被检索到的item（TP+FN）"的比例
# 一般来说，Precision就是检索出来的条目（比如：文档、网页等）有多少是准确的，Recall就是所有准确的条目有多少被检索出来了
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
# P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
F_measure= 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 
# AUC
pred <- prediction(pred_MCP[,2], pred_MCP[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "MCP\\table.csv")     

## 输出结果
i <- 7
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv")



stargazer(result)

###################### l0ara 模拟 Logistic 回归 + L0 惩罚 #####################
##############################################################################
# install.packages("l0ara")
## 构建for循环，得到十次交叉验证预测的AUC值。并纪录取值最大的一组，作为最优的训练集与测试集划分。
library("l0ara")

coef_L0 <- matrix()   # 存储系数结果
pred_L0 <- matrix()   # 存储预测结果

for(i in 1:30){
  # i <- 1
  ## 设置种子
  set.seed(i)
  ## 交叉验证，得到模型参数
  # lam <- seq(1,0.05,-0.05)
  lam <- seq(1,0.5,-0.5)
  cv.l0 <- cv.l0ara(x, y, family="logit", lam, measure = "mse")
  lambda.min <- cv.l0$lam.min
  
  print("*** lambda.min ***")
  print(lambda.min)
  
  ## 拟合
  l0.model <- l0ara(x, y, family = "logit", lambda.min)
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
  tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
  step_L0 <- data.frame(cbind(coef, tcross))
  coef_L0 <- cbind(coef_L0, step_L0)   #temp按行和pred合并
  
  ## 预测
  p <- predict(l0.model, newx = x.test, type = "response")
  kcross <- rep(i, length(p)) 
  temp_L0 <- data.frame(cbind(y.test, p, kcross)) 
  pred_L0 <- cbind(pred_L0,temp_L0)   #temp按行和pred合并
  
  print(paste("第：",i)) 
}


## 存储 整个的 y.test、p 和 折数
# write.csv(coef_L0, file = "L0\\coef_L0.csv")
# write.csv(pred_L0, file = "L0\\pred_L0.csv")

# 求平均 
p_L0 <- read.table("L0\\pred_L0.csv", header=TRUE, sep = ',')

my_pred <- function(x){
  p_ridge1 <- x[,-2]
  p_ridge2 <- matrix(data=0, nrow = dim(p_ridge1)[1], ncol = 1, byrow = FALSE, dimnames=list(c(as.character(x[,1])),c("prob")))
  for (j in 1:10){
    p_ridge2 <- p_ridge2[] + p_ridge1[,3*j]
  }
  p_ridge3 <- p_ridge2/10
  return(p_ridge3)
}
pred_L0 <- cbind(y.test,my_pred(p_L0))
# pred_L0 <- cbind(y,my_pred(p_L0))
# write.csv(pred_L0, file = "L0\\pred_L0.csv")

## ROC curve
# jpeg(file = "L0\\pAUC_L01.jpg")
plot.roc(pred_L0[,1], pred_L0[,2], print.auc=T, main="pAUC")
# dev.off()

## 性能指标
predict = ifelse(pred_L0[,2] > 0.5, 1, 0)
predict_value = predict
true_value = pred_L0[,1]
error = predict_value-true_value

# 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
# P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
F_measure = 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))
# AUC
pred <- prediction(pred_L0[,2], pred_L0[,1])
auc <- performance(pred,'auc')
auc <- unlist(slot(auc,'y.values'))

## 混淆矩阵，显示结果依次为TP、FN、FP、TN
# table(true_value, predict_value) %>% as.matrix(table) %>% write.csv(file = "L0\\table.csv")     

## 输出结果
i <- 5
result[i,1] <- accuracy
result[i,2] <- precision
result[i,3] <- recall
result[i,4] <- F_measure
result[i,5] <- specificity
result[i,6] <- auc
# result %>% write.csv(file = "result\\result.csv")



########################################  基因合并 ###########################################################
###################################################################################################################

## clear
rm(list = ls())


## load data (coef of seven methods)
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37")
Coef_ridge <- read.table("ridge\\coef_ridge.csv", header=TRUE, sep = ',')
Coef_lasso <- read.table("lasso\\coef_lasso.csv", header=TRUE, sep = ',')
Coef_Elastic <- read.table("elastic_net\\coef_elastic.csv", header=TRUE, sep = ',')
Coef_Bridge <- read.table("l0.5\\coef_Bridge.csv", header=TRUE, sep = ',')
Coef_l0 <- read.table("L0\\coef_L0name.csv", header=TRUE, sep = ',')
Coef_SCAD <- read.table("SCAD\\coef_SCAD.csv", header=TRUE, sep = ',')
Coef_MCP <- read.table("MCP\\coef_MCP.csv", header=TRUE, sep = ',')


## extract coef of 30 runs
# my_cbind 提取30次的coef -----------------------------------------------------
my_cbind <- function(x){
  x1 <- matrix()
  x1 <- as.character(x[,1])
  for (i in 0:29){
    x1 <- cbind(x1, x[, 3+2*i])
  }
  return(x1)
}


## compute -----------------------------------------------------------------------
coef_ridge1 <- my_cbind(Coef_ridge)
coef_lasso1 <- my_cbind(Coef_lasso)
coef_Elastic1 <- my_cbind(Coef_Elastic) 
coef_Bridge1 <- my_cbind(Coef_Bridge)
coef_l01 <- my_cbind(Coef_l0)  
coef_SCAD1 <- my_cbind(Coef_SCAD) 
coef_MCP1 <- my_cbind(Coef_MCP)  


## extract genes
# my_union_coef 取30次非0系数的并集 --------------------------------------------------
my_union_coef <- function(x){
  x1 <- matrix(data=NA)
  j <- 1
  for (i in 2:dim(x)[1]){
    # i <- 2
    if (x[i,2] !=0 | x[i,3] !=0 | x[i,4] !=0 | x[i,5] !=0
        | x[i,6] !=0 | x[i,7] !=0 | x[i,8] !=0 | x[i,9] !=0
        | x[i,10] !=0 | x[i,11] !=0 | x[i,12] !=0 | x[i,13] !=0
        | x[i,14] !=0 | x[i,15] !=0 | x[i,16] !=0 | x[i,17] !=0
        | x[i,18] !=0 | x[i,19] !=0 | x[i,20] !=0 | x[i,21] !=0
        | x[i,22] !=0 | x[i,23] !=0 | x[i,24] !=0 | x[i,25] !=0
        | x[i,26] !=0 | x[i,27] !=0 | x[i,28] !=0 | x[i,29] !=0
        | x[i,30] !=0 | x[i,31] !=0) {
      x1[j] <- as.character(x[i,1])
      j <- j+1
    }else{
      print("Wrong!")
    }
  }
  return(as.matrix(x1))
}  


# Selected gene -----------------------------------------------------------
coef_ridge <- my_union_coef(coef_ridge1)
coef_lasso <- my_union_coef(coef_lasso1)
coef_Elastic <- my_union_coef(coef_Elastic1)
coef_Bridge <- my_union_coef(coef_Bridge1)
coef_l0 <- my_union_coef(coef_l01)
coef_SCAD <- my_union_coef(coef_SCAD1)
coef_MCP <- my_union_coef(coef_MCP1)


# compute 2-2 Overlap -----------------------------------------------------------

gene_rl <- intersect(coef_ridge, coef_lasso)
gene_re <- intersect(coef_ridge, coef_Elastic)
gene_r1 <- intersect(coef_ridge, coef_l0)
gene_rb <- intersect(coef_ridge, coef_Bridge)
gene_rs <- intersect(coef_ridge, coef_SCAD)
gene_rm <- intersect(coef_ridge, coef_MCP)

gene_le <- intersect(coef_lasso, coef_Elastic)
gene_l1 <- intersect(coef_lasso, coef_l0)
gene_lb <- intersect(coef_lasso, coef_Bridge)
gene_ls <- intersect(coef_lasso, coef_SCAD)
gene_lm <- intersect(coef_lasso, coef_MCP)

gene_e1 <- intersect(coef_Elastic, coef_l0)
gene_eb <- intersect(coef_Elastic, coef_Bridge)
gene_es <- intersect(coef_Elastic, coef_SCAD)
gene_em <- intersect(coef_Elastic, coef_MCP)

gene_1b <- intersect(coef_l0, coef_Bridge)
gene_1s <- intersect(coef_l0, coef_SCAD)
gene_1m <- intersect(coef_l0, coef_MCP)

gene_bs <- intersect(coef_Bridge, coef_SCAD)
gene_bm <- intersect(coef_Bridge, coef_MCP)

gene_sm <- intersect(coef_SCAD, coef_MCP)

# View(gene_rl)  
# View(gene_re)  
# View(gene_r1)  
# View(gene_rb)  
# View(gene_rs)
# View(gene_rm)  
# 
# View(gene_le)  
# View(gene_l1)  
# View(gene_lb)  
# View(gene_ls)  
# View(gene_lm)  
# 
# View(gene_e1)
# View(gene_eb) 
# View(gene_es) 
# View(gene_em) 
# 
# View(gene_1b)
# View(gene_1s)
# View(gene_1m)
# 
# View(gene_bs)
# View(gene_bm)
# 
# View(gene_sm)



# compute number of features  --------------------------------------------------------------------
sum(coef_ridge != 0)
sum(coef_lasso != 0)
sum(coef_Elastic != 0)
sum(coef_l0 != 0)
sum(coef_Bridge != 0)
sum(coef_SCAD != 0)
sum(coef_MCP != 0)


## save data
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37\\result")
# write.csv(coef_ridge, "coef_ridge1.csv")
# write.csv(coef_lasso, "coef_lasso1.csv")
# write.csv(coef_Elastic, "coef_Elastic1.csv")
# write.csv(coef_l0, "coef_l01.csv")
# write.csv(coef_Bridge, "coef_Bridge1.csv")
# write.csv(coef_SCAD, "coef_SCAD1.csv")
# write.csv(coef_MCP, "coef_MCP1.csv")

## save biomarlers
gene <- intersect(intersect(intersect(coef_Bridge, coef_lasso), coef_Elastic), coef_SCAD)
# write.csv(gene, "gene_overlap_scale20.csv")
stargazer(gene) #使用stargazer生成LaTeX语言下的描述性统计表格


##################################    绘制公共的 ROC 曲线    ######################################################
###################################################################################################################

## clear
rm(list = ls())


## load data
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37")
coef_ridge <- read.table("ridge\\pred_ridge0.csv", header=TRUE, sep = ',')
coef_lasso <- read.table("lasso\\pred_lasso0.csv", header=TRUE, sep = ',')
coef_Elastic <- read.table("elastic_net\\pred_elastic0.csv", header=TRUE, sep = ',')
coef_Bridge <- read.table("l0.5\\pred_Bridge0.csv", header=TRUE, sep = ',')
coef_l0 <- read.table("l0\\pred_L00.csv", header=TRUE, sep = ',')
coef_SCAD <- read.table("SCAD\\pred_SCAD0.csv", header=TRUE, sep = ',')
coef_MCP <- read.table("MCP\\pred_MCP0.csv", header=TRUE, sep = ',')


library(pROC)
# pdf(file = "result\\pAUCall7.pdf")
roc_Elastic <- plot.roc( coef_Elastic[,2], coef_Elastic[,3], col="Salmon" )
roc_lasso <- lines.roc( coef_lasso[,2], coef_lasso[,3], col="Aquamarine" )
roc_SCAD <- lines.roc( coef_SCAD[,2], coef_SCAD[,3], col="Magenta" )# VioletRed
roc_Bridge <- lines.roc( coef_Bridge[,2], coef_Bridge[,3], col="Green" )
roc_l0 <- lines.roc( coef_l0[,2], coef_l0[,3], col="Tan" )
roc_MCP <- lines.roc( coef_MCP[,2], coef_MCP[,3], col="Cyan" )
roc_ridge <- lines.roc( coef_ridge[,2], coef_ridge[,3], col="Orchid" )  # Purple
legend("bottomright", legend=c("Elastic Net", "Lasso", "SCAD", "L1/2", "L0", "MCP", "Ridge"), col=c("Salmon", "Aquamarine", "Magenta", "Green", "Tan", "Cyan", "Orchid"), lwd=2)
# dev.off()


## load data
bar <- read.table("result\\bar.csv", header=TRUE, sep = ',')
colnames(bar) <- c("AUC","Penalty")
bar$AUC <- round(bar$AUC,3)      # 保留三位小数
library(ggplot2)
# pdf(file = "result\\bar_fig.pdf")
ggplot(data = bar, aes(x = Penalty, y = AUC, fill = Penalty))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = AUC))
# dev.off()


###################################################################################################################
###################################################################################################################
