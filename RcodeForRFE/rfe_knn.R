# install.packages("lattice")
library(caret)
library(lattice)
library(ggplot2)

setwd('D:\\E\\²©Ê¿\\R_³ÌÐò\\SVM_RFE\\Data\\limma')

eset <- read.table("matrix_DE.txt",header = TRUE,sep = "\t")

Control <- rfeControl(functions = caretFuncs, method = "cv",
                      verbose = FALSE , returnResamp = "final")

trControl1 <- trainControl( method = "cv",
                            classProbs=TRUE,
                            summaryFunction = twoClassSummary)
#seeds = seeds )

disease <- as.factor(c(rep("SPTB",98), rep("Term",228)))
# disease <- as.factor(c(rep(1,98), rep(0,228)))
# View(disease)

# KNN¡ªRFE -----------------------------------------------------------------

rf2 <- rfe(t(eset), disease, sizes =  c(1:694),
           rfeControl = Control, trControl = trControl1,
           metric = "Accuracy",  
           method = "knn")

feature_sele <-rf2$optVariables
# View(feature_sele)
# write.table(feature_sele, file = "ranklist_new_knn.txt", quote=F, sep="\t")


# NN¡ªRFE ------------------------------------------------------------------
# install.packages("neuralnet")
library(neuralnet)
eset1 <- t(eset)
eset1 <- lapply(eset1, factor)

rf2 <- rfe(t(eset), disease, sizes =  c(1:200,210,230,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000),
           rfeControl = Control, trControl = trControl1,
           metric = "ROC",  
           method = "nnet")
feature_sele <-rf2$optVariables
# View(feature_sele)
# write.table(feature_sele, file = "ranklist_new_nnet.txt", quote=F, sep="\t")
