# install.packages("lattice")
library(caret)
library(lattice)
library(ggplot2)

setwd('D:\\E\\²©Ê¿\\R_³ÌÐò\\SVM_RFE\\Data\\limma')

eset <- read.table("matrix_DE_T.txt",header = TRUE,sep = "\t")

Control <- rfeControl(functions = caretFuncs, method = "cv",
                      verbose = FALSE , returnResamp = "final")

trControl1 <- trainControl( method = "cv",
                            classProbs=TRUE,
                            summaryFunction = twoClassSummary)
#seeds = seeds )

disease <- as.factor(c(rep("SPTB",98), rep("Term",228)))
# disease <- as.factor(c(rep(1,98), rep(0,228)))
# class(disease)

# KNN¡ªRFE -----------------------------------------------------------------

rf2 <- rfe(eset, disease, sizes =  c(1:200,210,230,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000),
           rfeControl = Control, trControl = trControl1,
           metric = "Accuracy",  
           method = "knn")

feature_sele <-rf2$optVariables
# View(feature_sele)
# write.table(feature_sele, file = "ranklist_knn.txt", quote=F, sep="\t")


# NN¡ªRFE ------------------------------------------------------------------

rf2 <- rfe(eset, disease, sizes = c(1:200,210,230,250,300,400,500,600,700,800,900), 
           metric = "Accuracy", 
           method="nnet",
           rfeControl = rfeControl(functions = caretFuncs) 
)



# eset1 <- lapply(eset1, factor)

rf2 <- rfe(eset, disease, sizes =  c(886),
           rfeControl = Control, trControl = trControl1, method = "nnet",
           tuneGrid = expand.grid(size = c(8), decay = c(0.1)),
           maxit = 30, MaxNWts = 100000
           )
feature_sele <-rf2$optVariables
# View(feature_sele)
# write.table(feature_sele, file = "ranklist_nnet.txt", quote=F, sep="\t")




# # varimp
# 
# control <- trainControl( method = "cv",
#                             classProbs=TRUE,
#                             summaryFunction = twoClassSummary)
# # control <- trainControl(method="repeatedcv", number=10, repeats=3)
# model <- train(disease~., data=eset, method="nnet",trControl=control)
# imp<-varImp(model)
# plot(imp)
# 
# # rfe
# control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# results <- rfe(train[,!colnames %in% c("y")],train$y, sizes=c(1:(ncol(train)-1), rfeControl=control)
#                print(results)
#                # features selected 
#                predictors(results)
#                plot(results, type=c("g", "o"))
               