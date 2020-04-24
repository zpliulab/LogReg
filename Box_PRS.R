library(ggplot2)

pred_log = read.table("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37\\DE\\pre59491_73685_4cvlog.txt", header = T, check.names = FALSE)
dim(pred_log)

term0 <- which(pred_log$Lable == 0)
pred_log[term0, 2] <- "Term"
term1 <- which(pred_log$Lable == 1)
pred_log[term1, 2] <- "SPTB"

colnames(pred_log) <- c("preterm_risk_score","lable")


###########################################   画图   ##############################################

pred_log$lable <- as.factor(pred_log$lable)  # dose 为列名，即表达值

e <- ggplot(pred_log, aes(x = lable, y = preterm_risk_score))  # x 有三类，就有3个箱子

# Basic box plot
e + geom_boxplot()
# Box plot with mean points
e + geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 4, color = "blue")
# Change box plot colors by groups
e + geom_boxplot(aes(fill = lable))

# mean points and color
e + geom_boxplot(aes(fill = lable))+
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 4, color = "blue")

## 保存
setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37\\DE")
# jpeg(file = "Box_59491_73685_4cvlog.jpg")
# pdf(file = "Box_59491_73685_4cvlog.pdf")
e + geom_boxplot(aes(fill = lable))+
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 4, color = "blue")
# dev.off()

##########################################################################################

