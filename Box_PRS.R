## clear
rm(list = ls())

## package
library(ggplot2)
library(ggpubr)
library(Rcpp)

## load data
setwd("D:\\E\\²©Ê¿\\R_³ÌÐò\\GSE59491_15\\Data37\\DE")
pred_log = read.table("pre59491_46510_7cvlog_zf_new.txt", header = T, check.names = FALSE)
dim(pred_log)
term0 <- which(pred_log$Lable == 0)
pred_log[term0, 2] <- "Term"
term1 <- which(pred_log$Lable == 1)
pred_log[term1, 2] <- "SPTB"
colnames(pred_log) <- c("preterm_risk_score","label")

## plot
pred_log$lable <- as.factor(pred_log$label)   
ee <- ggboxplot(pred_log, x = "label", y = "preterm_risk_score", color = "label",
                palette = "npg", add = "jitter") 
my_comparisons <- list( c("SPTB", "Term") )
# pdf(file = "Box_59491_73685_17cvlog_zf_new.pdf",width = 4,height = 4)  # 2020.6.29 later
ee + stat_compare_means(comparisons = my_comparisons) + 
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "blue") 
# dev.off()
