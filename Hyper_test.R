library(tidyverse)
library(plyr)
library(stargazer) #转换latex 
# 学习 https://blog.csdn.net/linkequa/article/details/86491665 --------------

# fisher.test(matrix(c(2,614,176,17002),nrow = 2), alternative = "less")
# > p-value = 0.05143
# phyper(2,178,17618,616, lower.tail=T)
# > p-value = 0.05145337
# fisher.test(matrix(c(2,614,176,17002),nrow = 2), alternative = "great")
# p-value = 0.9864
# phyper(2-1,178,17618,616, lower.tail=F)
# Pvalue =  0.986356
# lower.tail logical; if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
# 仅在 lower.tail=F ,做enrichment时取q-1. (相比仅以q计算, 这个P值实际要略微大点) 
# 但lower.tail=T时, 如果是做depletion(选择接受假设), 此时取q-1, 结果实际上少算了一个值. 
# (实际要算的是P[X ≤ x],但如果q取q-1, 结果就成了P[X ＜ x]的情况, 将当前观测的P值给排除了,所得的P值要偏小)


# 计算超几何分布检验 ---------------------------------------------------------------

## 初始化
## 例子
# N <- 20000
# M <- 2005  
# n <- 805
# k <- 265
## 真实
N <- 359  # 总数
M <- 359  # 方法一
n <- 67   # 方法二
k <- 67   # overlap


# my_phyper 函数 ---------------------------------------------------------------

my_phyper <- function(N,M,n,k){
  q <- k
  m <- M
  p <- N-M
  l <- n
  p_value <- phyper(q-1, m, p, l, lower.tail=F)
  return(p_value)
}


# 结果 ----------------------------------------------------------------------

p_rl <- my_phyper(359,359,67,67)
p_re <- my_phyper(359,359,78,78)
p_r1 <- my_phyper(359,359,31,31)
p_rb <- my_phyper(359,359,25,25)
p_rs <- my_phyper(359,359,49,49)
p_rm <- my_phyper(359,359,27,27)

p_le <- my_phyper(359,67,78,67)
p_l1 <- my_phyper(359,67,31,28)
p_lb <- my_phyper(359,67,25,21)
p_ls <- my_phyper(359,67,49,48)
p_lm <- my_phyper(359,67,27,27)

p_e1 <- my_phyper(359,78,31,29)
p_eb <- my_phyper(359,78,25,21)
p_es <- my_phyper(359,78,49,48)
p_em <- my_phyper(359,78,27,27)

p_1b <- my_phyper(359,31,25,13)
p_1s <- my_phyper(359,31,49,19)
p_1m <- my_phyper(359,31,27,24)

p_bs <- my_phyper(359,25,49,20)
p_bm <- my_phyper(359,25,27,13)

p_sm <- my_phyper(359,49,27,27)


# 输出结果 --------------------------------------------------------------------


result <- matrix(0, 6, 6)
colnames(result) <- c( "lasso", "elatic net", "L1/2", "L0", "SCAD", "MCP")
rownames(result) <- c("ridge", "lasso", "elatic net", "L1/2", "L0", "SCAD")

i <- 1
result[i,1] <- signif(p_rl, digits = 3)
result[i,2] <- signif(p_re, digits = 3)
result[i,3] <- signif(p_r1, digits = 3)
result[i,4] <- signif(p_rb, digits = 3)
result[i,5] <- signif(p_rs, digits = 3)
result[i,6] <- signif(p_rm, digits = 3)

result[i+1,2] <- signif(p_le, digits = 3)
result[i+1,3] <- signif(p_l1, digits = 3)
result[i+1,4] <- signif(p_lb, digits = 3)
result[i+1,5] <- signif(p_ls, digits = 3)
result[i+1,6] <- signif(p_lm, digits = 3)

result[i+2,3] <- signif(p_e1, digits = 3) # format(p_e1, scientific = T)
result[i+2,4] <- signif(p_eb, digits = 3) 
result[i+2,5] <- signif(p_es, digits = 3)
result[i+2,6] <- signif(p_em, digits = 3)

result[i+3,4] <- signif(p_1b, digits = 3)
result[i+3,5] <- signif(p_1s, digits = 3)
result[i+3,6] <- signif(p_1m, digits = 3)

result[i+4,5] <- signif(p_bs, digits = 3)
result[i+4,6] <- signif(p_bm, digits = 3)

result[i+5,6] <- signif(p_sm, digits = 3)

setwd("D:\\E\\博士\\R_程序\\GSE59491_15\\Data37\\result")
result %>% write.csv(file = "phyper.csv") 

# 表格 ----------------------------------------------------------------------
stargazer(result) 

# 百分比填色 -------------------------------------------------------------------
result <- as.matrix(result)

p <- matrix(data = NA, nrow = 21, ncol = 1,byrow = F, dimnames = NULL)

for (j in 1:6){
  for (i in j:6){
    p[6*(j-1)+i] <- result[j,i]
  }
}

View(p)


# 删除缺失值 -------------------------------------------------------------------
p <- data.frame(p)
p_new <- na.omit(p)
View(p_new)



# normalized --------------------------------------------------------------
p_matrix <- as.matrix(p_new)[7:21]
View(p_matrix)
p_sort <- sort(p_matrix)
View(p_sort)
p_scale <- scale(p_matrix)
View(p_scale)


# sort --------------------------------------------------------------------
sort1 <- matrix(data = NA, nrow = 15, ncol = 1,byrow = F, dimnames = NULL)

for (k in 1:15){
  sort1[k] <- p_sort[k]/p_sort[15]
}
View(sort1)
