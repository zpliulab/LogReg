## ���ڶ� series_marix ���ӱ�ǩ������ݣ����� T ����
## ȷ�� T ����� p ֵ������ p/pdf <0.05 ��gene������

# ����������0-1����ǩ�ڵ�1�� ----------------------------------------------------------------
## ȫΪ 1 ��
my_test1 <- function(x){
  arry <- x[1,]
  test1 <- x[, arry == 1]
  return(test1)
}
## ȫΪ 0 ��
my_test0 <- function(x){
  arry <- x[1,]
  test0 <- x[, arry == 0]
  return(test0)
}


# ����������0-1����ǩ�ڵ�1�� ---------------------------------------------------------
## ȫΪ 1 ��
my_test_1 <- function(x){
  arry <- x[,1]
  test1 <- x[arry == 1, ]
  return(test1)
}
## ȫΪ 0 ��
my_test_0 <- function(x){
  arry <- x[,1]
  test0 <- x[arry == 0, ]
  return(test0)
}



# T���飬��������ǩ�ڵ�1�� -----------------------------------------------------------
my_p <- function(x){
  p <- matrix(data=NA, nrow = dim(x)[2]-1, ncol = 1, byrow = FALSE, dimnames=list(c(colnames(x[,-1])),c("pvalue")))
  for(i in 2:ncol(x)){
    p[i-1] <- t.test(test0[,i], test1[, i])$p.value
  }
  return(p)
}



# BHУ������������ǩ�ڵ�1�� ----------------------------------------------------------
my_BH_fdr <- function(x){
  BH_fdr <- as.matrix(p.adjust(x, "BH"))
  p_BH_fdr <- matrix(data=BH_fdr, nrow = dim(x)[1], ncol = 1, byrow = FALSE, dimnames=list(c(as.matrix(rownames(x))),c("p_BH_fdr")))
  return(p_BH_fdr)
}


# FDRУ������������ǩ�ڵ�1�� -------------------------------------------------------
library(fdrtool)
my_fdr <- function(x,y){
  x <- as.vector(x)
  fdr = fdrtool(x, statistic="pvalue")
  p_fdr <- as.matrix(fdr$qval)                       # estimated Fdr values
  y <- as.matrix(y)
  p_FDR <- matrix(data=p_fdr, nrow = dim(y)[1], ncol = 1, byrow = FALSE, dimnames=list(c(as.matrix(rownames(y))),c("p_fdr")))
  return(p_FDR)
}


# BHУ����ȡ ------------------------------------------------------------------
my_data <- function(x){
  x_hat <- t(x)[-1,]  
  x_hat_Right <- x_hat[Right, ]
  dim(x_hat_Right)    # 987  64
  y <- x[,1]
  y <- as.matrix(y)
  colnames(y) <- 'Lable'
  x_lab <- cbind(y, t(x_hat_Right))
  return(x_lab)
}  


# my_scale (��ǩ�ڵ�1��)-------------------------------------------------------------------
my_scale <- function(x){
  x1 <- cbind(t(x[1,]), scale(t(x[-1,])))
  x2 <- t(x1)
  return(x2)
}

# my_p_data ����pֵ����data��ȡ�������д�ɺ��� ---------------------------------------------------
my_p_data <- function(p, data ){
  Right <- which(p <= 0.05)
  # Right <- which(p <= 0.01)
  x_hat <- t(data)[-1,]          
  x_hat_Right <- x_hat[Right, ]
  y <- data[,1]
  y <- as.matrix(y)
  colnames(y) <- 'Lable'
  x_lab <- cbind(y, t(x_hat_Right))
  return(x_lab)
}
