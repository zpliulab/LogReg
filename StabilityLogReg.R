################# 2022.1.24 add stability index of our feature selection subset 


## clear 
rm(list = ls())


## package
library(stabm)
library(ggdendro)


## load data (genes selected by dofferent methods)
setwd("D:\\E\\²©Ê¿\\R_³ÌÐò\\GSE59491_15")
Ridge <- as.matrix(read.csv("Data37\\result\\coef_ridge1.csv", header=TRUE)[,2])
Lasso <- as.matrix(read.csv("Data37\\result\\coef_lasso1.csv", header=TRUE)[,2])
Enet <- as.matrix(read.csv("Data37\\result\\coef_elastic1.csv", header=TRUE)[,2])
L0 <- as.matrix(read.csv("Data37\\result\\coef_l01new.csv", header=TRUE)[,2])
Bridge <- as.matrix(read.csv("Data37\\result\\coef_Bridge1.csv", header=TRUE)[,2])
SCAD <- as.matrix(read.csv("Data37\\result\\coef_SCAD1.csv", header=TRUE)[,2])
MCP <- as.matrix(read.csv("Data37\\result\\coef_MCP1.csv", header=TRUE)[,2])


## some gene names are "." not "-"
library(stringr)
Ridge <- str_replace_all(Ridge, "[.]", "-") 
Lasso <- str_replace_all(Lasso, "[.]", "-") 
Enet <- str_replace_all(Enet, "[.]", "-") 
L0 <- str_replace_all(L0, "[.]", "-") 
Bridge <- str_replace_all(Bridge, "[.]", "-") 
SCAD <- str_replace_all(SCAD, "[.]", "-") 
MCP <- str_replace_all(MCP, "[.]", "-") 

## add labels
features <- list(Ridge, Lasso, Enet, L0, Bridge, SCAD, MCP)
listname <- alist(Ridge, Lasso, Enet, L0, Bridge, SCAD, MCP)
names(features) <- c("Ridge", "Lasso", "Enet", "L0", "Bridge", "SCAD", "MCP")


listwithname <- function(features) {
  names(features) <- eval(substitute(listname))
  return(features)
}

listwithname(features)
all.feats <- unique(unlist(features, use.names = FALSE))
p <- length(all.feats)


# Selected in each method
Reduce(intersect, features)
# Sorted selection frequency 
sort(table(unlist(features)), decreasing = TRUE)
## The selection frequency can be visualized with the plotFeatures() function:
plotFeatures(features)


################# all possiable Combination including 4 subsets  ###############
n <- 7
k <- 4
choose(n, k)
combn(n,k)
A <- as.matrix(combn(n,k))
colnames(A) <- rep(1:35, 1)
A[,1]


## C^4_7 £º stability
var <- NULL
index <- NULL
for (i in 1:35) {
  p <- length(unique(unlist(features[A[,i]], use.names = FALSE)))
  stab <- stabilityHamming(features[A[,i]], p)
  var <- c(var, i)
  index <- c(index, stab) 
}
plot(var, index, type='l', xaxt="n") 
axis(side=1, at=1:35, tck = -0.02, 
     labels = var)

## result
A[, which(index == max(index))]
## stability max : 2 3 6 7
## our selected subset : 2 3 5 6, the 
b <- as.vector(c(2, 3, 5, 6))    # b is 24
order(index,decreasing = T)    # 24 orders 4th
sort(index,decreasing = T)    # 24 orders 4th
var0 <- order(index,decreasing = T) 
index0 <- index[var0]
names(index0) <- var0
barplot(index0, col=c(rep("gray",3),"red",rep("gray",31)), 
                      xlab="Combination",ylab="Stability")

## the red bar is the stability index of our selected feature subset (2 3 5 6, i.e., Lasso, Elastic net, L1/2, SCAD)
## 