#path = ""
#setwd(path)
mydata <- read.csv('./seizure.csv') 
source("../functions.R")
set.seed(110)
selection <- names(mydata)
group1 <- subset(mydata, y == 1, select = selection[2:179])
group2 <- subset(mydata, y == 2, select = selection[2:179])
group3 <- subset(mydata, y == 3, select = selection[2:179])
group4 <- subset(mydata, y == 4, select = selection[2:179])
group5 <- subset(mydata, y == 5, select = selection[2:179])
#split each groups into 4 subgroups and calculate autocorrlation
`%notin%` <- Negate(`%in%`)
GenCorMatirx <- function(group, k = 4){
  num <- nrow(group)
  rawid <- c(1:num)
  id <- matrix(NA, ncol = k, nrow = num/k)
  subgroup <- list()
  R <- list()
  for (i in 1:k){
    id[,i] <- sample(rawid, num/k)
    for (j in 1:i){
      rawid <- rawid[! rawid %in% id[,i]]
    }
    subgroup[[i]] <- group[id[,i],]
    X <- as.matrix(subgroup[[i]])
    Sigma <- t(X)%*%X/(num/k)
    R[[i]] <- abs(diag((diag(Sigma)^(-1/2)))%*%Sigma%*%diag((diag(Sigma)^(-1/2))))
    diag(R[[i]]) <- 0
  }
  return(R)
}
group <- list()
corVal <- c()
k <- 4
for (i in 1:5){
  group[[i]] <- subset(mydata, y == i, select = selection[2:179])
  CorM <- GenCorMatirx(group[[i]], k = k)
  for (j in 1:k){
    upper.tri.ind <- upper.tri(CorM[[j]])
    corVal <- c(corVal, CorM[[j]][upper.tri.ind])
  }
}
t <- quantile(corVal, prob = 0.9)
GenAdjMatrix <- function(group_target, k = 4, t){
  CorM <- GenCorMatirx(group_target, k = k)
  nCorM <- dim(CorM[[1]])[1]
  for (j in 1:k){
    upper.tri.ind <- upper.tri(CorM[[j]])
    c.upper <- CorM[[j]][upper.tri.ind]
    c.upper[which(c.upper >= t)] <- 1
    c.upper[which(c.upper < t)] <- 0
    CorM[[j]] <- matrix(0, nCorM, nCorM)
    CorM[[j]][upper.tri.ind] <- c.upper
    CorM[[j]] <-CorM[[j]] + t(CorM[[j]])
  }
  return(CorM)
}
A1.list <- GenAdjMatrix(group_target = group1, k = 4, t = t)
A2.list <- GenAdjMatrix(group_target = group2, k = 4, t = t)
A3.list <- GenAdjMatrix(group_target = group3, k = 4, t = t)
A4.list <- GenAdjMatrix(group_target = group4, k = 4, t = t)
A5.list <- GenAdjMatrix(group_target = group5, k = 4, t = t)

# hypothesis testing starts here
GenPvalue <- function(list1, list2, d, rep = 100){
  A_bar_1 <- (list1[[1]]+list1[[2]]+list1[[3]]+list1[[4]])/4
  A_bar_2 <- (list2[[1]]+list2[[2]]+list2[[3]]+list2[[4]])/4
  rho <- Spearman(A_bar_1, A_bar_2, num = d, res = FALSE)
  rho.hat <- c()
  combine.list <- c(list1, list2)
  for(i in 1:rep){
    relabel <- sample(1:8, 4, replace = FALSE)
    A.perm <- matrix(0, nrow = 178, ncol = 178)
    B.perm <- matrix(0, nrow = 178, ncol = 178)
    for (j in 1:8){
      if (j %in% relabel){
        A.perm <- A.perm + combine.list[[j]]
      }else{
        B.perm <- B.perm + combine.list[[j]]
      }
    }
    A.perm_bar <- A.perm/4
    B.perm_bar <- B.perm/4
    rst <- Spearman(A.perm_bar, B.perm_bar, num = d, res = FALSE)
    rho.hat <- c(rho.hat, rst)
  }
  pvalue <- sum(rho.hat < rho)/length(rho.hat)
  return(pvalue)
}

rowPvalue <- function(A, d){
  return(c(GenPvalue(A, A1.list, d = d),
           GenPvalue(A, A2.list, d = d),
           GenPvalue(A, A3.list, d = d),
           GenPvalue(A, A4.list, d = d),
           GenPvalue(A, A5.list, d = d))
         )
}

tbl <- matrix(NA, nrow = 7, ncol = 5)
for (d in seq(2, 8)){
  tbl[d - 1, ] <- rowPvalue(A1.list, d)
}
tbl
