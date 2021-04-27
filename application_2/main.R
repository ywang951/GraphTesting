#path = ""
mydata <- read.csv('./seizure.csv')
source("../functions.R")

selection <- names(mydata)
group1 <- subset(mydata, y == 1, select = selection[2:179])
group2 <- subset(mydata, y == 2, select = selection[2:179])
group3 <- subset(mydata, y == 3, select = selection[2:179])
group4 <- subset(mydata, y == 4, select = selection[2:179])
group5 <- subset(mydata, y == 5, select = selection[2:179])
#split each groups into 4 subgroups and calculate autocorrlation
`%notin%` <- Negate(`%in%`)
GenCorMatirx <- function(group, k = 4) {
  num <- nrow(group)
  rawid <- c(1:num)
  id <- matrix(NA, ncol = k, nrow = num/k)
  subgroup <- list()
  R <- list()
  for (i in 1:k) {
    id[, i] <- sample(rawid, num/k)
    for (j in 1:i){
      rawid <- rawid[! rawid %in% id[, i]]
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
set.seed(9)
for (i in 1:5){
  group[[i]] <- subset(mydata, y == i, select = selection[2:179])
  CorM <- GenCorMatirx(group[[i]], k = k)
  for (j in 1:k){
    upper.tri.ind <- upper.tri(CorM[[j]])
    corVal <- c(corVal, CorM[[j]][upper.tri.ind])
  }
}

t <- quantile(corVal, prob = 0.75)
GenAdjMatrix <- function(group_target, k = 4, t) {
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
GenPvalue <- function(list1, list2, d, rep = 100, n0) {
  A_bar_1 <- (list1[[1]]+list1[[2]]+list1[[3]]+list1[[4]])/4
  A_bar_2 <- (list2[[1]]+list2[[2]]+list2[[3]]+list2[[4]])/4
  n <- dim(A_bar_1)[1]
  rho <- Spearman(A_bar_1, A_bar_2, num = d, res = FALSE)
  rhovec1 <- rhovec2 <- c()
  for (i in 1:rep) {
    A.perm <- B.perm <- matrix(0, n, n)
    for (k in 1:4) {
      A.temp <- list1[[k]]
      rownames(A.temp) <- colnames(A.temp) <- seq(n)
      a <- seq(n)
      s <- sample(seq(n), n0)
      a[sort(s)] <- s
      A.perm <- A.perm + permutation(A.temp, a)

      B.temp <- list2[[k]]
      rownames(B.temp) <- colnames(B.temp) <- seq(n)
      B.perm <- B.perm + permutation(B.temp, a)
    }
    rho1 <- Spearman(A_bar_1, A.perm/4, num = d, res = FALSE)
    rho2 <- Spearman(A_bar_2, B.perm/4, num = d, res = FALSE)
    rhovec1 <- c(rhovec1, rho1)
    rhovec2 <- c(rhovec2, rho2)
  }
  pval1 <- sum(rhovec1 < rho) / length(rhovec1)
  pval2 <- sum(rhovec2 < rho) / length(rhovec2)
  pvalue <- max(pval1, pval2)
  return(pvalue)
}

rowPvalue <- function(A, d, n0, rep) {
  return(c(GenPvalue(A, A1.list, d = d, rep, n0),
           GenPvalue(A, A2.list, d = d, rep, n0),
           GenPvalue(A, A3.list, d = d, rep, n0),
           GenPvalue(A, A4.list, d = d, rep, n0),
           GenPvalue(A, A5.list, d = d, rep, n0))
  )
}

tbl <- matrix(NA, nrow = 6, ncol = 5)
for (d in seq(2, 7)) {
  tbl[d - 1, ] <- rowPvalue(A1.list, d, n0 = 2, rep = 1000)
  print(c(d, tbl[d - 1, ]))
}
tbl
