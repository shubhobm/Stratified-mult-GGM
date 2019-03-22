rm(list=ls())
setwd("C:/Study/Stratified-mult-GGM/Real data/")
library(data.table)
final_model = readRDS("final_model.rds")
# final_model = readRDS("model_list.rds")[[2]]

# load data
data = readRDS("C:/Study/Stratified-mult-GGM/Real data/processed_data.rds")
X.names = colnames(data$X.list1[[1]])
Y.names = colnames(data$Y.list1[[1]])

# non-zero values in B
K = 2
groups = c("ER+","ER-")
B.values = vector("list",K)
for(k in 1:K){
  non.zero.B = which(final_model$B.refit[,,k] != 0, arr.ind=T)
  B.df = data.table(SampleGroup = groups[k],
                    mRNA = X.names[non.zero.B[,1]],
                    RNAseq = Y.names[non.zero.B[,2]])
  invisible(B.df[, Value := 0])
  for(i in 1:nrow(B.df)){
    invisible(B.df[i, Value := final_model$B.refit[non.zero.B[i,1],non.zero.B[i,2],k]])
  }
  B.values[[k]] = B.df[order(abs(Value), decreasing=T)]
}

B.values
fwrite(rbindlist(B.values), file="B_values.csv", sep=",")

# non-zero values in Theta
K = 2
Th.values = vector("list",K)
for(k in 1:K){
  non.zero.Th = which(final_model$Theta_refit$Theta[[k]] != 0, arr.ind=T)
  non.zero.Th = non.zero.Th[non.zero.Th[,1] < non.zero.Th[,2],] # just take lower triangle
  Th.df = data.table(SampleGroup = groups[k],
                     RNASeq1 = X.names[non.zero.Th[,1]],
                     RNAseq2 = Y.names[non.zero.Th[,2]])
  invisible(Th.df[, Value := 0])
  for(i in 1:nrow(Th.df)){
    invisible(Th.df[i, Value := final_model$Theta_refit$Theta[[k]][non.zero.Th[i,1],non.zero.Th[i,2]]])
  }
  # invisible(Th.df[, Value := round(Value, 5)])
  Th.values[[k]] = Th.df[order(abs(Value), decreasing=T)]
}

Th.values
fwrite(rbindlist(Th.values), file="Theta_values.csv", sep=",")


## Tune JSEM model for X ***********************************************
# **********************************************************************
X.indices = X.layer$indices
Zeta.groups = X.layer$groups
gamma = sqrt(log(p)/n) * seq(1, 0.4, -0.1)
bic.jsem <- sel.lambda.jsem(do.call(rbind, X.list), do.call(rbind, X.list),
                            unlist(X.indices), unlist(X.indices),
                            Zeta.groups,lambda=gamma)
gamma.min = gamma[which.min(bic.jsem$BIC)]
jsem.model = JSEM(do.call(rbind, X.list), unlist(X.indices),
                  Zeta.groups, lambda=gamma.min)
Zeta_new.array = array(0, c(p,p,K))
for(k in 1:K){
  Zeta_new.array[,,k] = jsem.model$Theta[[k]]
}

## Get debiased estimates **********************************************
# **********************************************************************
B.hat.array = final_model$B.refit
C.hat.array = B.hat.array
M = matrix(0,p,K)
for(k in 1:K){
  X.k = X.list[[k]]
  E.k = Y.list[[k]] - X.k %*% B.hat.array[,,k]
  for(i in 1:p){
    R.ik = X.k[,i] - X.k[,-i] %*% Zeta_new.array[i,-i,k]
    t.ik = as.numeric(t(R.ik) %*% X.k[,i]/n)
    C.hat.array[i,,k] = B.hat.array[i,,k] + t(R.ik) %*% E.k/n/t.ik
    M[i,k] = sqrt(n)*t.ik/sqrt(sum(R.ik^2/n))
  }
}

## Get eigenvectors and eigenvalues of precision matrices
Theta1 = solve(jmmle.model$Theta_refit$Omega[[1]])
Theta2 = solve(jmmle.model$Theta_refit$Omega[[2]])

## Global test statistics for i-th X-variable
D = rep(0,p)
d = matrix(0,p,q)
for(i in 1:p){
  Pooled.Cov.i = Theta1/M[i,1]^2 + Theta2/M[i,2]^2
  Diff.i = C.hat.array[i,,1] - C.hat.array[i,,2]
  ## pairwise test statistics
  d[i,] = (Diff.i)^2/diag(Pooled.Cov.i)
  ## overall test statistic
  D[i] = t(Diff.i) %*% solve(Pooled.Cov.i) %*% Diff.i
}

## determine threshold for i-th test
alpha = .2
d.ind.mat = matrix(0,p,q)
tau = rep(20,p)
which.i.reject = which(D > qchisq(.95, q))
for(i in which.i.reject){
  tau.vec = seq(0, 20, length.out=1e2)
  thres.vec = as.numeric(lapply(tau.vec, function(x) alpha/q * max(sum(d[i,]>x),1)))
  which.less = which((1 - pchisq(tau.vec,1)) <= thres.vec)
  if(length(which.less)>0){
    tau[i] = tau.vec[which.less[1]] # set tau as minimizer only if there is at least one tau entry less
  }
  # tau[i] = tau.vec[which.min(abs(1 - pchisq(tau.vec,1) - thres.vec))]
  
  d.ind.mat[i,] = as.numeric(d[i,]>tau[i])
}
