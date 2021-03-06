rm(list=ls())
setwd('d:/Study/My projects/Stratified-mult-GGM/Codes')
source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

library(glasso)
library(parallel)

##### Generate data
group = matrix(c(1, 2), nrow=2, ncol=2, byrow=T)           # grouping pattern
subnetSize.E = c(30, 30)
subnetSize.X = c(15, 15)    # subnet size
n = 100
p = sum(subnetSize.X)
q = sum(subnetSize.E)
K = 2

set.seed(12182017)
X.layer = GenerateLayer(n, subnetSize.X, group, D=1)
E.layer = GenerateLayer(n, subnetSize.E, group, D=1)

## generate group structure for coef array
B0.group.array = array(0, c(p,q,K))
g = 1
for(i in 1:p){
  for(j in 1:q){
    B0.group.array[i,j,] = g
    g = g+1
  }
}
B0.array = CoefArray2(B0.group.array[,,1], D=1)
# B0.array = CoefArray(B0.group.array)
B0.array = B0.array[[1]]
Diff.mat = B0.array[,,1] - B0.array[,,2]
Theta0.array = array(0, c(q,q,K))
for(k in 1:K){
  Theta0.array[,,k] = with(E.layer,
                           diag(diag(Omega[[k]])^(-0.5)) %*% Omega[[k]] %*% diag(diag(Omega[[k]])^(-0.5)))
}

## make Y-layer
Y.layer = E.layer
for(k in 1:K){
  Y.layer$data[[k]] = X.layer$data[[k]] %*% B0.array[,,k] + E.layer$data[[k]]
}

##### Given: X.list, Y.list, B.groups, Theta.groups
Y.list = lapply(Y.layer$data, as.matrix)
Y.indices = Y.layer$indices
Theta.groups = Y.layer$groups
X.list = lapply(X.layer$data, as.matrix)

Theta.group.array = array(0, c(q,q,K))
for(j in 1:q){
  Theta.group.array[j,-j,] = Y.layer$groups[[j]]
}

## Obtain JMMLE fit ****************************************************
# **********************************************************************
## tune JMMLE model
lambda.vec = sqrt(log(p)/n) * seq(1, 0.4, -0.1)
model.list = vector("list", length(lambda.vec))
nlambda = length(lambda.vec)

## get all models
loopfun1 = function(m){
  jmmle.1step(Y.list, Y.indices, X.list, B.group.array=B0.group.array, Theta.groups=Theta.groups,
       lambda = lambda.vec[m],
       gamma = sqrt(log(q)/n) * seq(1, 0.4, -0.1),
       init.option=1, tol=1e-3)
}
system.time(
  model.list <- lapply(1:nlambda, loopfun1)
)

## calculate HBIC
hbic.vec = rep(0, nlambda)
for(m in 1:nlambda){
  jmle.model = model.list[[m]]
  SSE.vec = rep(0,K)
  hbic.pen.vec = rep(0,K)
  
  for(k in 1:K){
    nk = nrow(Y.list[[k]])
    Theta.k = jmle.model$Theta_refit$Theta[[k]]
    for(j in 1:q)
    {
      Theta.k[j,j] = 0
    }
    SSE.vec[k] = sum(diag(crossprod((Y.list[[k]] - X.list[[k]] %*%
                                       jmle.model$B.refit[,,k]) %*% (diag(1,q) - Theta.k))))/nk
    hbic.pen.vec[k] = log(log(nk))*log(q*(q-1)/2)/nk * sum(Theta.k != 0)/2 +
      log(log(nk))*log(p*q)/nk * sum(jmle.model$B.refit[,,k] != 0)
  }
  hbic.vec[m] = sum(SSE.vec) + sum(hbic.pen.vec)
}

## select best model
jmmle.model = model.list[[which.min(hbic.vec)]]

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
B.hat.array = jmmle.model$B.refit
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
e1 = eigen(jmmle.model$Theta_refit$Omega[[1]])
P1 = e1$vectors
L1 = e1$values
e2 = eigen(jmmle.model$Theta_refit$Omega[[2]])
P2 = e2$vectors
L2 = e2$values

## Global test statistics for i-th X-variable
Omega1.sqrt = P1 %*% diag(sqrt(L1)) %*% t(P1)
Omega2.sqrt = P2 %*% diag(sqrt(L2)) %*% t(P2)
# O1.c1 = Omega1.sqrt%*%(C.hat.array[,,1])
# O2.c2 = Omega2.sqrt%*%(C.hat.array[,,2])
# 
# D = rep(0,p)
# for(i in 1:p){
#   D[i] = sum((M[i,1]*O1.c1[i,] - M[i,2]*O2.c2[i,])^2)
# }
# 
D = rep(0,p)
for(i in 1:p){
  D[i] = sum((M[i,1]*Omega1.sqrt%*%C.hat.array[i,,1] - M[i,2]*Omega2.sqrt%*%C.hat.array[i,,2])^2)
}

which(D > qchisq(.95, 2*q)) # indices where global test is rejected

## pairwise test statistics
Omega1i.sqrt = P1 %*% diag(1/sqrt(L1)) %*% t(P1)
Omega2i.sqrt = P2 %*% diag(1/sqrt(L2)) %*% t(P2)
d = matrix(0,p,q)
for(i in 1:p){
  for(j in 1:q){
    d[i,j] = ((C.hat.array[i,j,1] - C.hat.array[i,j,2])/
                (Omega1i.sqrt[j,j]/M[i,1] + Omega2i.sqrt[j,j]/M[i,2]))^2
  }
}

## determine threshold for i-th test
alpha = .05
d.ind.mat = matrix(0,p,q)
tau = rep(NA,p)
for(i in which(D > qchisq(.95, 2*q))){
  tau.vec = seq(0, 20, length.out=1e2)
  thres.vec = lapply(tau.vec, function(x) alpha/q * max(sum(d[i,]>x),1))
  thres.vec = as.numeric(thres.vec)
  tau[i] = tau.vec[which.min(abs(1 - pchisq(tau.vec,1) - thres.vec))]
  d.ind.mat[i,] = as.numeric(d[i,]>tau[i])
}

# tau.vec = seq(0, 20, length.out=1e2)
# thres.vec = lapply(tau.vec, function(x) alpha/(p*q) * max(sum(d>x),1))
# thres.vec = as.numeric(thres.vec)
# tau = tau.vec[which.min(abs(1 - pchisq(tau.vec,2) - thres.vec))]
# d.ind.mat = matrix(as.numeric(d>tau), nrow=p, ncol=q, byrow=F)

pow = sum(d.ind.mat == 1 & Diff.mat != 0, na.rm=T)/sum(Diff.mat != 0)
size = 1 - sum(d.ind.mat == 0 & Diff.mat == 0, na.rm=T)/sum(Diff.mat == 0)
FDP = sum(d.ind.mat == 1 & Diff.mat == 0, na.rm=T)/sum(d.ind.mat == 1, na.rm=T)
c(pow,size,FDP)

Theta_new.array = array(0, c(q,q,K))
for(k in 1:K){
  Theta_new.array[,,k] = jmmle.model$Theta_refit$Theta[[k]]
}
round(c(sum(B0.array != 0 & jmmle.model$B.refit != 0)/sum(B0.array != 0),
        sum(B0.array == 0 & jmmle.model$B.refit == 0)/sum(B0.array == 0),
        sqrt(sum((B0.array - jmmle.model$B.refit)^2)),
        sum(Theta0.array != 0 & Theta_new.array != 0)/sum(Theta0.array != 0),
        sum(Theta0.array == 0 & Theta_new.array == 0)/sum(Theta0.array == 0),
        sqrt(sum((Theta0.array - Theta_new.array)^2))), 3)

