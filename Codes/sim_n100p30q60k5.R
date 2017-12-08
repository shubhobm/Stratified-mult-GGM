rm(list=ls())
# setwd('D:/Study/My projects/Stratified-mult-GGM/Codes')
source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

library(glasso)
library(parallel)

##### Generate data
group = rbind(
  c(1, 2),
  c(1, 4),
  c(3, 2),
  c(3, 4),
  c(5, 2)
)                         # grouping pattern
subnetSize.E = c(30, 30)
subnetSize.X = c(15, 15)    # subnet size
n = 100
p = sum(subnetSize.X)
q = sum(subnetSize.E)
K = 5
set.seed(11192017)

## looping function
loopfun = function(rep){
  set.seed(rep*11222017)
  X.layer = GenerateLayer(n, subnetSize.X, group)
  E.layer = GenerateLayer(n, subnetSize.E, group)
  
  ## generate group structure for coef array
  B0.group.array = array(0, c(p,q,K))
  g = 1
  for(i in 1:p){
    for(j in 1:q){
      B0.group.array[i,j,] = g
      g = g+1
    }
  }
  B0.array = CoefArray(B0.group.array)
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
  
  ## tune JMLE model
  lambda.vec = sqrt(log(p)/n) * seq(0.8, 0.1, -0.1)
  nlambda = length(lambda.vec)
  
  ## get all models
  loopfun1 = function(m){
    jmle(Y.list, Y.indices, X.list, B.group.array=B0.group.array, Theta.groups=Theta.groups,
         lambda = lambda.vec[m],
         gamma = sqrt(log(q)/n) * seq(0.8, 0.2, -0.2),
         init.option=1, tol=1e-3)
  }
  system.time(
    model.list <- mclapply(1:nlambda, loopfun1, mc.cores=min(detectCores(),nlambda))
  )
  # system.time(
  #   model.list <- lapply(1:nlambda, loopfun1)
  # )
  
  ## calculate BIC
  bic.vec = rep(0, nlambda)
  hbic.vec = rep(0, nlambda)
  
  for(m in 1:nlambda){
    jmle.model = model.list[[m]]
    SSE.vec = rep(0,K)
    bic.pen.vec = rep(0,K)
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
      bic.pen.vec[k] = log(nk)/nk * (sum(Theta.k != 0)/2 + sum(jmle.model$B.refit[,,k] != 0))
      hbic.pen.vec[k] = log(log(nk))*log(q*(q-1)/2)/nk * sum(Theta.k != 0)/2 +
        log(log(nk))*log(p*q)/nk * sum(jmle.model$B.refit[,,k] != 0)
    }
    bic.vec[m] = sum(SSE.vec) + sum(bic.pen.vec)
    hbic.vec[m] = sum(SSE.vec) + sum(hbic.pen.vec)
  }
  
  ## all model evaluations
  Theta_new.array = array(0, c(q,q,K))
  eval.mat = matrix(0, ncol=9, nrow=nlambda)
  eval.mat[,1] = round(lambda.vec,4)
  eval.mat[,2] = round(bic.vec,4)
  eval.mat[,3] = round(hbic.vec,4)
  
  for(m in 1:length(lambda.vec)){
    model.m = model.list[[m]]
    for(k in 1:K){
      Theta_new.array[,,k] = model.m$Theta_refit$Theta[[k]]
    }
    eval.mat[m,-(1:3)] = round(c(sum(B0.array != 0 & model.m$B.refit != 0)/sum(B0.array != 0),
                                 sum(B0.array == 0 & model.m$B.refit == 0)/sum(B0.array == 0),
                                 sqrt(sum((B0.array - model.m$B.refit)^2)/sum(B0.array^2)),
                                 sum(Theta0.array != 0 & Theta_new.array != 0)/sum(Theta0.array != 0),
                                 sum(Theta0.array == 0 & Theta_new.array == 0)/sum(Theta0.array == 0),
                                 sqrt(sum((Theta0.array - Theta_new.array)^2)/sum(Theta0.array^2))), 3)
  }
  eval.mat = data.frame(eval.mat)
  names(eval.mat) = c("lambda","BIC","HBIC","B.TP","B.TN","B.rel.Fnorm","Theta.TP","Theta.TN","Theta.rel.Fnorm")
  eval.mat
}

eval.list = lapply(1:1e2, loopfun)
save(eval.list, file="out_n100p30q60k5.Rda")