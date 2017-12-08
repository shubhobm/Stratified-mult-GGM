rm(list=ls())
# setwd('D:/Study/My projects/Stratified-mult-GGM/Codes/')
setwd('JMLR_revision_code')
source('l1ML_Main.R')
source('../Generator.R')
source('../jsem.R')
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
  
  ## separate estimation
  B_sep.array = array(0, c(p,q,K))
  Theta_sep.array = array(0, c(q,q,K))
  lambda.vec = sqrt(log(p)/n)*seq(0.8,0.2,-0.2)
  rho.vec = sqrt(log(q)/n)*seq(0.8,0.2,-0.2)
  nlambda = length(lambda.vec)
  nrho = length(rho.vec)
  eval.mat = matrix(0, ncol=9, nrow=nlambda*nrho)
  m3 = 0
  
  ## loop over tuning parameter grid
  loopfun1 = function(m3){
    m1 = floor((m3-1)/nrho) + 1
    m2 = m3 - (m1-1)*nrho
    
    ## get estimates
    for(k in 1:K){
      model.k = l1ML_Main(Y.list[[k]], X.list[[k]],
                          lambda=lambda.vec[m1], rho=rho.vec[m2],
                          initializer="Lasso", StabilizeTheta=F, VERBOSE=F)
      B_sep.array[,,k] = model.k$B.est
      Theta_sep.array[,,k] = model.k$Theta.est
    }
    
    ## calculate BIC
    bic.vec = rep(0, K)
    for(k in 1:K){
      nk = nrow(Y.list[[k]])
      B.k = B_sep.array[,,k]
      Theta.k = Theta_sep.array[,,k]
      bic.vec[k] = -log(det(Theta.k)) +
        sum(diag(crossprod(Y.list[[k]] - X.list[[k]] %*% B.k) %*% Theta.k))/nk +
        log(nk)/nk * ((sum(Theta.k != 0)-q)/2 + sum(B.k != 0))
    }
    
    ## return
    cat("Done- lambda=",lambda.vec[m1],"rho=",rho.vec[m2],"\n")
    c(round(c(lambda.vec[m1], rho.vec[m2], sum(bic.vec)),4),
      round(c(sum(B0.array != 0 & B_sep.array != 0)/sum(B0.array != 0),
              sum(B0.array == 0 & B_sep.array == 0)/sum(B0.array == 0),
              sqrt(sum((B0.array - B_sep.array)^2)/sum(B0.array^2)),
              sum(Theta0.array != 0 & Theta_sep.array != 0)/sum(Theta0.array != 0),
              sum(Theta0.array == 0 & Theta_sep.array == 0)/sum(Theta0.array == 0),
              sqrt(sum((Theta0.array - Theta_sep.array)^2)/sum(Theta0.array^2))), 3))
  }
  
  eval.mat1 = mclapply(1:(nlambda*nrho), loopfun1, mc.cores=8)
  
  ## remove error entries
  eval.mat = list()
  m4 = 0
  for(m3 in 1:(nlambda*nrho)){
    if(class(eval.mat1[[m3]])=="numeric"){
      m4 = m4 + 1
      eval.mat[[m4]] = eval.mat1[[m3]]
    }
  }
  eval.mat = matrix(unlist(eval.mat), ncol=9, byrow=T)
  cat("Replication",rep,"done!\n")
  # eval.mat = matrix(unlist(eval.mat1), nrow=nlambda*nrho, byrow=T)
  # eval.mat = data.frame(eval.mat)
  # names(eval.mat) = c("lambda","rho","BIC","B.TP","B.TN","B.rel.Fnorm","Theta.TP","Theta.TN","Theta.rel.Fnorm")
  eval.mat
}

eval.list = lapply(1:2, loopfun)
save(eval.list, file="../out_n100p30q60k5_sep.Rda")
