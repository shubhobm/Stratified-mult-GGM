rm(list=ls())
# setwd('d:/Study/My projects/Stratified-mult-GGM/Codes')
source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

library(glasso)
library(parallel)

##### Function to calculate evaluation metrics
evaluate = function(TP,TN,FP,FN){
  SEN = TP/(TP+FN)
  SPE = TN/(TN+FP)
  MCC = (TP*TN - FP*FN)/sqrt(TP+FP)/sqrt(TP+FN)/sqrt(TN+FP)/sqrt(TN+FN)
  F1 = 2*TP/(2*TP+FP+FN)
  c(SEN,SPE,MCC,F1)
}

##### Common wrapper function
get.outputs = function(n=100, subnetSize.X=rep(10,2), subnetSize.E=rep(10,2),
                       sparsity.B=5, sparsity.Theta=5, K=5, nrep=50, filename=NULL){
  
  ## Set up some quantities
  group = rbind(
    c(1, 2),
    c(1, 4),
    c(3, 2),
    c(3, 4),
    c(5, 2)
  )  # grouping pattern
  p = sum(subnetSize.X)
  q = sum(subnetSize.E)
  
  loopfun = function(rep){
    set.seed(1e3*rep)
    
    ## Generate data *******************************************************
    # **********************************************************************
    X.layer = GenerateLayer(n, subnetSize.X, group, D=1, sparsity=sparsity.Theta/p)
    E.layer = GenerateLayer(n, subnetSize.E, group, D=1, sparsity=sparsity.Theta/q)
    
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
    
    ## Tune JSEM model *****************************************************
    # **********************************************************************
    for(k in 1:K){
      Y.list[[k]] = scale(Y.list[[k]], center=T, scale=F)
    }
    
    gamma = sqrt(log(p)/n) * seq(1, 0.4, -0.1)
    bic.jsem <- sel.lambda.jsem(do.call(rbind, Y.list), do.call(rbind, Y.list),
                                unlist(Y.indices), unlist(Y.indices),
                                Theta.groups,lambda=gamma)
    gamma.min = gamma[which.min(bic.jsem$BIC)]
    jsem.model = JSEM(do.call(rbind, Y.list), unlist(Y.indices),
                      Theta.groups, lambda=gamma.min)
    Theta_new.array = array(0, c(q,q,K))
    for(k in 1:K){
      Theta_new.array[,,k] = jsem.model$Theta[[k]]
    }
    
    TP.Theta = sum(Theta0.array != 0 & Theta_new.array != 0)
    TN.Theta = sum(Theta0.array == 0 & Theta_new.array == 0)
    FP.Theta = sum(Theta0.array != 0) - TP.Theta
    FN.Theta = sum(Theta0.array == 0) - TN.Theta
    cat("=============\nReplication",rep,"done!\n=============\n")
    rbind(c(evaluate(TP.Theta,TN.Theta,FP.Theta,FN.Theta),
            sqrt(sum((Theta0.array - Theta_new.array)^2)/sum(Theta0.array^2)))
    )
  }
  
  out.mat = mclapply(1:nrep, loopfun, mc.cores=10)
  # out.mat = lapply(1:nrep, loopfun)
  if(is.null(filename)){
    filename = paste0("jsem_n",n,"p",p,"q",q,".Rda")
  }
  save(out.mat, file=filename)
}

##### Generate data
# get.outputs(n = 100, subnetSize.X = c(30, 30), subnetSize.E = c(15, 15))
# get.outputs(n = 100, subnetSize.X = c(15, 15), subnetSize.E = c(30, 30))
get.outputs(n = 150, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100))
get.outputs(n = 150, subnetSize.X = c(150, 150), subnetSize.E = c(150, 150))
get.outputs(n = 200, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="jsem_n200p200q200modelB.Rda")
get.outputs(n = 100, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="jsem_n100p200q200modelB.Rda")

# 
# 
# out.mat = matrix(unlist(out.mat), ncol=4, byrow=T)
# rbind(round(apply(out.mat,2,mean),3),
#       round(apply(out.mat,2,sd),3))