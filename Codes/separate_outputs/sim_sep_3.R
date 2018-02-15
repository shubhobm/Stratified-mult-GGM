rm(list=ls())
# setwd('D:/Study/My projects/Stratified-mult-GGM/Codes/')
setwd('JMLR_revision_code/')
source('l1ML_Main.R')
setwd('../')
source('Generator.R')
source('jsem.R')

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
                       sparsity.B=5, sparsity.Theta=5, K=5, seed.vec=1:10, filename=NULL){
  
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
    
    ## Separate estimation *************************************************
    # **********************************************************************
    lambda.vec = sqrt(log(p)/n)*seq(1.8,0.4,-0.4)
    rho.vec = sqrt(log(q)/n)*seq(1,0.4,-0.2)
    nlambda = length(lambda.vec)
    nrho = length(rho.vec)
    eval.mat = matrix(0, ncol=9, nrow=nlambda*nrho)
    m3 = 0
    
    ## loop over tuning parameter grid
    loopfun1 = function(m3){
      m1 = floor((m3-1)/nrho) + 1
      m2 = m3 - (m1-1)*nrho
      
      ## get estimates
      model.m3.list = list()
      for(k in 1:K){
        model.m3.list[[k]] = l1ML_Main(Y.list[[k]], X.list[[k]],
                            lambda=lambda.vec[m1], rho=rho.vec[m2],
                            initializer="Lasso", StabilizeTheta=F, VERBOSE=F)
      }
      cat("Done- lambda=",lambda.vec[m1],"rho=",rho.vec[m2],"\n")
      model.m3.list
    }
    
    model.list0 = mclapply(1:(nlambda*nrho), loopfun1, mc.cores=8)
    # model.list0 = lapply(1:(nlambda*nrho), loopfun1)
    
    ## remove error entries
    model.list = list()
    m4 = 0
    for(m3 in 1:(nlambda*nrho)){
      error = FALSE
      for(k in 1:K){
        if(class(model.list0[[m3]])!="list"){
          error = TRUE
        }
      }
      if(!error){
        m4 = m4 + 1
        model.list[[m4]] = model.list0[[m3]]
      }
    }
    rm(model.list0)
    
    ## choose best model by BIC ********************************************
    # **********************************************************************
    ntuning = length(model.list)
    bic.vec = rep(NA, ntuning)
    
    for(m in 1:ntuning){
      sep.model = model.list[[m]]
      
      if(class(sep.model)=="list"){ ## get BIC if no error in training the model
        bic.part.vec = rep(0, K)
        for(k in 1:K){
          model.k = sep.model[[k]]
          nk = nrow(Y.list[[k]])
          B.k = model.k$B.est
          Theta.k = model.k$Theta.est
          bic.part.vec[k] = -log(det(Theta.k)) +
            sum(diag(crossprod(Y.list[[k]] - X.list[[k]] %*% B.k) %*% Theta.k))/nk +
            log(nk)/nk * ((sum(Theta.k != 0)-q)/2 + sum(B.k != 0))
        }
        bic.vec[m] = sum(bic.part.vec)
      }
    }
    
    best.model = model.list[[which.min(bic.vec)]]
    
    ## return
    B_sep.array = array(0, c(p,q,K))
    Theta_sep.array = array(0, c(q,q,K))
    for(k in 1:K){
      B_sep.array[,,k] = model.k$B.est
      Theta_sep.array[,,k] = model.k$Theta.est
    }
    
    TP.B = sum(B0.array != 0 & B_sep.array != 0)
    TN.B = sum(B0.array == 0 & B_sep.array == 0)
    FN.B = sum(B0.array != 0) - TP.B
    FP.B = sum(B0.array == 0) - TN.B
    TP.Theta = sum(Theta0.array != 0 & Theta_sep.array != 0)
    TN.Theta = sum(Theta0.array == 0 & Theta_sep.array == 0)
    FP.Theta = sum(Theta0.array != 0) - TP.Theta
    FN.Theta = sum(Theta0.array == 0) - TN.Theta
    
    cat("=============\nReplication",rep,"done!\n=============\n")
    rbind(c(evaluate(TP.B,TN.B,FP.B,FN.B),sqrt(sum((B0.array - B_sep.array)^2)/sum(B0.array^2))),
          c(evaluate(TP.Theta,TN.Theta,FP.Theta,FN.Theta),
            sqrt(sum((Theta0.array - Theta_sep.array)^2)/sum(Theta0.array^2))))
  }
  
  # out.mat = mclapply(1:nrep, loopfun, mc.cores=8)
  out.mat = lapply(20+seed.vec, loopfun)
  if(is.null(filename)){
    filename = paste0("estsep_n",n,"p",p,"q",q,"_3.Rda")
  }
  save(out.mat, file=filename)
}

##### Generate data
get.outputs(n = 100, subnetSize.X = c(30, 30), subnetSize.E = c(15, 15))
get.outputs(n = 100, subnetSize.X = c(15, 15), subnetSize.E = c(30, 30))
# get.outputs(n = 150, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100))
# get.outputs(n = 150, subnetSize.X = c(150, 150), subnetSize.E = c(150, 150))
# get.outputs(n = 200, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
#             sparsity.B=30, sparsity.Theta=30, filename="estsep_n200p200q200modelB_3.Rda")
# get.outputs(n = 100, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
#             sparsity.B=30, sparsity.Theta=30, filename="estsep_n100p200q200modelB_3.Rda")
# 
