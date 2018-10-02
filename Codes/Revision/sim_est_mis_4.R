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
                       sparsity.B=5, sparsity.Theta=5, K=5, seed.vec=1:30, filename=NULL){
  
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
    B0.array = CoefArray1(B0.group.array, missing.prob=0.2)
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
    lambda.vec = sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
    model.list = vector("list", length(lambda.vec))
    nlambda = length(lambda.vec)
    
    ## get all models
    loopfun1 = function(m){
      model = jmmle.1step(Y.list, Y.indices, X.list, B.group.array=B0.group.array, Theta.groups=Theta.groups,
                          lambda = lambda.vec[m],
                          gamma = sqrt(log(q)/n) * seq(1, 0.4, -0.1),
                          init.option=1, tol=1e-3)
      
      ## post-process to tackle misspecification ***************************
      
      ## Tune JSEM model for X
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
      
      ## calculate debiased coefficients and scaling factors
      B.hat.array = model$B.refit
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
      
      Omega.hat.array = array(0, c(q,q,K))
      for(k in 1:K){
        Omega.hat.array[,,k] = model$Theta_refit$Omega[[k]]
      }
      
      D.array = C.hat.array
      for(i in 1:p){
        for(j in 1:q){
          D.array[i,j,] = Omega.hat.array[j,j,]* (M[i,]* C.hat.array[i,j,])^2
        }
      }
      
      ## run the FDR procedure
      B.thr.array = B.hat.array
      tau.vec = seq(0, 20, length.out=1e2)
      alpha = 0.2
      for(k in 1:K){
        for(i in 1:p){
          thres.vec = as.numeric(lapply(tau.vec, function(x) alpha/q * max(sum(D.array[i,,k]>x),1)))
          which.less = which((1 - pchisq(tau.vec,1)) <= thres.vec)
          if(length(which.less)>0){
            tau.i = tau.vec[which.less[1]] # set tau as minimizer only if there is at least one tau entry less
          }
          B.thr.array[i,which(D.array[i,,k] < tau.i),k] = 0 # set elements below threshold to zero
        }
      }
      
      ## return new coefficients
      model$B.refit = B.thr.array
      model
    }
    model.list <- mclapply(1:nlambda, loopfun1, mc.cores=nlambda)
    # model.list <- lapply(1:nlambda, loopfun1)
    
    ## calculate HBIC
    hbic.vec = rep(NA, nlambda)
    for(m in 1:nlambda){
      jmle.model = model.list[[m]]
      
      if(class(jmle.model)=="list"){ ## if no error in training the model
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
    }
    
    ## select best model
    jmmle.model = model.list[[which.min(hbic.vec)]]
    
    ## Calculate metrics ***************************************************
    # **********************************************************************
    Theta_new.array = array(0, c(q,q,K))
    for(k in 1:K){
      Theta_new.array[,,k] = jmmle.model$Theta_refit$Theta[[k]]
    }
    
    TP.B = sum(B0.array != 0 & jmmle.model$B.refit != 0)
    TN.B = sum(B0.array == 0 & jmmle.model$B.refit == 0)
    FN.B = sum(B0.array != 0) - TP.B
    FP.B = sum(B0.array == 0) - TN.B
    FDR.B = sum(B0.array == 0 & jmmle.model$B.refit != 0)/max(sum(jmmle.model$B.refit != 0),1)
    TP.Theta = sum(Theta0.array != 0 & Theta_new.array != 0)
    TN.Theta = sum(Theta0.array == 0 & Theta_new.array == 0)
    FP.Theta = sum(Theta0.array != 0) - TP.Theta
    FN.Theta = sum(Theta0.array == 0) - TN.Theta
    FDR.Theta = 0
    # for(m in 1:length(lambda.vec)){
    #   model.m = model.list[[m]]
    # 
    #   eval.mat[m,-(1:3)] = round(c(sum(B0.array != 0 & model.m$B.refit != 0)/sum(B0.array != 0),
    #                                sum(B0.array == 0 & model.m$B.refit == 0)/sum(B0.array == 0),
    #                                sqrt(sum((B0.array - model.m$B.refit)^2)/sum(B0.array^2)),
    #                                sum(Theta0.array != 0 & Theta_new.array != 0)/sum(Theta0.array != 0),
    #                                sum(Theta0.array == 0 & Theta_new.array == 0)/sum(Theta0.array == 0),
    #                                sqrt(sum((Theta0.array - Theta_new.array)^2)/sum(Theta0.array^2))), 3)
    # }
    # eval.mat = data.frame(eval.mat)
    cat("=============\nReplication",rep,"done!\n=============\n")
    rbind(c(evaluate(TP.B,TN.B,FP.B,FN.B),
            sqrt(sum((B0.array - jmmle.model$B.refit)^2)/sum(B0.array^2)), FDR.B),
          c(evaluate(TP.Theta,TN.Theta,FP.Theta,FN.Theta),
            sqrt(sum((Theta0.array - Theta_new.array)^2)/sum(Theta0.array^2)), FDR.Theta)
    )
  }
  
  # out.mat = mclapply(1:nrep, loopfun, mc.cores=8)
  out.mat = lapply(90+seed.vec, loopfun)
  if(is.null(filename)){
    filename = paste0("estmis_n",n,"p",p,"q",q,"_4.Rda")
  }
  save(out.mat, file=filename)
}

##### Generate data
get.outputs(n = 100, subnetSize.X = c(30, 30), subnetSize.E = c(15, 15))
get.outputs(n = 100, subnetSize.X = c(15, 15), subnetSize.E = c(30, 30))
get.outputs(n = 150, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100))
get.outputs(n = 150, subnetSize.X = c(150, 150), subnetSize.E = c(150, 150))
get.outputs(n = 200, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="estmis_n200p200q200modelB_4.Rda")
get.outputs(n = 100, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="estmis_n100p200q200modelB_4.Rda")
