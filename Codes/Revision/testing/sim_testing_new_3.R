rm(list=ls())
# setwd('d:/Study/My projects/Stratified-mult-GGM/Codes')
source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

library(glasso)
library(parallel)

##### Common wrapper function
get.outputs = function(n=100, subnetSize.X=rep(10,2), subnetSize.E=rep(10,2),
                       sparsity.B=5, sparsity.Theta=5, K=2, seed.vec=1:30, filename=NULL){
  
  ## Set up some quantities
  group = matrix(c(1, 2), nrow=2, ncol=2, byrow=T)           # grouping pattern
  p = sum(subnetSize.X)
  q = sum(subnetSize.E)
  
  loopfun = function(seed){
    set.seed(1e3*seed)
    
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
    B0.array = CoefArray2(B0.group.array[,,1], D=1, sparsity=sparsity.B/p)
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
    lambda.vec = sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
    model.list = vector("list", length(lambda.vec))
    nlambda = length(lambda.vec)
    
    ## get all models
    loopfun1 = function(m){
      jmmle.1step(Y.list, Y.indices, X.list, B.group.array=B0.group.array, Theta.groups=Theta.groups,
                  lambda = lambda.vec[m],
                  gamma = sqrt(log(q)/n) * seq(1, 0.4, -0.1),
                  init.option=1, tol=1e-3)
    }
    model.list <- mclapply(1:nlambda, loopfun1, mc.cores=nlambda)
    
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
    
    pow.simul = length(which.i.reject)/p
    pow = sum(d.ind.mat == 1 & Diff.mat != 0, na.rm=T)/sum(Diff.mat != 0)
    size = 1 - sum(d.ind.mat == 0 & Diff.mat == 0, na.rm=T)/sum(Diff.mat == 0)
    FDP = sum(d.ind.mat == 1 & Diff.mat == 0, na.rm=T)/max(sum(d.ind.mat == 1, na.rm=T),1)
    cat("=============\nReplication",seed,"done!\n=============\n")
    c(pow.simul,pow,size,FDP)
  }
  
  # out.mat = mclapply(1:nrep, loopfun, mc.cores=8)
  out.mat = lapply(60+seed.vec, loopfun)
  if(is.null(filename)){
    filename = paste0("testnew_n",n,"p",p,"q",q,"_3.Rda")
  }
  save(out.mat, file=filename)
}

##### Generate data
get.outputs(n = 100, subnetSize.X = c(30, 30), subnetSize.E = c(15, 15))
get.outputs(n = 200, subnetSize.X = c(30, 30), subnetSize.E = c(15, 15))

get.outputs(n = 100, subnetSize.X = c(15, 15), subnetSize.E = c(30, 30))
get.outputs(n = 200, subnetSize.X = c(15, 15), subnetSize.E = c(30, 30))

get.outputs(n = 150, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100))

get.outputs(n = 150, subnetSize.X = c(150, 150), subnetSize.E = c(150, 150))
get.outputs(n = 300, subnetSize.X = c(150, 150), subnetSize.E = c(150, 150))

get.outputs(n = 100, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="testnew_n100p200q200modelB_3.Rda")
get.outputs(n = 200, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="testnew_n200p200q200modelB_3.Rda")
get.outputs(n = 300, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="testnew_n300p200q200modelB_3.Rda")
