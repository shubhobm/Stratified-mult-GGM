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
                       sparsity.B=5, sparsity.Theta=5, K=2, seed.vec=1:10, filename=NULL){
  
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
    
    ## Obtain separate fits ************************************************
    # **********************************************************************
    lambda = sqrt(log(p)/n) * seq(1, 0.4, -0.1)
    
    skeleton.hat = array(1, c(p,q)); # assuming all directed edges possibly exist
    Existing.edges = diag(1:p) %*% skeleton.hat
    B.array = array(0, dim=c(p,q,K))
    bic.vec = rep(0, length(lambda))
    
    # For each choice of lambda, obtain B
    loopfun.B = function(nlam){
      B.array = array(0, c(p,q,K))
      lam = lambda[nlam]
      for(k in 1:K){
        LS = l1LS_Main(Y.list[[k]], X.list[[k]], skeleton.hat=skeleton.hat,
                       lambda=lam,initializer="Lasso")
        B.array[,,k] = LS$B0
        bic.vec[nlam] = bic.vec[nlam] +
          sum((Y.list[[k]] - X.list[[k]] %*% LS$B0)^2) +
          log(n)*sum(B.array[,,k]!=0)
      }
      
      # store B estimate
      B.array
    }
    
    B_list = mclapply(1:length(lambda), loopfun.B, mc.cores=min(length(lambda),32))
    B.hat.array = B_list[[which.min(bic.vec)]]
    
    ## Get nbhd coefs and precision matrices for Y *************************
    # **********************************************************************
    gamma = sqrt(log(p)/n) * seq(1, 0.4, -0.1)
    bic.vec = rep(0, length(gamma))
    
    # calculate covariance matrices
    empcov = list()
    Ehat.list = list()
    for(k in 1:K){
      Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% B.hat.array[,,k]
      empcov[[k]] = cov(Ehat.list[[k]]) 
      while (kappa(empcov[[k]]) > 1e+2){
        empcov[[k]] = empcov[[k]] + 0.05 * diag(p)      
      }
    }
    
    loopfun.Theta = function(ngam){
      Theta.array = array(0, c(q,q,K))
      Omega.y.array = array(0, c(q,q,K))
      gam = gamma[ngam]
      for(k in 1:K){
        for(i in 1:q){
          lasso.model = glmnet(Ehat.list[[k]][,-i], Ehat.list[[k]][,i], lambda=gam)
          Theta.array[-i,i,k] = as.numeric(lasso.model$beta)
        }
        
        # get glasso fit
        Ahat = Theta.array[,,k]
        Ahat[which(abs(Ahat)>0, arr.ind=T)] = 1
        if(sum(abs(Ahat))>0){
          Omega.y.array[,,k] = glasso(empcov[[k]],rho = 0.1*log(q)/n, zero = Ahat)$w
        } else{
          Omega.y.array[,,k] = diag(1/diag(empcov[[k]]))
        }
        bic.vec[ngam] = bic.vec[ngam] +
          matTr(empcov[[k]] %*% Omega.y.array[,,k]) -log(det(Omega.y.array[,,k])) +
          log(n) * sum(Ahat)/(2*n)
      }
      
      list(Theta.array, Omega.y.array)
    }
    
    Theta_list = mclapply(1:length(gamma), loopfun.Theta,
                          mc.cores=min(length(lambda),32))
    Theta_new.array = Theta_list[[which.min(bic.vec)]][[1]]
    Omega.y = Theta_list[[which.min(bic.vec)]][[2]]
    
    ## Get nbhd coefs and precision matrices for X *************************
    # **********************************************************************
    gamma = sqrt(log(p)/n) * seq(1, 0.4, -0.1)
    Zeta_list = list()
    bic.vec = rep(0, length(gamma))
    
    # calculate covariance matrices
    empcov = list()
    for(k in 1:K){
      empcov[[k]] = cov(X.list[[k]]) 
      while (kappa(empcov[[k]]) > 1e+2){
        empcov[[k]] = empcov[[k]] + 0.05 * diag(p)      
      }
    }
    
    loopfun.Zeta = function(ngam){
      Zeta.array = array(0, c(p,p,K))
      Omega.array = array(0, c(p,p,K))
      gam = gamma[ngam]
      for(k in 1:K){
        for(i in 1:p){
          lasso.model = glmnet(X.list[[k]][,-i], X.list[[k]][,i], lambda=gam)
          Zeta.array[-i,i,k] = as.numeric(lasso.model$beta)
        }
        
        # get glasso fit
        Ahat = Zeta.array[,,k]
        Ahat[which(abs(Ahat)>0, arr.ind=T)] = 1
        if(sum(abs(Ahat))>0){
          Omega.array[,,k] = glasso(empcov[[k]],rho = 0.1*log(p)/n, zero = Ahat)$w
        } else{
          Omega.array[,,k] = diag(1/diag(empcov[[k]]))
        }
        bic.vec[ngam] = bic.vec[ngam] +
          matTr(empcov[[k]] %*% Omega.array[,,k]) -log(det(Omega.array[,,k])) +
          log(n) * sum(Ahat)/(2*n)
      }
      
      list(Zeta.array, Omega.array)
    }
    
    Zeta_list = mclapply(1:length(gamma), loopfun.Zeta,
                          mc.cores=min(length(lambda),32))
    Zeta_new.array = Zeta_list[[which.min(bic.vec)]][[1]]

    ## Get debiased estimates **********************************************
    # **********************************************************************
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
    Theta1 = solve(Omega.y[,,1])
    Theta2 = solve(Omega.y[,,2])
    
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
  out.mat = lapply(seed.vec, loopfun)
  if(is.null(filename)){
    filename = paste0("testsep_n",n,"p",p,"q",q,"_1.Rda")
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
            sparsity.B=30, sparsity.Theta=30, filename="testsep_n100p200q200modelB_1.Rda")
get.outputs(n = 200, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="testsep_n200p200q200modelB_1.Rda")
get.outputs(n = 300, subnetSize.X = c(100, 100), subnetSize.E = c(100, 100),
            sparsity.B=30, sparsity.Theta=30, filename="testsep_n300p200q200modelB_1.Rda")
