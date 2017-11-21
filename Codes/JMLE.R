## The Joint Multi-level Estimation (JMLE) Method
jmle = function(Y.list, Y.indices=NULL, X.list,
                B.group.array, Theta.groups,
                lambda=NULL, gamma=NULL,
                init.option=1, refit.B=TRUE, tol=1e-4, maxit=20){
  
  ## initialize quantities
  arraydims = dim(B.group.array)
  p = arraydims[1]; q = arraydims[2]; K = arraydims[3]
  n = nrow(Y.list[[1]])
  Theta.group.array = array(0, c(q,q,K))
  for(j in 1:q){
    Theta.group.array[j,-j,] = Y.layer$groups[[j]]
  }
  
  ## default values if NULL
  if(is.null(Y.indices)){
    Y.indices = list()
    for(k in 1:K){
      Y.indices[[k]] = rep(k, nrow(Y.list[[k]]))
    }
  }
  if(is.null(lambda)){
    lambda = .5 * sqrt(log(p*q)/n)
  }
  if(is.null(gamma)){
    gamma = sqrt(log(q)/n) * seq(1, 0.1, -0.1)
  }

  skeleton.hat = array(1, c(p,q)); # assuming all directed edges possibly exist
  Existing.edges = diag(1:p) %*% skeleton.hat
  B_init.array = array(0, dim=c(p,q,K))
  Ehat.list = list()
  
  ## Two options for intialization
  ## Option 1: initialize B, then run JSEM to initialize Theta
  if(init.option==1){
    ## initialize B
    cat("Initializing B array\n")
    skeleton.hat = array(1, c(p,q)); # assuming all directed edges possibly exist
    Existing.edges = diag(1:p) %*% skeleton.hat
    B_init.array = array(0, dim=c(p,q,K))
    Ehat.list = list()

    for(k in 1:K){
      LS = l1LS_Main(Y.list[[k]], X.list[[k]], skeleton.hat=skeleton.hat,
                     lambda=lambda,initializer="Lasso")
      B_init.array[,,k] = LS$B0
      Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% LS$B0
    }


    ## initialize Theta
    cat("Initializing Theta array: ")
    init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                                     unlist(Y.indices), unlist(Y.indices),
                                     Theta.groups, lambda=gamma)
    init.jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                           Theta.groups, lambda=gamma[which.min(init.bic.jsem$BIC)])
    Theta_init.array = array(0, c(q,q,K))
    for(k in 1:K){
      Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
    }
  }

  ## Option 2: initialize Theta only by running JSEM once
  if(init.option==2){
    ## initialize Theta
    cat("Initializing Theta array: ")
    init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Y.list), do.call(rbind, Y.list),
                                     unlist(Y.indices), unlist(Y.indices),
                                     Theta.groups, lambda=gamma)
    init.jsem.model = JSEM(do.call(rbind, Y.list), unlist(Y.indices),
                           Theta.groups, lambda=gamma[which.min(init.bic.jsem$BIC)])
    Theta_init.array = array(0, c(q,q,K))
    for(k in 1:K){
      Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
    }
    B_init.array = array(0, dim=c(p,q,K))
    Ehat.list = Y.list
  }

  ## make long X and Y matrices
  Y = do.call(rbind, Y.list)
  X = as.matrix(do.call(bdiag, X.list))
  
  # initialize
  Normfunc = c(); iter = 0; CONVERGE=FALSE; refit.B=TRUE; update.counter=0;
  updateTheta = FALSE; # we don't update Theta until B is stabilized a bit
  B_new.array = B_init.array
  Theta_new.array = Theta_init.array
  
  # # store bic values from JSEM models
  # bic.mat = init.bic.jsem$BIC
  
  ## Here we start with the alternating procedure
  cat('-----\n')
  while(!CONVERGE){
    iter = iter + 1;
    B_old.array = B_new.array
    Theta_old.array = Theta_new.array
    
    cat("Iteration ", iter, ":\n")
    # Updating B
    cat("Updating B array\n")
    for(j in 1:q){
      
      # make long vector or errors for j-th column for all k
      Et.j.list = list()
      for(k in 1:K){
        Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_old.array[-j,j,k]
      }
      Ehat.theta.j = unlist(Et.j.list)
      rm(Et.j.list)
      
      # build model
      temp = grpreg(X, Y[,j] + Ehat.theta.j, unlist(as.numeric(B0.group.array[,j,])),
                    family="gaussian", penalty="grLasso", lambda=lambda)
      B_new.array[,j,] = matrix(temp$beta[-1], ncol=K, byrow=F)
      
      ## refit if necessary
      if(refit.B){
        # make long vector or errors for j-th column for all k
        temp1 = temp$beta[-1]
        B.j.support = which(abs(temp1)>1e-6)
        if (length(B.j.support)>0){
          # build model
          temp2 = lm(Y[,j] + Ehat.theta.j~X[,B.j.support]+0)
          
          # save non-zero coefs
          temp1[B.j.support] = temp2$coef
          B_new.array[,j,] = matrix(temp1, ncol=K, byrow=F)
        }
      }
      
      # now update j-th column of all K matrices in B_new, and Ehat
      for(k in 1:K){
        Ehat.list[[k]][,j] = Y.list[[k]][,j] - X.list[[k]] %*% as.matrix(B_new.array[,j,k], ncol=1)
      }
    }
    
    # # update array of E
    # for(k in 1:K){
    #   Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% B_new.array[,,k]
    # }
    
    # Updating Theta
    if (iter >=10 | sqrt(sum(B_new.array - B_old.array)^2)<0.1 | updateTheta == TRUE){
      updateTheta = TRUE; # once we start updating Theta, we just start from now on;
      cat("Updating Theta array: ")
      bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                                  unlist(Y.indices), unlist(Y.indices),
                                  Theta.groups, lambda=gamma)
      jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                        Theta.groups, lambda=gamma[which.min(bic.jsem$BIC)])
      Theta_new.array = array(0, c(q,q,K))
      for(k in 1:K){
        Theta_new.array[,,k] = jsem.model$Theta[[k]]
      }
      update.counter = update.counter + 1
      # bic.mat = rbind(bic.mat, bic.jsem$BIC)
    }
    
    # check convergence
    # Objfunc[iter] = Obj(Y.list, X.list, Theta_new.array, B_new.array,
    #                     Theta.group.array, B0.group.array,
    #                     lambda=sqrt(log(p)/n)*.5, gamma=lambda.jsem)
    Normfunc[iter] = sqrt(sum(B_new.array - B_old.array)^2)/sqrt(sum(B_new.array^2)) +
      sqrt(sum(Theta_new.array - Theta_old.array)^2)/sqrt(sum(Theta_new.array^2))
    if (iter == 1){
      Norm_diff = Normfunc[1]
    }
    else{
      Norm_diff = Normfunc[iter] - Normfunc[iter-1]
    }
    cat("Norm_diff =",round(abs(Norm_diff),4),'\n-----\n')
    CONVERGE = (abs(Norm_diff)<tol)
    
    if (iter == maxit){
      cat("Max iterations reached",'\n')
      break;
    }
  }
  
  if(CONVERGE){
    cat("Converged after",iter,"iterations.\n")
  } else{
    warning("algorithm didn't converge.\nRequired epsilon for convergence is ", tol,
            " while current value is ", round(abs(Norm_diff), 4))
  }
  
  ## If refitting of B matrices hasn't been done inside the loop then refit in the end
  B_refit.array = B_new.array
  if(!refit.B){ 
    for(j in 1:q){

            # make long vector or errors for j-th column for all k
      B.j.support = which(abs(B_new.array[,j,])>1e-6)
      if (length(B.j.support)>0){
        Et.j.list = list()
        for(k in 1:K){
          Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_old.array[-j,j,k]
        }
        Ehat.theta.j = unlist(Et.j.list)
        rm(Et.j.list)
        
        # build model
        temp = lm(Y[,j] + Ehat.theta.j~X[,B.j.support]+0)
        
        # save non-zero coefs
        temp1 = B_refit.array[,j,]
        temp1[B.j.support] = temp$coef
        B_refit.array[,j,] = temp1
      }
    }
  }
  
  ## Refit Theta
  Ahat <- jsem.model$Ahat
  Info = list()
  for (k in 1:K){
    Info[[k]] = zeroInd(Ahat[[k]], 1)$zeroArr
  }
  Theta_refit = multi.glasso(do.call(rbind, Ehat.list), unlist(Y.indices),
                             gamma[which.min(bic.jsem$BIC)], Info)
  
  ## return
  return(list(B.refit=B_refit.array, Theta_refit=Theta_refit))
}



##################################
# EOF
##################################