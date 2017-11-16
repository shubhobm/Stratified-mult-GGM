##### Given: X, Y, B0, Theta0: all K-length lists of matrices
## initialization of the residual matrix
Ehat = list();
for(k in 1:K){
  Ehat[[k]] = (Y[[k]] - X[[k]] %*% B0[[k]])
}

## record the result from every iteration
## for algorithm testing/debugging purpose only
if (TRACE){
  TRACE.est =  vector("list")
}

## make long X and Y matrices
Y.long = lapply(Y, rbind)
X.long = lapply(X, rbind)

## Here we start with the alternating procedure
B_new = B0; Theta_new = Theta0
while(!CONVERGE){
  iter = iter + 1;
  B_old = B_new ; Theta_old = Theta_new; 
  
  # Updating B
  output.list_B = foreach (j=1:q)%dopar% {
    
    # make long vector or errors for j-th column for all k
    Et.j.list = list()
    for(k in 1:K){
      Et.j.list[[k]] = Ehat[[k]][,-j] %*% Theta_old[[k]][-j,j]
    }
    Et.j.long = unlist(Et.j.list)
    rm(Et.j.list)
    
    # if (length(which(Existing.edges[,j]!=0))==0)
    #   Bj = Bj;
    # if (length(which(Existing.edges[,j]!=0))==1){
    #   coef1 = lm((Y[,j]+R[,j])~X[,which(Existing.edges[,j]!=0)]+0)$coef;
    #   ## this is the least square solution, but we need to apply a soft-thresholding
    #   coef1_st = max(0,coef1-lambda/(2*Theta_old[j,j]))*sign(coef1);
    #   B_j[which(Existing.edges[,j]!=0)] = coef1_st;
    # }				
    # if (length(which(Existing.edges[,j]!=0))>=2){
    #   temp = glmnet(X[,which(Existing.edges[,j]!=0)],(Y[,j]+R[,j]),intercept=FALSE);
    #   B_j[which(Existing.edges[,j]!=0)] = predict(temp,s=lambda/(2*Theta_old[j,j]),type="coefficients")[-1];
    # }
    # B_j;
    temp = grpreg(X.long, Y.long[,j] - Et.j.long, B.group, # B.group tbd
                  family="gaussian", penalty="grLasso", lambda=lambda)
    B_new.j = matrix(temp$beta[-1], ncol=K, byrow=F)
    E.j = matrix(Y.long[,j] - X.long %*% temp$beta[-1], ncol=K, byrow=F)
    
    # now update j-th column of all K matrices in B_new, and Ehat
    for(k in 1:K){
      B_new[[k]][,j] = B.new.j[,k]
      Ehat[[k]][,j] = E.j[,k]
    }
  }
  
  # Updating Theta
  if (iter >=10 | norm(B_new-B_old,"F")<0.1 | updateTheta == TRUE){
    updateTheta = TRUE; # once we start updating Theta, we just start from now on;
    jsem.model = JSEM(Ehat, Y.indices, index, gamma=gamma) ## Y.indices, index tbd
    Theta_new = jsem.model$Ahat
    update.counter = update.counter + 1;
  }
  
  # check convergence
}