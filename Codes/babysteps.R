rm(list=ls())
setwd('D:/Study/My projects/Stratified-mult-GGM/Codes')
source('code_jsem.R')


##### Given: X.list, Y.list, B.groups, Theta.groups
Objfunc = c(); iter = 0; CONVERGE=FALSE; update.counter=0;
updateTheta = FALSE; # we don't update Theta until B is stabilized a bit

## initialize
init.jsem.model = JSEM(Y.list, Y.indices, Theta.groups, gamma=gamma) ## Y.indices, Theta.groups.list tbd
Theta0.array = list()
for(k in 1:K){
  Theta0.array[,,k] = init.jsem.model$Ahat[[k]]
}

## initialization of the residual matrix
Ehat.list = Y.list
# for(k in 1:K){
#   Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% B0.array[,,k]
# }

## record the result from every iteration
## for algorithm testing/debugging purpose only
if (TRACE){
  TRACE.est =  vector("list")
}
library(Matrix)
## make long X and Y matrices
Y = lapply(Y, rbind)
X = bdiag(X)

## Here we start with the alternating procedure
B0.array = array(0, dim=c(p,q,K))
B_new.array = B0.array; Theta_new.array = Theta0.array
while(!CONVERGE){
  iter = iter + 1;
  B_old.array = B_new.array
  Theta_old.array = Theta_new.array
  
  # Updating B
  foreach (j=1:q)%dopar% {
    
    # make long vector or errors for j-th column for all k
    Et.j.list = list()
    for(k in 1:K){
      Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_old.array[-j,j,k]
    }
    Ehat.theta.j = unlist(Et.j.list)
    rm(Et.j.list)
    
    # build model
    temp = grpreg(X, Y[,j] - Ehat.theta.j, B.group, # B.group tbd
                  family="gaussian", penalty="grLasso", lambda=lambda)
    B_new.array[,j,] = matrix(temp$beta[-1], ncol=K, byrow=F)

    # now update j-th column of all K matrices in B_new, and Ehat
    for(k in 1:K){
      Ehat.list[[k]][,j] = Y.list[[k]] - X.list[[k]] %*% B_new.array[,j,k]
    }
  }
  
  # Updating Theta
  if (iter >=10 | norm(B_new-B_old,"F")<0.1 | updateTheta == TRUE){
    updateTheta = TRUE; # once we start updating Theta, we just start from now on;
    jsem.model = JSEM(Ehat.list, Y.indices, Theta.groups.list, gamma=gamma) ## Y.indices, index tbd
    Theta_new.array = list()
    for(k in 1:K){
      Theta_new.array[,,k] = jsem.model$Ahat[[k]]
    }
    update.counter = update.counter + 1
  }
  
  # check convergence
  Objfunc[iter] = Obj(Y.list, X.list, Theta_new.array, B_new.array, Theta.groups, B.groups, lambda, gamma)
  if (iter == 1){
    Obj_diff = Objfunc[1]
  }
  else{
    Obj_diff = Objfunc[iter] - Objfunc[iter-1]
  }
  CONVERGE = (abs(Obj_diff)<1e-4)
  
  if (iter == 100){
    cat("Exceeds the maximum number of iterations allowed",'\n')
    break;
  }
  
}