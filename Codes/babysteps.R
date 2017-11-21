rm(list=ls())
setwd('D:/Study/My projects/Stratified-mult-GGM/Codes')
source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')

library(glasso)

##### Generate data
group = rbind(
  c(1, 2),
  c(1, 4),
  c(3, 2),
  c(3, 4),
  c(5, 2),
  c(5, 4),
  c(6, 2),
  c(6, 4),
  c(7, 2),
  c(7, 4)
)                         # grouping pattern
subnetSize = c(10, 10)    # subnet size
n = 100
p = 20
q = 20
K = 10

## generate the two layers
set.seed(11192017)
X.layer = GenerateLayer(n, subnetSize, group)
E.layer = GenerateLayer(n, subnetSize, group)

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

## record the result from every iteration
## for algorithm testing/debugging purpose only
if (TRACE){
  TRACE.est =  vector("list")
}

## initialize B
skeleton.hat = array(1, c(p,q)); # assuming all directed edges possibly exist
Existing.edges = diag(1:p) %*% skeleton.hat
B_init.array = array(0, dim=c(p,q,K))
Ehat.list = list()

for(k in 1:K){
  LS = l1LS_Main(Y.list[[k]], X.list[[k]], skeleton.hat=skeleton.hat,
                 lambda=0.02,initializer="Lasso")
  B_init.array[,,k] = LS$B0
  Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% LS$B0
}


## initialize Theta
jsem.grid <- sqrt(log(q)/n) * seq(1, 0.1, -0.1)
init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                            unlist(Y.indices), unlist(Y.indices),
                            Theta.groups, lambda=jsem.grid)
lambda.jsem <- jsem.grid[which.min(init.bic.jsem$BIC)]
init.jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                       Theta.groups, lambda=lambda.jsem) ## Y.indices, Theta.groups.list tbd
Theta_init.array = array(0, c(q,q,K))
for(k in 1:K){
  Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
}

## initialize Theta
# jsem.grid <- sqrt(log(p)/n) * seq(1, 0.1, -0.1)
# init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Y.list), do.call(rbind, Y.list),
#                             unlist(Y.indices), unlist(Y.indices),
#                             Theta.groups, lambda=jsem.grid)
# lambda.jsem <- jsem.grid[which.min(init.bic.jsem$BIC)]
# init.jsem.model = JSEM(do.call(rbind, Y.list), unlist(Y.indices),
#                        Theta.groups, lambda=lambda.jsem) ## Y.indices, Theta.groups.list tbd
# Theta_init.array = array(0, c(q,q,K))
# for(k in 1:K){
#   Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
# }
# B_init.array = array(0, dim=c(p,q,K))

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
                  family="gaussian", penalty="grLasso", lambda=sqrt(log(p)/n)*.5)
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
  cat("Updating Theta array\n")
  if (iter >=10 | sqrt(sum(B_new.array - B_old.array)^2)<0.1 | updateTheta == TRUE){
    updateTheta = TRUE; # once we start updating Theta, we just start from now on;
    bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                                unlist(Y.indices), unlist(Y.indices),
                                Theta.groups, lambda=jsem.grid)
    lambda.jsem <- jsem.grid[which.min(bic.jsem$BIC)]
    jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                           Theta.groups, lambda=lambda.jsem) ## Y.indices, Theta.groups.list tbd
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
  cat("Norm_diff =",round(abs(Norm_diff),4),"\n")
  CONVERGE = (abs(Norm_diff)<1e-4)
  
  if (iter == 1e2){
    cat("Max iterations reached",'\n')
    break;
  }
}

if(CONVERGE){
  cat("Converged.\n")
} else{
  warning("algorithm didn't converge.\nRequired epsilon for convergence is ", 1e-4,
      " while current value is ", format(abs(Norm_diff), scientific=T))
}

##### Refit if converged
## if converge, we start the refitting procedure
if(CONVERGE){ 

  ## Refit B
  B_refit.array = array(0,c(p,q,K));
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
  
  ## Refit Theta
  Ahat <- jsem.model$Ahat
  Info = list()
  for (k in 1:K){
    Info[[k]] = zeroInd(Ahat[[k]], 1)$zeroArr
  }
  Theta_refit = multi.glasso(do.call(rbind, Ehat.list), unlist(Y.indices),
                                 lambda.jsem, Info)
}

c(sum(B0.array != 0 & B_init.array != 0)/sum(B0.array != 0),
  sum(B0.array == 0 & B_init.array == 0)/sum(B0.array == 0),
  sqrt(sum((B0.array - B_init.array)^2)))

c(sum(B0.array != 0 & B_refit.array != 0)/sum(B0.array != 0),
  sum(B0.array == 0 & B_refit.array == 0)/sum(B0.array == 0),
  sqrt(sum((B0.array - B_new.array)^2)),
  sqrt(sum((B0.array - B_refit.array)^2)))

c(sum(Theta0.array != 0 & Theta_init.array != 0)/sum(Theta0.array != 0),
  sum(Theta0.array == 0 & Theta_init.array == 0)/sum(Theta0.array == 0),
  sqrt(sum((Theta0.array - Theta_init.array)^2)))

c(sum(Theta0.array != 0 & Theta_new.array != 0)/sum(Theta0.array != 0),
  sum(Theta0.array == 0 & Theta_new.array == 0)/sum(Theta0.array == 0),
  sqrt(sum((Theta0.array - Theta_new.array)^2)))

c(sum(Theta0.array != 0 & Theta_refit$Omega != 0)/sum(Theta0.array != 0),
  sum(Theta0.array == 0 & Theta_refit$Omega != 0)/sum(Theta0.array == 0),
  sqrt(sum((Theta0.array - Theta_refit$Omega)^2)))

  
## plot bic values
plot(bic.mat[2,], type='l', lwd=2, col=adjustcolor("darkblue", alpha.f = .5+1/nrow(bic.mat)/2), ylim=c(50,250))
for(it in 1:nrow(bic.mat)){
  lines(bic.mat[it,], type='l', lwd=2, col=adjustcolor("darkblue", alpha.f = .5+it/nrow(bic.mat)/2))
}
lines(bic.mat[1,], type='l', lwd=2, col="black")

## separate estimation
setwd('JMLR_revision_code')
source('l1ML_Main.R')
B_lin.array = array(0, c(p,q,K))
Theta_lin.array = array(0, c(q,q,K))

for(k in 1:K){
  model.k = l1ML_Main(Y.list[[k]], X.list[[k]], lambda=0.02,initializer="Lasso",
                      rho=sqrt(log(q)/n))
  B_lin.array[,,k] = model.k$B.est
  Theta_lin.array[,,k] = model.k$Theta.est
}

c(sum(B0.array != 0 & B_lin.array != 0)/sum(B0.array != 0),
  sum(B0.array == 0 & B_lin.array == 0)/sum(B0.array == 0),
  sqrt(sum((B0.array - B_lin.array)^2)))

c(sum(Theta0.array != 0 & Theta_lin.array != 0)/sum(Theta0.array != 0),
  sum(Theta0.array == 0 & Theta_lin.array == 0)/sum(Theta0.array == 0),
  sqrt(sum((Theta0.array - Theta_lin.array)^2)))
