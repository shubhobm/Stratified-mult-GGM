library(gdata)
squaredError = function(Y.list, X.list, Theta.array, B.array){
  # this function calculates the total squared error for all list entries
  
  err.list = list()
  K = dim(B.array)[3]
  for(k in 1:K){
    nk = nrow(X.list[[k]])
    qk = ncol(Y.list[[k]])
    Ek = Y.list[[k]] - X.list[[k]] %*% B.array[,,k]
    Tk = diag(1,qk) - Theta.array[,,k]
    err.list[[k]] = sum(diag(crossprod(Ek %*% Tk)))/ nk
  }
  do.call(sum, err.list)
}

Obj = function(Y.list, X.list, Theta.array, B.array,
               Theta.group.array, B.group.array, lambda, gamma){
  # this function calculates the objective function
  
  # calculate penalties
  # For each group index in the group array, collect corresponding elements in the main array ...
  # and sum their l2 norms
  unique.Theta.groups = unique(paste(Theta.group.array))
  Theta.norm = sapply(unique.Theta.groups, function(g)
    sum(Theta.array[which(Theta.group.array==g, arr.ind=T)]^2))
  Theta.norm = sum(sqrt(Theta.norm))
  
  unique.B.groups = unique(paste(B.group.array))
  B.norm = sapply(unique.B.groups, function(g)
    sum(B.array[which(B.group.array==g, arr.ind=T)]^2))
  B.norm = sum(sqrt(B.norm))
  
  squaredError(Y.list, X.list, Theta.array, B.array) + lambda*Theta.norm + gamma*B.norm
}

BICfunc = function(Y,X,Theta,B)
{	#BIC = -log(det(Theta)) + tr(S%*%Theta) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n + (log(n)/n)*(nonzeros in upperTriangle(Theta) + nonzeros in B)
	#And we prefer smaller BIC value
	n = nrow(X)
	PENs = log(n)/n*(sum(abs(upperTriangle(Theta))>1e-6) + sum(abs(B)>1e-6))
	FIT = sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n - log(det(Theta)) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n 
	BIC = FIT + PENs
	return(BIC)
}



		
		
				

