library(gdata)
squaredError = function(Y.list, X.list, Theta.list, B.list){
  # this function calculates the total squared error for all list entries
  
  err.list = list()
  for(k in 1:K){
    nk = nrow(X.list[[k]])
    qk = ncol(Y.list[[k]])
    Ek = Y.list[[k]] - X.list[[k]] %*% B.list[[k]]
    Tk = diag(1,qk) - Theta.list[[k]]
    err.list[[k]] = sum(diag(crossprod(Ek %*% Tk)))/ nk
  }
  lapply(err.list, sum)
}

Obj = function(Y.list, X.list, Theta.array, B.array,
               Theta.group.array, B.group.array, lambda, gamma){
  # this function calculates the objective function
  
  # calculate penalties
  # For each group index in the group array, collect corresponding elements in the main array ...
  # and sum their l2 norms
  unique.Theta.groups = unique(Theta.group.array)
  Theta.norm = 0
  for(ng in 1:unique.Theta.groups){
    Theta.norm = Theta.norm + sqrt(sum(Theta.group.array[which(Theta.group.array==ng, arr.ind=T)]^2))
  }
  
  unique.B.groups = unique(B.group.array)
  B.norm = 0
  for(ng in 1:unique.B.groups){
    B.norm = B.norm + sqrt(sum(B.group.array[which(B.group.array==ng, arr.ind=T)]^2))
  }
  
  squaredError(Y.list, X.list, Theta.list, B.list) + lambda*Theta.pen + gamma*B.norm
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



		
		
				

