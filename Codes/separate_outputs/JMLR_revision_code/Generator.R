## this file stores the functions that are used to generate data 
require(matrixcalc); require(MASS); require(gdata);

# RegMat() generates the regression matrix
RegMat = function(p,q,sparsity=NULL,SNR=NULL){
	# p regressors, q regressions
	# default sparsity is 5/p (model 1 in Cai, model 2 in Cai has a sparsity of 30/p)
	if (is.null(sparsity))
		sparsity = 5/p;
	
	B = array(0,c(p,q));
	if (is.null(SNR)){
		for (i in 1:p)
			for (j in 1:q)
				B[i,j] = rbinom(1,1,sparsity)*sample(c(-1,1),1)*runif(1,0.5,1);
	}
	else{
		for (i in 1:p)
			for (j in 1:q)
				B[i,j] = rbinom(1,1,sparsity)*sample(c(-1,1),1)*runif(1,0,SNR);
	}
	return(B)
}

# ThetaMat() generates the precision matrix
ThetaMat = function(q,type="random",CN=NULL,Theta1=NULL,Theta2=NULL){
	# by default, generates a random graph with sparsity level being 5/q;
	# if conditional number is not specified, the default will be the dimension of the matrix (as in Cai)
	# alternatively, type could be "band", with the off-diagonals provided
	
	Theta = array(0,c(q,q))
	if (type=="random"){
		diag(Theta) = 0;
		if (is.null(CN))
			CN = q;
		# for off-diagonals:
		for (i in 1:(q-1)){
			for (j in (i+1):q){
				Theta[j,i] = rbinom(1,1,5/q)*sample(c(-1,1),1)*runif(1,0.5,1);
				Theta[i,j] = Theta[j,i];
			}
		}
		# now bump-up the diagonal to get the desired condition number:
		Theta = Theta + diag(q)*(0.001+abs(min(eigen(Theta)$values))); 
		egval = eigen(Theta)$values; 
		CN_Theta = max(egval)/min(egval)
		
		while(CN_Theta>CN){
			Theta = Theta + 0.01*diag(q);
			egval = eigen(Theta)$values;
			CN_Theta = max(egval)/min(egval);
		}
		return(list(Theta=Theta,CN_Theta=CN_Theta))
	}
	if (type=="band"){
		diag(Theta) = 1;
		if (is.null(Theta2)){ # tridiagonal, inverse matrix of AR(1)
			if (abs(Theta1)>=1){
				stop("Wrong Input Value for Theta1, Theta needs to be diagonally dominant.\n")
			}
			for (i in 1:(q-1)){
				for (j in (i+1):q){
					Theta[i,j] = ifelse(j-i==1,Theta1,0);
					Theta[j,i] = Theta[i,j];
				}
			}	
		}
		else{	# two bands 
			if ((abs(Theta1)+abs(Theta2))>1)
				stop("Wrong Input Value for Off-Diagonals, Theta needs to be diagonally dominant.\n")
			for (i in 1:(q-1)){
				for (j in (i+1):q){
					Theta[i,j] = ifelse(j-i==1,sample(c(-1,1),1)*Theta1,ifelse(j-i==2,sample(c(-1,1),1)*Theta2,0));
					Theta[j,i] = Theta[i,j];
				}
			}
		}
		egval = eigen(Theta)$values;
		CN_Theta = max(egval)/min(egval);
		return(list(Theta=Theta,CN_Theta=CN_Theta))
	}
}

gen.2layer = function(n,B,Theta,standardize=TRUE,scale=FALSE){
	## standardize: if we standardize the covariance matrix, default=TRUE
	## scale: if we scale the regression matrix to impose the exact SNR, default=FALSE
	if (standardize){
		Sigma = cov2cor(solve(Theta));
	}
	else{
		Sigma = solve(Theta);
	}
	
	if (scale)
		B = scale(B,center=FALSE,scale=1/sqrt(diag(Sigma)))
	
	p1 = nrow(B); p2 = nrow(Theta);
	
	## Assuming each entry in X comes from iid standard Gaussian, centered
	X = array(rnorm(n*p1),c(n,p1));
	X = scale(X,center=TRUE,scale=FALSE);
	
	# generate error terms in layer 2
	E = mvrnorm(n,mu=rep(0,p2),Sigma=Sigma);
	E = scale(E,center=TRUE,scale=FALSE);
	
	# generate Y
	Y = X %*% B + E;
	
	FullData = cbind(X,Y)
	return(list(XY=FullData,X=X,Y=Y))
}




























