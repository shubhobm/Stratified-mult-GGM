# initial work
rm(list=ls())
setwd("/n/subho-data/JMMLE-outputs/real-data")
Required.Packages <- c("data.table", "glmnet","glasso","parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')
# source('./Codes/lasso_inference.r')

# load data
data = readRDS("processed_data.rds")

## tune JMLE model
p = ncol(data$X.list1[[1]])
q = ncol(data$Y.list1[[1]])
K = length(data$X.list1)

## generate splits
set.seed(08062020)
nrep = 100
train1 = vector("list", nrep)
n1 = nrow(data$X.list1[[1]])
train2 = train1
n2 = nrow(data$X.list1[[2]])
for(i in 1:nrep){
	train1[[i]] = sample(1:n1, ceiling(.8*n1), replace=F) 
	train2[[i]] = sample(1:n2, ceiling(.8*n2), replace=F) 
}
saveRDS(list(train1=train1, train2=train2), file="/n/subho-data/JMMLE-outputs/real-data/cv_splits.rds")

# initialize and run
best_model = vector("list", nrep)
hbic = rep(0, nrep)
mse = hbic

for(i in 1:nrep){

X.list1.i = list(data$X.list1[[1]][train1[[i]],], data$X.list1[[2]][train2[[i]],])
Y.list1.i = list(data$Y.list1[[1]][train1[[i]],], data$Y.list1[[2]][train2[[i]],])
n = max(sapply(X.list1.i,nrow))
lambda.vec = sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
nlambda = length(lambda.vec)

## get all models
loopfun1 = function(m){
	jmmle.1step(
		Y.list=Y.list1.i, X.list=X.list1.i,
		B.group.array=data$B.group.array, Theta.groups=data$Theta.groups,
		lambda = lambda.vec[m],
		gamma = sqrt(log(q)/n) * seq(0.8, 0.1, -0.1),
		init.option=1, tol=1e-3
	)
}
system.time(
  model.list <- mclapply(1:nlambda, loopfun1, mc.cores=min(detectCores(),nlambda))
)

## calculate HBIC
hbic.vec = rep(0, nlambda)
for(m in which(sapply(model.list,length)==2)){
  jmle.model = model.list[[m]]
  SSE.vec = rep(0,K)
  hbic.pen.vec = rep(0,K)
  
  for(k in 1:K){
    nk = nrow(Y.list1.i[[k]])
    Theta.k = jmle.model$Theta_refit$Theta[[k]]
    for(j in 1:q)
    {
      Theta.k[j,j] = 0
    }
    SSE.vec[k] = sum(diag(crossprod((Y.list1.i[[k]] - X.list1.i[[k]] %*% jmle.model$B.refit[,,k]) %*% (diag(1,q) - Theta.k))))/nk
    hbic.pen.vec[k] = log(log(nk))*log(q*(q-1)/2)/nk * sum(Theta.k != 0)/2 +
      log(log(nk))*log(p*q)/nk * sum(jmle.model$B.refit[,,k] != 0)
  }
  hbic.vec[m] = sum(SSE.vec) + sum(hbic.pen.vec)
}

# find minimum, save model and hbic
min_ind = which.min(hbic.vec[hbic.vec>0])
if(length(min_ind)>0){
	best_model[[i]] = model.list[[min_ind[1]]]
	hbic[i] = hbic.vec[min_ind[1]]
}

## calculate RMSE
#mse.part.vec = rep(0,K)
#for(k in 1:K){
#	nk = nrow(Y.list1.i[[k]])
#	Theta.k = jmle.model$Theta_refit$Theta[[k]]
#	for(j in 1:q)
#	{
#		Theta.k[j,j] = 0
#    	}
#	mse.part.vec[k] = sum(diag(crossprod(
#		data$Y.list1[[k]][-train1[[i]],] - data$X.list1[[k]][-train1[[i]],] %*% jmle.model$B.refit[,,k])))/nk
#}
#mse[i] = sum(mse.part.vec)
cat("====================\nSplit",i,"done!!\n====================\n")
}

saveRDS(list(best_model=best_model, hbic=hbic, mse=mse), file="/n/subho-data/JMMLE-outputs/real-data/cv_results.rds")
