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
n = max(sapply(data$X.list1,nrow))
p = ncol(data$X.list1[[1]])
q = ncol(data$Y.list1[[1]])
K = length(data$X.list1)
lambda.vec = sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
nlambda = length(lambda.vec)

## get all models
loopfun1 = function(m){
  with(data, jmmle.1step(
    Y.list=Y.list1, X.list=X.list1,
    B.group.array=B.group.array, Theta.groups=Theta.groups,
    lambda = lambda.vec[m],
    gamma = sqrt(log(q)/n) * seq(0.8, 0.1, -0.1),
    init.option=1, tol=1e-3))
}
system.time(
  model.list <- mclapply(1:nlambda, loopfun1, mc.cores=min(detectCores(),nlambda))
)
saveRDS(model.list, "model_list.rds")

## calculate HBIC
hbic.vec = rep(0, nlambda)
for(m in which(sapply(model.list,length)==2)){
  jmle.model = model.list[[m]]
  SSE.vec = rep(0,K)
  hbic.pen.vec = rep(0,K)
  
  for(k in 1:K){
    nk = nrow(data$Y.list1[[k]])
    Theta.k = jmle.model$Theta_refit$Theta[[k]]
    for(j in 1:q)
    {
      Theta.k[j,j] = 0
    }
    SSE.vec[k] = with(data, sum(diag(crossprod((Y.list1[[k]] - X.list1[[k]] %*% jmle.model$B.refit[,,k]) %*% (diag(1,q) - Theta.k))))/nk)
    hbic.pen.vec[k] = log(log(nk))*log(q*(q-1)/2)/nk * sum(Theta.k != 0)/2 +
      log(log(nk))*log(p*q)/nk * sum(jmle.model$B.refit[,,k] != 0)
  }
  hbic.vec[m] = sum(SSE.vec) + sum(hbic.pen.vec)
}

final.model = model.list[[which.min(hbic.vec[hbic.vec>0])]]
saveRDS(final.model, "final_model.rds")
