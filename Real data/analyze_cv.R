# initial work
rm(list=ls())
# setwd("/n/subho-data/JMMLE-outputs/real-data")
setwd("C:/Study/Stratified-mult-GGM/Real data/")
Required.Packages <- c("data.table", "parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

## load data ***********************************************************
# **********************************************************************
data = readRDS("processed_data.rds")
cv_results = readRDS("cv_results.rds")
cv_splits = readRDS("cv_splits.rds")
X = data$X.list1
Y = data$Y.list1
n = sapply(X, nrow)
p = ncol(X[[1]])
q = ncol(Y[[1]])
K = 2
nrep = 100

## compute MSSPE = scaled prediction error
spe.vec = rep(0, nrep)
Bprop.vec = rep(0, nrep)
Oprop.vec = rep(0, nrep)

which_good = which(cv_results$hbic>0)
for(i in which_good){
    
    iBhat = cv_results$best_model[[i]]$B.refit
    iOmega = cv_results$best_model[[i]]$Theta_refit$Omega
    
    for(k in 1:K){
        train.ik = cv_splits[[k]][[i]]
        spe.vec[i] = spe.vec[i] + sum(diag(
            crossprod(Y[[k]][-train.ik,] - X[[k]][-train.ik,] %*% iBhat[,,k]) %*% iOmega[[k]]
        ))/n[k]
    }
    Bprop.vec[i] = mean(iBhat!=0)
    
    # non-zero coef proportion in Omega_y
    iTheta = array(0, c(q,q,K))
    for(k in 1:K){
        iTheta[,,k] = iOmega[[k]]
    }
    for(j in 1:q){
        iTheta[j,j,] = 0
    }
    Oprop.vec[i] = mean(iTheta!=0)
}
summarize = function(x) c(mean(x), sd(x))
summarize(sqrt(spe.vec[which_good]))
summarize(Bprop.vec[which_good])
summarize(Oprop.vec[which_good])

## compute MSSPE = scaled prediction error
spe0.vec = rep(0, nrep)
for(i in which(cv_results$hbic>0)){
    
    iBhat = cv_results$best_model[[i]]$B.refit
    iOmega = cv_results$best_model[[1]]$Theta_refit$Omega
    
    for(k in 1:K){
        train.ik = cv_splits[[k]][[i]]
        spe0.vec[i] = spe.vec[i] + sum(diag(
            crossprod(Y[[k]][-train.ik,]) %*% iOmega[[k]]
        ))/n[k]
    }
}
summarize(sqrt(spe0.vec))
