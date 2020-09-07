# initial work
rm(list=ls())
# setwd("/n/subho-data/JMMLE-outputs/real-data")
setwd("C:/Study/Stratified-mult-GGM/Real data/")
Required.Packages <- c("data.table", "parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

## load data ***********************************************************
# **********************************************************************
data = readRDS("processed_data.rds")
cv_results = readRDS("cv_results_jsem.rds")
cv_splits = readRDS("cv_splits.rds")
X = data$X.list1
Y = data$Y.list1
n = sapply(X, nrow)
p = ncol(X[[1]])
q = ncol(Y[[1]])
K = 2
nrep = length(cv_results)

## compute MSSPE = scaled prediction error
Oprop.vec = rep(0, nrep)
for(i in 1:nrep){
    
    # non-zero coef proportion in Omega_y
    iTh = cv_results[[i]]$Theta
    iTheta = array(0, c(q,q,K))
    for(k in 1:K){
        iTheta[,,k] = iTh[[k]]
    }
    for(j in 1:q){
        iTheta[j,j,] = 0
    }
    Oprop.vec[i] = mean(iTheta!=0)
}
summarize = function(x) c(mean(x), sd(x))
summarize(Oprop.vec)
