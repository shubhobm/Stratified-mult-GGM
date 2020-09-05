# initial work
rm(list=ls())
# setwd("/n/subho-data/JMMLE-outputs/real-data")
setwd("C:/Study/Stratified-mult-GGM/Real data/")
Required.Packages <- c("data.table", "parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

## load data ***********************************************************
# **********************************************************************
data = readRDS("processed_data.rds")

# load cv results
chunk_inds = c(0:3,6:9)
nchunk = length(chunk_inds)
cv_results = vector("list",nchunk)
for(i in 1:nchunk){
    cv_results[[i]] = readRDS(paste0("cv_results_sep_",chunk_inds[i],".rds"))
}
best_model = do.call(c, lapply(cv_results, function(x) x$best_model))
bic = do.call(c, lapply(cv_results, function(x) x$bic))
cv_results = list(best_model=best_model, bic=bic)

# load cv splits
cv_splits = readRDS("cv_splits.rds")
cv_splits[[1]] = cv_splits[[1]][c(1:40,61:100)]
cv_splits[[2]] = cv_splits[[2]][c(1:40,61:100)]

X = data$X.list1
Y = data$Y.list1
n = sapply(X, nrow)
p = ncol(X[[1]])
q = ncol(Y[[1]])
K = 2
nrep = 80

## compute MSSPE = scaled prediction error
spe.vec = rep(0, nrep)
Bprop.vec = rep(0, nrep)
Oprop.vec = rep(0, nrep)

which_good = which(cv_results$bic>0)
for(i in which_good){
    
    iBhat = array(0,c(p,q,K))
    for(k in 1:K){
        iBhat[,,K] = cv_results$best_model[[i]][[k]]$B.est
    }
    iOmega = array(0,c(q,q,K))
    for(k in 1:K){
        iOmega[,,K] = cv_results$best_model[[i]][[k]]$Theta.est
    }

    for(k in 1:K){
        train.ik = cv_splits[[k]][[i]]
        spe.vec[i] = spe.vec[i] + sum(diag(
            crossprod(Y[[k]][-train.ik,] - X[[k]][-train.ik,] %*% iBhat[,,k]) %*% iOmega[,,k]
            # crossprod(Y[[k]][-train.ik,]) %*% iOmega[,,k]
        ))/n[k]
    }
    Bprop.vec[i] = mean(iBhat!=0)
    
    # non-zero coef proportion in Omega_y
    iTheta = array(0, c(q,q,K))
    for(k in 1:K){
        iTheta[,,k] = iOmega[,,k]
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
