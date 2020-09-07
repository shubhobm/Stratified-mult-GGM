## prelims *************************************************************
# **********************************************************************
rm(list=ls())
setwd("/n/subho-data/JMMLE-outputs/real-data")
source('jsem.R')

Required.Packages <- c("data.table", "glmnet","glasso","parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

## load data and initialize ********************************************
# **********************************************************************
data = readRDS("processed_data.rds")
n = sapply(data$X.list1, nrow)
p = ncol(data$X.list1[[1]])
q = ncol(data$Y.list1[[1]])
K = 2
rho = sqrt(log(q)/min(n))*10^seq(1,-2,by=-.5)

## indices and groups **************************************************
# **********************************************************************
# Theta groups: group structure in X
Theta.groups = vector("list", q)
Yg = data$Yg
for (j in 1:q){
    Theta.groups[[j]] = matrix(0, K, q)
    for (k in 1:K){
        Theta.groups[[j]][k,] = match(Yg, unique(Yg))
    }
    Theta.groups[[j]] = Theta.groups[[j]][,-j]
}

## load splits
nrep = 1e2
splits = readRDS("cv_splits.rds")
train1 = splits$train1
train2 = splits$train2

# set looping function and run
repfun = function(rep){

    # set samples
    Y.list = list(data$Y.list1[[1]][train1[[rep]],], data$Y.list1[[2]][train2[[rep]],])
    
    # indices of Y
    Y.indices <- list(rep(1, length(train1[[rep]])), rep(2, length(train2[[rep]])))
    
    # choose bic
    bic.jsem = sel.lambda.jsem(do.call(rbind, Y.list), do.call(rbind, Y.list),
                               unlist(Y.indices), unlist(Y.indices),
                               Theta.groups,lambda=rho)
    # select best model
    which_good = which(bic.jsem$BIC!=0)
    gamma.min = rho[which_good][which.min(bic.jsem$BIC[which_good])]
    # retrain best model and return
    JSEM(do.call(rbind, Y.list), unlist(Y.indices), Theta.groups, lambda=gamma.min)
}
model.list = mclapply(1:nrep, function(x) try(repfun(x)), mc.cores=50)

saveRDS(model.list, file="cv_results_jsem.rds")
