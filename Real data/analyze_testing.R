# initial work
rm(list=ls())
# setwd("/n/subho-data/JMMLE-outputs/real-data")
setwd("C:/Study/Stratified-mult-GGM/Real data/")
Required.Packages <- c("data.table", "glmnet","glasso","parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

## load data ***********************************************************
# **********************************************************************
data = readRDS("processed_data.rds")
final_model = readRDS("final_model.rds")
final_model_sep = readRDS("final_model_sep.rds")
model_list_sep = readRDS("model_list_sep.rds")
n = sapply(data$X.list1, nrow)
p = ncol(data$X.list1[[1]])
q = ncol(data$Y.list1[[1]])
K = 2

## Tune JSEM model for X ***********************************************
# **********************************************************************
# indices of X
X.indices <- vector("list", K)
for (k in 1:K){
    X.indices[[k]] <- rep(k, n[[k]])  
}

# Zeta groups: group structure in X
Zeta.groups = vector("list", p)
Xg = data$Xg
for (i in 1:p){
    Zeta.groups[[i]] = matrix(0, K, p)
    for (k in 1:K){
        Zeta.groups[[i]][k,] = match(Xg, unique(Xg))
    }
    Zeta.groups[[i]] = Zeta.groups[[i]][,-i]
}
gamma = sqrt(log(p)/min(n)) * seq(1, 0.4, -0.1)
X.list = data$X.list1
bic.jsem <- sel.lambda.jsem(do.call(rbind, X.list), do.call(rbind, X.list),
                            unlist(X.indices), unlist(X.indices),
                            Zeta.groups,lambda=gamma)
gamma.min = gamma[which.min(bic.jsem$BIC)]
jsem.model = JSEM(do.call(rbind, X.list), unlist(X.indices),
                  Zeta.groups, lambda=gamma.min)
Zeta_new.array = array(0, c(p,p,K))
for(k in 1:K){
    Zeta_new.array[,,k] = jsem.model$Theta[[k]]
}

## Refit to get Omega
if(is.null(jsem.model)){
    Ahat = list()
    for(k in 1:K){
        Ahat[[k]] = matrix(0, q, q)
        Ahat[[k]][which(abs(Theta_new.array[,,k])>eps, arr.ind=T)] = 1
        diag(Ahat[[k]]) = 0
    }
} else{
    Ahat = jsem.model$Ahat
}
Info = list()
for (k in 1:K){
    Info[[k]] = zeroInd(Ahat[[k]], 1)$zeroArr
}
Theta_refit = multi.glasso(do.call(rbind, X.list), unlist(X.indices), jsem.model$lambda, Info)

# non-zero values in Omega-x
Om.values = vector("list",K)
groups = c("ER+","ER-")
X.names = colnames(X.list[[1]])
for(k in 1:K){
    non.zero.Om = which(Theta_refit$Omega[[k]] != 0, arr.ind=T)
    non.zero.Om = non.zero.Om[non.zero.Om[,1] < non.zero.Om[,2],] # just take lower triangle
    Om.df = data.table(SampleGroup = groups[k],
                       mRNA1 = X.names[non.zero.Om[,1]],
                       mRNA2 = X.names[non.zero.Om[,2]])
    invisible(Om.df[, Value := 0])
    for(i in 1:nrow(Om.df)){
        Om.df[i, Value := Theta_refit$Omega[[k]][non.zero.Om[i,1],non.zero.Om[i,2]]]
    }
    Om.values[[k]] = Om.df[order(abs(Value), decreasing=T)]
}

Om.values
fwrite(rbindlist(Om.values), file="Omegax_values.csv", sep=",")

## Get debiased estimates **********************************************
# **********************************************************************
B.hat.array = final_model$B.refit
B.hat.array = final_model_sep$B_sep_array
C.hat.array = B.hat.array
Y.list = data$Y.list1
M = matrix(0,p,K)
for(k in 1:K){
    X.k = X.list[[k]]
    E.k = Y.list[[k]] - X.k %*% B.hat.array[,,k]
    for(i in 1:p){
        R.ik = X.k[,i] - X.k[,-i] %*% Zeta_new.array[i,-i,k]
        t.ik = as.numeric(t(R.ik) %*% X.k[,i]/n[k])
        C.hat.array[i,,k] = B.hat.array[i,,k] + t(R.ik) %*% E.k/n[k]/t.ik
        M[i,k] = sqrt(n[k])*t.ik/sqrt(sum(R.ik^2/n[k]))
    }
}

## Get eigenvectors and eigenvalues of precision matrices
Theta1 = solve(final_model$Theta_refit$Omega[[1]])
Theta2 = solve(final_model$Theta_refit$Omega[[2]])
alpha = .05

## Global test statistics for i-th X-variable
D = rep(0,p)
d = matrix(0,p,q)
for(i in 1:p){
    Pooled.Cov.i = Theta1/M[i,1]^2 + Theta2/M[i,2]^2
    Diff.i = C.hat.array[i,,1] - C.hat.array[i,,2]
    ## pairwise test statistics
    d[i,] = (Diff.i)^2/diag(Pooled.Cov.i)
    ## overall test statistic
    D[i] = t(Diff.i) %*% solve(Pooled.Cov.i) %*% Diff.i
}
mean(D > qchisq(1-alpha, q) & (B.hat.array[,,1]!=0 | B.hat.array[,,2]!=0))
mean((B.hat.array[,,1]!=0 | B.hat.array[,,2]!=0))

## determine threshold for i-th test
alpha = .2
d.ind.mat = matrix(0,p,q)
tau = rep(20,p)
which.i.reject = which(D > qchisq(.95, q))
for(i in which.i.reject){
    tau.vec = seq(0, 20, length.out=1e2)
    thres.vec = as.numeric(lapply(tau.vec, function(x) alpha/q * max(sum(d[i,]>x),1)))
    which.less = which((1 - pchisq(tau.vec,1)) <= thres.vec)
    if(length(which.less)>0){
        tau[i] = tau.vec[which.less[1]] # set tau as minimizer only if there is at least one tau entry less
    }
    # tau[i] = tau.vec[which.min(abs(1 - pchisq(tau.vec,1) - thres.vec))]
    
    d.ind.mat[i,] = as.numeric(d[i,]>tau[i])
}
