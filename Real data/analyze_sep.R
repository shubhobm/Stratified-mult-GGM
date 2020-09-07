## prelims*************************************************************
# **********************************************************************
rm(list=ls())
setwd("/n/subho-data/JMMLE-outputs/real-data")
source('l1ML_Main.R')
source('Generator.R')
source('jsem.R')

Required.Packages <- c("data.table", "glmnet","glasso","parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

## load data ***********************************************************
# **********************************************************************
data = readRDS("processed_data.rds")
n = sapply(data$X.list1, nrow)
p = ncol(data$X.list1[[1]])
q = ncol(data$Y.list1[[1]])
K = 2

## Separate estimation *************************************************
# **********************************************************************
#lambda.vec = sqrt(log(p))*seq(1.8,0.4,-0.4)
#rho.vec = sqrt(log(q))*seq(1,0.4,-0.2)
lambda.vec = sqrt(log(p))*10^seq(-1,-3,by=-.5)
rho.vec = sqrt(log(q))*10^seq(1,-2,by=-.5)
nlambda = length(lambda.vec)
nrho = length(rho.vec)
ntune = nlambda*nrho
eval.mat = matrix(0, ncol=9, nrow=ntune)

## loop over tuning parameter grid
X.list = data$X.list1
Y.list = data$Y.list1
loopfun1 = function(m3){
    m1 = floor((m3-1)/nrho) + 1
    m2 = m3 - (m1-1)*nrho
    
    ## get estimates
    model.m3.list = list()
    for(k in 1:K){
        model.m3.list[[k]] = l1ML_Main(Y.list[[k]], X.list[[k]],
                                       lambda=lambda.vec[m1]/sqrt(n[k]), rho=rho.vec[m2]/sqrt(n[k]),
                                       initializer="Lasso", StabilizeTheta=F, VERBOSE=F)
    }
    cat("Done- lambda=",lambda.vec[m1],"rho=",rho.vec[m2],"\n")
    model.m3.list
}
model.list0 = mclapply(1:ntune, loopfun1, mc.cores=min(detectCores(), ntune))

## remove error entries
model.list = list()
m4 = 0
for(m3 in 1:(nlambda*nrho)){
    error = FALSE
    for(k in 1:K){
        if(class(model.list0[[m3]])!="list"){
            error = TRUE
        }
    }
    if(!error){
        m4 = m4 + 1
        model.list[[m4]] = model.list0[[m3]]
    }
}
rm(model.list0)
# save outputs
saveRDS(model.list, file="model_list_sep.rds")

## choose best model by BIC ********************************************
# **********************************************************************
ntuning = length(model.list)
bic.vec = rep(NA, ntuning)

for(m in 1:ntuning){
    sep.model = model.list[[m]]
    
    if(class(sep.model)=="list"){ ## get BIC if no error in training the model
        bic.part.vec = rep(0, K)
        for(k in 1:K){
            model.k = sep.model[[k]]
            B.k = model.k$B.est
            Theta.k = model.k$Theta.est
            bic.part.vec[k] = -log(det(Theta.k)) +
                sum(diag(crossprod(Y.list[[k]] - X.list[[k]] %*% B.k) %*% Theta.k))/n[k] +
                log(n[k])/n[k] * ((sum(Theta.k != 0)-q)/2 + sum(B.k != 0))
        }
        bic.vec[m] = sum(bic.part.vec)
    }
}

## get best model coefs and save
best.model = model.list[[which.min(bic.vec)]]
B_sep.array = array(0, c(p,q,K))
Theta_sep.array = array(0, c(q,q,K))
for(k in 1:K){
    model.k = best.model[[k]]
    B_sep.array[,,k] = model.k$B.est
    Theta_sep.array[,,k] = model.k$Theta.est
}
saveRDS(list(B_sep_array=B_sep.array, Theta_sep_array=Theta_sep.array), file="final_model_sep.rds")