## prelims *************************************************************
# **********************************************************************
rm(list=ls())
setwd("/n/subho-data/JMMLE-outputs/real-data")
source('l1ML_Main.R')
source('Generator.R')
source('jsem.R')

Required.Packages <- c("data.table", "glmnet","glasso","parallel","optparse")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

## parse arguments******************************************************
# **********************************************************************
option_list = list(
	make_option(c("-o", "--offset"), type="double", default=0, help="offset for simulation replication [default = %default]", metavar="double")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
offset = opt$o

## load data ***********************************************************
# **********************************************************************
data = readRDS("processed_data.rds")
n = sapply(data$X.list1, nrow)
p = ncol(data$X.list1[[1]])
q = ncol(data$Y.list1[[1]])
K = 2
lambda.vec = sqrt(log(p))*10^seq(-1,-3,by=-.5)
rho.vec = sqrt(log(q))*10^seq(1,-2,by=-.5)
nlambda = length(lambda.vec)
nrho = length(rho.vec)
ntune = nlambda*nrho

## load splits
nrep = 10
splits = readRDS("/n/subho-data/JMMLE-outputs/real-data/cv_splits.rds")
train1 = splits$train1
train2 = splits$train2

# initialize and run
best_model = vector("list", nrep)
bic = rep(0, nrep)
for(i in 1:nrep){

oi = 10*offset + i
X.list = list(data$X.list1[[1]][train1[[oi]],], data$X.list1[[2]][train2[[oi]],])
Y.list = list(data$Y.list1[[1]][train1[[oi]],], data$Y.list1[[2]][train2[[oi]],])
eval.mat = matrix(0, ncol=9, nrow=ntune)

## get all models
loopfun1 = function(m3){
    m1 = floor((m3-1)/nrho) + 1
    m2 = m3 - (m1-1)*nrho
    
    ## get estimates
    model.m3.list = list()
    for(k in 1:K){
        model.m3.list[[k]] = l1ML_Main(Y.list[[k]], X.list[[k]],
                                       lambda=lambda.vec[m1]/sqrt(n[k]), rho=rho.vec[m2]/sqrt(n[k]),
                                       initializer="Lasso", StabilizeTheta=F, VERBOSE=F, maxit=50)
    }
    cat("Done- lambda=",lambda.vec[m1],"rho=",rho.vec[m2],"\n")
    model.m3.list
}
model.list = mclapply(1:ntune, loopfun1, mc.cores=min(detectCores(), ntune))

## choose best model by BIC ********************************************
# **********************************************************************
ntuning = length(model.list)
bic.vec = rep(0, ntuning)
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

# find minimum, save model and hbic
min_ind = which.min(bic.vec[bic.vec>0])
if(length(min_ind)>0){
	best_model[[i]] = model.list[[min_ind[1]]]
	bic[i] = bic.vec[min_ind[1]]
}

cat("====================\nSplit",oi,"done!!\n====================\n")
}

saveRDS(list(best_model=best_model, bic=bic),
	file=paste0("/n/subho-data/JMMLE-outputs/real-data/cv_results_sep_",offset,".rds"))
