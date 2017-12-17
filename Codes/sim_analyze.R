setwd("D:/Study/My projects/Stratified-mult-GGM/Codes")

## *******************************************
## Setting 1: n=100, p=60, q=30, k=5
## *******************************************
## for jmle estimators
load('out_n100p60q30k5.Rda')
nrep = length(eval.list)
eval.list = lapply(eval.list, function(x) x[which(x$BIC>0),])
HBIC.mat = matrix(0, nrow=nrep, ncol=6)
for(i in 1:nrep) HBIC.mat[i,] = as.numeric(eval.list[[i]][which.min(eval.list[[i]][,3]),-(1:3)])
BIC.mat = matrix(0, nrow=nrep, ncol=6)
for(i in 1:nrep) BIC.mat[i,] = as.numeric(eval.list[[i]][which.min(eval.list[[i]][,2]),-(1:3)])
z = data.frame(rbind(round(apply(HBIC.mat,2,mean),3),
                     round(apply(HBIC.mat,2,sd),3),
                     round(apply(BIC.mat,2,mean),3),
                     round(apply(BIC.mat,2,sd),3)))
names(z) = c("B.TP","B.TN","B.relFnorm","Theta.TP","Theta.TN","Theta.relFnorm")
row.names(z) = c("mean.HBIC","sd.HBIC","mean.BIC","sd.BIC")
z

## for separate estimators
## output not there yet *******************************************
load('out_n100p60q30k5_sep.Rda')
BIC.mat = matrix(0, nrow=10, ncol=6)
for(i in 1:10) BIC.mat[i,] = as.numeric(eval.list[[i]][which.min(eval.list[[i]][,3]),-(1:3)])
zsep = data.frame(rbind(round(apply(BIC.mat,2,mean),3),
                        round(apply(BIC.mat,2,sd),3)))
names(zsep) = c("B.TP","B.TN","B.relFnorm","Theta.TP","Theta.TN","Theta.relFnorm")
row.names(zsep) = c("mean.BIC","sd.BIC")
zsep

## *******************************************
## Setting 2: n=100, p=30, q=60, k=5
## *******************************************
## for jmle estimators
load('out_n100p30q60k5_biggrid.Rda')
nrep = length(eval.list)
HBIC.mat = matrix(0, nrow=nrep, ncol=6)
for(i in 1:nrep) HBIC.mat[i,] = as.numeric(eval.list[[i]][which.min(eval.list[[i]][,3]),-(1:3)])
BIC.mat = matrix(0, nrow=nrep, ncol=6)
for(i in 1:nrep) BIC.mat[i,] = as.numeric(eval.list[[i]][which.min(eval.list[[i]][,2]),-(1:3)])
z = data.frame(rbind(round(apply(HBIC.mat,2,mean),3),
                     round(apply(HBIC.mat,2,sd),3),
                     round(apply(BIC.mat,2,mean),3),
                     round(apply(BIC.mat,2,sd),3)))
names(z) = c("B.TP","B.TN","B.relFnorm","Theta.TP","Theta.TN","Theta.relFnorm")
row.names(z) = c("mean.HBIC","sd.HBIC","mean.BIC","sd.BIC")
z

## for separate estimators
load('out_n100p30q60k5_sep.Rda')
BIC.mat = matrix(0, nrow=100, ncol=6)
for(i in 1:100) BIC.mat[i,] = as.numeric(eval.list[[i]][which.min(eval.list[[i]][,3]),-(1:3)])
zsep = data.frame(rbind(round(apply(BIC.mat,2,mean),3),
                        round(apply(BIC.mat,2,sd),3)))
names(zsep) = c("B.TP","B.TN","B.relFnorm","Theta.TP","Theta.TN","Theta.relFnorm")
row.names(zsep) = c("mean.BIC","sd.BIC")
zsep

