# initial work
# setwd("D:/Study/My projects/Stratified-mult-GGM")
setwd("C:/Study/Stratified-mult-GGM/")
Required.Packages <- c("data.table", "glmnet")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

# Sample groups
sample.groups = fread("Data/samplegroup.csv")

# X and Y data
# The top layer is the mRNA and the bottom RNAseq
X.data = fread("Data/BRCA_mRNA_Array.csv")
Y.data = fread("Data/BRCA_RNASeq.csv")

## check all samples are in X and Y data
#sum(names(X.data)[-1] %in% sample.groups$CASE_ID)
#sum(names(Y.data)[-1] %in% sample.groups$CASE_ID)
# they are.

## assign groups to each row
groups = fread("Data/pathway.csv")[-1]
X.data[, X.groups := "0"]
X.genes = gsub("\\_.*","",X.data[,V1])
for(i in 1:ncol(groups)){
  which.i = which(X.genes %in% sapply(groups[,..i], paste))
  X.data$X.groups[which.i] = paste(X.data$X.groups[which.i], i, sep=".")
}

Y.data[, Y.groups := "0"]
Y.genes = gsub("\\_.*","",Y.data[,V1])
for(i in 1:ncol(groups)){
  which.i = which(Y.genes %in% sapply(groups[,..i], paste))
  Y.data$Y.groups[which.i] = paste(Y.data$Y.groups[which.i], i, sep=".")
}

## overlapping groups :(
# choose genes in only 1 group
length(which(sapply(X.data$X.groups,nchar) %in% 3:4)) # 2424 such genes
X.data1 = X.data[sapply(X.groups,nchar) %in% 3:4]
# X.data1[1:100, X.groups]

length(which(sapply(Y.data$Y.groups,nchar) %in% 3:4)) # 2775 such genes
Y.data1 = Y.data[sapply(Y.groups,nchar) %in% 3:4]
# Y.data1[1:100, Y.groups]

# now make lists
plus.cases = sample.groups[ER=="Positive", CASE_ID]
minus.cases = sample.groups[ER=="Negative", CASE_ID]

X.list = vector("list", 2)
X.list[[1]] = X.data1[,.SD, .SDcols=c(1,which(names(X.data1) %in% plus.cases))]
X.list[[2]] = X.data1[,.SD, .SDcols=c(1,which(names(X.data1) %in% minus.cases))]

Y.list = vector("list", 2)
Y.list[[1]] = Y.data1[,.SD, .SDcols=c(1,which(names(Y.data1) %in% plus.cases))]
Y.list[[2]] = Y.data1[,.SD, .SDcols=c(1,which(names(Y.data1) %in% minus.cases))]

# now do marginal selection in X and Y to reduce dimensions
Y1 = t(log(Y.list[[1]][,-1]+1))
colnames(Y1) = Y.list[[1]][,V1]
Y2 = t(log(Y.list[[2]][,-1]+1))
colnames(Y2) = Y.list[[2]][,V1]
# 
# # take only top 100 for each Y, then their union
# top.inds = union(order(apply(Y1,2,mad))[1:100], order(apply(Y2,2,mad))[1:100])
# table(Y.data1[top.inds, Y.groups])
# Y.list1 = list(Y1[, top.inds], Y2[,top.inds])
# # selects a lot of 76 group

# use l1_LS
Rsq1.vec = rep(0, ncol(Y1))
lambda = .01
X = t(X.list[[1]][,-1])
Xsd = apply(X,2,sd)
pb = txtProgressBar(0, ncol(Y1))
B1 = matrix(0, ncol(X), ncol(Y1))

for(j in 1:ncol(Y1)){
  Yj = Y1[,j]
  if(var(Yj)>0){
    temp = glmnet(X,Yj,intercept=FALSE)
    B1[,j] = predict(temp, s=lambda, type="coefficients")[-1]
    Rsq1.vec[j] = 1 - var(Yj - X %*% B1[,j])/var(Yj)
  }
  B1[,j] = B1[,j]*Xsd # standardize
  setTxtProgressBar(pb,j)
}
close(pb)

Rsq2.vec = rep(0, ncol(Y2))
X = t(X.list[[2]][,-1])
Xsd = apply(X,2,sd)
pb = txtProgressBar(0, ncol(Y2))
B2 = matrix(0, ncol(X), ncol(Y2))

for(j in 1:ncol(Y2)){
  Yj = Y2[,j]
  if(var(Yj)>0){
    temp = glmnet(X,Yj,intercept=FALSE)
    B2[,j] = predict(temp, s=lambda, type="coefficients")[-1]
    Rsq2.vec[j] = 1 - var(Yj - X %*% B2[,j])/var(Yj)
  }
  B2[,j] = B2[,j]*Xsd # standardize
  setTxtProgressBar(pb,j)
}
close(pb)

# take only top 100 for each Y, then their union
top.inds.Y = union(order(Rsq1.vec)[1:100], order(Rsq2.vec)[1:100])
table(Y.data1[top.inds.Y, Y.groups])
# more balanced
Y.list1 = list(Y1[, top.inds.Y], Y2[,top.inds.Y])

# change X matrix structures
X1 = t(X.list[[1]][,-1])
colnames(X1) = X.list[[1]][,V1]
X2 = t(X.list[[2]][,-1])
colnames(X2) = X.list[[2]][,V1]

# take top 200 for each X with highest median coefficients across top.inds of Y, then their union
top.inds = union(order(apply(X1,2,mad))[1:300], order(apply(X2,2,mad))[1:300])
table(X.data1[top.inds, X.groups])
# top.inds.X = union(order(apply(B1[,top.inds.Y],1,function(x) mean(abs(x))))[1:200],
#                    order(apply(B2[,top.inds.Y],1,function(x) mean(abs(x))))[1:200])
top.inds.X = union(order(apply(B1,1,function(x) mean(abs(x))))[1:200],
                   order(apply(B2,1,function(x) mean(abs(x))))[1:200])
table(X.data1[top.inds.X, X.groups])
X.list1 = list(X1[, top.inds.X], X2[,top.inds.X])

# final groups
n = sapply(X.list1, nrow)
p = ncol(X.list1[[1]])
q = ncol(Y.list1[[1]])
K = 2

# Beta groups
B.group.array = array("", c(p,q,2))
Xg = X.data1[top.inds.X, X.groups]
Yg = Y.data1[top.inds.Y, Y.groups]

# if both genes are in same group, then that group, else paste
for(i in 1:p){
  for(j in 1:q){
    if(Xg[i]==Yg[j])
    {
      B.group.array[i,j,] = Xg[i]
    } else{
      B.group.array[i,j,] = paste(Xg[i], Yg[j], sep=".")
    }
  }
}

# Theta groups
Theta.groups = vector("list", q)
for (i in 1:q){
  Theta.groups[[i]] = matrix(0, K, q)
  for (k in 1:K){
    Theta.groups[[i]][k,] = match(Yg, unique(Yg))
  }
  Theta.groups[[i]] = Theta.groups[[i]][,-i]
}

saveRDS(list(X.list1=X.list1, Y.list1=Y.list1,
             Xg=Xg, Yg=Yg,
             B.group.array=B.group.array, Theta.groups=Theta.groups),
        file="./Data/processed_data.rds")
