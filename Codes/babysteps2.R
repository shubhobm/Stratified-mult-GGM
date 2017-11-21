rm(list=ls())
setwd('D:/Study/My projects/Stratified-mult-GGM/Codes')
source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

library(glasso)

##### Generate data
group = rbind(
  c(1, 2),
  c(1, 4),
  c(3, 2),
  c(3, 4),
  c(5, 2),
  c(5, 4),
  c(6, 2),
  c(6, 4),
  c(7, 2),
  c(7, 4)
)                         # grouping pattern
subnetSize = c(10, 10)    # subnet size
n = 100
p = 20
q = 20
K = 10

## generate the two layers
set.seed(11192017)
X.layer = GenerateLayer(n, subnetSize, group)
E.layer = GenerateLayer(n, subnetSize, group)

## generate group structure for coef array
B0.group.array = array(0, c(p,q,K))
g = 1
for(i in 1:p){
  for(j in 1:q){
    B0.group.array[i,j,] = g
    g = g+1
  }
}
B0.array = CoefArray(B0.group.array)
Theta0.array = array(0, c(q,q,K))
for(k in 1:K){
  Theta0.array[,,k] = with(E.layer,
                           diag(diag(Omega[[k]])^(-0.5)) %*% Omega[[k]] %*% diag(diag(Omega[[k]])^(-0.5)))
}

## make Y-layer
Y.layer = E.layer
for(k in 1:K){
  Y.layer$data[[k]] = X.layer$data[[k]] %*% B0.array[,,k] + E.layer$data[[k]]
}

##### Given: X.list, Y.list, B.groups, Theta.groups
Y.list = lapply(Y.layer$data, as.matrix)
Y.indices = Y.layer$indices
Theta.groups = Y.layer$groups
X.list = lapply(X.layer$data, as.matrix)

Theta.group.array = array(0, c(q,q,K))
for(j in 1:q){
  Theta.group.array[j,-j,] = Y.layer$groups[[j]]
}

## tune JMLE model
source('JMLE.R')
lambda.vec = sqrt(log(p*q)/n) * seq(1, 0.1, -0.1)
bic.vec = rep(0, length(lambda.vec))
model.list = list()

system.time(
  for(m in 1:length(lambda.vec)){
    jmle.model = jmle(Y.list, Y.indices, X.list, B.group.array=B0.group.array, Theta.groups=Theta.groups,
                      lambda = lambda.vec[m],
                      gamma = sqrt(log(q)/n) * seq(1, 0.1, -0.1),
                      init.option=1, tol=1e-3)
    
    ## calculate BIC
    jmle.bic.vec = rep(0,K)
    for(k in 1:K){
      nk = nrow(Y.list[[k]])
      Theta.k = jmle.model$Theta_refit$Theta[[k]]
      for(j in 1:q)
      {
        Theta.k[j,j] = 0
      }
      jmle.bic.vec[k] = sum(diag(crossprod((Y.list[[k]] - X.list[[k]] %*%
                                              jmle.model$B.refit[,,k]) %*% (diag(1,q) - Theta.k))))/nk +
        log(nk)/nk * (sum(Theta.k != 0)/2 + sum(jmle.model$B.refit[,,k] != 0))
    }
    bic.vec[m] = sum(jmle.bic.vec)
    model.list[[m]] = jmle.model
  }
)

best.model = model.list[[which.min(bic.vec)]]

c(sum(B0.array != 0 & best.model$B.refit != 0)/sum(B0.array != 0),
  sum(B0.array == 0 & best.model$B.refit == 0)/sum(B0.array == 0),
  sqrt(sum((B0.array - best.model$B.refit)^2)))

Theta_new.array = array(0, c(q,q,K))
for(k in 1:K){
  Theta_new.array[,,k] = best.model$Theta_refit$Theta[[k]]
}
c(sum(Theta0.array != 0 & Theta_new.array != 0)/sum(Theta0.array != 0),
  sum(Theta0.array == 0 & Theta_new.array == 0)/sum(Theta0.array == 0),
  sqrt(sum((Theta0.array - Theta_new.array)^2)))

## separate estimation
setwd('JMLR_revision_code')
source('l1ML_Main.R')
B_lin.array = array(0, c(p,q,K))
Theta_lin.array = array(0, c(q,q,K))

for(k in 1:K){
  model.k = l1ML_Main(Y.list[[k]], X.list[[k]],
                      lambda=.5*sqrt(log(p*q)/n),
                      rho=.5*sqrt(log(q)/n),
                      initializer="Lasso", StabilizeTheta=FALSE)
  B_lin.array[,,k] = model.k$B.est
  Theta_lin.array[,,k] = model.k$Theta.est
}

c(sum(B0.array != 0 & B_lin.array != 0)/sum(B0.array != 0),
  sum(B0.array == 0 & B_lin.array == 0)/sum(B0.array == 0),
  sqrt(sum((B0.array - B_lin.array)^2)))

c(sum(Theta0.array != 0 & Theta_lin.array != 0)/sum(Theta0.array != 0),
  sum(Theta0.array == 0 & Theta_lin.array == 0)/sum(Theta0.array == 0),
  sqrt(sum((Theta0.array - Theta_lin.array)^2)))
