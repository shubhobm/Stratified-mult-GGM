###############################################
# Example
## K = 10, n = 50, p = 50

#Set the seed
set.seed(2015)

#Load packages
library(MASS)
library(Matrix)
library(igraph)
library(glasso)
library(gdata)
library(RANN)
library(grpreg)


#general parameters
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
p =  sum(subnetSize)      # number of variables
rho = 0                   # misspecification ratio
K = dim(group)[1]         # number of models/networks 
n = 100                   # sample size for each model.
m = 2                     # m determines the density of the scale-free network.
nreps = 50                # number of replication
rho.joint = 0.01          # tuning parameter for 2nd setp in JGrpL

##-----------------------------------------
#1-1) Generate the sparsity pattern for all variables
##----------------------------------------- 
ix = vector("list", p)
for (i in 1:p){
  ix[[i]] = matrix(0, K, p)
  for (j in 1:p){
    if (i <= subnetSize[1] || j <= subnetSize[1]) {
      colj = match(group[,1], unique(group[,1]))
    } else {
      colj = match(group[,2], unique(group[,2]))
    }
    ix[[i]][, j] = colj
  }
  ix[[i]] = ix[[i]][, -i]
}

##-----------------------------------------
#1-2) Generate subnetworks with different structures
##-----------------------------------------
nSet = length(unique(c(group)))

subnet.adj = vector("list", nSet)
for (i in unique(group[,1])){
  subnet.adj[[i]] = sf.net(p,  m)$A
  subnet = graph.adjacency(subnet.adj[[i]], mode="undirected")
  neworder = rank(-degree(subnet), ties.method = "first")
  subnet.adj[[i]] = subnet.adj[[i]][neworder, neworder]
}
for (i in unique(group[, 2])){
  subnet.adj[[i]] = sf.net(subnetSize[2], m)$A
  subnet = graph.adjacency(subnet.adj[[i]], mode="undirected")
  neworder = rank(-degree(subnet), ties.method = "first")
  subnet.adj[[i]] = subnet.adj[[i]][neworder, neworder]
}

Amat = vector("list", K)
for (k in 1:K){
  Amat[[k]] = subnet.adj[[group[k, 1]]]
  Amat[[k]][(subnetSize[1] + 1):p, (subnetSize[1] + 1):p] = subnet.adj[[group[k, 2]]]
}

##-----------------------------------------
#(3) Generate edge weights
Omega <- vector("list", K)
Sigma <- vector("list", K)
Ip <- diag(1,p)
for (k in 1:K){
  Amat.g <- graph.adjacency(Amat[[k]], mode = "undirected")
  Amat.el <- get.edgelist(Amat.g)
  
  ## Find the complement of Amat[[k]]
  Amat.c <- matrix(1, p, p) - Amat[[k]] - Ip
  
  Amat.c.g <- graph.adjacency(Amat.c, mode = "undirected")
  Amat.c.el <- get.edgelist(Amat.c.g)
  
  ## sample rho percent of the variables from the complementary graph
  add.el <- Amat.c.el[sample(seq(1, dim(Amat.c.el)[1]), rho*dim(Amat.el)[1]) , ]
  tmp.g <- graph.edgelist(rbind(add.el, Amat.el), directed = F)
  
  Amat[[k]] <- as.matrix(get.adjacency(tmp.g, type="both"))
  
  weights <- matrix(0, p, p)
  upperTriangle(weights, diag = F) <- runif((p*(p - 1))/2, 0.5, 1)*(2*rbinom((p*(p - 1))/2, 1, 0.5) - 1)
  weights <- weights + t(weights)
  
  pdobj <- pd(weights * Amat[[k]])
  Omega[[k]] <- pdobj$A
  Sigma[[k]] <- pdobj$Ainv
}

##---------------------------------------
# y is the class indicator
##-----------------------------------------
y <- vector("list", K)
for (k in 1:K){
  y[[k]] <- rep(k, n)  
}
trainY <- unlist(y)        

x <- vector("list", K)
for (k in 1:K){
  x[[k]] <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma[[k]])
  x[[k]] <- scale(x[[k]], center = T, scale = T)
}
trainX <- do.call(rbind, x) # training data

jsem.grid <- sqrt(log(p)/n) * seq(1, 0.1, -0.1)
bic.jsem <- sel.lambda.jsem(trainX, trainX, trainY, trainY, ix, lambda = jsem.grid)
lambda.jsem <- jsem.grid[which.min(bic.jsem$BIC)]
fit.joint <- JSEM(trainX, trainY, ix, lambda = lambda.jsem)


