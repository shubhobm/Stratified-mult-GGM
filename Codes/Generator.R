## this file stores the functions that are used to generate data 
require(matrixcalc)
require(MASS)
require(gdata)
library(Matrix)
library(igraph)
library(RANN)
library(grpreg)

#***********************************************************#
# Generates the p times q times K array of regression coefficients
#***********************************************************#
CoefArray = function(B.group.array, sparsity=NULL, SNR=NULL){

  # default sparsity is 5/p (model 1 in Cai, model 2 in Cai has a sparsity of 30/p)
  arraydims = dim(B.group.array)
	if (is.null(sparsity))
		sparsity = 5/arraydims[1]
	
	# default SNR is 1
	if (is.null(SNR))
	  SNR = 1
	
  unique.elems = unique(as.numeric(B.group.array))
  signal.groups = unique.elems[which(rbinom(length(unique.elems), 1, sparsity)==1)] # randomly assign non-zero groups
  B.array = array(0, arraydims)
  for(h in signal.groups){ # generate entries in non-zero groups
    h.indices = which(B.group.array==h, arr.ind=T)
    B.array[h.indices] = sample(c(-1,1),nrow(h.indices), replace=T)*runif(nrow(h.indices),0.5,SNR)
  }
  
  return(B.array)
}

#***********************************************************#
# generates a single layer of data: a K-length list of matrices
#***********************************************************#
GenerateLayer = function(n, subnetSize, group, m=2, rho=0, rho.joint=0.01){
  
  p =  sum(subnetSize)      # number of variables
  rho = 0                   # misspecification ratio
  K = dim(group)[1]         # number of models/networks 
  
  ## Generate the sparsity pattern for all variables
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
  
  ## Generate subnetworks with different structures
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
  
  ## Generate edge weights
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
  
  ## y is the class indicator
  y <- vector("list", K)
  for (k in 1:K){
    y[[k]] <- rep(k, n)  
  }
  # trainY <- unlist(y)        
  
  x <- vector("list", K)
  for (k in 1:K){
    x[[k]] <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma[[k]])
    x[[k]] <- as.matrix(scale(x[[k]], center = T, scale = T))
  }
  # trainX <- do.call(rbind, x) # training data
  return(list(data=x, indices=y, groups=ix, Omega=Omega))
}