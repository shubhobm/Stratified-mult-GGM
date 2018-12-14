rm(list=ls())
library(doMC)
registerDoMC()
options(cores=2)

setwd("Desktop/Sparsity_SIM")
source("Generator_new.R")
source("mixMRF_OLS.R")
source("Eval.R")
library(MASS)

p1 = 10; p2 = 10; n = 3000; 

#************************
#
## Generate parameters:
#
#************************
standardizeE = TRUE;
GenerateData = FALSE; 

if (GenerateData){
## generate B matrix 
B = RegMat(p1,p2,sparsity=2/p1,SNR=0.3); 

## generate ThetaX and ThetaY
 
ThetaX.info = ThetaMat(p1,"band",Theta1=0.2); ThetaY.info = ThetaMat(p2,"band",Theta1=0.3); 
ThetaX = ThetaX.info$Theta; ThetaY = ThetaY.info$Theta;

## solve for Sigma_X matrix (marginally), get Theta_XY and Theta_giant
par_list = solve_for_pars(B,ThetaX,ThetaY);
SigmaX = par_list$SigmaX; Theta_XY = par_list$Theta_XY; ThetaY = par_list$ThetaY
Theta_joint = par_list$Theta_giant;

# check the condition number
max(eigen(ThetaX)$values)/min(eigen(ThetaX)$values);
max(eigen(ThetaY)$values)/min(eigen(ThetaY)$values);
max(eigen(Theta_joint)$values)/min(eigen(Theta_joint)$values);
(sum(!!Theta_joint)-(p1+p2))/((p1+p2)*(p1+p2-1));
}
#paras = list(B=B,ThetaX=ThetaX,ThetaY=ThetaY,Theta_joint=Theta_joint,SigmaX=SigmaX)
#save(paras,file="goodboy.Rda")

load("goodboy.Rda"); # set of parameters that have good properties
B = paras$B; SigmaX = paras$SigmaX; ThetaY = paras$ThetaY; 
Theta_XY = paras$Theta_joint[1:p1,(p1+1):(p1+p2)]

## generate the data based on the conditional model:
Fulldata = gen.2layer(n,B,ThetaY,SigmaX,scale=FALSE); 
X = Fulldata$X; Y = Fulldata$Y;

# Theta estimate through inverting the sample covariance matrix:
Theta.inv = solve(cov(cbind(X,Y))); 
Theta.inv = ifelse(abs(Theta.inv)>1e-5,Theta.inv,0);
ThetaXY.inv = Theta.inv[1:p1,(p1+1):(p1+p2)]; ThetaY.inv = Theta.inv[(p1+1):(p1+p2),(p1+1):(p1+p2)];
B_inv = -ThetaXY.inv %*% solve(ThetaY.inv);
B_inv_cheated = ifelse(abs(B_inv)>0.1,B_inv,0)   # this is a cheated estimate of B

Eval(B,B_inv,directed=T)
Eval(B,B_inv_cheated,directed=T)

# Theta estimate through OLS:
B_OLS = mixMRF_OLS(X,Y,ht=T);
B_OLS_noht = mixMRF_OLS(X,Y,ht=F);
Eval(B,B_OLS$B_calculated,directed=T)
Eval(B,B_OLS_noht$B_calculated,directed=T)

norm(B_OLS_noht$B_calculated-B_inv,"f")/norm(B_OLS$B_calculated,"f")
norm(B_OLS$B_calculated-B_inv_cheated,"f")/norm(B_OLS$B_calculated,"f")






