library(doMC)
registerDoMC()
options(cores=2)

setwd("Desktop/Sparsity_SIM")
source("Generator_new.R")
source("l1ML_Main.R")
source("Eval.R")


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

load("goodboy.Rda");
B = paras$B; SigmaX = paras$SigmaX; ThetaY = paras$ThetaY; 
Theta_XY = paras$Theta_joint[1:p1,(p1+1):(p1+p2)]

## generate the data based on the conditional model:
Fulldata = gen.2layer(n,B,ThetaY,SigmaX,scale=FALSE); 
X = Fulldata$X; Y = Fulldata$Y;

# how glasso works for this 
GL = function(X,Y,rho.GL){
	# rho.GL is the glasso penalty
	p1 = ncol(X); p2 = ncol(Y); n = nrow(X);
	fulldata = cbind(X,Y);
	GL.est = huge(fulldata,lambda=rho.GL,method="glasso",verbose=FALSE)$icov[[1]];
	Theta_XX = GL.est[1:p1,1:p1];
	Theta_XY = GL.est[1:p1,(p1+1):(p1+p2)];
	Theta_YY = GL.est[(p1+1):(p1+p2),(p1+1):(p1+p2)];
	return(list(Theta_XX = Theta_XX,Theta_XY=Theta_XY,Theta_YY=Theta_YY,Theta_full=GL.est));
}

## first need to check how GL works 
rho.GL.seq = seq(from=0.2, to=1.5,length=10)*sqrt(log(p1+p2)/n);

## this array records the evaluation result for graphical Lasso applied to data generated conditionally
EvalGL = array(NA,c(length(rho.GL.seq),3,4)); 
dimnames(EvalGL)[[2]] = c("SEN","SPC","Rel.Ferror");
dimnames(EvalGL)[[3]] = c("ThetaX","ThetaXY","ThetaY","GiantTheta")

for (i in 1:length(rho.GL.seq)){
	
	## using the data that is generated based on the conditional mode:
	myGL = GL(X,Y,rho.GL.seq[i]);
	Theta_XX.est = myGL$Theta_XX; Theta_XY.est = myGL$Theta_XY; Theta_YY.est = myGL$Theta_YY; 
	Theta.full.est = myGL$Theta_full;
	
	EvalGL[i,,1] = as.numeric(Eval(ThetaX,Theta_XX.est,directed=FALSE)[c(1,2,6)]);
	EvalGL[i,,2] = as.numeric(Eval(Theta_XY,Theta_XY.est,directed=TRUE)[c(1,2,6)]);
	EvalGL[i,,3] = as.numeric(Eval(ThetaY,Theta_YY.est,directed=FALSE)[c(1,2,6)]);
	EvalGL[i,,4] = as.numeric(Eval(Theta_joint,Theta.full.est,directed=FALSE)[c(1,2,6)])
}
EvalGL

# this is using our estimation procedure
lambda.seq = seq(from=2,to=4,length=5)*sqrt(log(p1)/n);
rho.seq = seq(from=1.5,to=3,length=5)*sqrt(log(p2)/n)
EvalB = array(NA,c(length(lambda.seq),length(rho.seq),3)); 
dimnames(EvalB)[[3]] = c("SEN","SPC","Rel.Ferror");
EvalTheta = array(NA,c(length(lambda.seq),length(rho.seq),3));
dimnames(EvalTheta)[[3]] = c("SEN","SPC","Rel.Ferror");

for (i in 1:length(lambda.seq)){
	for (j in 1:length(rho.seq)){
		BTheta = l1ML_Main(Y,X,initializer="Lasso",lambda=lambda.seq[i],rho=rho.seq[j],Screening=FALSE,StabilizeTheta=FALSE);	
		EvalB[i,j,] = as.numeric(Eval(B,BTheta$B.est,directed=TRUE)[c(1,2,6)]);
		EvalTheta[i,j,] = as.numeric(Eval(ThetaY,BTheta$Theta.est,directed=FALSE)[c(1,2,6)]);
	}
}
EvalB; 
EvalTheta;
 