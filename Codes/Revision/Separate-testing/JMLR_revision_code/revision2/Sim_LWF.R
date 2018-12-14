library(doMC)
registerDoMC()
options(cores=2)

setwd("Desktop/Sparsity_SIM")
source("Generator_new2.R")
source("mixMRF.R")
source("l1ML_Main.R")
source("Eval.R")

p1 = 30; p2 = 60; n = 100; 

ThetaX.info = ThetaMat(p1,"band",Theta1=0.2); ThetaX = ThetaX.info$Theta; 
ThetaY.info = ThetaMat(p2,"band",Theta1=0.3); ThetaY = ThetaY.info$Theta;
ThetaXY = ThetaXYMat(p1,p2);

B = -ThetaXY %*% solve(ThetaY)

SigmaGiant.info = SigmaMatGiant(ThetaX,ThetaXY,ThetaY);
SigmaGiant = SigmaGiant.info$SigmaGiant; 
eigen(SigmaGiant)$values


Fulldata = gen.2layer(n,SigmaGiant,p1,p2);
X = Fulldata$X; Y = Fulldata$Y;

# =================================================================================
lambda1.seq = seq(from=0.1,to=1,length=5)*sqrt(log(p1)/n);
lambda2.seq = seq(from=0.001,to=0.005,length=5)*sqrt(log(p2)/n);

EvalThetaXY = array(NA,c(length(lambda1.seq),length(lambda2.seq),4)); 
dimnames(EvalThetaXY)[[3]] = c("SEN","SPC","MCC","FDR");
EvalThetaY = array(NA,c(length(lambda1.seq),length(lambda2.seq),4));
dimnames(EvalThetaY)[[3]] = c("SEN","SPC","MCC","FDR");

for (i in 1:length(lambda1.seq)){
	for (j in 1:length(lambda2.seq)){
		cat("i=",i,"j=",j,"\n")
		result = mixMRF(X,Y,lambda1=lambda1.seq[i],lambda2=lambda2.seq[j],scaledLasso=FALSE,rule="or")
		EvalB[i,j,] = as.numeric(Eval(ThetaXY,result$ThetaXY,directed=TRUE)[c(1,2,3,4)]);
		EvalThetaY[i,j,] = as.numeric(Eval(ThetaY,result$ThetaY,directed=FALSE)[c(1,2,3,4)]);
	}
}
EvalThetaXY;
EvalThetaY;

# ===================================================================================
# this is using our estimation procedure
lambda.seq = seq(from=2,to=4,length=5)*sqrt(log(p1)/n);
rho.seq = seq(from=1.5,to=3,length=5)*sqrt(log(p2)/n)
EvalThetaXY2 = array(NA,c(length(lambda.seq),length(rho.seq),4)); 
dimnames(EvalB)[[3]] = c("SEN","SPC","MCC","FDR");
EvalTheta2 = array(NA,c(length(lambda.seq),length(rho.seq),4));
dimnames(EvalTheta)[[3]] = c("SEN","SPC","MCC","FDR");

for (i in 1:length(lambda.seq)){
	for (j in 1:length(rho.seq)){
		BTheta = l1ML_Main(Y,X,initializer="Lasso",lambda=lambda.seq[i],rho=rho.seq[j],Screening=T,StabilizeTheta=T);	
		B.est = BTheta$B.est;
		Theta.est = BTheta$Theta.est;
		ThetaXY.est = -B %*% Theta.est;
		EvalThetaXY2[i,j,] = as.numeric(Eval(ThetaXY,ThetaXY.est,directed=TRUE)[c(1,2,3,4)]);
		EvalTheta2[i,j,] = as.numeric(Eval(ThetaY,Theta.est,directed=FALSE)[c(1,2,3,4)]);
	}
}
EvalThetaXY2; 
EvalTheta2;