DEBUG = TRUE
library(doMC)
registerDoMC()
options(cores=2)

setwd("Desktop/Sparsity_SIM")
source("Generator.R")
source("mixMRF.R")
source("Eval.R")

set.seed(612);

p1 = 200; p2 = 200; n = 150; 

## generate Theta matrix and B matrix 
B = RegMat(p1,p2); 
ThetaY.info = ThetaMat(p2); ThetaY = ThetaY.info$Theta; ThetaY_CN = ThetaY.info$CN_Theta;

Fulldata = gen.2layer(n,B,ThetaY);
X = Fulldata$X; Y = Fulldata$Y;
		
lambda1.seq = seq(from=1,to=3,length=5)*sqrt(log(p1)/n);
lambda2.seq = seq(from=0.1,to=1,length=5)*sqrt(log(p2)/n);

EvalB = array(NA,c(length(lambda1.seq),length(lambda2.seq),3)); 
dimnames(EvalB)[[3]] = c("SEN","SPC","MCC");
EvalThetaY = array(NA,c(length(lambda1.seq),length(lambda2.seq),3));
dimnames(EvalThetaY)[[3]] = c("SEN","SPC","MCC");

for (i in 1:length(lambda1.seq)){
	for (j in 1:length(lambda2.seq)){
		cat("i=",i,"j=",j,"\n")
		result = mixMRF(X,Y,lambda1=lambda1.seq[i],lambda2=lambda2.seq[j],scaledLasso=FALSE,rule="or")
		EvalB[i,j,] = as.numeric(Eval(B,result$B_calculated,directed=TRUE)[c(1,2,3)]);
		EvalThetaY[i,j,] = as.numeric(Eval(ThetaY,result$ThetaY,directed=FALSE)[c(1,2,3)]);
	}
}
EvalB;
# EvalThetaY;