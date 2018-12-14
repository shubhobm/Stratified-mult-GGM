library(doMC)
registerDoMC()
options(cores=2)

setwd("Desktop/Sparsity_SIM")
source("Generator_new.R")
source("l1ML_Main.R")
source("GL.R")
source("Eval.R")


p1 = 10; p2 = 10; n = 2000; n.iter = 5;

#************************
#
## Generate parameters:
#
#************************
standardizeE = TRUE;
GenerateData = FALSE; 

# ******************************************
#
#  generated model parameters / load generated parameters
#
# ******************************************

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
#paras = list(B=B,ThetaX=ThetaX,ThetaY=ThetaY,Theta_joint=Theta_joint,SigmaX=SigmaX)
#save(paras,file="goodboy.Rda")
}
if (!GenerateData){
	load("goodboy.Rda") 
	# paras is loaded
	B = paras$B; SigmaX = paras$SigmaX; ThetaY = paras$ThetaY; 
	Theta_XY = paras$Theta_joint[1:p1,(p1+1):(p1+p2)];
}


# ********************************
#	set tuning parameters
# ********************************

## tuning parameters for Glasso:
rho.GL.seq = seq(from=0.5, to=3,length=10)*sqrt(log(p1+p2)/n);
# tuning parameters for our procedure
lambda.seq = seq(from=2,to=4,length=5)*sqrt(log(p1)/n);
rho.seq = seq(from=1.5,to=3,length=5)*sqrt(log(p2)/n)

# **************************************
#	set evaluation recording arrays
# **************************************
## this array records the evaluation result for graphical Lasso applied to data generated conditionally
EvalGL = array(NA,c(length(rho.GL.seq),3,4,n.iter)); 
dimnames(EvalGL)[[2]] = c("SEN","SPC","Rel.Ferror");
dimnames(EvalGL)[[3]] = c("ThetaXY","ThetaY","CalculatedB","ThresholdedB")

EvalB = array(NA,c(length(lambda.seq),length(rho.seq),3,n.iter)); 
dimnames(EvalB)[[3]] = c("SEN","SPC","Rel.Ferror");
EvalTheta = array(NA,c(length(lambda.seq),length(rho.seq),3,n.iter));
dimnames(EvalTheta)[[3]] = c("SEN","SPC","Rel.Ferror");



for (iter in 1:n.iter){
	cat("***********************ITER =",iter,"***********************\n")
	## generate the data based on the conditional model:
	Fulldata = gen.2layer(n,B,ThetaY,SigmaX,scale=FALSE); 
	X = Fulldata$X; Y = Fulldata$Y;


	for (k in 1:length(rho.GL.seq)){
	
		## using the data that is generated based on the conditional mode:
		myGL = GL(X,Y,rho.GL.seq[k]);
		Theta_XY.est = myGL$Theta_XY; Theta_YY.est = myGL$Theta_YY; 
		B_calculated = myGL$B_calculated; 
		B_calculated_thr = ifelse(abs(as.matrix(B_calculated))>5e-3,as.matrix(B_calculated),0)
	
		EvalGL[k,,1,iter] = as.numeric(Eval(Theta_XY,Theta_XY.est,directed=TRUE)[c(1,2,6)]);
		EvalGL[k,,2,iter] = as.numeric(Eval(ThetaY,Theta_YY.est,directed=FALSE)[c(1,2,6)]);
		EvalGL[k,,3,iter] = as.numeric(Eval(B,B_calculated,directed=TRUE)[c(1,2,6)]);
		EvalGL[k,,4,iter] = as.numeric(Eval(B,B_calculated_thr,directed=TRUE)[c(1,2,6)]);
	}

	for (i in 1:length(lambda.seq)){
		for (j in 1:length(rho.seq)){
			BTheta = l1ML_Main(Y,X,initializer="Lasso",lambda=lambda.seq[i],rho=rho.seq[j],Screening=FALSE,StabilizeTheta=FALSE);	
			EvalB[i,j,,iter] = as.numeric(Eval(B,BTheta$B.est,directed=TRUE)[c(1,2,6)]);
			EvalTheta[i,j,,iter] = as.numeric(Eval(ThetaY,BTheta$Theta.est,directed=FALSE)[c(1,2,6)]);
		}
	}
}

EvalB.avg = apply(EvalB,c(1,2,3),mean);
EvalTheta.avg = apply(EvalTheta,c(1,2,3),mean);
EvalGL.avg = apply(EvalGL,c(1,2,3),mean)

## consider fixed rho when rho = rho.seq[2]
plot(1-EvalB.avg[,2,"SPC"],EvalB.avg[,2,"SEN"],type="b",cex=0.5,xlim=c(0,1),ylim=c(0,1),xlab="1-specificity",ylab="sensitivity",main="ROC curve",cex.lab=0.7,cex.axis=0.7,cex.main=0.7);
lines(1-EvalTheta.avg[,2,"SPC"],EvalTheta.avg[,2,"SEN"],type="b",cex=0.5,col=2)

lines(1-EvalGL.avg[,2,"CalculatedB"],EvalGL.avg[,1,"CalculatedB"],type="b",cex=0.5,col=3);
lines(1-EvalGL.avg[,2,"ThetaY"],EvalGL.avg[,1,"ThetaY"],type="b",cex=0.5,col=4);

lines(1-EvalGL.avg[,2,"ThetaXY"],EvalGL.avg[,1,"ThetaXY"],type="b",cex=0.5,col=5)















 