library(doMC)
registerDoMC()
options(cores=2)

setwd("Desktop/Sparsity_SIM")
source("Generator.R")
source("l1ML_Main.R")
source("Eval.R")


p1 = 20; p2 = 20; n = 500; 

#************************
#
## Generate parameters:
#
#************************

## generate a giant Theta matrix
Theta_giant.info = ThetaMat(p1+p2,type="band",Theta1=0.2,Theta2=0.1); Theta_joint= Theta_giant.info$Theta;

# extract the corresponding block
Theta_XY = Theta_joint[1:p1,(p1+1):(p1+p2)];
ThetaY = Theta_joint[(p1+1):(p1+p2),(p1+1):(p1+p2)];

# get the B matrix:
B = - Theta_XY %*% solve(ThetaY)


# check the condition number
max(eigen(ThetaX)$values)/min(eigen(ThetaX)$values);
max(eigen(ThetaY)$values)/min(eigen(ThetaY)$values);
max(eigen(Theta_joint)$values)/min(eigen(Theta_joint)$values);
(sum(!!Theta_joint)-(p1+p2))/((p1+p2)*(p1+p2-1));

## generate the data based on the conditional model:
Fulldata = mvrnorm(n,rep(0,p1+p2),solve(Theta_joint))
X = Fulldata[,1:p1] Y = Fulldata[,(p1+1):(p1+p2)];

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

## this array records the evaluation result for graphical Lasso applied to data generated jointly
EvalGLjoint = array(NA,c(length(rho.GL.seq),3,4)); 
dimnames(EvalGLjoint)[[2]] = c("SEN","SPC","Rel.Ferror");
dimnames(EvalGLjoint)[[3]] = c("ThetaX","ThetaXY","ThetaY","GiantTheta")

for (i in 1:length(rho.GL.seq)){
	
	## using the data that is generated based on the conditional mode:
	myGL = GL(X,Y,rho.GL.seq[i]);
	Theta_XX.est = myGL$Theta_XX; Theta_XY.est = myGL$Theta_XY; Theta_YY.est = myGL$Theta_YY; 
	Theta.full.est = myGL$Theta_full;
	
	EvalGL[i,,1] = as.numeric(Eval(ThetaX,Theta_XX.est,directed=FALSE)[c(1,2,6)]);
	EvalGL[i,,2] = as.numeric(Eval(Theta_XY,Theta_XY.est,directed=TRUE)[c(1,2,6)]);
	EvalGL[i,,3] = as.numeric(Eval(ThetaY,Theta_YY.est,directed=FALSE)[c(1,2,6)]);
	EvalGL[i,,4] = as.numeric(Eval(Theta_joint,Theta.full.est,directed=FALSE)[c(1,2,6)])
	
	## using the joint data:
	GL_on_joint = huge(XYjoint,lambda=rho.GL.seq[i],method='glasso',verbose=F)$icov[[1]]
	
	EvalGLjoint[i,,1] = as.numeric(Eval(ThetaX,GL_on_joint[1:p1,1:p1],directed=FALSE)[c(1,2,6)]);
	EvalGLjoint[i,,2] = as.numeric(Eval(Theta_XY,GL_on_joint[1:p1,(p1+1):(p1+p2)],directed=TRUE)[c(1,2,6)]);
	EvalGLjoint[i,,3] = as.numeric(Eval(ThetaY,GL_on_joint[(p1+1):(p1+p2),(p1+1):(p1+p2)],directed=FALSE)[c(1,2,6)]);
	EvalGLjoint[i,,4] = as.numeric(Eval(Theta_joint,GL_on_joint,directed=FALSE)[c(1,2,6)])
}
EvalGL
EvalGLjoint


## the following lines of code just to let us see if we can "trust" glasso
E = Y - X %*% B;
evalE = array(NA,c(10,3)); rho.E = seq(from=0.5,to=2,length=10)*sqrt(log(p2)/n);
for (k in 1:10){
	E.est = huge(E,lambda=rho.E[k],method="glasso",verbose=FALSE)$icov[[1]]
	evalE[k,] = as.numeric(Eval(ThetaY,E.est,directed=FALSE)[c(1,2,6)])
}
evalE



		

# this is the 		
# result = l1ML_Main(Y,X,initializer="Lasso",lambda=lambda.opt,rho=rho.opt);