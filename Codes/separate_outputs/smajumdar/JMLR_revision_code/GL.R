# how glasso works for this 
GL = function(X,Y,rho.GL){
	# rho.GL is the glasso penalty
	p1 = ncol(X); p2 = ncol(Y); n = nrow(X);
	fulldata = cbind(X,Y);
	GL.est = huge(fulldata,lambda=rho.GL,method="glasso",verbose=FALSE)$icov[[1]];
	Theta_XX = GL.est[1:p1,1:p1];
	Theta_XY = GL.est[1:p1,(p1+1):(p1+p2)];
	Theta_YY = as.matrix(GL.est[(p1+1):(p1+p2),(p1+1):(p1+p2)]);
	Theta_YY = ifelse(abs(Theta_YY)>1e-5,Theta_YY,0)
	B_calculated = - Theta_XY %*% solve(Theta_YY);
	return(list(Theta_XX = Theta_XX,Theta_XY=Theta_XY,Theta_YY=Theta_YY,Theta_full=GL.est,B_calculated=B_calculated));
}