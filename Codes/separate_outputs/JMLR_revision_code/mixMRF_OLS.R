mixMRF_OLS = function(X,Y,ht=TRUE,alpha=0.1,rule="or"){
	p1 = ncol(X); p2 = ncol(Y); n = nrow(X);
	# regress each Y[,j] on cbind(X,Y[,-j]);
	ThetaXYY = array(NA,c(p1+p2,p2)); 
	if (!ht){	
		for (j in 1:p2){
			response = Y[,j];
			regressor = cbind(X,Y[,-j])
			OLS = lm(response~regressor+0);
			coef.est = OLS$coefficients;
			noise.est = sum(OLS$residuals^2)/(n-(p1+p2-1));
		
			# diagonal
			omega_jj = 1/noise.est^2; 
			# off diagonal
			Omega_j_but_j = -coef.est * omega_jj;
			
			ThetaXYY[-(p1+j),j] = Omega_j_but_j; ThetaXYY[(p1+j),j] = omega_jj
		}
	}
	if (ht){
		for (j in 1:p2){
			response = Y[,j];
			regressor = cbind(X,Y[,-j]);
			OLS = lm(response~regressor+0);
			coef.est = OLS$coefficients;
			coef.skeleton = which(summary(OLS)$coefficients[,4]<alpha);
			
			if (length(coef.skeleton)==0){
				coef.est = 0;
				noise.est = sum(response^2)/n
			}
			else{
			# redo OLS:
			redo.result = lm(response~regressor[,coef.skeleton]+0);
			noise.est = sum(redo.result$residuals^2)/(n-length(coef.skeleton));
			coef.est[coef.skeleton] = redo.result$coefficients;
			coef.est[-coef.skeleton] = 0;
			}
			
			# diagonal
			omega_jj = 1/noise.est^2; 
			# off diagonal
			Omega_j_but_j = -coef.est * omega_jj;
			
			ThetaXYY[-(p1+j),j] = Omega_j_but_j; ThetaXYY[(p1+j),j] = omega_jj
			
		}
	}
	
	ThetaY = ThetaXYY[(p1+1):(p1+p2),]; ThetaXY = ThetaXYY[1:p1,]
	if (rule=="or"){
		for (i in 1:(p2-1)){
			for (j in (i+1):p2){
				if (abs(ThetaY[i,j]) == 0 | abs(ThetaY[j,i])==0)
					ThetaY[i,j] = ThetaY[j,i] = 0;  
			}
		}
	}
	B_calculated = -ThetaXY %*% solve(ThetaY);
	
	return(list(ThetaXY = ThetaXY, ThetaY = ThetaY, B_calculated = B_calculated))
}