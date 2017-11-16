Eval=function(true, est, directed=FALSE){
        if (directed==FALSE){
                true.v=upperTriangle(!!true, diag=FALSE)
                est.v=upperTriangle(!!est, diag=FALSE)
        }
        if (directed==TRUE){
                true.v=as.vector(!!true)
                est.v=as.vector(!!est)
        }

        TP = sum(est.v & true.v)
        TN = sum((!est.v) & (!true.v))
        FP = sum(est.v & !true.v)
        FN = sum(!est.v & true.v)

        FDR = FP/max(1,sum(est.v))

        TTP = sum(true.v)
        TTN = sum(!true.v)

        SEN = TP/TTP
        SPE = TN/TTN

	DE_MCC = sqrt((TP+FP))*sqrt((TP+FN))*sqrt((TN+FP))*sqrt((TN+FN))
	MCC = (TP*TN-FP*FN)/DE_MCC
	
	Abs.Error = norm(est-true,"F")
        Rel.Error = Abs.Error/norm(true,"F")

        result = data.frame(SEN,SPE,MCC,FDR,Abs.Error,Rel.Error)

        return(result)
}

