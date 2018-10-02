rm(list=ls())
setwd('c:/Study/Stratified-mult-GGM/Codes/Revision/')

analyze = function(list, broken=FALSE, array=NULL, range=NULL){
  matrix.list = list()

  if(broken){
    output = list()
    index = 1
    for(i in array){
      load(paste0(list,"_",i,".Rda"))
      output[[index]] = rbind(apply(simplify2array(out.mat), 1:2, mean),
                          apply(simplify2array(out.mat)/sqrt(5), 1:2, sd))
      index = index+1
    }
    output = round(apply(simplify2array(output), 1:2, mean),3)
  } else{
    load(paste0(list,".Rda"))
    if(is.null(range)){
      range = 1:length(out.mat)
    }
    output =  rbind(round(apply(simplify2array(out.mat)[,,range], 1:2, mean),3),
                     round(apply(simplify2array(out.mat)[,,range], 1:2, sd)/sqrt(5),3))
  }
  
  output
}

analyze("est_n100p60q30", broken=T, array=1:5)
analyze("estfull_n100p60q30")
analyze("est_n100p30q60", broken=T, array=1:5)
analyze("estfull_n100p30q60")
analyze("est_n150p200q200", broken=T, array=1:5)
analyze("est_n150p300q300", broken=T, array=1)
analyze("est_n100p200q200modelB", broken=T, array=1:5)
analyze("est_n200p200q200modelB", broken=T, array=1:5)

analyze("estmis_n100p60q30")
analyze("estmis_n100p30q60")
analyze("estmis_n150p200q200", broken=T, array=1:4)
analyze("estmis_n150p300q300", broken=T, array=1:4)
analyze("estmis_n100p200q200modelB", broken=T, array=1:4)
analyze("estmis_n200p200q200modelB", broken=T, array=1:4)
