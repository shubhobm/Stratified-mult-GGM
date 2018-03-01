rm(list=ls())
setwd('d:/Study/My projects/Stratified-mult-GGM/Codes/jsem_outputs/')

analyze = function(filename, broken=FALSE, array=NULL, range=NULL){
  load(paste0(filename,".Rda"))
  out.mat1 = matrix(unlist(out.mat), ncol=5, byrow=T)
  rbind(round(apply(out.mat1, 2, mean),3),
        round(apply(out.mat1, 2, sd)/sqrt(5),3))
}

analyze("jsem_n100p60q30")
analyze("jsem_n100p30q60")
analyze("jsem_n150p200q200")
analyze("jsem_n150p300q300")
analyze("jsem_n100p200q200modelB")
analyze("jsem_n200p200q200modelB")
