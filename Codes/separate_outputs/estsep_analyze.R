rm(list=ls())
setwd('d:/Study/My projects/Stratified-mult-GGM/Codes/separate_outputs/')

analyze = function(list, broken=FALSE, array=NULL){
  matrix.list = list()

  if(broken){
    output = list()
    index = 1
    for(i in array){
      load(paste0(list,"_",i,".Rda"))
      output[[index]] = rbind(apply(simplify2array(out.mat), 1:2, mean),
                          apply(simplify2array(out.mat), 1:2, sd))
      index = index+1
    }
    output = round(apply(simplify2array(output), 1:2, mean),3)
  } else{
    load(paste0(list,".Rda"))
    matrix = matrix(unlist(out.mat), ncol=4, byrow=T)
    output =  rbind(round(apply(simplify2array(out.mat), 1:2, mean),3),
                     round(apply(simplify2array(out.mat), 1:2, sd),3))
  }
  
  output
}

analyze.size = function(list, broken=FALSE){
  vec.list = list()
  
  if(broken){
    for(i in 1:5){
      load(paste0(list,"_",i,".Rda"))
      vec.list[[i]] = as.numeric(unlist(out.mat))
      vector = as.numeric(unlist(vec.list))
    }
  } else{
    load(paste0(list,".Rda"))
    vector = as.numeric(unlist(out.mat))
  }
  c(round(mean(vector),3), round(sd(vector),3))
}

analyze("estsep_n100p30q60", broken=T, array=2:3)
analyze("estsep_n100p60q30", broken=T, array=2:3)
analyze("estsep_n150p200q200", broken=T, array=2)
analyze("est_n150p300q300", broken=T, array=1:5)
analyze("est_n100p200q200modelB", broken=T, array=1:5)
analyze("est_n200p200q200modelB", broken=T, array=1:5)

analyze("estmis_n100p60q30")
analyze("estmis_n100p30q60")
analyze("estmis_n150p200q200", broken=T, array=1:4)
analyze("estmis_n150p300q300", broken=T, array=1:4)
analyze("estmis_n100p200q200modelB", broken=T, array=1:4)
analyze("estmis_n200p200q200modelB", broken=T, array=1:4)
