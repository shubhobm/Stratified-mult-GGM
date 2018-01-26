rm(list=ls())
setwd('d:/Study/My projects/Stratified-mult-GGM/Codes/testing_outputs/')

analyze = function(list, broken=FALSE){
  matrix.list = list()

  if(broken){
    for(i in 1:4){
      load(paste0(list,"_",i,".Rda"))
      matrix.list[[i]] = matrix(unlist(out.mat), ncol=4, byrow=T)
      matrix = lapply(matrix.list, rbind)
    }
  } else{
    load(paste0(list,".Rda"))
    matrix = matrix(unlist(out.mat), ncol=4, byrow=T)
  }
  rbind(round(apply(matrix,2,mean),3),
        round(apply(matrix,2,sd),3))
}

analyze.size = function(list, broken=FALSE){
  vec.list = list()
  
  if(broken){
    for(i in 1:4){
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

analyze("testnew_n100p60q30")
analyze("testnew_n100p30q60")
analyze("testnew_n150p200q200")
analyze("testnew_n150p300q300")
analyze("testnew_n100p200q200modelB")
analyze("testnew_n200p200q200modelB")

analyze.size("outtestsizenew_n100p60q30")
analyze.size("outtestsizenew_n100p30q60")
analyze.size("outtestsizenew_n150p200q200")
analyze.size("outtestsizenew_n150p300q300", broken=T)
analyze.size("outtestsizenew_n200p200q200modelB", broken=T)
