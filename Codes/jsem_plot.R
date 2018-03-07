##
rm(list=ls())
setwd('D:/Study/My projects/Stratified-mult-GGM/Codes/')
library(Matrix)

make.symm = function(elements,k){
  mat = matrix(0,k,k)
  nelem = 1
  for(i in 2:k){
    for(j in 1:(i-1)){
      mat[i,j] = elements[nelem]
      nelem = nelem + 1
    }
  }
  mat = mat + t(mat)
  mat
}

set.seed(03052018)
p = 10
p1 = p/2
b11 = make.symm(rbinom(p1*(p1-1)/2,1,prob=0.5), k=p1)
b12 = make.symm(rbinom(p1*(p1-1)/2,1,prob=0.5), k=p1)
b21 = make.symm(rbinom(p1*(p1-1)/2,1,prob=0.5), k=p1)
b22 = make.symm(rbinom(p1*(p1-1)/2,1,prob=0.5), k=p1)

corrplot(a1, outline=T, cl.pos="n", tl.pos="n", method="square", col="black",
         main="1", mar=c(0,1,3,1), addgrid.col=NA)

pdf("adj1.pdf",width=4,height=4)
image(bdiag(b11,b21),ylab="",sub="",xlab="", axes=F, main="1", breaks=1)
dev.off()
pdf("adj2.pdf",width=4,height=4)
image(bdiag(b11,b22),ylab="",sub="",xlab="", axes=F, main="2")
dev.off()
pdf("adj3.pdf",width=4,height=4)
image(bdiag(b12,b21),ylab="",sub="",xlab="", axes=F, main="3")
dev.off()
pdf("adj4.pdf",width=4,height=4)
image(bdiag(b12,b22),ylab="",sub="",xlab="", axes=F, main="4")
dev.off()

# cplot = function(mat) corrplot(mat, method="square",col="black", tl.pos="n", cl.pos="n", addgrid.col=NA,
#                                pch.col="white")
# par(mfrow=c(1,4))
# cplot(a1)
# cplot(a2)
# cplot(a3)
# cplot(a4)
# par(mfrow=c(1,1))
