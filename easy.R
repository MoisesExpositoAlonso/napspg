set.seed(1)

library(dplyr)
library(moiR)
library(data.table)
devtools::load_all("~/napspg")

# bedfile<-"/home/mexpositoalonso/nap/napML/databig/515g.bed"
# mapfile<-"/home/mexpositoalonso/napML/napML/databig/515g.map"
# famfile<-"/home/mexpositoalonso/napML/nap/databig/515g.fam"

setwd('~/napspg')
bedfile<-"easy.bed"
mapfile<-"easy.map"
famfile<-"easy.fam"

map<-fread(mapfile)
fam<-fread(famfile)
bed<-readbed(bedfile, N = nrow(fam), p = nrow(map), myrows = 1:nrow(fam),mycols=1:nrow(map))

y<-fam[,6]
N<-myrows<-nrow(fam)
p<-mycols<-nrow(map)


res<-nap(
bedfile,famfile,mapfile,
myrows=NULL,mycols=NULL,
s=NULL,mod=1,epi=1,iter=20,k=1
)


# par<-
# epi<-
# res<-napspgC(bedfile, N, p, myrows, mycols, y, par, epi, mod, maxit = 20L, verbose = 1L)  
#   
# library(microbenchmark)
# nsnps=100000
# sstart<-rnorm(nsnps,0,0.01)
# 
# # X<-readbed(bedfile,N,p , 1:N,1:nsnps)
# # G<-attach.big.matrix("example.desc")
# 
# 
# r<-napSPG_R(bedfile,N,p,1:N,1:nsnps,fn(y),sstart,mod=1,epi=1,iter=10)
# r2<-nap(bedfile,famfile,mapfile,
#         1:N,1:nsnps,
#         sstart,mod=1,epi=1,iter=100,k=5)
# r2
# 
# som<-sample(1:length(r2$w),0.1*length(r2$w))
# 
# strue=ssimC(1:nsnps,0.01)
# plot(sinf,sstart)
# plot(sinf,strue)
# plot(r2$s,strue)
# 
# hist(fn(y))
# 
# microbenchmark(times = 3,
#   memory=napSPG_(fn(y),A = readbed(bedfile,N,p , 1:N,1:nsnps),s = sstart,iter=50) ,
#   filebacked=napSPG(fn(y),A = G,h_ = 1:nrow(G),m_=1:nsnps,s = sstart,iter=50),
#   onlyC=nap(bedfile,N,p,1:N,1:nsnps,fn(y),sstart,mod=1,epi=1,iter=50)
# )
