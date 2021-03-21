## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=FALSE,message=FALSE----------------------------------------------
devtools::install_local("Bayestrat_0.1.0.tar.gz")

## ----warning=F,message=F------------------------------------------------------
library(Bayestrat)

## -----------------------------------------------------------------------------
data(data.test)
dim(data.test)
data.test[1:5,1:7]

## -----------------------------------------------------------------------------
data(pc)
dim(pc)
pc[1:5,1:5]

## -----------------------------------------------------------------------------
data(data.null)
dim(data.null)
data.null[1:5,1:5]

## ----warning=F,message=F------------------------------------------------------
setwd("~/Downloads/plink1.9")
out1<-bayestrat(data.test,data.null=NULL,data.pc=NULL,race=F,plink=T,
								data.null.name="data.null",n.pc=40,n.cov=2,priorType="Laplace",
								priorError=c(2,1),priorLaplace=list(type="random",G=c(1.01,0.01),
																										startlambda=NULL),
								priorNormal=list(type="random",IG=c(1,1),startvarPC=NULL),
								nIter=5000,burnIn=500,minAbsBeta=1e-09,alpha=c(0.95,0.999,0.9999),
								n.core=1,save.pc=F,checkConvergence=F,plot.interval=NULL,n.chains=1,
								check.whichPC=c(1,2,3),check.whichSNP=1,seed=30)

## -----------------------------------------------------------------------------
out2<-bayestrat(data.test,data.null,race=F,plink=F,n.pc=40,n.cov=2,
								priorLaplace=list(type="random",G=c(1.01,0.01),startlambda=NULL),
								nIter=5000,burnIn=500,save.pc=F,seed=30)

## -----------------------------------------------------------------------------
out3<-bayestrat(data.test,data.pc=pc,race=F,n.pc=40,n.cov=2,priorType="Normal",
								nIter=5000,burnIn=500,seed=30)

## -----------------------------------------------------------------------------
out4<-bayestrat(data.test,race=T,n.pc=0,n.cov=3,nIter=5000,burnIn=500,seed=30)

## -----------------------------------------------------------------------------
out5<-bayestrat(data.test[,-4],race=F,n.pc=0,n.cov=2,nIter=5000,burnIn=500,seed=30)

## -----------------------------------------------------------------------------
names(out2)

## -----------------------------------------------------------------------------
result<-bayestratSummary(out2) #40PCs
names(result)
result$snp

## -----------------------------------------------------------------------------
bayestratSummary(out5)$snp #noAdj
bayestratSummary(out4)$snp #race

## ----message=F,warning=F------------------------------------------------------
out6<-bayestrat(data.test,data.null,plink=F,n.pc=40,n.cov=2,
								priorLaplace=list(type="random",G=c(1.01,0.01),startlambda=NULL),
								nIter=5000,burnIn=500,save.pc=F,n.core=2,seed=30)

## ----message=F,warning=F------------------------------------------------------
setwd("~/Documents/Genetics/bayestratDraft")
bayestrat(data.test,data.pc=pc,n.pc=40,n.cov=2,nIter=5000,burnIn=500,checkConvergence=T,
					n.chains=2,check.whichPC=c(1:4,seq(5,40,5)),check.whichSNP=1,seed=30)

## ----message=F,warning=F------------------------------------------------------
data(pc)
out6<-pcCluster(data.test,data.pc=pc,race=T,cluster=T,cluster.whichPC=c(1,2,3),
								n.cluster=4,iter.max=10,algorithm="Hartigan-Wong",plotPC=F,
								plot.whichPC=c(1,2,3),colorby="clusters",colors=NULL)
out6$confusion.matrix

## ----message=F,warning=F------------------------------------------------------
out7<-pcCluster(data.test,data.pc=pc,race=T,cluster=T,n.cluster=4,plotPC=T,
								plot.whichPC=c(1,2,3),colorby="clusters")

## -----------------------------------------------------------------------------
data.null<-data.null[,1:100]
r<-nrow(data.null)
c<-ncol(data.null)
ped<-matrix(nrow=r,ncol=c*2+6)
ped<-as.data.frame(ped)
ped[,1]<-1:r
ped[,2]<-rep(1,r)
ped[,3]<-rep(0,r)
ped[,4]<-rep(0,r)
ped[,5]<-rep(1,r)
ped[,6]<-data.test$Y
for(i in 1:r){
	for(j in 1:c){
		if(data.null[i,j]=="0"){
			ped[i,2*j-1+6]<-"C" #C major allele, A minor allele
			ped[i,2*j+6]<-"C"
		}else if(data.null[i,j]=="1"){
			ped[i,2*j-1+6]<-"A"
			ped[i,2*j+6]<-"C"
		}else if(data.null[i,j]=="2"){
			ped[i,2*j-1+6]<-"A"
			ped[i,2*j+6]<-"A"
		}
	}
}
ped[1:5,1:10]
write.table(ped,file="data.null.ped",row.names=F,col.names = F,quote=F)

## -----------------------------------------------------------------------------
map<-matrix(0,nrow=c,ncol=4)
map<-as.data.frame(map)
map[,1]<-rep("1",c)
name<-c()
for(i in 1:c){
	name<-c(name,paste("SNP",i,sep=""))
}
map[,2]<-name
map[,3]<-rep(0,c)
map[,4]<-seq(5000,5990,10)
head(map)
write.table(map,file="data.null.map",row.names=F,col.names = F,quote=F)

## ----comment='',echo=F--------------------------------------------------------
cat(readLines("BayestratDiagnosticNumericalSummary.txt"), sep = '\n')

