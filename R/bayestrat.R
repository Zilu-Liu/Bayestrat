#' bayestrat
#'
#' \code{bayestrat} conducts genetic association tests for single nucleotide polymorphisms (SNPs) while accounting for population stratification using principal components (PCs) under the Bayesian framework.
#'
#' @param data.test A data frame consisting of a continuous phenotype, (non-genetic) covariates and genotypes for candidate SNPs as columns, individuals as rows. The column name for race information has to be \code{Race}, if applicable.
#' @param data.null (optional) A data frame or matrix of genotypes for null SNPs. Rows represent individuals corresponding to \code{data.test}. Necessary if \code{plink=FALSE}.
#' @param data.pc (optional) A data frame or matrix of PC scores. Rows represent individuals corresponding to \code{data.test}. Columns represent scores for a particular PC. PCs to be included in the regression model have to be contained in this dataset. Default is NULL, in which case Bayestrat will conduct PC computation.
#' @param race TRUE or FALSE, whether to include race information as a covariate.
#' @param plink TRUE or FALSE, whether to call PLINK1.9 (installed under the current working directory) for PC calculation, suitable for large null data. If FALSE, PCs will be calculated within R.
#' @param data.null.name (optional) A character of the common filename prefix of the PED and MAP files (.ped, .fam) for the null data. Necessary if \code{plink=TRUE}.
#' @param n.pc The total number of PCs to be included into model. Set to 0 if no PC to be added.
#' @param n.cov The number of covariates (other than PCs) added to model. Set to 0 if there is no covariates.
#' @param priorType "Laplace" or "Normal", whether to assign PC effects with Lapalce or normal priors.
#' @param priorError A vector of values of \eqn{a_e} and \eqn{b_e} for the Inverse-Gamma distributed random error term, i.e., \eqn{\sigma_e^2} follows \eqn{IG(a_e,b_e)}.
#' @param priorLaplace A list for \code{priorType="Laplace"}, consisitng of:
#' \describe{
#' \item{type}{"fixed" or "random", whether to assign \eqn{\lambda} with a fixed constant or not. The prior variance of PC effects is 2/\eqn{\lambda^2}.}
#' \item{G}{If \code{type="random"}, values of the shape parameter \eqn{a_\lambda} and rate parameter \eqn{b_\lambda} for the Gamma-distributed \eqn{\lambda^2}, i.e. \eqn{\lambda^2} follows \eqn{G(a_\lambda,b_\lambda)}.}
#' \item{startlambda}{The start value of \eqn{\lambda} in the MCMC iterations. If \code{type="fixed"}, this is the constant value assigned to \eqn{\lambda}.}
#'}
#'
#' @param priorNormal A list for \code{priorType="Normal"}, consisitng of:
#' \describe{
#' \item{type}{"fixed" or "random", whether to assign \eqn{\sigma^2_\gamma} (prior variance of PC effects) with a fixed constant or not.}
#' \item{IG}{If \code{type="random"}, values of \eqn{a_\gamma} and \eqn{b_\gamma} for the Inverse-Gamma distributed \eqn{\sigma^2_\gamma}, i.e. \eqn{\sigma^2_\gamma} follows \eqn{IG(a_\gamma,b_\gamma)}.}
#' \item{startvarPC}{The start value of \eqn{\sigma^2_\gamma} in the MCMC iterations. If \code{type="fixed"}, this is the constant value assigned to \eqn{\sigma^2_\gamma}.}
#'}
#'
#' @param nIter The number of total iterations. Default is 10000.
#' @param burnIn The number of burn-in. Default is 1000.
#' @param minAbsBeta If \code{priorType="Laplace"}, this is the minimum absolute value of the PC effects to avoid numeric problems. Default is 1e-9.
#' @param alpha A vector of length three consisting of confidence levels. Credible intervals for each level will be calculated.
#' @param n.core The number of cores to be used in testing SNPs. Parallelization computing is realized by setting n.core greater than 1. Default is 1.
#' @param save.pc TRUE or FALSE, whether to save the calculated PC data to file.
#' @param checkConvergence TRUE or FALSE, whether to conduct convergence diagnostics.
#' @param plot.interval The interval between two posterior samples that are extracted to draw trace plots.
#' @param n.chains The number of chains to run in convergence diagnostics. Default if 1.
#' @param check.whichPC On which PCs the convergence diagnostics will be performed. Set to NULL if no diagnostics for PCs.
#' @param check.whichSNP On which SNP the convergence diagnostics will be performed.
#' @param seed Seed to be used for the implemented MCMC. Default is a random seed.
#'
#' @return A list with the following components:
#' \describe{
#' \item{n}{A vector of sample sizes for testing each SNP.}
#' \item{mu}{The estimated value (posterior mean) of \eqn{\mu}.}
#' \item{varE}{The estimated value (posterior mean) of \eqn{\sigma_e^2}.}
#' \item{CIsnp}{A list of three data frames of credible intervals for SNP effects at each confidence level.}
#' \item{SDsnp}{The estimated standard errors for SNP effects.}
#' \item{CIcov}{A list of three data frames of credible intervals for covariate effects at each confidence level.}
#' \item{SDcov}{The estimated standard errors for covariate effects.}
#' \item{CIpc}{A list of three data frames of credible intervals for PC effects at each confidence level.}
#' \item{SDpc}{The estimated standard errors for PC effects.}
#' \item{lambda}{If \code{priorType="Laplace"}, the estimated value of \eqn{\lambda}.}
#' \item{SDlambda}{If \code{priorType="Laplace"}, the estimated standard error of \eqn{\lambda}.}
#' \item{varPC}{If \code{priorType="Normal"}, the estimated value of \eqn{\sigma_\gamma^2}.}
#' \item{SDvarPC}{If \code{priorType="Normal"}, the estimated standard error of \eqn{\sigma_\gamma^2}.}
#' }
#'
#' @export
#' @importFrom data.table fread
#' @importFrom coda mcmc mcmc.list raftery.diag rejectionRate gelman.diag gelman.plot
#' @importFrom flashpcaR flashpca
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline lines par plot
#' @importFrom stats cmdscale dist dnorm kmeans model.matrix na.pass prcomp quantile rchisq rgamma rnorm runif sd var
#' @importFrom utils browseURL install.packages
#'
#' @examples
#' data(data.test)
#' data(data.null)
#' out<-bayestrat(data.test,data.null,race=FALSE,plink=FALSE,n.pc=2,n.cov=2,save.pc=FALSE,seed=30)
#'
#' @useDynLib Bayestrat sample_beta


bayestrat<-function(data.test,data.null=NULL,data.pc=NULL,race=F,plink=F,
										data.null.name=NULL,n.pc,n.cov,priorType="Laplace",priorError=c(2,1),
										priorLaplace=list(type="random",G=c(1.4,0.4),startlambda=NULL),
										priorNormal=list(type="random",IG=c(1,1),startvarPC=NULL),
										nIter=10000,burnIn=1000,minAbsBeta=1e-09,alpha=c(0.95,0.999,0.9999),
										n.core=1,save.pc=T,checkConvergence=F,plot.interval=NULL,n.chains=1,
										check.whichPC=c(1,2,3),check.whichSNP=1,seed=NULL)
{
	if(n.core>1){
		if(!require("foreach")){
			install.packages("foreach")
			library(foreach)
		}
		if(!require("doParallel")){
			install.packages("doParallel")
			library(doParallel)
		}
	}

	if(is.null(seed)==F){
		set.seed(seed)
	}

	if(race==F){
		if(is.null(data.test$Race)==F){ #have race data but do not add in
			data.test<-subset(data.test,select=-Race)
			cat("Column Race deleted\n")
		}
	}
	m<-ncol(data.test)-1-n.cov

	n<-nrow(data.test)

	#### calculate PC ####
	if(n.pc>0){
		if(!is.null(data.pc)){
			pc<-data.pc
		}else{
			if(plink){ #Use plink
				if(is.null(data.null.name)){
					stop("Null data name has to be provided\n")
				}
				cat("Start PC calculation using PLINK\n")
				map<-fread(paste(data.null.name,".map",sep=""))
				system(paste("./plink --file",data.null.name,"--pca",min(n,nrow(map)),"--out pca")) #plink exe in the same folder
				eigenvec<-fread("pca.eigenvec")
				eigenvec<-as.data.frame(eigenvec[,-c(1:2)])
				eigenval<-fread("pca.eigenval")
				eigenval<-as.vector(as.matrix(eigenval))
				if(any(eigenval<=0)){
					cat("Negative eigen values exist, will be deleted in PC calculation\n")
					eigenval<-eigenval[-which(eigenval<=0)]
				}
				pc<-matrix(0,nrow=n,ncol=length(eigenval))
				for(i in 1:ncol(pc)){
					pc[,i]<-eigenvec[,i]*sqrt(eigenval[i])
				}
			}else{ #prcomp or flashpca
				#check null data
				if(is.null(data.null)){
					stop("Null dataset has to be provided\n")
				}
				if(nrow(data.null)!=n){
					stop("data.null and data.test must have the same number of rows\n")
				}
				tmp<-apply(data.null,2,sd)==0
				if(any(tmp)){
					cat("WARNING: Null SNPs with no sequence variance will be deleted for PC calculation\n")
					data.null<-data.null[,-which(tmp==T)]
				}
				if(any(is.na(data.null))){
					tmp<-unlist(apply(data.null,2,FUN=function(v){if(any(is.na(v))){return(TRUE)}else(return(FALSE))}))
					cat("WARNING: Null SNPs with missing value will be deleted for PC calculation\n")
					data.null<-data.null[,-which(tmp==T)]
				}

				k<-ncol(data.null)

				if(n.pc<min(n,k)/2){ #Use flashpca: fast but only returns up to n/2-1 PCs
					cat("Start PC calculation using flashpca\n")
					pc<-flashpca(as.matrix(data.null),ndim=n.pc,stand="binom2")
					pc<-pc$projection
				}else{ #Use prcomp or plink
					if(n*k<=10^7){ #small dataset use prcomp
						cat("Start PC calculation using prcomp\n")
						pc<-prcomp(data.null,center=T,scale. = T)
						pc<-pc$x
						#pc<-FastPCA(scale(as.matrix(data.null)),PCs=n.pc) #best for calculating under 30 PCs
						#pc<-pc$scores
					}else{ #large data set needs plink
						stop("PC calculation by PLINK is required for large null data")
					}
				}
			}
			if(save.pc==T){
				pc<-as.data.frame(pc)
				namesPC<-rep(0,ncol(pc))
				for(i in 1:ncol(pc)){
					namesPC[i]<-paste("PC",i,sep="")
				}
				colnames(pc)<-namesPC
				save(pc,file="pc.RData")
			}
			cat("End PC calculation\n")
		}
	}else{ #n.pc=0
		pc<-NULL
		cat("No PC added to model\n")
	}

	#### dummy covariates and organize cov data ####
	if(n.cov>0){
		cov<-as.data.frame(data.test[,c(2:(1+n.cov))])
		colnames(cov)<-colnames(data.test)[c(2:(1+n.cov))]
	}else{
		cov<-NULL
	}

	if(n.cov>0){
		cov.dum<-as.data.frame(matrix(0,nrow=n,ncol=1))
		cov.c<-as.data.frame(matrix(0,nrow=n,ncol=1))
		n.ctnCov<-0
		for(i in 1:n.cov){
			if(class(cov[1,i])=="factor"){
				l<-table(cov[,i])
				options(na.action=na.pass)
				tmp<-as.data.frame(model.matrix(~cov[,i]-1)[,-which.max(l)]) #the value of cov with max individuals is the baseline.
				colnames(tmp)<-names(l)[-which.max(l)]
				cov.dum<-cbind(cov.dum,tmp)
			}else if(class(cov[1,i])=="numeric" | class(cov[1,i])=="integer"){
				cov.c<-cbind(cov.c,cov[,i])
				colnames(cov.c)[i+1]<-colnames(cov)[i]
				n.ctnCov<-n.ctnCov+1
			}else{
				cat("Warning: Covariates in data.test have to be factor (categorical variables) or numeric/integer (continuous variables)")
			}
		}
		data.cov<-as.data.frame(cbind(cov.c[,-1],cov.dum[,-1]))
		colnames(data.cov)<-c(colnames(cov.c)[-1],colnames(cov.dum)[-1])
		n.cov.new<-ncol(data.cov)
	}else{ #no cov
		data.cov<-NULL
		n.cov.new<-0
		n.ctnCov<-0
	}

	MAF<-apply(data.test[,-c(1:(n.cov+1))],2,mean,na.rm=T)/2
	data.test[,1]<-scale(data.test[,1])

	#### chenck convergence ####
	if(checkConvergence){

		if(is.null(plot.interval)==T){
			plot.interval<-(nIter-burnIn)%/%300
		}

		if(is.null(check.whichPC)==F){
			PCname<-rep(NA,length(check.whichPC))
			for(i in 1:length(check.whichPC)){
				PCname[i]<-paste("PC",check.whichPC[i],sep="")
			}
		}

		my.draws<-vector("list",n.chains)

		snp.vector<-matrix(numeric(n.chains*(nIter-burnIn)),ncol=n.chains)
		cov.matrix<-vector("list",n.chains)
		pc.matrix<-vector("list",n.chains)
		lambdaVarPC.vector<-snp.vector
		if(priorType=="Laplace"){
			xname<-"lambda"
		}else{
			xname<-"varPC"
		}

		for(i in 1:n.chains){
			startbF<-runif(1,-1,1)
			if(n.pc>0){
				if(priorType=="Laplace"){
					start.lambda<-rgamma(1,shape=priorLaplace$G[1],rate=priorLaplace$G[2])
					start.lambda<-sqrt(start.lambda)
					priorLaplace$startlambda<-start.lambda
				}else if(priorType=="Normal"){
					start.varPC<-rgamma(1,shape=priorNormal$IG[1],rate=priorNormal$IG[2])
					start.varPC<-1/start.varPC
					priorNormal$startvarPC<-start.varPC
				}
			}
			out<-testSNP(check.whichSNP,data.test,data.cov,n.ctnCov,pc,MAF,n.pc,n.cov,
									 priorType,priorError,priorLaplace,priorNormal,nIter,burnIn,
									 minAbsBeta,alpha,startbF)

			#extract posterior samples
			if(is.null(check.whichPC)==F){
				if(priorType=="Laplace"){
					draw<-cbind(out$post_bF_sample,out$post_lambda_sample,
											out$post_bL_sample[,check.whichPC])
					lambdaVarPC.vector[,i]<-out$post_lambda_sample
					pc.matrix[[i]]<-out$post_bL_sample[,check.whichPC]
				}else if(priorType=="Normal"){
					draw<-cbind(out$post_bF_sample,out$post_varPC_sample,
											out$post_bR_sample[,check.whichPC])
					lambdaVarPC.vector[,i]<-out$post_varPC_sample
					pc.matrix[[i]]<-out$post_bR_sample[,check.whichPC]
				}
				colnames(draw)<-c("snp",colnames(data.cov),xname,PCname)
				my.draws[[i]]<-mcmc(draw)
				snp.vector[,i]<-out$post_bF_sample[,1]
				cov.matrix[[i]]<-out$post_bF_sample[,-1]
			}else{ #check.whichPC=NULL
				if(n.pc>0){
					if(priorType=="Laplace"){
						draw<-cbind(out$post_bF_sample,out$post_lambda_sample)
						lambdaVarPC.vector[,i]<-out$post_lambda_sample
					}else if(priorType=="Normal"){
						draw<-cbind(out$post_bF_sample,out$post_varPC_sample)
						lambdaVarPC.vector[,i]<-out$post_varPC_sample
					}
					colnames(draw)<-c("snp",colnames(data.cov),xname)
					my.draws[[i]]<-mcmc(draw)
					snp.vector[,i]<-out$post_bF_sample[,1]
					cov.matrix[[i]]<-out$post_bF_sample[,-1]
				}else{ #n.pc=0
					draw<-out$post_bF_sample
					colnames(draw)<-c("snp",colnames(data.cov))
					my.draws[[i]]<-mcmc(draw)
					snp.vector[,i]<-out$post_bF_sample[,1]
					cov.matrix[[i]]<-out$post_bF_sample[,-1]
				}
			}
		}

		mh.list<-mcmc.list(my.draws)
		if(n.chains==1){
			D_Raf<-raftery.diag(my.draws[[1]],q=0.025,r=0.005,s=0.95)
			acceptance.rate<-1-rejectionRate(my.draws[[1]])

			#Convergence Numeric Summary
			sink("BayestratDiagnosticNumericalSummary.txt")
			cat("Gelman Diagnostic needs at least 2 chains.\n\n\n")
			cat("Raftery Diagnostic \n")
			print(D_Raf)
			cat("Acceptance Rate \n")
			cat(acceptance.rate,"\n")
			sink()

			#Convergence Graphical Summary
			pdf(file="BayestratDiagnosticPlot.pdf")

			seq.get<-seq(1,length(snp.vector[,1]),plot.interval)
			par(mfrow=c(3,2))

			ymin<-min(snp.vector[seq.get,1])-0.1
			ymax<-max(snp.vector[seq.get,1])+0.1
			plot(seq.get,snp.vector[seq.get,1],type="l",col=1,ylim=c(ymin,ymax),
					 ylab="snp",xlab="iteration")
			abline(h=0,col="red")

			if(n.cov.new>0){
				for (i in 1:ncol(cov.matrix[[1]])){
					ymin=min(cov.matrix[[1]][seq.get,i])-0.1
					ymax=max(cov.matrix[[1]][seq.get,i])+0.1
					plot(seq.get,cov.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(ymin,ymax),
							 ylab=colnames(data.cov)[i],xlab="iteration")
					abline(h=0,col="red")
				}
			}

			if(n.pc>0){
				ymin<-min(lambdaVarPC.vector[seq.get,1])-0.1
				ymax<-max(lambdaVarPC.vector[seq.get,1])+0.1
				plot(seq.get,lambdaVarPC.vector[seq.get,1],type="l",col=1,ylim=c(ymin,ymax),
						 ylab=xname,xlab="iteration")
				abline(h=0,col="red")
			}

			if(is.null(check.whichPC)==F){
				for (i in 1:length(check.whichPC)){
					ymin=min(pc.matrix[[1]][seq.get,i])-0.1
					ymax=max(pc.matrix[[1]][seq.get,i])+0.1
					plot(seq.get,pc.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(ymin,ymax),
							 ylab=PCname[i],xlab="iteration")
					abline(h=0,col="red")
				}
			}

			dev.off()

		}else if(n.chains>1){
			D_Gel<-gelman.diag(mh.list) #ideal reduction factor<1.1
			D_Raf<-vector("list",n.chains)
			acceptance.rate<-vector("list",n.chains)

			for (i in 1:n.chains){
				D_Raf[[i]]<-raftery.diag(my.draws[[i]],q=0.025,r=0.005,s=0.95)
				acceptance.rate[[i]]<-1-rejectionRate(my.draws[[i]])   #acceptance rate of one MCMC
			}

			sink("BayestratDiagnosticNumericalSummary.txt")
			cat("Gelman Diagnostic \n")
			print(D_Gel)
			cat("\n\n\n")

			for (i in 1:n.chains){
				cat("Result for chain ",i,"\n",sep="")
				#cat("Chain ",i,"'s Result \n",sep="")
				cat("Raftery Diagnostic \n")
				print(D_Raf[[i]])
				cat("Acceptance Rate \n")
				print(acceptance.rate[[i]])
				cat("\n\n\n")
			}
			sink()


			pdf(file="BayestratDiagnosticPlot.pdf")

			gelman.plot(mh.list)

			seq.get<-seq(1,length(snp.vector[,1]),plot.interval)
			par(mfrow=c(3,2))

			ymin<-min(snp.vector[seq.get,1])-0.1
			ymax<-max(snp.vector[seq.get,1])+0.1
			plot(seq.get,snp.vector[seq.get,1],type="l",col=1,ylim=c(ymin,ymax),
					 ylab="snp",xlab="iteration")
			for (i in 2: n.chains){
				lines(seq.get,snp.vector[seq.get,i],col=i)
			}
			abline(h=0,col=(n.chains+1))

			if(n.cov.new>0){
				for (i in 1:ncol(cov.matrix[[1]])){
					ymin<-min(cov.matrix[[1]][seq.get,i])-0.1
					ymax<-max(cov.matrix[[1]][seq.get,i])+0.1
					plot(seq.get,cov.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(ymin,ymax),
							 ylab=colnames(data.cov)[i],xlab="iteration")
					for (j in 2:n.chains){
						lines(seq.get,cov.matrix[[j]][seq.get,i],col=j)
					}
					abline(h=0,col=(n.chains+1))
				}
			}

			if(n.pc>0){
				ymin<-min(lambdaVarPC.vector[seq.get,1])-0.1
				ymax<-max(lambdaVarPC.vector[seq.get,1])+0.1
				plot(seq.get,lambdaVarPC.vector[seq.get,1],type="l",col=1,ylim=c(ymin,ymax),
						 ylab=xname,xlab="iteration")
				for (i in 2: n.chains){
					lines(seq.get,lambdaVarPC.vector[seq.get,i],col=i)
				}
				abline(h=0,col=(n.chains+1))
			}

			if(is.null(check.whichPC)==F){
				for (i in 1:length(check.whichPC)){
					ymin<-min(pc.matrix[[1]][seq.get,i])-0.1
					ymax<-max(pc.matrix[[1]][seq.get,i])+0.1
					plot(seq.get,pc.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(ymin,ymax),
							 ylab=PCname[i],xlab="iteration")
					for (j in 2:n.chains){
						lines(seq.get,pc.matrix[[j]][seq.get,i],col=j)
					}
					abline(h=0,col=(n.chains+1))
				}
			}

			tmp<-dev.off()
			cat("Diagnostic results saved\n")
		}

	}else{ #checkConvergence=F
		out<-list(n=rep(0,m))
		out$mu<-0
		out$varE<-0
		out$CIsnp<-list()
		out$CIsnp[[1]]<-data.frame()
		out$CIsnp[[2]]<-data.frame()
		out$CIsnp[[3]]<-data.frame()
		out$SDsnp<-rep(0,m)

		if(n.cov.new>0){
			out$CIcov<-list()
			out$CIcov[[1]]<-as.data.frame(matrix(0,nrow=n.cov.new,ncol=2))
			out$CIcov[[2]]<-out$CIcov[[1]]
			out$CIcov[[3]]<-out$CIcov[[1]]
			out$SDcov<-rep(0,n.cov.new)
		}
		if(n.pc>0){
			out$CIpc<-list()
			out$CIpc[[1]]<-as.data.frame(matrix(0,nrow=n.pc,ncol=2))
			out$CIpc[[2]]<-out$CIpc[[1]]
			out$CIpc[[3]]<-out$CIpc[[1]]
			out$SDpc<-rep(0,n.pc)
			if(priorType=="Laplace"){
				out$lambda<-0
				out$SDlambda<-0
			}else if(priorType=="Normal"){
				out$varPC<-0
				out$SDvarPC<-0
			}
		}

		if(n.core==1){
			cat("Start testing\n")
			time<-proc.time()[3]
			for(i in 1:m){
				tmp<-testSNP(i,data.test,data.cov,n.ctnCov,pc,MAF,n.pc,n.cov,priorType,
										 priorError,priorLaplace,priorNormal,nIter,burnIn,minAbsBeta,
										 alpha,startbF=0)

				out$n[i]<-tmp$n
				out$mu<-out$mu+tmp$mu/m
				out$varE<-out$varE+tmp$varE/m
				for(j in 1:3){
					out$CIsnp[[j]]<-rbind(out$CIsnp[[j]],tmp$CI.bF[[j]][1+n.cov.new,])
				}
				out$SDsnp[i]<-tmp$SD.bF[1+n.cov.new]

				if(n.cov.new>0){
					for(j in 1:3){
						out$CIcov[[j]]<-out$CIcov[[j]]+tmp$CI.bF[[j]][1:n.cov.new,]/m
					}
					out$SDcov<-out$SDcov+tmp$SD.bF[1:n.cov.new]/m
				}

				if(n.pc>0){
					if(priorType=="Laplace"){
						for(j in 1:3){
							out$CIpc[[j]]<-out$CIpc[[j]]+tmp$CI.bL[[j]]/m
						}
						out$SDpc<-out$SDpc+tmp$SD.bL/m
						out$lambda<-out$lambda+tmp$lambda/m
						out$SDlambda<-out$SDlambda+tmp$SD.lambda/m
					}else if(priorType=="Normal"){
						for(j in 1:3){
							out$CIpc[[j]]<-out$CIpc[[j]]+tmp$CI.bR[[j]]/m
						}
						out$SDpc<-out$SDpc+tmp$SD.bR/m
						out$varPC<-out$varPC+tmp$varPC/m
						out$SDvarPC<-out$SDvarPC+tmp$SD.varPC/m
					}
				}
			}
			cat("End testing\n")
			cat(paste("Times elapsed",round(proc.time()[3]-time,3)))
			cat("\n")

		}else if(n.core>1){
			cat("Start testing\n")
			time<-proc.time()[3]

			registerDoParallel(n.core)
			tmp<-foreach(i=1:m) %dopar% {
				testSNP(i,data.test,data.cov,n.ctnCov,pc,MAF,n.pc,n.cov,priorType,
								priorError,priorLaplace,priorNormal,nIter,burnIn,minAbsBeta,
								alpha,startbF=0)
			}
			stopImplicitCluster()

			cat("End testing\n")
			cat(paste("Times elapsed",round(proc.time()[3]-time,3)))
			cat("\n")
			cat("Start summarizing results\n")

			for(i in 1:m){
				out$n[i]<-tmp[[i]]$n
				out$mu<-out$mu+tmp[[i]]$mu/m
				out$varE<-out$varE+tmp[[i]]$varE/m
				for(j in 1:3){
					out$CIsnp[[j]]<-rbind(out$CIsnp[[j]],tmp[[i]]$CI.bF[[j]][1+n.cov.new,]) #rows: SNP
				}
				out$SDsnp[i]<-tmp[[i]]$SD.bF[1+n.cov.new]

				if(n.cov.new>0){
					for(j in 1:3){
						out$CIcov[[j]]<-out$CIcov[[j]]+tmp[[i]]$CI.bF[[j]][1:n.cov.new,]/m  #averaged CI for cov and PCs
					}
					out$SDcov<-out$SDcov+tmp[[i]]$SD.bF[1:n.cov.new]/m
				}

				if(n.pc>0){
					if(priorType=="Laplace"){
						for(j in 1:3){
							out$CIpc[[j]]<-out$CIpc[[j]]+tmp[[i]]$CI.bL[[j]]/m
						}
						out$SDpc<-out$SDpc+tmp[[i]]$SD.bL/m
						out$lambda<-out$lambda+tmp[[i]]$lambda/m
						out$SDlambda<-out$SDlambda+tmp[[i]]$SD.lambda/m
					}else if(priorType=="Normal"){
						for(j in 1:3){
							out$CIpc[[j]]<-out$CIpc[[j]]+tmp[[i]]$CI.bR[[j]]/m
						}
						out$SDpc<-out$SDpc+tmp[[i]]$SD.bR/m
						out$varPC<-out$varPC+tmp[[i]]$varPC/m
						out$SDvarPC<-out$SDvarPC+tmp[[i]]$SD.varPC/m
					}
				}
			}
			cat("End summarizing results\n")

		}#end if n.core>1


		SNPname<-colnames(data.test)[-c(1:(1+n.cov))]
		alphaname<-c()
		for(j in 1:3){
			alphaname<-c(alphaname,paste(alpha[j]*100,"%",sep=""))
		}
		Covname<-colnames(data.cov)
		CIname<-c("lower","upper")
		if(n.pc>0){
			if(is.null(colnames(pc))){
				PCname<-c(1:n.pc)
			}else{
				PCname<-colnames(pc)[c(1:n.pc)]
			}
		}

		names(out$n)<-SNPname
		names(out$CIsnp)<-alphaname
		for(j in 1:3){
			colnames(out$CIsnp[[j]])<-CIname
			rownames(out$CIsnp[[j]])<-SNPname
		}
		names(out$SDsnp)<-SNPname
		if(n.cov.new>0){
			names(out$CIcov)<-alphaname
			for(j in 1:3){
				colnames(out$CIcov[[j]])<-CIname
				rownames(out$CIcov[[j]])<-Covname
			}
			names(out$SDcov)<-Covname
		}

		if(n.pc>0){
			names(out$CIpc)<-alphaname
			for(j in 1:3){
				colnames(out$CIpc[[j]])<-CIname
				rownames(out$CIpc[[j]])<-PCname
			}
			names(out$SDpc)<-PCname
		}

		return(out)
		rm(tmp)
	}
}


testSNP<-function(which.snp,data.test,data.cov,n.ctnCov,pc,MAF,n.pc,n.cov,priorType,
									priorError,priorLaplace,priorNormal,nIter,burnIn,minAbsBeta,alpha,
									startbF){

	#delete missing y, cov and genotype
	if(is.null(data.cov)){ #no covariates
		data<-cbind(data.test[,1],data.test[,1+which.snp])
	}else{
		data<-cbind(data.test[,1],data.cov,data.test[,1+n.cov+which.snp])
	}
	if(any(is.na(data))){
		mis<-unlist(apply(data,2,FUN=function(v){if(any(is.na(v))) return(which(is.na(v)))}))
		mis<-sort(unique(mis))
		data<-data[-mis,]
		pc<-pc[-mis,]
	}

	#standardize pc and the continuous covariates
	if(n.pc>0){
		pc.s<-scale(pc)*sqrt(MAF[which.snp]*(1-MAF[which.snp]))
		colnames(pc.s)<-colnames(pc)
	}
	if(n.ctnCov>0){#standardize continuous cov
		data<-cbind(data[,1],scale(data[,2:(1+n.ctnCov)])*sqrt(MAF[which.snp]*(1-MAF[which.snp])),data[,-c(1:(1+n.ctnCov))])
	}

	#run the analysis
	y<-as.numeric(data[,1])
	XF<-as.matrix(data[,-1])
	XL<-NULL
	XR<-NULL
	n<-length(y)

	if(n.pc>0){
		if(priorType=="Laplace"){
			XL<-as.matrix(pc.s[,1:n.pc])
		}else if(priorType=="Normal"){
			XR<-as.matrix(pc.s[,1:n.pc])
		}
	}

	nSums <- 0
	mu<-mean(y)
	e<-(y-mu)
	varE<-var(e, na.rm=TRUE)/2

	post_mu<-0
	post_varE<-0
	post_logLik<-0
	post_yHat<-rep(0, n)
	post_yHat2<-rep(0, n)
	Normal<-!is.null(XR)
	Laplace<-!is.null(XL)

	SVD.XF<-svd(XF) #XF=u%*%diag(d)%*%v'
	SVD.XF$Vt<-t(SVD.XF$v)
	SVD.XF<-SVD.XF[-3]
	pF0<-length(SVD.XF$d)
	pF<-ncol(XF)
	bF0<-rep(0, pF0)
	bF<-rep(startbF, pF)
	post_bF<-bF
	post_bF2<-bF
	post_bF_sample<-matrix(0,nrow=nIter-burnIn,ncol=pF)
	rm(XF)

	if(Laplace){
		pL<-ncol(XL)
		xL2<-colSums(XL*XL)
		bL<-rep(0, pL)
		tmp<-1/2/sum(xL2/n)
		tau2<-rep(tmp, pL)
		if(is.null(priorLaplace$startlambda)){
			lambda<-1
		}else{
			lambda<-priorLaplace$startlambda
		}
		lambda2<-lambda^2
		post_lambda_sample<-numeric(nIter-burnIn)
		post_lambda<-0
		post_bL<-rep(0, pL)
		post_bL2<-post_bL
		post_bL_sample<-matrix(0,nrow=nIter-burnIn,ncol=pL)
		post_tau2<-rep(0, pL)
		XLstacked<-as.vector(XL)
		rm(XL)
	}else {
		lambda=NA
	}
	if(Normal){
		pR<-ncol(XR)
		xR2<-colSums(XR * XR)
		bR<-rep(0, pR)
		if(is.null(priorNormal$startvarPC)){
			varPC<-varE/2/sum(xR2/n)
		}else{
			varPC<-priorNormal$startvarPC
		}
		post_bR<-rep(0, pR)
		post_bR2<-post_bR
		post_bR_sample<-matrix(0,nrow=nIter-burnIn,ncol=pR)
		post_varPC_sample<-numeric(nIter-burnIn)
		post_varPC<-0
		XRstacked<-as.vector(XR)
		rm(XR)
	}

	#dyn.load("sample_betas.so")
	for(i in 1:nIter){
		sol<-(crossprod(SVD.XF$u, e)+bF0) #t(u)%*%e+bF0
		tmp<-sol+rnorm(n=pF0, sd=sqrt(varE))
		bF<-crossprod(SVD.XF$Vt, tmp/SVD.XF$d)
		e<-e+SVD.XF$u %*% (bF0 - tmp)
		bF0<-tmp

		if(Normal){
			ans<- .Call("sample_beta", n, pR, XRstacked, xR2,
									bR, e, rep(varPC, pR), varE, minAbsBeta)
			bR<-ans[[1]]
			e<-ans[[2]]
			SS<-crossprod(bR)+2*priorNormal$IG[2]
			df<-pR+2*priorNormal$IG[1]
			if(priorNormal$type=="random"){
				varPC<-SS/rchisq(df=df, n=1)
			}
		}
		if(Laplace){
			varBj<-tau2 * varE
			ans<- .Call("sample_beta", n, pL, XLstacked, xL2,
									bL, e, varBj, varE, minAbsBeta)
			bL<-ans[[1]]
			e<-ans[[2]]
			nu<-sqrt(varE)*lambda/abs(bL)
			tmp<-NULL
			try(tmp<-rinvGauss(n = pL, mu = nu, lambda = lambda2))
			if(!is.null(tmp))
			{
				if(!any(is.na(sqrt(tmp))))
				{
					tau2<-1/tmp
				}else{
					cat("WARNING: tau2 was not updated due to numeric problems with beta\n");
				}
			}else{
				cat("WARNING: tau2 was not updated due to numeric problems with beta\n");
			}
			if(priorLaplace$type=="random"){
				rate<-sum(tau2)/2 + priorLaplace$G[2] #G(shape,rate)
				shape<-pL + priorLaplace$G[1]
				lambda2<-rgamma(rate = rate, shape = shape,n = 1)
				if(!is.na(lambda2))
				{
					lambda<-sqrt(lambda2)
				}else{
					cat("WARNING: lambda was not updated due to numeric problems with beta\n");
				}
			}#else lambda not updated
		} #end if Laplace

		e<-e+mu
		rhs<-sum(e)/varE
		C<-n/varE
		sol<-rhs/C
		mu<-rnorm(n = 1, sd = sqrt(1/C)) + sol
		e<-e-mu
		SS<-crossprod(e)+2*priorError[2]
		df<-n+2*priorError[1]
		if(Laplace){
			if(!any(is.na(sqrt(tau2))))
			{
				SS<-SS+as.numeric(crossprod(bL/sqrt(tau2)))
			}else{
				cat("WARNING: SS was not updated due to numeric problems with beta\n");
			}
			df<-df + pL
		}
		varE<-as.numeric(SS)/rchisq(n = 1, df = df)
		sdE<-sqrt(varE)
		yHat<-(y-e)

		if(i>burnIn){
			nSums<-nSums+1
			k<-(nSums-1)/(nSums)
			tmpE<-e
			tmpSD<-sqrt(varE)
			logLik <- sum(dnorm(tmpE, sd = tmpSD, log = TRUE))
			post_logLik <- post_logLik * k + logLik/nSums
			post_mu <- post_mu * k + mu/nSums
			post_varE <- post_varE * k + varE/nSums # calculating the mean
			post_yHat <- post_yHat * k + yHat/nSums
			post_yHat2 <- post_yHat2 * k + (yHat^2)/nSums
			post_bF <- post_bF * k + bF/nSums
			post_bF2 <- post_bF2 * k + (bF^2)/nSums
			post_bF_sample[i-burnIn,]<-bF
			if(Laplace){
				post_lambda_sample[i-burnIn]<-lambda
				post_lambda <- post_lambda * k + lambda/nSums
				post_bL <- post_bL * k + bL/nSums
				post_bL2 <- post_bL2 * k + (bL^2)/nSums
				post_bL_sample[i-burnIn,]<-bL
				post_tau2 <- post_tau2 * k + tau2/nSums
			}
			if(Normal){
				post_bR <- post_bR * k + bR/nSums
				post_bR2 <- post_bR2 * k + (bR^2)/nSums
				post_bR_sample[i-burnIn,]<-bR
				post_varPC_sample[i-burnIn]<-varPC
				post_varPC <- post_varPC * k + varPC/nSums
			}
		}

	} #end for i in 1:nIter

	tmp<-sqrt(post_yHat2 - (post_yHat^2))

	CI.bF<-list()
	CI.bF[[1]]<-as.data.frame(matrix(0,nrow=pF,ncol=2))
	CI.bF[[2]]<-CI.bF[[1]]
	CI.bF[[3]]<-CI.bF[[1]]
	for (i in 1:pF){
		for(j in 1:3){
			CI.bF[[j]][i,]<-c(quantile(post_bF_sample[,i],probs=(1-alpha[j])/2,na.rm=T),quantile(post_bF_sample[,i],probs=(1+alpha[j])/2,na.rm=T))
		}
	}

	if(Laplace){
		CI.bL<-list()
		CI.bL[[1]]<-as.data.frame(matrix(0,nrow=pL,ncol=2))
		CI.bL[[2]]<-CI.bL[[1]]
		CI.bL[[3]]<-CI.bL[[1]]
		for (i in 1:pL){
			for(j in 1:3){
				CI.bL[[j]][i,]<-c(quantile(post_bL_sample[,i],probs=(1-alpha[1])/2,na.rm=T),quantile(post_bL_sample[,i],probs=(1+alpha[j])/2,na.rm=T))
			}
		}
		#CI.lambda<-c(quantile(post_lambda_sample,probs=(1-alpha[1])/2,na.rm=T),quantile(post_lambda_sample,probs=(1+alpha[j])/2,na.rm=T))
	}

	if(Normal){
		CI.bR<-list()
		CI.bR[[1]]<-as.data.frame(matrix(0,nrow=pR,ncol=2))
		CI.bR[[2]]<-CI.bR[[1]]
		CI.bR[[3]]<-CI.bR[[1]]
		for (i in 1:pR){
			for(j in 1:3){
				CI.bR[[j]][i,]<-c(quantile(post_bR_sample[,i],probs=(1-alpha[1])/2,na.rm=T),quantile(post_bR_sample[,i],probs=(1+alpha[j])/2,na.rm=T))
			}
		}
		#CI.varPC<-c(quantile(post_varPC_sample,probs=(1-alpha[1])/2,na.rm=T),quantile(post_varPC_sample,probs=(1+alpha[j])/2,na.rm=T))
	}

	#output
	out <- list(n=n, mu = post_mu, varE = post_varE)
	tmpE <- (y - post_yHat)
	tmpSD <- sqrt(post_varE)
	out$bF <- as.vector(post_bF)
	out$SD.bF <- as.vector(sqrt(post_bF2 - post_bF^2))
	out$CI.bF<-CI.bF
	out$post_bF_sample<-as.data.frame(post_bF_sample)
	if (Laplace) {
		out$bL <- as.vector(post_bL)
		tmp <- as.vector(sqrt(post_bL2 - (post_bL^2)))
		out$SD.bL <- tmp
		out$CI.bL<-CI.bL
		#out$tau2 <- post_tau2
		out$post_bL_sample<-as.data.frame(post_bL_sample)
		out$lambda <- post_lambda
		out$SD.lambda <- sd(post_lambda_sample)
		#out$CI.lambda<-CI.lambda
		out$post_lambda_sample<-post_lambda_sample
	}
	if (Normal) {
		out$bR <- as.vector(post_bR)
		tmp <- as.vector(sqrt(post_bR2 - (post_bR^2)))
		out$SD.bR <- tmp
		out$CI.bR<-CI.bR
		out$post_bR_sample<-as.data.frame(post_bR_sample)
		out$varPC <- post_varPC
		out$SD.varPC <- sd(post_varPC_sample)
		#out$CI.varPC<-CI.varPC
		out$post_varPC_sample<-post_varPC_sample
	}
	return(out)
	rm(post_bF_sample)
	rm(post_bR_sample)
	rm(post_bL_sample)
	rm(post_lambda_sample)
	rm(post_varPC_sample)

}


rinvGauss=function(n,mu,lambda)
{
	#As in the case of normal distribution, check lengths
	if(length(n)>1) n<-length(n)

	#Check that mu and lambda are positive
	if(any(mu<=0)) stop("mu must be positive")
	if(any(lambda<0)) stop("lambda must be positive")

	#Check lengths and adjust them recycling
	if(length(mu)>1 && length(mu)!=n) mu<- rep(mu,length=n)
	if(length(lambda)>1 && length(lambda)!=n) lambda = rep(lambda,length=n)

	#Generate random sample from standard normal
	g<-rnorm(n,mean=0,sd=1)

	#Transform to  a sample from chi-squared with 1 df
	v<-g*g

	#Compute roots, equation 5 in reference paper
	#see Fortran code below equation (6)
	w<-mu*v
	cte<-mu/(2*lambda)
	sol1<-mu+cte*(w-sqrt(w*(4*lambda+w)))
	sol2<-mu*mu/sol1

	#Uniform random numbers (0,1)
	u<-runif(n)

	ifelse(u<mu/(mu+sol1),sol1,sol2)
}

