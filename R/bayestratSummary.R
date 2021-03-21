#' bayestratSummary
#'
#'\code{bayestratSummary} summarizes the output from \code{bayestrat}.
#'
#' @param output The output from \code{bayestrat}.
#'
#' @return A list with the following components:
#' \describe{
#' \item{snp}{A data frame of the scaled effect (effect size/standard error), effect size, standard error and confidence level for significant SNPs.}
#' \item{cov}{A similar data frame for significant covariates.}
#' \item{pc}{A similar data frame for significant PCs.}
#' }
#' @export
#'

bayestratSummary<-function(output){
	alpha<-as.numeric(sub("%","",names(output$CIsnp)))/100
	snp<-data.frame()
	c<-NULL
	for(j in 3:1){
		SD<-output$SDsnp
		CI<-output$CIsnp[[j]]
		if(is.null(c)==F){
			tmp<-match(c,names(SD))
			SD<-SD[-tmp]
			CI<-CI[-tmp,]
		}
		signs<-sign(CI)
		signs<-apply(signs,1,prod)
		if(any(signs>=0)){ #signs=0 or 1 means snp is significant
			s<-which(signs>=0)
			tmp<-cbind(apply(CI[s,],1,mean),SD[s],alpha[j])
			tmp<-cbind(tmp[,1]/tmp[,2],tmp)
			orders<-order(abs(tmp[,1]),decreasing=T)
			if(length(s)>1){
				tmp<-tmp[orders,]
			}
			snp<-rbind(snp,tmp)
			c<-c(c,names(s))
		}
	}
	if(nrow(snp)>0){
		colnames(snp)<-c("ScaledEffect","EffectSize","SD","alpha")
	}

	if(is.null(output$SDcov)==F){
		cov<-data.frame()
		c<-NULL
		for(j in 3:1){
			SD<-output$SDcov
			CI<-output$CIcov[[j]]
			if(is.null(c)==F){
				SD<-SD[-c]
				CI<-CI[-c,]
			}
			signs<-sign(CI)
			signs<-apply(signs,1,prod)
			if(any(signs>=0)){ #signs=0 or 1 means cov is significant
				s<-which(signs>=0)
				tmp<-cbind(apply(CI[s,],1,mean),SD[s],alpha[j])
				tmp<-cbind(tmp[,1]/tmp[,2],tmp)
				orders<-order(abs(tmp[,1]),decreasing=T)
				if(length(s)>1){
					tmp<-tmp[orders,]
				}
				cov<-rbind(cov,tmp)
				c<-c(c,s)
			}
		}
		if(nrow(cov)>0){
			colnames(cov)<-c("ScaledEffect","EffectSize","SD","alpha")
		}
	}

	if(is.null(output$SDpc)==F){
		pc<-data.frame()
		c<-NULL
		for(j in 3:1){
			SD<-output$SDpc
			CI<-output$CIpc[[j]]
			if(is.null(c)==F){
				SD<-SD[-c]
				CI<-CI[-c,]
			}
			signs<-sign(CI)
			signs<-apply(signs,1,prod)
			if(any(signs>=0)){ #signs=0 or 1 means pc is significant
				s<-which(signs>=0)
				tmp<-cbind(apply(CI[s,],1,mean),SD[s],alpha[j])
				tmp<-cbind(tmp[,1]/tmp[,2],tmp)
				orders<-order(abs(tmp[,1]),decreasing=T)
				if(length(s)>1){
					tmp<-tmp[orders,]
				}
				pc<-rbind(pc,tmp)
				c<-c(c,s)
			}
		}
		if(nrow(pc)>0){
			colnames(pc)<-c("ScaledEffect","EffectSize","SD","alpha")
		}
	}

	result<-list(snp=snp)
	if(is.null(output$SDcov)==F){
		result$cov<-cov
	}
	if(is.null(output$SDpc)==F){
		result$pc<-pc
	}

	return(result)

}
