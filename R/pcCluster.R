#' pcCluster
#'
#' \code{pcCluster} performs individual clustering based on PC scores and generates a 3D plot in a pop-up window.
#'
#' @param data.test (opional) A data frame consisting of a continuous phenotype, covariates and genotypes for candidate SNPs as columns, individuals as rows. If \code{race=TRUE}, this is necessary and has to include a column named "Race".
#' @param data.pc A data frame or matrix of PC scores. Rows represent individuals. Columns represent PCs.
#' @param race TRUE or FALSE, whether to comapre clustering results to race categories.
#' @param cluster TRUE or FALSE, whether to cluster individuals by the k-means method.
#' @param cluster.whichPC A numerical vector, on which PC(s) the clustering is based. Necessary if \code{cluster=TRUE}.
#' @param n.cluster The number of clusters for the k-means method.
#' @param iter.max The maximum number of iterations allowed in the k-means algorithm.
#' @param algorithm The algorithm used in k-means clustering. One of "Hartigan-Wong", "Lloyd", "Forgy" and "MacQueen". Default is "Hartigan-Wong".
#' @param plotPC TRUE or FALSE, whether to generate a 3D PC plot.
#' @param plot.whichPC  A vector of length three containing which columns of \code{data.pc} to be plotted. Default is (1,2,3).
#' @param colorby One of "clusters", "race" or "none", whether to color the 3D PC plot by the clustering result, the race information or no color.
#' @param colors A vector of color names to be used for plots.
#'
#' @return If \code{cluster=T}, returns the object from \code{kmeans}, and a confusion matrix is also included if \code{race=T}.
#' @export
#' @import rgl
#'
#' @examples
#'data(data.test)
#'data(pc)
#'out<-pcCluster(data.test,data.pc=pc,race=TRUE,cluster=TRUE,n.cluster=4,
#'plotPC=FALSE,plot.whichPC=c(1,2,3),colorby="clusters")
#'
#'

pcCluster<-function(data.test,data.pc,race=T,cluster=T,cluster.whichPC=c(1,2,3),n.cluster,
										iter.max=10,algorithm="Hartigan-Wong",plotPC=T,plot.whichPC=c(1,2,3),
										colorby="clusters",colors=NULL){

	if(race==T){
		Race<-data.test$Race
	}else{
		Race<-NULL
	}

	if(cluster==T){
		pc.sub<-data.pc[,cluster.whichPC]
		result<-kmeans(pc.sub,centers=n.cluster,iter.max=iter.max,algorithm=algorithm)
		if(is.null(Race)==F){
			#confusion matrix
			clusters<-as.factor(result$cluster)
			Race<-as.factor(Race)
			lc<-levels(clusters)
			lr<-levels(Race)
			tab<-table(Race,clusters)
			if(length(lr)<=length(lc)){  #change the order of columns/rows of tab so that max values in each row/column is on the diagonal line
				for(i in 1:nrow(tab)){ #change order of columns
					m<-which.max(tab[i,])
					if(m!=i){
						tmp<-tab[,m]
						tab[,m]<-tab[,i]
						tab[,i]<-tmp
						tmp<-lc[m]
						lc[m]<-lc[i]
						lc[i]<-tmp
					}
				}
				colnames(tab)<-lc
			}else{ #change order of rows
				for(i in 1:ncol(tab)){
					m<-which.max(tab[,i])
					if(m!=i){
						tmp<-tab[m,]
						tab[m,]<-tab[i,]
						tab[i,]<-tmp
						tmp<-lr[m]
						lr[m]<-lr[i]
						lr[i]<-tmp
					}
				}
				rownames(tab)<-lr
			}
			tab<-cbind(tab,apply(tab,1,sum))
			colnames(tab)[ncol(tab)]<-"total"
			tab<-rbind(tab,apply(tab,2,sum))
			rownames(tab)[nrow(tab)]<-"total"
			result$confusion.matrix<-tab
		}
	}
	if(plotPC==T){
		pc.sub<-data.pc[,plot.whichPC]
		xlim<-c(min(pc.sub[,1])-0.1,max(pc.sub[,1])+0.1)
		ylim<-c(min(pc.sub[,2])-0.1,max(pc.sub[,2])+0.1)
		zlim<-c(min(pc.sub[,3])-0.1,max(pc.sub[,3])+0.1)
		if(is.null(colnames(pc.sub))){
			pc.name<-c(paste("PC",plot.whichPC[1],sep=""),paste("PC",plot.whichPC[2],sep=""),paste("PC",plot.whichPC[3],sep=""))
		}else{
			pc.name<-colnames(pc.sub)
		}
		options(rgl.printRglwidget = TRUE)
		if(colorby=="none"){
			plot3d(pc.sub[,1],pc.sub[,2],pc.sub[,3],xlim=xlim,ylim=ylim,zlim=zlim,
						 xlab=pc.name[1],ylab=pc.name[2],zlab=pc.name[3],main="")
		}else if(colorby=="clusters"){ #color by clustering results from k-means
			if(is.null(colors)==T){
				colors=c(1:n.cluster)
			}
			pop.name<-as.character(c(1:n.cluster))
			pop<-list()
			for(i in 1:n.cluster){
				pop[[i]]<-which(result$cluster==i)
			}
			plot3d(pc.sub[pop[[1]],1],pc.sub[pop[[1]],2],pc.sub[pop[[1]],3],
						 xlim=xlim,ylim=ylim,zlim=zlim,xlab=pc.name[1],ylab=pc.name[2],
						 zlab=pc.name[3],col=colors[1],main="")
			for(i in 2:n.cluster){
				points3d(pc.sub[pop[[i]],1],pc.sub[pop[[i]],2],pc.sub[pop[[i]],3],col=colors[i])
			}
			legend3d("topright",col=c(1:n.cluster),pch=rep(20,n.cluster),legend=pop.name)
		}else if(colorby=="race"){ #color by race info
			pop.name<-names(table(Race))
			n.pop<-length(pop.name)
			pop<-list()
			for(i in 1:n.pop){
				pop[[i]]<-which(Race==pop.name[i])
			}
			if(is.null(colors)==T){
				colors=c(1:n.pop)
			}
			plot3d(pc.sub[pop[[1]],1],pc.sub[pop[[1]],2],pc.sub[pop[[1]],3],
						 xlim=xlim,ylim=ylim,zlim=zlim,xlab=pc.name[1],ylab=pc.name[2],
						 zlab=pc.name[3],col=colors[1],main="")
			for(i in 2:n.pop){
				points3d(pc.sub[pop[[i]],1],pc.sub[pop[[i]],2],pc.sub[pop[[i]],3],col=colors[i])
			}
			legend3d("topright",col=c(1:n.pop),pch=rep(20,n.pop),legend=pop.name)
		}
		browseURL(paste("file://", writeWebGL(dir=file.path(tempdir(), "webGL"), width=500), sep=""))
	}

	return(result)
}
