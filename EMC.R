# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <yenchic@andrew.cmu.edu>
# Reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "Enhanced mode clustering." arXiv preprint arXiv:1406.1780 (2014).
# Date: 10/13/2015
#' @import RANN
library(RANN)
#' @import MASS
library(MASS)
#' @import Rcpp
library(Rcpp)

cppFunction('NumericVector MSCpp(NumericMatrix data, NumericMatrix query, double h, int max_iterations , double eps){
	int n = data.nrow(), d = data.ncol(), m = query.nrow();
	double dist_tmp; 
	NumericVector K(n);
	double K_tot;
	
	NumericVector MS(d);
	NumericVector newPt(d);
	NumericVector oldPt(d);
	int iter_now;
	double err_now;
	
	NumericVector result(d*m);
	

	for(int w=0; w<m; w++){
		// for each point
		for(int j = 0; j<d;j++){
			newPt(j) = query(w,j);
		}
		err_now = 1e14;
		iter_now = 0; 
		
		while((iter_now < max_iterations)&&(err_now > eps)){	// ignore those nearly not change
		
			for(int j =0; j<d; j++){
				MS(j) =0;
				oldPt(j) = newPt(j);
			}
			
			K_tot = 0;
			for(int i = 0; i<n; i++){
				dist_tmp = 0;
				for(int j =0; j<d; j++){
					dist_tmp += (data(i,j) - newPt(j))*(data(i,j) - newPt(j))/h/h;
				}
				K(i) = exp(-1*dist_tmp/2);
				K_tot += K(i);
				for(int j =0; j<d; j++){
					MS(j) += data(i,j)*K(i);
				}
			}
				
			// updates & errors
			err_now = 0;
			for(int j =0; j<d; j++){
				newPt(j) = MS(j)/K_tot;
				err_now =+ (newPt(j) - oldPt(j)) * (newPt(j) - oldPt(j));
			}
			err_now = sqrt(err_now);
			iter_now++;
		}
		
		for(int j =0; j<d; j++){
			//result(w*d+j) = newPt(j);
			result(w+m*j) = newPt(j);
		}
	}
	
	return result;

}')


#' Mean shift algorithm using Rcpp
#' 
#' @param data Input data matrix.
#' @param query The mesh points that you want to apply mean shift to.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @return The mesh points after mean shift.
#' @examples
#' x = matrix(rnorm(1000), ncol=2)
#' x_MS = MS(x, x, 0.5)
#'  # mean shift
#' 
#' plot(x, cex=0.5)
#' points(x_MS, col="red", pch=20, cex=2)
#'  # plot the shifted result versus original case
#' 
#' @export
MS = function(data, query, h, max.iterations=200, eps= 1e-15){
	tmp = MSCpp(data=data, query=query, h=h, max_iterations= max.iterations, eps=eps)
	return(matrix(tmp, ncol=ncol(query)))
}


#' Hitting probability for soft mode clustering.
#' @param data Input data points.
#' @param h Smoothing parameters.
#' @param modes Local modes.
#' @return A matrix about the hitting probabaility from each data point to each local modes.
#' @export
HP.soft = function(data,h,modes){
	n.modes = nrow(modes)
	X.m = rbind(modes, data)
		#adding modes to the original data

	Dist.X = as.matrix(dist(X.m))
		#distance matrix
	
	P0.X = exp(-1*Dist.X^2/2/h^2)
	P.X = P0.X/rowSums(P0.X)
		#transition prob.

	T.X = P.X[(n.modes+1): nrow(X.m),(n.modes+1): nrow(X.m)]
	S.X = P.X[(n.modes+1): nrow(X.m),1:n.modes]
	H.X = ginv(diag(,nrow=nrow(T.X),ncol=ncol(T.X))-T.X)%*% S.X
	return(H.X)
}


#' Fast mean shift using heirachical clustering.
#' 
#' @param data Input data matrix.
#' @param query The mesh points that you want to apply mean shift to.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @param cut The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
#' @return The mode clustering result; a list consisting of
#'  \item{label}{the cluster labels for query points.}
#'  \item{modes}{the local modes corresponding to each label.}
#'  
#' @export
fms = function(data, query, h, eps=1.0e-8, max.iterations=100, cut = 0.1){

	result = list()

	pt_ms = MS(data, query, h=h, eps=eps, max.iterations= max.iterations)
	D1 = dist(pt_ms)
	D1_hclust = hclust(D1)
	cluster_lab= cutree(D1_hclust, h=cut*h)
		# find the cluster lable
		
	modes = matrix(NA, nrow=max(cluster_lab), ncol=ncol(data))
	for(i in 1:max(cluster_lab)){
		w_tmp = which(cluster_lab==i)
		if(length(w_tmp)>1){
			modes[i,] = colSums(pt_ms[w_tmp,])/length(w_tmp)
		}
		if(length(w_tmp)==1){
			modes[i,] = pt_ms[w_tmp,]
		}
	}
	result$label = cluster_lab
	result$modes = modes
	
	return(result)
}
	## Note: Here we use hierachical clustering to speed up. This is because at high dimension, it takes a long time for mean shift to completely stop. We can early stop the mean shift. For some clusters, points are closed but haven't arrived local modes. Then we do hierachical clustering which will put these points together since they're close enough.

#' Visualization for mode clustering.
#' 
#' @param data Input data matrix.
#' @param labels Cluster labels for input data.
#' @param modes Local modes for corresonding labels.
#' @param rho The contrast parameter for visualization. Default is to use the method given in Chen et al. (2014).
#' 
#' @return Visualization coordinates. A list consisting:
#'  \item{data}{the visualization coordinates for data points.}
#'  \item{modes}{the visualization coordinates for local modes.}
#'  \item{rho}{the contrast paramter used for visualization.}
#' @export
vis.cluster = function(data, labels, modes, rho=NULL){
	## check if data number is the same as labels and modes is also the same as labels
	result = list()
	
	mds_list = list()
	m_r0 = rep(NA, nrow(modes))
	idx_tmp = NULL
	
	### Stage 1: MDS on modes
	m_mds = cmdscale(dist(modes))
	
	### Stage 2: MDS within clusters
	for(i_c in 1:nrow(modes)){
		w_tmp = which(labels ==i_c)
		data_c = data[w_tmp,]
		data_c_m = rbind(modes[i_c,], data_c)
			# Add mode to each cluster
		
		mds_tmp = cmdscale(dist(data_c_m))
			# MDS within clusters
	
		mds_list[[i_c]] = mds_tmp
		m_r0[i_c] = max(dist(mds_tmp)[1:nrow(mds_tmp[-1,])])
			# size of cluster
		idx_tmp = c(idx_tmp, w_tmp)
	}
	if(!is.null(rho)){
		m_mds_adj = m_mds*rho
	}
	
	if(is.null(rho)){
		rho = 2*max(m_r0)/quantile(dist(m_mds),0.05)
			# Adjustment: to make clusters look separated from each others
		m_mds_adj = m_mds*rho
			# Adjusted locations of modes	
	}
		
	### Placing points around clusters	
	data_tmp = NULL
	for(i_c in 1:nrow(modes)){
		data_c = t(t(mds_list[[i_c]]) - mds_list[[i_c]][1,]+m_mds_adj[i_c,])
		data_c = cbind(data_c, rep(i_c, nrow(data_c)))
		data_tmp = rbind(data_tmp, data_c[-1,])
	}
	result$data = data_tmp[order(idx_tmp),]
	result$modes = m_mds_adj
	result$rho
	return(result)	
}

#' Old method for Enhanced mode clustering (does not create a class).
#' 
#' @param data Input data matrix.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @param n0 The thresholding size for tiny clusters. Default is to use the method given in Chen et al. (2014).
#' @param rho The contrast parameter for visualization. Default is to use the method given in Chen et al. (2014).
#' @param noisy True or False. To desplay noisy clusters (without thresholding). Default is False.
#' @param T_denoise Maximal number of denoising. If tiny clusters presence, we will remove them and redo mode clustering. This is the maximal number of redoing mean shift clustering.
#' @param cut The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
#' 
#' @return Summary informations using enhanced mode clustering. A list consisting:
#'  \item{label}{the cluster labels for query points.}
#'  \item{modes}{the local modes corresponding to each label.}
#'  \item{c.matrix}{the connectivity matrix.}
#'  \item{vis.data}{the visualization coordinates for data points.}
#'  \item{vis.modes}{the visualization coordinates for local modes.}
#'  \item{SC.plot}{the size of ordered clusters before denoising.}
#'  \item{size.threshold}{the size threshold for denoising tiny clusters.}
#'  \item{bandwidth}{the smoothing bandwidth.}
#'  \item{rho}{the contrast paramter used for visualization.}
#'  \item{noisy.label}{the cluster labels for query points before denoising.}
#'  \item{noisy.modes}{the local modes corresponding to each label before denoising.}
#' @export
EMC_old = function(data, h=NULL, max.iterations=1000, eps=1.0e-8, n0= NULL, rho=NULL, cut=0.1, noisy = F, T_denoise =5){
	result = list()
	
	cat("Sample size: ", nrow(data),"\n")
	cat("Data dimension: ", ncol(data),"\n")
	
	if(is.null(h)){
		d = ncol(data)
		n = nrow(data)
		h = (4/(d+4))^(1/(d+6))/n^(1/(d+6))*mean(apply(data,2,sd))
	}
	cat("Smooth bandwidth: ", h,"\n")
	
	### Step 1: original mode clustering
	cat("Step 1: Mode clustering... ")
	cluster = fms(X,X,h=h, cut=cut, eps= eps, max.iterations= max.iterations)
	cluster_lab= cluster$label
	modes = cluster$modes
	cat(" Done.\n")

	### Step 2: denoise tiny clusters
	cat("Step 2: Denoise tiny clusters...\n ")
	if(is.null(n0)){
		n0 = (n*log(n)/20)^((d/(d+6)))
	}
	cluster_lab_tmp = cluster_lab
	cluster_tmp = cluster
	idx_tiny = NULL
	idx_tiny_mode = which(table(factor(cluster_lab_tmp))<n0)
	
	i_count =0
	while(length(idx_tiny_mode)>0){	
		i_count = i_count +1
		for(idx_tmp in idx_tiny_mode){
			w_tmp = which(cluster_lab_tmp==idx_tmp)
			idx_tiny = unique(c(idx_tiny, w_tmp))
		}
		
		cluster_tmp = fms(X[-idx_tiny,],X, h)
		cluster_lab_tmp = cluster_tmp$label
		idx_tiny_mode = which(table(factor(cluster_lab_tmp))<n0)
		cat("Iteration to denoise: ")
		cat(i_count)
		cat("\n")
		if(i_count == T_denoise){
			cat("WARNING: there are still tiny clusters! \n")
			break
		}
	}
	
	cluster_sig = cluster_tmp
	cluster_lab_sig = cluster_sig$label
	modes_sig = cluster_sig$modes
	n_modes_sig = nrow(modes_sig)
	cat(" Done.\n")
	cat("Size  threshold: ", n0,"\n")
	cat("Number of significant clusters: ", n_modes_sig,"\n")

	### Step 3: measuring connectivity 
	cat("Step 3: Measuring connectivity... ")
	soft_hp = HP.soft(X, h, modes_sig)
	Cn_matrix = matrix(NA, n_modes_sig, n_modes_sig)
	for(i.c in 1:n_modes_sig){
		tmp = colSums(soft_hp[which(cluster_lab_sig ==i.c),])/nrow(soft_hp[which(cluster_lab_sig ==i.c),])
		Cn_matrix[i.c,] = tmp
	}
	Cn_matrix = (Cn_matrix+t(Cn_matrix))/2
	cat(" Done.\n")

	### Step 4: visualization coordinates
	cat("Step 4: Dimension reducting... ")
	vis_data = vis.cluster(X, cluster_lab_sig, modes_sig, rho=rho)
	cat(" Done.\n")
	
	result$labels = cluster_lab_sig
	result$modes = modes_sig	
	result$c.matrix = Cn_matrix
	result$vis.data = vis_data$data
	result$vis.modes = vis_data$modes

	result$SC.plot = sort(table(factor(cluster_lab)),decreasing=T)
	result$size.threshold = n0
	result$bandwidth = h
	result$rho = vis_data$rho

	if(noisy){
		result$noisy.labels = cluster_lab
		result$noisy.modes = modes
	}
	
	return(result)
}

#' Creating object EMC.
#' 
#' @param data Input data matrix.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @param ... See \emph{EMC.default}.
#' 
#' @return An S4 object \emph{EMC}; see \emph{EMC.default} for more details.
#' @examples
#' ## Getting olive oils data.
#' library(freqparcoord)
#' data(oliveoils)
#' D=oliveoils
#' 
#' X0= D[,3:10]
#' X = scale(X0)
#' Y = as.numeric(D[,1])
#' 
#' ## EMC
#' X_emc = EMC(X)
#' 
#' ## Visualization plot
#' par(mar=c(0.5,0.5,2,0.5))
#' plot(X_emc, pch=20, cex=0.7,xlab="", ylab="", main="Color by Cluster", xaxt="n", yaxt="n", txt.pos=4)
#'
#' ## Colored by produced area
#' col_sel = rainbow(max(Y))
#' col_sel[4] = "brown"
#' col_sel[3] = "limegreen"
#' par(mar=c(0.5,0.5,2,0.5))
#' plot(X_emc, pch=20, cex=0.7,xlab="", ylab="", main="Color by Produce Area (Connectivity: hitting prob.)", xaxt="n", yaxt="n", col=col_sel[Y], txt.pos=4)
#' legend("topleft", levels(D[,1])[1:4], col=col_sel[1:4],pch=c(19,19,19,19),cex=1.3)
#' legend("bottomleft", levels(D[,1])[5:9], col=col_sel[5:9],pch=c(19,19,19,19,19),cex=1.3)
#' 
#' ## summary statistics
#' X_emc
#' summary(X_emc)
#' 
#' ## SC_plot
#' plot(X_emc$SC.plot , xlim=c(1,20), ylim=c(0,300), ylab="", main="SC-plot (Olive Oil)", xlab="Index of ordered cluster", cex.lab=2, cex.main=2, cex.axis=1.5, pch=19)
#' mtext("Size of cluster", side=2, line=2.2, cex=2)
#' abline(h=X_emc$size.threshold, lwd=2, col="purple")
#' legend("topright",expression((n*log(n)/20)^{frac(d,d+6)}), col="purple", lwd=8, cex=2)
#' @export
EMC = function(data, h=NULL, eps=1.0e-8, max.iterations=100, ...) UseMethod("EMC")

#' Enhanced mode clustering.
#' 
#' @param data Input data matrix.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @param n0 The thresholding size for tiny clusters. Default is to use the method given in Chen et al. (2014).
#' @param rho The contrast parameter for visualization. Default is to use the method given in Chen et al. (2014).
#' @param noisy True or False. To desplay noisy clusters (without thresholding). Default is False.
#' @param T_denoise Maximal number of denoising. If tiny clusters presence, we will remove them and redo mode clustering. This is the maximal number of redoing mean shift clustering.
#' @param cut The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
#' 
#' @return An S4 object about summary informations using enhanced mode clustering. A list consisting:
#'  \item{label}{the cluster labels for query points.}
#'  \item{modes}{the local modes corresponding to each label.}
#'  \item{c.matrix}{the connectivity matrix.}
#'  \item{vis.data}{the visualization coordinates for data points.}
#'  \item{vis.modes}{the visualization coordinates for local modes.}
#'  \item{SC.plot}{the size of ordered clusters before denoising.}
#'  \item{size.threshold}{the size threshold for denoising tiny clusters.}
#'  \item{bandwidth}{the smoothing bandwidth.}
#'  \item{rho}{the contrast paramter used for visualization.}
#'  \item{noisy.label}{the cluster labels for query points before denoising.}
#'  \item{noisy.modes}{the local modes corresponding to each label before denoising.}
#'  @export
EMC.default = function(data, h=NULL, eps=1.0e-8, max.iterations=100, n0= NULL, rho=NULL, cut=0.1, noisy = F, T_denoise =5){
	result = list()
	
	cat("Sample size: ", nrow(data),"\n")
	cat("Data dimension: ", ncol(data),"\n")
	
	if(is.null(h)){
		d = ncol(data)
		n = nrow(data)
		h = (4/(d+4))^(1/(d+6))/n^(1/(d+6))*mean(apply(data,2,sd))
	}
	cat("Smooth bandwidth: ", h,"\n")
	
	### Step 1: original mode clustering
	cat("Step 1: Mode clustering... ")
	cluster = fms(X,X,h=h, cut=cut, eps= eps, max.iterations= max.iterations)
	cluster_lab= cluster$label
	modes = cluster$modes
	cat(" Done.\n")

	### Step 2: denoise tiny clusters
	cat("Step 2: Denoise tiny clusters...\n ")
	if(is.null(n0)){
		n0 = (n*log(n)/20)^((d/(d+6)))
	}
	cluster_lab_tmp = cluster_lab
	cluster_tmp = cluster
	idx_tiny = NULL
	idx_tiny_mode = which(table(factor(cluster_lab_tmp))<n0)
	
	i_count =0
	while(length(idx_tiny_mode)>0){	
		i_count = i_count +1
		for(idx_tmp in idx_tiny_mode){
			w_tmp = which(cluster_lab_tmp==idx_tmp)
			idx_tiny = unique(c(idx_tiny, w_tmp))
		}
		
		cluster_tmp = fms(X[-idx_tiny,],X, h)
		cluster_lab_tmp = cluster_tmp$label
		idx_tiny_mode = which(table(factor(cluster_lab_tmp))<n0)
		cat("Iteration to denoise: ")
		cat(i_count)
		cat("\n")
		if(i_count == T_denoise){
			cat("WARNING: there are still tiny clusters! \n")
			break
		}
	}
	
	cluster_sig = cluster_tmp
	cluster_lab_sig = cluster_sig$label
	modes_sig = cluster_sig$modes
	n_modes_sig = nrow(modes_sig)
	cat(" Done.\n")
	cat("Size  threshold: ", n0,"\n")
	cat("Number of significant clusters: ", n_modes_sig,"\n")

	### Step 3: measuring connectivity 
	cat("Step 3: Measuring connectivity... ")
	soft_hp = HP.soft(X, h, modes_sig)
	Cn_matrix = matrix(NA, n_modes_sig, n_modes_sig)
	for(i.c in 1:n_modes_sig){
		tmp = colSums(soft_hp[which(cluster_lab_sig ==i.c),])/nrow(soft_hp[which(cluster_lab_sig ==i.c),])
		Cn_matrix[i.c,] = tmp
	}
	Cn_matrix = (Cn_matrix+t(Cn_matrix))/2
	cat(" Done.\n")

	### Step 4: visualization coordinates
	cat("Step 4: Dimension reducting... ")
	vis_data = vis.cluster(X, cluster_lab_sig, modes_sig, rho=rho)
	cat(" Done.\n")
	
	result$labels = cluster_lab_sig
	result$modes = modes_sig	
	result$c.matrix = Cn_matrix
	result$vis.data = vis_data$data
	result$vis.modes = vis_data$modes

	result$SC.plot = sort(table(factor(cluster_lab)),decreasing=T)
	result$size.threshold = n0
	result$bandwidth = h
	result$rho = vis_data$rho

	if(noisy){
		result$noisy.labels = cluster_lab
		result$noisy.modes = modes
	}
	
	class(result) = "EMC"
	return(result)
}

#' Print EMC object.
#' @param x An 'EMC' object.
print.EMC = function(x, ... ){
	cat("Clustering Labels:\n")
	print(x$labels)
	cat("\n")
	cat("Local Modes:\n")
	print(x$modes)
	cat("\n")
	cat("Connectivity Matrix:\n")
	print(x$c.matrix)
	cat("\n")	
}
#' Summary of an EMC object.
#' @param x An 'EMC' object.
summary.EMC = function(x, ...){
	cat("Clusters Size:")
	print(table(x$label))
	cat("\n")
	cat("SC-plot (ordered cluster size before denoising):")
	print(x$SC.plot)
	cat("\n")
	cat("Connectivity Matrix:\n")
	print(x$c.matrix)
	cat("\n")	
	cat("Smoothing Bandwidth:\n")
	print(x$bandwidth)
	cat("\n")	
	cat("Threshold for Denoising:\n")
	print(x$size.threshold)
	cat("\n")	
}


#' Visualization plot for EMC object.
#' @param x An 'EMC' object.
#' @param Cn_lv Connectivity level to be displayed. Default is to automatically select by number of clusters.
#' @param col Color for each data point.
#' @param c.lwd 'lwd' for the lines connecting clusters.
#' @param mode.pch 'pch' for modes.
#' @param mode.cex 'cex' for modes.
#' @param mode.col 'col' for modes.
#' @param txt.pos 'pos' for showing clusters ('pos' is from function \emph{text}).
#' @param txt.cex 'cex' for showing clusters ('cex' is from function \emph{text}).
#' @param txt.offset 'offset' for showing clusters ('offset' is from function \emph{text}).
plot.EMC = function(x, Cn_lv = NULL, col = NULL, c.lwd = 20, mode.pch = 19, mode.cex = 0.6, mode.col = "gray", txt.pos = 2, txt.cex =2.2, txt.offset = 1.5, ... ){
		if(is.null(col))
			col = rainbow(nrow(x$modes))[x$vis.data[,3]]
			
		if(is.null(Cn_lv))
			Cn_lv = 1/nrow(x$modes)/2
			
	plot(x$vis.data[,1:2], col=col, ...)
		# Connectivity measures
		for(i in 1: (nrow(x$modes)-1)){
			for(j in (i+1): nrow(x$modes)){
				if(x$c.matrix[i,j]>Cn_lv){
					segments(x0=x$vis.modes[i,1], y0 = x$vis.modes[i,2], x1=x$vis.modes[j,1], y1=x$vis.modes[j,2], lwd=x$c.matrix[i,j]*c.lwd)
				}
			}
		}	
		points(x$vis.modes, pch=mode.pch, cex=mode.cex, col=mode.col)
		for(i_c in 1: nrow(x$modes)){
			text(x= x$vis.modes[i_c,1], y= x$vis.modes[i_c,2], labels=i_c, pos= txt.pos, cex=txt.cex, offset=txt.offset)
		}
}

