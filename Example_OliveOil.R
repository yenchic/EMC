# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <yenchic@andrew.cmu.edu>
# Reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "Enhanced mode clustering." arXiv preprint arXiv:1406.1780 (2014).
# Date: 04/07/2015
source("EMC.R")


##### Oliveoils data
library(freqparcoord)
data(oliveoils)
D=oliveoils

head(D)

X0= D[,3:10]
X = scale(X0)
Y = as.numeric(D[,1])


## (1) Enhanced mode clustering
X_emc = EMC(X)




## (2) Visualization
par(mar=c(0.5,0.5,2,0.5))
plot(X_emc, pch=20, cex=0.7,xlab="", ylab="", main="Color by Cluster", xaxt="n", yaxt="n", txt.pos=4)





## (3) plot for produce area
col_sel = rainbow(max(Y))
col_sel[4] = "brown"
col_sel[3] = "limegreen"

par(mar=c(0.5,0.5,2,0.5))
plot(X_emc, pch=20, cex=0.7,xlab="", ylab="", main="Color by Produce Area (Connectivity: hitting prob.)", xaxt="n", yaxt="n", col=col_sel[Y], txt.pos=4)
legend("topright", levels(D[,1])[1:4], col=col_sel[1:4],pch=c(19,19,19,19),cex=1.3)
legend("bottomright", levels(D[,1])[5:9], col=col_sel[5:9],pch=c(19,19,19,19,19),cex=1.3)





## (4) Summary statistics
names(X_emc)

X_emc

summary(X_emc)

	# Significant modes
X_emc$modes

	# Cluster label
X_emc$labels

	# Connectivity matrix
X_emc$c.matrix

	# Visualization coordinates
X_emc$vis.data

	# Modes in visualization coordinates
X_emc$vis.modes





## (5) SC_plot
plot(X_emc$SC.plot , xlim=c(1,20), ylim=c(0,300), ylab="", main="SC-plot (Olive Oil)", xlab="Index of ordered cluster", cex.lab=2, cex.main=2, cex.axis=1.5, pch=19)
mtext("Size of cluster", side=2, line=2.2, cex=2)
abline(h=X_emc$size.threshold, lwd=2, col="purple")
legend("topright",expression((n*log(n)/20)^{frac(d,d+6)}), col="purple", lwd=8, cex=2)
	# This leaves 7 clusters







## (6) Confusion matrix
CC_matrix = table(Y, X_emc$labels)
row.names(CC_matrix)<-levels(D[,1])
colnames(CC_matrix)<-c(1: nrow(X_emc$modes))
CC_matrix


#### Remark: One can use the above as a template to analyze other dataset. This gives a procedure of high dimensional mode clustering.
