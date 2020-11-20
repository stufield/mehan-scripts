


library(MASS)
library(lattice)

plot_3d_hist <- function(xdata, ydata, xbins=20, ybins=20,  screen = list(z = 250, x=-60), ...) {
	data=data.frame(cbind(ydata,xdata))
	h=hist(data[,2], plot=F, xbins)
	xbreaks = h$breaks
	h=hist(data[,1], plot=F, ybins)
	ybreaks = h$breaks
	#plot_data = matrix(0, nrow=length(xbreaks)-1, ncol=length(ybreaks)-1, dimnames=list(xbreaks[2:length(xbreaks)], ybreaks[2:length(ybreaks)]))
	plot_data = matrix(0, nrow=0, ncol=3, dimnames=list(c(), c("X", "Y", "Z")))
	for (i in seq(2,length(xbreaks))) {
		cur_data = data[data[,2]<xbreaks[i],]
		cur_data = cur_data[cur_data[,2]>xbreaks[i-1],]
		h=hist(cur_data[,1], plot=F, breaks=ybreaks)
		for (j in seq(2,length(ybreaks))) {
			plot_data = rbind(plot_data, c(xbreaks[i], ybreaks[j], h$counts[j-1]))
			#plot_data[i-1,j-1] = h$counts[j-1]
		}
		
	}
	wireframe(Z~ X* Y, data=data.frame(plot_data), drape=T, pretty=T, scales=list(arrows=FALSE), col.regions=(rainbow(1000, start=.1, end=.05)), ylab="Y", xlab="X", zlab="Frequency", screen=screen,...)
	#return(data.frame(plot_data))
}
	
