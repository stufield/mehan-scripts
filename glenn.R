library(splines)

load.assay = function(filename) {
	data = read.table(filename, row.names=T, header=T, sep=",")
}

get.assay.cols = function(assay.data) {
	cols = names(assay.data)
	cols = cols[cols!="Aptamer"]
	cols = cols[cols!="Target"]
	cols = cols[cols!="Probe.Type"]
	cols = cols[cols!="Diluting.Factor"]
	cols
}

my.plot.matrix = function(mat, plot.type="Scatter", x.axis=NULL, labels = NULL, cols=NULL, fit.type="None", ylim=NULL, xlim=NULL, legend.pos="topright", xgrid=NULL, ygrid=NULL, ...) {
	mat[mat==-Inf] = NA
	mat = as.data.frame(mat)
	if (is.null(x.axis)) 
		x.axis = 1:ncol(mat)
	if (is.null(cols)) 
		cols = 2:(nrow(mat)+1)
	
	if (is.null(labels)) {
		#labels = names(mat)
		if (is.null(labels)) {
			labels = as.character(x.axis)
		}
	}
	if (is.null(xlim)) {
		xlim = range(x.axis)
	}
	else {
		if (xlim[1] == -1234)
			xlim[1] = min(x.axis)
		if (xlim[2] == -1234)
			xlim[2] = max(x.axis)
	}
	if (is.null(ylim)) {
		ylim = get.lim(mat)
		ylim[1] = ylim[1]*.95
		ylim[2] = ylim[2]*1.05
	}
	else {
		if (ylim[1] == -1234)
			ylim[1] = min(mat)
		if (ylim[2] == -1234)
			ylim[2] = max(mat)
	}
	print("X Grid")
	print(xgrid)
	if (!is.null(xgrid)) {
		print(c(xlim, 1/(xgrid/100)+1))
		print(seq(xlim[1], xlim[2], length.out=1/(xgrid/100)))
	}
	print("Y Grid")
	print(ygrid)
	print("hmm")
	if (!is.null(ygrid)) {
		print(c(ylim, 1/(ygrid/100)+1))
		print(seq(ylim[1], ylim[2], length.out=1/(ygrid/100)))
	}
	if (plot.type == "Scatter") {
		par(mar=c(5,4,4,2), las=1)
		for (i in 1:nrow(mat)) {
			all.y = as.numeric(mat[i,])
			y.vals = all.y[is.finite(all.y)]
			x.vals = x.axis[is.finite(all.y)]
			if (i == 1) {
				plot.default(-100000, -1000000, axes=F, frame.plot=T, pch=19, type='p', ylim=ylim, xlim=xlim, ...)
			}
			if (!is.null(xgrid))  {
					
				abline(v=seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])*xgrid/100), col=8,lty=3)
			}
			if (!is.null(ygrid))  {
			
				abline(h=seq(ylim[1], ylim[2], by=(ylim[2]-ylim[1])*ygrid/100), col=8,lty=3)
			}

			points(x.vals, y.vals, col=cols[i], type='p', pch=19)
			if (fit.type == "Interpolate (Straight)") {
				#cat("Performing Interpolated Fit\n")
				points(x.vals, y.vals, col=cols[i], lty=2, type='l')
			}
			if (fit.type == "Interpolate (Smooth)") {
				cat("Performing Interpolated Fit\n")
				lines(predict(interpSpline(x.vals, y.vals)), col=cols[i], lty=2)
			}
			else if (fit.type == "Fit (LOESS)") {
				#cat("Performing Smooth Fit\n")
				smooth1 = try(points(loess.smooth(x.vals, y.vals), type='l', col=cols[i]), silent=T)
			}
			#else 
				#cat("No fitting selected\n")
		}
		axis(1, at=x.axis, labels = labels)
		axis(2)
		if (nrow(mat) > 1)
			legend(legend.pos, row.names(mat), col=cols, pch=19, lty=1, lwd=1, bg='white')
	}
	else if (plot.type == "Barplot") {
		par(mar=c(10,4,4,2), las=2)
		max.height = max(apply(mat, 2, sum))
		barplot(as.matrix(mat), col=cols,  ylim=c(0, max.height*1.3), names.arg=labels, ...)
		if (nrow(mat) > 1)
			legend(legend.pos, row.names(mat), col=cols, lwd=4)
	}
			
}


plot.titrations.adat = function(adat.data, sample.ids.list, x.axis=NULL, labels=NULL, preview=F, path="./", xlog=F, ylog=F, width=5, height=5, par.list=list(), ...) {
	apts = get.aptamers(adat.data)
	
	if (substr(path, nchar(path), nchar(path)) != "/")
		path=paste(path, "/", sep="")

	if (is.null(labels)) 
		labels = as.character(sample.ids.list[[1]])
	if (is.null(names(sample.ids.list)))
		stop("Must have names in sample.ids.list")
	if (!is.null(x.axis) & xlog)
		x.axis = log(x.axis, base=10)
	for (apt in apts) {
		mat = NULL
		apt.data = adat.data[,c("SampleId", apt)]
		for (series in names(sample.ids.list)) {
			sample.ids = sample.ids.list[[series]]
			mat = rbind(mat, sapply(sample.ids, function(sample.id) apt.data[apt.data$SampleId == sample.id,apt]))
		}
		#mat = t(mat)
		dimnames(mat) = list(names(sample.ids.list), labels)

		filename = paste(path, apt, ".png", sep="")
		if (!preview)
			plot.to.file(file=filename, height, width);
		#par(mar=c(15,4,4,2))
		par(par.list)
		log.str = ""
		if (ylog)
			log.str = paste(log.str,"y", sep="")
			#mat = log(mat, base=10)
		my.plot.matrix(mat, x.axis=x.axis, labels=labels, main=apt, ...)
		if (preview)
			break
		else
			dev.off()
	}
}
plot.titrations = function(data, column.list, x.axis=NULL, labels=NULL, preview=F, path="./", xlog=F, ylog=F, ...) {
	
	if (substr(path, nchar(path), nchar(path)) != "/")
		path=paste(path, "/", sep="")

	if (is.null(labels)) 
		labels = as.character(column.list[[1]])
	if (is.null(names(column.list)))
		stop("Must have names in column.list")
	if (!is.null(x.axis) & xlog)
		x.axis = log(x.axis, base=10)
	for (i in row.names(data)) {
		mat = matrix(ncol=0, nrow=length(labels))
		for (col.name in names(column.list)) {
			columns = column.list[[col.name]]
			mat = cbind(mat, as.numeric(data[i,columns]))
		}
		mat = t(mat)
		dimnames(mat) = list(names(column.list), labels)

		if ("Aptamer" %in% names(data)) 
			apt = data[i,"Aptamer"]
		else
			apt = as.character(i)
		filename = paste(path, gsub("[+/]", ".", data[i,"Target"]), "_",  apt, ".png", sep="")
		if (!preview)
			png(file=filename, width=640, height=640)
		#par(mar=c(15,4,4,2))
		title = paste(data[i,"Target"], " (", apt, ")", sep="")
		if (ylog)
			mat = log(mat, base=10)
		print(title)
		my.plot.matrix(mat, x.axis=x.axis, labels=labels, main=title, ...)
		if (preview)
			break
		else
			dev.off()
	}
}

	

plot.titration = function(data1, names, data2=NULL, ...) {
	cols =  get.assay.cols(data1)
	print(cols)
	if (!is.null(data2)) {
		cols2 = get.assay.cols(data2)
		if (! all(cols == cols2))
			stop("Data1 and Data2 don't match")
	}
	labels = lapply(cols, kill.X)
	cat = F
	if (all(labels==cols)) # Not numbers
		cat = T
	else
		x.axis = as.numeric(labels)
	for (i in row.names(data1)) {
		all.data = as.numeric(data1[i,cols])
		if (!is.null(data2)) 
			all.data = c(all.data, as.numeric(data2[i,cols]))
		png(file=paste(data1[i,"Target"], "_",  i, ".png", sep=""), width=640, height=640)
		par(mar=c(15,4,4,2))
		title = paste(data1[i,"Target"], " (", i, ")", sep="")
		# Ex1
		if (cat) {
			barplot(as.matrix(data1[i,cols]), main=title, ylab="RFU", col='blue', las=2)
		}
		else {
			plot.default(x.axis, as.numeric(data1[i,cols]), col='red', axes=F, frame.plot=T, main=title, ylab="RFU", pch=19, type='p', ylim=get.lim(all.data), ...)
			smooth1 = try(points(loess.smooth(x.axis, data1[i,cols]), type='l', col='red'), silent=T)
			if (!is.null(smooth1)){
				lines(predict(interpSpline(x.axis, data1[i,cols])), col='red', lty=2)
				print(sprintf("No loess for %s on data2", i))
			}

			axis(1, at=x.axis, labels = labels)
			axis(2)
			if (!is.null(data2))  {
				points(x.axis, data2[i,cols], col='blue', type='p', pch=19)
				smooth2 = try(points(loess.smooth(x.axis, data2[i,cols]), type='l', col='blue'), silent=T)
				if (!is.null(smooth2)) {
					lines(predict(interpSpline(x.axis, data2[i,cols])), col='blue', lty=2)
					print(sprintf("No loess for %s on data2", i))
				}
			}
			legend("topright", names, col=c('red', 'blue'), pch=19, lty=1, lwd=1)
		}

		dev.off()
	}


}

kill.X = function(name) {
	if (substr(name, 1, 1) == "X") 
		return(substring(name, 2))
	name
}


plot.case.control = function(filename) {
	shear.time.spin = read.csv(filename, header=T, row.names=1)
	shear = 3:10
	control = 11:18
	labels = c("0", "0.5", "1", "2", "4", "8", "12", "20")
	x.axis = as.numeric(labels)
	for (i in seq(1, dim(shear.time.spin)[1])) {
		#png(file=paste("shear_spin_plots/",sub("/", "_", shear.time.spin[i,"Target"]), "_",  row.names(shear.time.spin)[i], ".png", sep=""), width=640, height=640)
		shear.data = as.numeric(shear.time.spin[i,shear])
		control.data = as.numeric(shear.time.spin[i,control])
		title = paste(shear.time.spin[i,"Target"], " (", row.names(shear.time.spin)[i], ")", sep="")
		# Ex1
		plot.default(x.axis, shear.data, col='red', axes=F, frame.plot=T, main=title, xlab="Time (hours)", ylab="RFU", pch=19, type='o', ylim=get.ylim(shear.data, control.data))
		axis(1, at=x.axis, labels = labels)
		axis(2)
		points(x.axis, control.data, col='blue', pch=19, type='o')

		# Ex1


		legend("topright", c("Shear", "Control"), col=c('red', 'blue'), pch=19, lty=1, lwd=1)
		#dev.off()
		break
		
	}
}

plot.rates = function(filename, x.vals, cols1, path, cols2=NULL, do.plot=T...) {
	data = read.csv(filename, header=T, row.names=1, sep="\t")
	k.values = c()
	k.values2 = c()
	labels = as.character(x.vals)
	for (i in seq(1, dim(data)[1])) {
		png(file=paste(path, "/",sub("/", "_", data[i,"Target"]), "_",  row.names(data)[i], ".png", sep=""), width=640, height=640)

		# Plot
		plot.data = data.frame(x=x.vals, y=as.numeric(data[i,cols1]))
		title = paste(data[i,"Target"], " (", row.names(data)[i], ")", sep="")
		if (do.plot) {
			plot.default(plot.data, col='red', axes=F, frame.plot=T, main=title, ylab="RFU", pch=19, ...)
			axis(1, at=x.vals, labels = labels)
			axis(2)
		}

		# Model
		cat(noquote(paste(i,"Fitting model for", title, "\n")))
		data.nls = NULL
		data.nls = try(nls( y ~ (max-min)*(1-exp(-k*x))+min, data=plot.data, start=list(max=max(plot.data$y), min=min(plot.data$y), k=3), control=nls.control(maxiter=200, minFactor=0, warnOnly=T)))
		if (class(data.nls) != "try-error") {
			max = coef(data.nls)[1]
			min = coef(data.nls)[2]
			k = coef(data.nls)[3]
			k.values = c(k.values, k)
			if (do.plot) {
				lines(plot.data$x, (max-min)*(1-exp(-k*plot.data$x))+min, col='red')
				
				text(max(plot.data$x), min(plot.data$y), labels=paste("k=", sprintf("%0.2f", k)), cex=1.5, col='red', adj=c(1.2,-1.5))
			}
		}

		if (!is.null(cols2)) {
			print("Support for cols2 not implemented")
			return(NULL)
		}

		#legend("topleft", c("Shear", "Control"), col=c('red', 'blue'), pch=19, lty=1, lwd=1)

		dev.off()
		
	}
	k.values
}


plot.single = function(filename, path) {
	data = read.csv(filename, header=T, row.names=1, sep="\t")
	shear = 4:9
	labels = c("0.25", "0.5", "0.75", "1", "1.25", "1.5")
	x.axis = as.numeric(labels)
	for (i in seq(1, dim(data)[1])) {
		png(file=paste(path, "/",sub("/", "_", data[i,"Target"]), "_",  row.names(data)[i], ".png", sep=""), width=640, height=640)
		plot.data = as.numeric(data[i,shear])
		title = paste(data[i,"Target"], " (", row.names(data)[i], ")", sep="")
		# Ex1
		plot.default(x.axis, plot.data, col='red', axes=F, frame.plot=T, main=title, xlab="[NHS-biotin], (mM)", ylab="RFU", pch=19, type='o')
		axis(1, at=x.axis, labels = labels)
		axis(2)

		# Ex1


		#legend("topleft", c("Shear", "Control"), col=c('red', 'blue'), pch=19, lty=1, lwd=1)
		dev.off()
		
	}
}


set.dev.size = function(width=5, height=5) {
	if (dev.cur() > 1)
		dev.off()
	x11(width=width, height=height)
}




