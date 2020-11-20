# pca() - various functions to prove PCA functionality
#
# Revision History
#---------------------
# $Id: $: $Date: $ 
# 

library(MASS)
library(caTools)
library(Matrix)


pca_wrapper = function(training.data, plot=T, num.pcs=3, dims=1:2,
                       test.data=NULL, max.disease.auc=1.1, do.log=T,
                       layout=matrix(1:9, nrow=3, byrow=F), apt.color=NULL,
                       center=T, scale=T, ...) {

	par(mar=c(4,2.5, 2.5,2))
	if (plot)
		layout(layout)
	if (!"Response" %in% names(training.data))
		response = factor(1)
	# Save some variables for later use
	else
		response = factor(training.data$Response) # Response
	control = response[1]
	cat("Peeling", as.character(control), "\n")
	disease = response[nrow(training.data)]
	which1 = response == control 
	which2 = response == disease


	# Do initial PCA on full dataset for reference
	apts = get.aptamers(training.data)
	if (length(apts) == 0)
		apts = 1:(ncol(training.data)-1)
	train.norm = center.scale.data(training.data[,apts], center=center, scale=scale, do.log=do.log)
	train.pca = prcomp(train.norm, center=F, scale.=F)


	# Perform PCA on controls for peeling	
	train.control.data = training.data[which1,apts]
	train.control.norm = center.scale.data(train.control.data, center=center, scale=scale, do.log=do.log)
	train.control.pca = prcomp(train.control.norm, center=F, scale.=F)
	train.control.pca.var = train.control.pca$sdev^2/sum(train.control.pca$sdev^2)
	# If number of PCs to peel not provided, do heuristic approach
	if (is.null(num.pcs)) {
		for (i in 1:(length(train.control.pca.var)-1)) {
			if (train.control.pca.var[i]/train.control.pca.var[i+1] > 2.5)
				break
		}
		num.pcs = i	
	}
	# Project the whole dataset, transformed by the controls, into the space of the controls
	train.peel = apply.center.scale(train.control.data, training.data, center=center, scale=scale, do.log=do.log) # The dataset we will eventually peel
	train.peel = train.peel[,get.aptamers(train.peel)]
	if ("Response" %in% names(train.peel))
		train.peel = train.peel[,1:(ncol(train.peel)-1)]
	train.control.pca.proj =  as.matrix(train.peel) %*% train.control.pca$rotation

	# Perform Peeling
	peeled.dims = c() # Used to track which PCs are peeled 
	peel.aucs = c() # List of aucs associated with PCs
	for (peel.col in 1:ncol(train.peel)) {
		if (length(peeled.dims) >= num.pcs) 
			break
		if (length(levels(training.data$Response)) > 1) {
			auc = colAUC(train.control.pca.proj[,peel.col], training.data$Response)
		}
		else
			auc=.5
		peel.aucs = c(peel.aucs, auc)
		if (auc < max.disease.auc) {
			train.peel =  t(apply(train.peel, 1, pca.peel, train.control.pca$rotation[,peel.col]))
			peeled.dims = c(peeled.dims, peel.col)
		}
		else {
			cat(sprintf("Skipped Component %i with auc of %0.3f\n", peel.col, auc))
		}

	}
	train.peel = as.data.frame(train.peel)
	train.peel.pca = prcomp(train.peel, scale.=F, center=F)
	train.peel.data = as.data.frame(apply.center.scale(train.control.data, train.peel, undo=T, do.log=do.log, center=center, scale=scale))
	train.peel.data = cbind(training.data[,get.meta(training.data)], train.peel.data, Response=training.data$Response)
	aucs = c()
	for (i in 1:min(100, ncol(train.peel.pca$x))) {
		auc = pca.auc(training.data, train.peel.pca, i)
		aucs = c(aucs, auc)
	}

	#test.data = training.data
	test.peel.pca = NULL
	if (!is.null(test.data)) {
		test.peel = as.matrix(apply.center.scale(train.control.data, test.data[,get.aptamers(test.data)], center=center, scale=scale, do.log=do.log))
		for (peel.col in peeled.dims) {
			test.peel =  t(apply(test.peel, 1, pca.peel, train.control.pca$rotation[,peel.col]))
		}
		test.peel.pca = prcomp(test.peel, center=F, scale.=F, do.log=do.log)
		test.proj = test.peel %*% train.peel.pca$rotation
		test.peel.pca$x = test.proj
		test.peel = as.data.frame(test.peel)
		test.peel$Response = test.data$Response
		
	}
	#train.peel.disease.pca = prcomp(train.control.peel[which2,])
	if (plot) {
		if (!is.null(test.data) & plot)
			layout(matrix(1:12, nrow=3, byrow=F))

		if (is.null(apt.color) & "Response" %in% names(training.data) & length(levels(response)) > 1) {
			training.ks = calc.ks(training.data, fdr=F)
			ks.dists = training.ks[get.aptamers(training.data),1]
			apt.color = map.color(ks.dists)
		}
		# Plot initial
		plot.pca.wrapper(train.pca, dims=dims, training.data=training.data, skip.layout=T, ...)
		#screeplot.auc(train.pca, auc.train=training.data, main=sprintf("Full Dataset",num.pcs) )
		#plot.rotation(train.pca, dims=dims, col=apt.color, ...)
		#plot.projection(train.pca, dims=dims, classes=response)

		# Plot Peel
		#screeplot.auc(train.control.pca, auc.train=training.data, auc.proj=train.control.pca.proj, main=sprintf("Peeling Dataset (Using %i PCs)",num.pcs) )
		#plot.rotation(train.control.pca, dims=dims, col=apt.color)
		train.control.pca$x = train.control.pca.proj 
		#blah = list(x=train.control.pca.proj)
		#plot.projection(train.control.pca, dims=dims, classes=response)
		#plot.projection(train.control.pca, dims=dims, col='red', add=T)
		#screeplot.auc(training.data, train.peel.disease.pca, npcs=20, main=sprintf("Disease After Peeling"))
		plot.pca.wrapper(train.control.pca, dims=dims, training.data=training.data, skip.layout=T, ...)
		
		#Plot Final
		#screeplot.auc(train.peel.pca, auc.train=training.data, main=sprintf("After Peeling"))
		#plot.rotation(train.peel.pca, dims=dims, col=apt.color)
		#plot.projection(train.peel.pca, dims=dims, classes=response)
		plot.pca.wrapper(train.peel.pca, dims=dims, training.data=training.data, skip.layout=T, ...)
		xlim = par('usr')[1:2]
		ylim = par('usr')[3:4]


		# Plot Test 
		if (!is.null(test.data)) {
			screeplot.auc(test.peel.pca, auc.train = test.data, main=sprintf("Test Dataset"))
			plot.rotation(test.peel.pca, dims=dims, col=apt.color)
			max1 = 0
			max1.i = NULL
			max2 = 0
			max2.i = NULL
			for (i in 1:20) {
				auc = pca.auc(training.data, train.peel.pca, i)
				if (auc > max1) {
					if (max1 > max2) {
						max2 = max1
						max2.i = max1.i
					}
					max1 = auc
					max1.i = i
				}
				else if (auc > max2) {
					max2 = auc
					max2.i = i
				}
			}

			max1.i = 1
			max2.i = 2

			#plot(test.proj[,max1.i], test.proj[,max2.i], main="Test Projected onto Train", col=as.numeric(test.data$Response)+1, xlab=sprintf("Component %i (Train AUC=%0.2f)", max1.i, max1), ylab=sprintf("Component %i (Train AUC=%0.2f)", max2.i, max2), xlim=xlim, ylim=ylim)
			plot.projection(test.peel.pca, dims=dims, classes=test.data$Response)
			legend("bottomright", levels(test.data$Response), col=unique(as.numeric(test.data$Response)+1))
		}

		#plot.rotation(train.peel.disease.pca, dims=dims, col=ks.color)

		
		#plot.projection(train.peel.disease.pca, dims=dims, classes=response[which2])
		#proj = train.norm %*% train.peel.disease.pca$rotation
		#print(length(proj[which1,dims[1]]))
		#points(proj[which1,dims[1]], proj[which1,dims[2]], col=2)

	}

	cluster1 = row.names(train.norm)[which(train.peel.pca$x[,1]>0)]
	cluster2 = row.names(train.norm)[which(train.peel.pca$x[,1]<0)]
	disease = row.names(training.data[which2,])
	score1 = inter.union(cluster1, disease)
	score2 = inter.union(cluster2, disease)
	invisible(list(final=train.peel.pca, peel=train.control.pca, orig=train.pca, test=test.peel.pca, aucs=aucs, train.peel=train.peel, peeled.data=train.peel.data))


}

pca.norm = function(data, center=T, scale=F, num.pcs=1, do.log=T) {
	train.peel = data[,get.aptamers(data)]
	train.peel = center.scale.data(train.peel, center=center, scale=scale, do.log=do.log)
	train.pca = prcomp(train.peel, scale.=F)
	for (i in 1:num.pcs) {
		train.peel =  t(apply(train.peel, 1, pca.peel, train.pca$rotation[,i]))
	}
	peeled.data = as.data.frame(apply.center.scale(data[,get.aptamers(data)], train.peel, undo=T, do.log=do.log, center=center, scale=scale))
	if ("Response" %in% names(data))
		peeled.data = cbind(peeled.data, Response=data$Response)
	else if (length(get.meta(data) > 1))
		peeled.data = cbind(data[,get.meta(data)], peeled.data)
	peeled.data

}
create.mean.sd.data = function(data) {
	apply(data, 2, function(x) { c(mean(x), sd(x)) })
}

center.scale.data = function(data, do.log=T, ...) {
	apts = get.aptamers(data)
	new.data = data[,apts]
	if (do.log) 
		new.data = log(new.data)
	new.data = scale(new.data, ...)
	#new.data$Response = data$Response
	new.data
}

apply.center.scale = function(orig.data, new.data, center=T, scale=T, do.log=T, undo=F) {
	orig.apts = get.aptamers.list(names(as.data.frame(orig.data)))
	new.apts = get.aptamers.list(names(as.data.frame(new.data)))
	if (is.null(new.apts) | !(length(get.aptamers(new.data)) >= 0))
		stop("Bad Names on new.data in apply.center.scale")
	if (! all(orig.apts == new.apts))
		stop("Different apt lists")

	if (!undo) {
		if (do.log) {
			orig.data[,orig.apts] = log(orig.data[,orig.apts])
			new.data[,new.apts] = log(new.data[,new.apts])
		}
		if (center & scale) {

			orig.apt.data = orig.data[,orig.apts]
			new.apt.data = new.data[,orig.apts]
			nrows = nrow(new.apt.data)
			new.apt.data = rbind(new.apt.data, apply(orig.apt.data, 2, mean), apply(orig.apt.data, 2, sd))
			new.data = as.data.frame(apply(new.apt.data, 2, function(x) (x[1:nrows] - x[nrows+1])/x[nrows+2]))
			names(new.data) = orig.apts

			#for (apt in orig.apts) 
				#new.data[,apt] = (new.data[,apt]-mean(orig.data[,apt]))/sd(orig.data[,apt])
		}
		else if (center) {
			for (apt in orig.apts) 
				new.data[,apt] = (new.data[,apt]-mean(orig.data[,apt]))
		}
		else if (scale)
			stop("Just scale not supported")

	}
	else{
		if (do.log) {
			orig.data[,orig.apts] = log(orig.data[,orig.apts])
		}
		if (center & scale) {
			orig.apt.data = orig.data[,orig.apts]
			new.apt.data = new.data[,orig.apts]
			nrows = nrow(new.apt.data)
			new.apt.data = rbind(new.apt.data, apply(orig.apt.data, 2, mean), apply(orig.apt.data, 2, sd))
			new.data = as.data.frame(apply(new.apt.data, 2, function(x) (x[1:nrows]*x[nrows+2] + x[nrows+1])))
			names(new.data) = orig.apts
		}
		else if (center & !scale){ 
			orig.apt.data = orig.data[,orig.apts]
			new.apt.data = new.data[,orig.apts]
			nrows = nrow(new.apt.data)
			new.apt.data = rbind(new.apt.data, apply(orig.apt.data, 2, mean))
			new.data = as.data.frame(apply(new.apt.data, 2, function(x) (x[1:nrows] + x[nrows+1])))
			names(new.data) = orig.apts
		}

		else if (scale)
			stop("Just scale not supported")
		if (do.log)
			new.data[,new.apts] = exp(new.data[,new.apts])
	}
	new.data
}

pca.scale.factors = function(raw.train, scaled.train) {
	common.apts = intersect(get.aptamers(raw.train), get.aptamers(scaled.train))
	med.scale.factors = rep(NULL, nrow(scaled.train))
	for (i in 1:nrow(scaled.train))  
		med.scale.factors = c(med.scale.factors, median(as.numeric(scaled.train[i,common.apts])/as.numeric(raw.train[i,common.apts])))
	med.scale.factors
}
pca.norm.old = function(raw.data, scale.factors.df=NULL, plot=T, file=NULL, ...) {
	if ("Response" %in% names(raw.data))
		stop("Must be raw data")
	# Save some variables for later use
	buffers = raw.data[raw.data$ClassName == 'Buffer',]
	
	train.data = create.training.data(raw.data[raw.data$ClassName != 'Buffer',], class1="Norm Sample", class2="Remaining Samples", ...)
	if (!is.null(file))
		png(file=file, width=10, height=10, units='in', res=300)
	pca.wrap = pca_wrapper(train.data, layout=matrix(1:12, byrow=F, nrow=3), num.pcs=1)
	scale.factors = pca.scale.factors (train.data, pca.wrap$peeled.data)
	hist(scale.factors, breaks=50)
	subarray.boxplot(train.data)
	subarray.boxplot(pca.wrap$peeled.data)
	if (!is.null(file))
		dev.off()
	pca.wrap$scale.factors = as.matrix(scale.factors)
	row.names(pca.wrap$scale.factors) = row.names(train.data)
	return(pca.wrap)
	# Do initial PCA on full dataset
	train.norm = center.scale.data(training.data[,get.aptamers(training.data)])
	train.pca = prcomp(train.norm, scale.=F, center=F)
	train.peel =  t(apply(train.norm, 1, pca.peel, train.pca$rotation[,1]))
	train.peel.pca =  prcomp(train.peel, scale.=F, center=T)
	training.data.peeled = as.data.frame(apply.center.scale(training.data, train.peel, undo=T))
	training.data.peeled = cbind(training.data[,get.meta(training.data)], training.data.peeled)



	# Plotting
	par(mfrow=c(3,3))
	dims=1:2
	training.ks = calc.ks(training.data, fdr=F)
	ks.dists = training.ks[get.aptamers(training.data),1]
	ks.dists.int = sapply(ks.dists, function(x) { max(1, floor(x*100)) })
	scale.factor = 100/(max(ks.dists.int)-min(ks.dists.int))
	ks.dists.int2 = sapply(ks.dists.int,  function(x) { max(1, floor((x-min(ks.dists.int))*scale.factor))})
	ks.color = rev(topo.colors(100))[ks.dists.int2]

	screeplot.auc(train.pca, auc.train=training.data, main=sprintf("Full Dataset",1) )
	screeplot.auc(train.peel.pca, auc.train=training.data, main=sprintf("After Peeling"))
	subarray.boxplot(training.data, main="Raw Data", do.log=T)

	plot.rotation(train.pca, dims=dims, col=ks.color)
	plot.rotation(train.peel.pca, dims=dims, col=ks.color)
	subarray.boxplot(training.data.peeled, main="PCA Normalized Data", do.log=T)


	plot.projection(train.pca, dims=dims, classes=response)
	plot.projection(train.peel.pca, dims=dims, classes=response)
	if (!is.null(scale.factors.df))
		plot(log(scale.factors.df[row.names(training.data),1], base=2), train.pca$x[row.names(training.data),1], xlab="Scale Factors (log base 2)", ylab="PCA Scale Factors", main="Scale Factor Comparison")


	list(final=train.peel.pca, orig=train.pca, peeled.data=training.data.peeled)
}

	

pca.control = function(training.data, plot=T, num.pcs=NULL) {
	if (!"Response" %in% names(training.data))
		stop("Must be training data")
	par(mfrow=c(3,3))
	response = factor(training.data$Response)
	which1 = response == levels(response)[1] # Controls
	which2 = response == levels(response)[2] # Disease
	train.norm = center.scale.data(training.data[,get.aptamers(training.data)])
	train.pca = prcomp(train.norm)
	prcomp(train.norm[which1,])
}


plot.apts = function(train.data, pca.data, dim, cutoff, flip=F) {
	if (!flip) {
		if (cutoff < 0) 
			bools = pca.data$rotation[,dim]  <= cutoff
		else
			bools = pca.data$rotation[,dim]  >= cutoff
	}
	else {
		if (cutoff > 0) 
			bools = pca.data$rotation[,dim]  <= cutoff
		else
			bools = pca.data$rotation[,dim]  >= cutoff
	}
	blah = as.data.frame(cbind(row.names(pca.data$rotation), pca.data$rotation[,dim]))
	blah.bool = blah[bools,]
	blah.bool = blah.bool[order(blah.bool[,2], decreasing=(cutoff < 0)),]
	apts = blah.bool[,1]
	plot.matrix(t(train.data[,apts]), cex=.3, force.smooth=F, do.legend=F)
}
	

score.pca.lists = function(dataset, lists, ...) {
	par(mfrow=get.row.col(length(lists)))
	spcas = lapply(names(lists), function(x) {cat(sprintf("Performing Supervised PCA for %s\n", x)); supervised.peel(dataset, row.names(lists[[x]]), do.plot=F)})
	names(spcas) = names(lists)
	sapply(names(lists), function(list.type) plot.projection(spcas[[list.type]]$unweighted, main=list.type, ...))
	invisible(spcas)
}

peel.lists = function(dataset, lists, classes=NULL, ...) {
	dataset.ks = calc.ks(dataset)[get.aptamers(dataset),1]
	layout(matrix(c(1,0,2:(length(lists)*2+1)), byrow=F, nrow=2))
	cur.dat = dataset
	first = T
	for (i in 1:length(lists)) {
		apt.bools = as.numeric(get.aptamers(dataset) %in% match.seq.ids(row.names(lists[[i]]), get.aptamers(dataset)))
		pch = c(21, 19)[apt.bools+1] # 21 if False, 20 if True
		spca = supervised.peel(cur.dat, row.names(lists[[i]]), do.plot=F)
		if (first) {
			plot.projection(spca$orig, main="Original Data", classes=classes)
			first=F
		}
		plot.rotation(spca$unweighted, scores=dataset.ks, pch=pch)
		plot.projection(spca$peeled, main=sprintf("Peeled %s", names(lists)[i]), classes=classes)
		cur.dat = spca$peeled.data
	}
	#plot.projection(spca$orig, main="Original Data", classes=classes)

	invisible(cur.dat )
}



multi.apt.pch = function(apt.list, aptamers, aptamers2=NULL, aptamers3=NULL, aptamers4=NULL, aptamers5=NULL) {
	if (!is.null(aptamers2))
		aptamers = c(aptamers, aptamers2)
	if (!is.null(aptamers3))
		aptamers = c(aptamers, aptamers3)
	if (!is.null(aptamers4))
		aptamers = c(aptamers, aptamers4)
	if (!is.null(aptamers5))
		aptamers = c(aptamers, aptamers5)
	weight.mask = as.numeric(apt.list %in% match.seq.ids(aptamers, apt.list))
	weight.mat = diag(weight.mask)
	pch = c(21, 19)[weight.mask+1] # 21 if False, 20 if True
	cex = rep(1, length(pch))
	if (!is.null(aptamers2))
		pch[apt.list %in% match.seq.ids(aptamers2, apt.list)] = 17 # 21 if False, 20 if True
	if (!is.null(aptamers3)) {
		pch[apt.list %in% match.seq.ids(aptamers3, apt.list)] = 18 # 21 if False, 20 if True
		cex[apt.list %in% match.seq.ids(aptamers3, apt.list)] = 1.5 # 21 if False, 20 if True
	}
	if (!is.null(aptamers4)) {
		pch[apt.list %in% match.seq.ids(aptamers4, apt.list)] = 15 # 21 if False, 20 if True
		cex[apt.list %in% match.seq.ids(aptamers4, apt.list)] = 1.5 # 21 if False, 20 if True
	}
	if (!is.null(aptamers5))
		pch[apt.list %in% match.seq.ids(aptamers5, apt.list)] = 8 # 21 if False, 20 if True

	list(weight.mask=weight.mask, weight.mat=weight.mat, pch=pch, aptamers=unique(aptamers), cex=cex)
}



supervised.peel = function(data, aptamers, do.plot=T, sparse=F, dims=1:2, center=T, scale=F, do.log=T, num.pcs=1, aptamers2=NULL, aptamers3=NULL, aptamers4=NULL, aptamers5=NULL, samples=NULL, ...) {
	#if (!("Response" %in% names(data)))
		#stop("Must be training data for supervised peel")
	mcall = match.call()
	#print(is.null(aptamers))
	#print(class(aptamers))
	#print(length(aptamers))
	if (is.null(aptamers) | class(aptamers) != "character" | length(aptamers) == 0) {
		print(aptamers)
		stop("Problem with aptamers")
	}
	# Create Weight mask from aptamers
	ret.list = multi.apt.pch(get.aptamers(data), aptamers, aptamers2, aptamers3, aptamers4, aptamers5)
	weight.mat = ret.list$weight.mat
	weight.mask = ret.list$weight.mask
	pch = ret.list$pch
	aptamers = ret.list$aptamers

	# Perform initial PCA
	orig.data = center.scale.data(data, center=center, scale=scale, do.log=do.log) 
	orig.pca = prcomp(orig.data, center=F, scale.=F)

	# Perform weighted PCA
	weighted.data = orig.data %*% weight.mat
	colnames(weighted.data) = get.aptamers(data)
	weighted.pca = prcomp(weighted.data, center=F, scale.=F)

	# Apply weighted basis to unweighted data
	unweighted.pca = weighted.pca
	#unweighted.pca$rotation = ginv(as.matrix(orig.data)) %*% weighted.pca$basis %*% diag(weighted.pca$single.vals)
	unweighted.pca$rotation = t(ginv(diag(weighted.pca$single.vals)) %*% t(weighted.pca$basis) %*% as.matrix(orig.data))
	unweighted.pca$x = orig.data %*% unweighted.pca$rotation

	# Peel away first dimension
	if (!do.log)print("warning implelent do.log=false")
	peeled.data = pca.peel.dim(orig.data, unweighted.pca, 1:num.pcs)
	peeled.pca = prcomp(peeled.data, scale.=scale)
	


	peeled.data.undo = as.data.frame(apply.center.scale(data, peeled.data, undo=T, do.log=T))
	if ("Response" %in% names(data))
		peeled.data.undo$Response = data$Response
	else
		peeled.data.undo = cbind(data[,get.meta(data)], peeled.data.undo)
	ret.list = list(orig=orig.pca, weighted=weighted.pca, unweighted=unweighted.pca, peeled=peeled.pca, peeled.data=peeled.data.undo)
	if (do.plot) {
		training.data = NULL
		if ("Response" %in% names(data))
			training.data = data
		plot.supervised.peel(ret.list, dims=dims, training.data=training.data, aptamers=aptamers, aptamers2=aptamers2, aptamers3=aptamers3, aptamers4=aptamers4, aptamers5=aptamers5, samples=samples, ...)
	}
	invisible(ret.list)
}
plot.supervised.peel = function(sp.data, aptamers=NULL, dims=1:2, training.data=NULL, apt.pch=21, aptamers2=NULL, aptamers3=NULL, aptamers4=NULL, aptamers5=NULL,...) {
	par(mar=c(4,4,3.5,1.5))
	layout(matrix(1:12, nrow=3, byrow=F))
	ret.list = multi.apt.pch(row.names(sp.data$orig$rotation), aptamers, aptamers2, aptamers3, aptamers4, aptamers5)
	plot.pca.wrapper(sp.data$orig, training.data=training.data, apt.pch = ret.list$pch, skip.layout=T, main="Original: Dataset", dims=dims, ...)
	plot.pca.wrapper(sp.data$weighted, training.data=training.data, apt.pch = ret.list$pch, skip.layout=T, main="Weighted Dataset", dims=dims, ...)
	plot.pca.wrapper(sp.data$unweighted, training.data=training.data, apt.pch = ret.list$pch, skip.layout=T, main="Unweighted Dataset", dims=dims, ...)
	plot.pca.wrapper(sp.data$peeled, training.data=training.data, apt.pch = ret.list$pch, skip.layout=T, main="Peeled Dataset", dims=dims, ...)
}

training.pca.wrapper = function(training.data, ...) {
	training.ks = calc.ks(training.data, fdr=F)
	ks.dists = training.ks[apts,1]
	ks.color = map.color(ks.dists)
	pca.wrapper(training.data, apt.color=ks.color)

}

pca.who.am.i = function(pca.data) {
	for (object in ls(1)) {
		if (object == 'pca.data')
			next	
		if (!is.null(row.names(get(object)))) {
			if (nrow(get(object)) == nrow(pca.data$orig$x)) {
				if(all(row.names(get(object)) == row.names(pca.data$orig$x))) {
					if(all(get.aptamers(get(object)) == row.names(pca.data$orig$rotation))) {
						return(object)
						break
					}
				}
			}
		}
	}
	print("not found")
}

prcomp.default <-
    function(x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, ...)
{
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    s <- svd(x)
	D = s$d
    s$d <- s$d / sqrt(max(1, nrow(x) - 1))
    if (!is.null(tol)) {
        ## we get rank at least one even for a 0 matrix.
        rank <- sum(s$d > (s$d[1L]*tol))
        if (rank < ncol(x)) {
            s$v <- s$v[, 1L:rank, drop = FALSE]
            s$d <- s$d[1L:rank]
        }
    }
    dimnames(s$v) <-
        list(colnames(x), paste("PC", seq_len(ncol(s$v)), sep = ""))
    r <- list(sdev = s$d, rotation = s$v, basis=s$u, single.vals=D,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc)
    if (retx) r$x <- x %*% s$v
    class(r) <- "prcomp"
    r
}
screeplot.auc = function(pca.data, auc.train, auc.proj=NULL, num.pcs = 20, scores=NULL, ...) {
	if (is.null(auc.proj))
		auc.proj = pca.data$x
	num.pcs = min(ncol(auc.proj), num.pcs)
	aucs = c()
	for (i in 1:num.pcs) {
		if (length(levels(auc.train$Response)) > 1)
			auc = colAUC(auc.proj[,i], auc.train$Response)
		else 
			auc = .5
		aucs = c(aucs, auc[1])
	}
	col.indices = sapply(aucs, function(x) {max(1, floor((x-.5)*200))})

	cols = rev(topo.colors(100))[col.indices]
	screeplot(pca.data, npcs=20, col=cols, ...)
}

pca.vote.data = function(pca.data, comp, training.data, rev=F) {
	vd = data.frame(vote=pca.data$x[,comp], class=training.data[row.names(pca.data$x), "Response"])
	if (rev)
		vd[,1] = -vd[,1]
	vd
}



pca.cv = function(training.data, num.pcs=NULL, cv.num=10, train.dims=1:20, max.disease.auc=1, plot=F, min.train.auc=0, ...) {
	all.rf.vote.data = NULL
	all.knn.vote.data = NULL
	all.first.pc.vote.data = NULL
	for (i in 1:cv.num) {
		#cat(sprintf("CV bin %i\n", i))
		test.indices = seq(i, nrow(training.data), by=cv.num)
		training.subset = training.data[-test.indices,]
		testing.subset = training.data[test.indices,]
		if (plot)
			x11()
		training.pca_wrap = pca_wrapper(training.subset, test=testing.subset, num.pcs=num.pcs, plot=plot, max.disease.auc=max.disease.auc, ...)

		final.train.dims = c()
	
		for (i in train.dims) {
			if (training.pca_wrap$aucs[i] > min.train.auc) {
				final.train.dims = c(final.train.dims, i)
			}
		}
		cat(sprintf("Training using %i PCs\n", length(final.train.dims)))

		train.data = cbind(as.data.frame(training.pca_wrap$final$x[,final.train.dims]), Response=training.subset$Response)
		names(train.data) = c(paste("PC", final.train.dims, sep=""), "Response")
		test.data = cbind(as.data.frame(training.pca_wrap$test$x[,final.train.dims]), Response=testing.subset$Response)
		names(test.data) = c(paste("PC", final.train.dims, sep=""), "Response")
		rf.model = randomForest(Response ~ ., data=train.data)
		rf.vote.data = create.vote.data(rf.model, test=test.data)

		knn.vd = kknn.vote.data(train.data, 15, test=test.data)

		first.pc.vote.data = data.frame(vote = training.pca_wrap$test$x[,final.train.dims[1]], class = testing.subset$Response)

		all.knn.vote.data = rbind(all.knn.vote.data, knn.vd)
		all.first.pc.vote.data = rbind(all.first.pc.vote.data, first.pc.vote.data)
		all.rf.vote.data = rbind(all.rf.vote.data, rf.vote.data)
	}
	true.roc(vote.data=all.rf.vote.data, pos.class=levels(training.data$Response)[2], col='blue')
	true.roc(vote.data=all.knn.vote.data, pos.class=levels(training.data$Response)[2], add=T, col='red')
	if (mean(all.first.pc.vote.data[all.first.pc.vote.data$class==levels(training.data$Response)[2],1]) < 0)
		all.first.pc.vote.data[,1] = -all.first.pc.vote.data[,1]
	true.roc(vote.data=all.first.pc.vote.data, pos.class=levels(training.data$Response)[2], add=2, col='green')

	legend('bottomright', c("RF", "KNN", "PC1"), col=c('blue', 'red', 'green'), lty=1)
}





pca.rocs = function(training.data, test.data, wrapper.output, dims=1:20) {
	model = randomForest(wrapper.output$final$x[,1:20], training.data$Response)
	true.roc(model=model, col='blue')
	true.roc(model=model, test=cbind(as.data.frame(wrapper.output$test$x[,dims]), Response=test.data$Response), col='red', add=T)
	legend("bottomright", c("Train", "Test"), col=c('blue', 'red'), lty=1)
}


pca.peel.dim = function(data.to.peel, prcomp.data, peel.dims) {
	new.peel = data.to.peel
	for (peel.dim in peel.dims) {
		new.peel =  t(apply(new.peel, 1, pca.peel, prcomp.data$rotation[,peel.dim]))
	}
	new.peel
}
apply.peel.pca = function(orig.data, new.data, num.pcs=NULL, do.plot=T) {
	response = factor(orig.data$Response)
	which1 = response == levels(response)[1] # Controls
	which2 = response == levels(response)[2] # Disease
	orig.norm = center.scale.data(orig.data[,get.aptamers(orig.data)])
	orig.pca = prcomp(orig.norm)
	orig.control.pca = prcomp(orig.norm[which1,])
	pca.var = orig.control.pca$sdev^2/sum(orig.control.pca$sdev^2)
	if (is.null(num.pcs)) {
		for (i in 1:(length(pca.var)-1)) {
			if (pca.var[i]/pca.var[i+1] > 2.5)
				break
		}
		num.pcs = i	
	}
	new.peel = apply.center.scale(orig.data, new.data)
	for (peel.col in 1:num.pcs) {
		new.peel =  t(apply(new.peel, 1, pca.peel, orig.control.pca$rotation[,peel.col]))
	}
	new.peel.pca = prcomp(new.peel, scale.=T)
	proj = new.peel %*% orig.control.pca$rotation
	if (do.plot) {
		par(mfrow=c(3,2))
		orig.ks = calc.ks(orig.data, fdr=F)
		ks.dists = orig.ks[get.aptamers(orig.data),1]
		ks.dists.int = sapply(ks.dists, function(x) { max(1, floor(x*100)) })
		scale.factor = 100/(max(ks.dists.int)-min(ks.dists.int))
		ks.dists.int2 = sapply(ks.dists.int,  function(x) { max(1, floor((x-min(ks.dists.int))*scale.factor))})
		ks.color = rev(topo.colors(100))[ks.dists.int2]

		plot(orig.pca, npcs=20, main=sprintf("Original Dataset",num.pcs) )
		plot(orig.control.pca, npcs=20, main=sprintf("Peeling Dataset (Using %i PCs)",num.pcs) )
		plot(orig.peel.disease.pca, npcs=20, main=sprintf("Disease After Peeling"))
		plot(orig.peel.pca, npcs=20, main=sprintf("After Peeling"))

		plot.rotation(train.pca, dims=dims, col=ks.color)
		plot.rotation(train.control.pca, dims=dims, col=ks.color)
		plot.rotation(train.peel.disease.pca, dims=dims, col=ks.color)
		plot.rotation(train.peel.pca, dims=dims, col=ks.color)

		plot.projection(train.pca, dims=dims, classes=response)
		plot.projection(train.control.pca, dims=dims, classes=response[which1])
		proj = train.norm %*% train.control.pca$rotation
		points(proj[which2,dims[1]], proj[which2,dims[2]], col=3)
		
		plot.projection(train.peel.disease.pca, dims=dims, classes=response[which2])
		proj = train.norm %*% train.peel.disease.pca$rotation
		points(proj[which1,dims[1]], proj[which1,dims[2]], col=2)

		plot.projection(train.peel.pca, dims=dims, classes=response)

	}
	list(proj=proj, pca=new.peel.pca)
}


plot.pca.auc = function(training.data, pca.data, dims=1:20, ...) {
	aucs = c()
	for (i in dims) 
		aucs = c(aucs, pca.auc(training.data, pca.data, i))
	plot(dims, aucs, type='l', ...)


}

pca.auc = function(training.data, pca.data, dim) {
	colAUC((pca.data$x[,dim]), training.data$Response)
}

inter.union = function(x,y) {
	length(intersect(x,y)) / length(union(x,y))
}

peel.term = function(data.row, eigen.vec) {
	#print(data.row[1:5])
	#print(eigen.vec[1:5])
	dot.p = sum(data.row * eigen.vec)
	#print(dot.p)
	norm.f = sum(eigen.vec * eigen.vec)
	#print(norm.f)
	proj.term = (dot.p/norm.f) 
	proj.term
}


pca.peel = function(data.row, eigen.vec) {
	#print(data.row[1:5])
	#print(eigen.vec[1:5])
	#print("Data")
	#print(data.row[1:5])
	#print(eigen.vec[1:5])
	dot.p = sum(data.row * eigen.vec)
	#print(dot.p[1:5])
	norm.f = sum(eigen.vec * eigen.vec)
	#print(norm.f)
	proj.term = (dot.p/norm.f) * eigen.vec
#	print((dot.p/norm.f))
#	print(proj.term[1:5])

	data.row-proj.term
}


col.string = c('green', 'red', 'blue', 'purple', 'cyan', 'orange', 'black', 'grey', '#990066', '#006600')

plot.pca.dims = function(data.prcomp, dims, value, classes, scores=NULL, main = sprintf("PCA Plot (%s)",ifelse(value=="x", "projection", value)), col=NULL, add=F, do.ident=F, ident.labels=NULL, pch = 21, xlab=NULL, ylab=NULL, do.legend=T, legend.pos="bottomright", xlim=NULL, ylim=NULL, bg=NULL, class.colors=col.string, ...) {
	if (is.null(col)) {
		col = 1
		if (!is.null(classes)) {
			classes = factor(classes)
			col = class.colors[(as.numeric(classes)-1)%%length(class.colors)+1]
		}
		if (!is.null(scores)) {
			#scores[scores == 'None'] = NA
			col = map.color(scores)
		}
	}
	if (add) {
		fun = points
	}
	else {
		fun = plot
	}
	if (is.null(xlab))
		xlab = sprintf("Component %i (%0.2f %%)", dims[1], (data.prcomp$sdev[dims[1]]^2)/sum( data.prcomp$sdev^2)*100)
	if (is.null(ylab))
		ylab=sprintf("Component %i (%0.2f %%)", dims[2], (data.prcomp$sdev[dims[2]]^2)/sum( data.prcomp$sdev^2)*100 )
	fun(data.prcomp[[value]][,dims[1]], data.prcomp[[value]][,dims[2]], col=col, main=main, xlab=xlab, ylab=ylab, pch=pch, xlim=xlim, ylim=ylim, bg=bg, ...)
	if (!is.null(classes) & !add) {
		if (length(pch) > 1) {
			pch.data = unique(cbind(pch, as.character(classes)))
			row.names(pch.data) = pch.data[,2]
			pch = as.numeric(pch.data[as.character(levels(classes)),1])
		}
		if (do.legend){
			legend(legend.pos, levels((classes)), col=class.colors[(seq(length(levels(classes)))-1)%%length(class.colors)+1],  pch=pch, pt.bg=bg, ...)
		}
	}
	if (do.ident) {
		if (value == 'x') {
			if (is.null(ident.labels))
				ident.labels = row.names(data.prcomp[[value]])
			return(row.names(data.prcomp[[value]])[identify(data.prcomp[[value]][,dims[1]], data.prcomp[[value]][,dims[2]], ident.labels)])
		}
		else if (value == 'rotation'){
			ident = identify(data.prcomp[[value]][,dims[1]], data.prcomp[[value]][,dims[2]], sapply(row.names(data.prcomp[[value]]), remove.seq_id))
			ident.apts = row.names(data.prcomp[[value]])[ident]
			return(ident.apts)
		}
	}

}


plot.rotation = function(data.prcomp, dims=1:2, classes=NULL, scores=NULL, col=NULL, training.data=NULL, aptamers=NULL, aptamers2=NULL, aptamers3=NULL, aptamers4=NULL, aptamers5=NULL, pch=NULL, cex=NULL, auto.ident=NULL, pos=4, auto.cex=.8, auto.top10=T, ...) {
	if (class(data.prcomp) == "prcomp")
		value = "rotation"
	if (class(data.prcomp) == "spca")
		value = "loadings"
	value='rotation'
	if (!is.null(aptamers)){
		ret.list = multi.apt.pch(row.names(data.prcomp$rotation), aptamers, aptamers2, aptamers3, aptamers4, aptamers5)
		pch = ret.list$pch
		cex = ret.list$cex
	}
	if (!is.null(training.data) & is.null(scores) & is.null(classes))
		scores = calc.ks(training.data, fdr=F)[get.aptamers(training.data),1]
	ret.val = plot.pca.dims(data.prcomp, dims, value, classes, scores, col=col, pch=pch, cex=cex, ...)
	if (is.null(auto.ident) & auto.top10) {
		rnames1 = row.names(data.prcomp$rotation[order(data.prcomp$rotation[,dims[1]]),])
		rnames2 = row.names(data.prcomp$rotation[order(data.prcomp$rotation[,dims[2]]),])
		pos4 = c(rnames1[1:5], rnames2[1:5], rnames2[(length(rnames2)-4):length(rnames2)])
		pos2 = c(rnames1[(length(rnames1)-4):length(rnames1)])
		pos4 = setdiff(pos4, pos2)
		auto.identify(data.prcomp$rotation[,dims[1]], data.prcomp$rotation[,dims[2]], labels=row.names(data.prcomp$rotation), auto=pos4, FUN=remove.seq_id, pos=4, cex=auto.cex)
		auto.identify(data.prcomp$rotation[,dims[1]], data.prcomp$rotation[,dims[2]], labels=row.names(data.prcomp$rotation), auto=pos2, FUN=remove.seq_id, pos=2, cex=auto.cex)
	}
	if (!is.null(auto.ident))
		auto.identify(data.prcomp$rotation[,dims[1]], data.prcomp$rotation[,dims[2]], labels=row.names(data.prcomp$rotation), auto=auto.ident, FUN=remove.seq_id, pos=pos, cex=auto.cex)
		

	ret.val 
}
plot.projection = function(data.prcomp, dims=1:2, classes=NULL, scores=NULL, col=NULL, samples=NULL, ...) {
	if (!is.null(samples))
		classes = c('Rest', "Samples")[as.numeric(row.names(data.prcomp$x) %in% samples)+1]
	plot.pca.dims(data.prcomp, dims, "x", classes, scores, col=col, ...)
}

plot.pca.wrapper = function(data.prcomp, dims=1:2, sample.bg=NULL, sample.classes=NULL, sample.scores=NULL, apt.scores=NULL, apt.pch=NULL, sample.pch=21, apt.col=NULL, sample.col=NULL, training.data=NULL, skip.layout=F, main="", aptamers=NULL, aptamers2=NULL, aptamers3=NULL, aptamers4=NULL, aptamers5=NULL, samples=NULL, sample.xlim=NULL, sample.ylim=NULL, apt.bg=NULL, apt.classes=NULL, apt.auto.ident=NULL, apt.xlim=NULL, apt.ylim=NULL, ...) {
	if (!skip.layout)
		par(mfrow=c(3,1))
	if (!is.null(training.data)) {
		apt.scores = calc.ks(training.data, fdr=F)[get.aptamers(training.data),1]
		if (is.null(sample.classes)) 
			sample.classes = training.data$Response
	}
	if (!is.null(data.prcomp$sdev))
		screeplot.auc(data.prcomp, auc.train=training.data, main=main) 
	
	apt.cex=1
	if (is.null(apt.pch)) {
		ret.list = multi.apt.pch(row.names(data.prcomp$rotation), aptamers, aptamers2, aptamers3, aptamers4, aptamers5)
		apt.pch = ret.list$pch
		apt.cex = ret.list$cex
	}

	plot.rotation(data.prcomp, dims=dims, scores=apt.scores, col=apt.col, pch=apt.pch, cex=apt.cex, bg=apt.bg, classes=apt.classes, auto.ident=apt.auto.ident, xlim=apt.xlim, ylim=apt.ylim, ...)
	plot.projection(data.prcomp, dims=dims, classes=sample.classes, bg=sample.bg, scores=sample.scores, col=sample.col, pch=sample.pch, samples=samples, xlim=sample.xlim, ylim=sample.ylim, ...)
}

get.pca.names = function(data.prcomp, type, dim, value) {
	if (type == 'r')
		type = "rotation"
	else if (type == 'p')
		type = 'x'
	else 
		stop("type must be r (rotation) or p (projection)")
	temp.dat = data.prcomp[[type]]
	temp.dat = temp.dat[order(temp.dat[,dim], decreasing=value>0),]
	if (value < 0) 
		return(row.names(temp.dat)[temp.dat[,dim] < value])
	else {
		return(row.names(temp.dat)[temp.dat[,dim] > value])
	}
}


ica.wrapper = function(data, n.comp, ...) {
	library(fastICA)
	apts = get.aptamers(data)
	data.ica = fastICA(data[,apts], n.comp=n.comp, ...)
	rot = t(data.ica$A)
	proj = data.ica$S
	comps = paste("IC", 1:n.comp)
	dimnames(rot) = list(apts, comps)
	dimnames(proj) = list(row.names(data), comps)
	list(rotation=rot, x=proj, sdev=rep(1, min(nrow(rot), nrow(proj))))
}


#pca.pairs = function(data.prcomp, dims=1:5, apt.scores, sample.classes, sample.scores, ...) {
pca.sample.pairs = function(data.prcomp, dims=1:5, training.data=NULL, sample.classes1=NULL, sample.scores1=NULL, sample.pch1=NULL, sample.classes2=NULL, sample.scores2=NULL, sample.pch2=NULL, sample.col=NULL, ...) {

	par(mfrow=c(length(dims), length(dims)), mar=c(2,2,0.5,.5))
	do.legend=F
	for (i in 1:length(dims)) {
		for (j in 1:length(dims)) {
			if (i == j) {
				exp.var = (data.prcomp$sdev[dims[i]]^2)/sum( data.prcomp$sdev^2)*100
				plot (-10, -10, axes=F, xlab="", ylab="", main="", frame.plot=T, cex=0)
				my.text(.5, .5, sprintf("PC %i\n(%0.2f %%)", dims[i],exp.var))
			}
			else if (i < j) {
				#if (j==length(dims) & i == j-1)
					#do.legend=T
				if (is.null(sample.classes1) & !is.null(training.data))
					sample.classes1=training.data$Response
				plot.projection(data.prcomp, dims=c(dims[j], dims[i]), classes=sample.classes1, scores=sample.scores1, xlab="", ylab="", main="", do.legend=do.legend, pch=sample.pch1, col=sample.col, ...)
			}
			else {
				if (is.null(sample.classes2) & !is.null(training.data))
					sample.classes2=training.data$Response
				plot.projection(data.prcomp, dims=c(dims[j], dims[i]), classes=sample.classes2, scores=sample.scores2, xlab="", ylab="", main="", do.legend=do.legend, pch=sample.pch2, ...)
			}
		}
	}
}
pca.pairs = function(data.prcomp, dims=1:5, training.data=NULL, apt.scores=NULL, sample.classes=NULL, sample.scores=NULL, apt.col=NULL, aptamers=NULL, aptamers2=NULL, aptamers3=NULL, aptamers4=NULL, aptamers5=NULL, apt.classes=NULL, sample.pch=NULL, apt.auto.ident=NULL, ...) {

	ret.list = multi.apt.pch(row.names(data.prcomp$rotation), aptamers, aptamers2, aptamers3, aptamers4, aptamers5)

	par(mfrow=c(length(dims), length(dims)), mar=c(2,2,0.5,.5))
	do.legend=F
	for (i in 1:length(dims)) {
		for (j in 1:length(dims)) {
			if (i == j) {
				exp.var = (data.prcomp$sdev[dims[i]]^2)/sum( data.prcomp$sdev^2)*100
				plot (-10, -10, axes=F, xlab="", ylab="", main="", frame.plot=T, cex=0)
				my.text(.5, .5, sprintf("PC %i\n(%0.2f %%)", dims[i],exp.var))
			}
			else if (i < j) {
				#if (j==length(dims) & i == j-1)
					#do.legend=T
				if (is.null(sample.classes) & !is.null(training.data))
					sample.classes=training.data$Response
				plot.projection(data.prcomp, dims=c(dims[j], dims[i]), classes=sample.classes, scores=sample.scores, xlab="", ylab="", main="", do.legend=do.legend, pch=sample.pch, ...)
			}
			else
				plot.rotation(data.prcomp, dims=c(dims[j], dims[i]), scores=apt.scores, xlab="", ylab="", main="", training.data=training.data, col=apt.col, pch = ret.list$pch, classes=apt.classes, cex=ret.list$cex, do.legend=F, auto.ident=apt.auto.ident, ...)
		}
	}
}

quick.proj = function(data.row, eigen.vec) {
	dot.p = sum(data.row * eigen.vec)
	#print(dot.p[1:5])
	norm.f = sum(eigen.vec * eigen.vec)
	#print(norm.f)
	proj = (dot.p/norm.f) * eigen.vec
	sign(proj[1])*sqrt(proj^2)
}
pca.get.slope.scores = function(x12, y12, data) {
	loc.points.lm = lm(y ~ x, data=data.frame(x=x12, y=y12))
	abline(loc.points.lm$coefficients, lty=2, col='red')
	if (ncol(data) != 2)
		stop("Data to pca.get.slope.bools must have 2 cols")
	apply(data, 1, function(x) quick.proj(c(x[1], x[2]), c(x[1], x[1]*loc.points.lm$coefficients[2] + loc.points.lm$coefficients[1]))[2])
}
pca.get.slope.bools = function(x12, y12, data) {
	loc.points.lm = lm(y ~ x, data=data.frame(x=x12, y=y12))
	abline(loc.points.lm$coefficients, lty=2, col='red')
	if (ncol(data) != 2)
		stop("Data to pca.get.slope.bools must have 2 cols")
	apply(data, 1, function(x) x[2] > x[1]*loc.points.lm$coefficients[2] + loc.points.lm$coefficients[1])
}

pca.get.min.slope.list = function(bools, row_names) {
	above = row_names[bools]
	below = row_names[!bools]
	if (length(above) < length(below))
		return(above)
	else
		return(below)
}

pca.get.common.slope.list = function(cur.apts, bools, names) {
	above = row.names(data.prcomp$rotation)[bools]
	below = row.names(data.prcomp$rotation)[!bools]
	if (intersect(above, cur.apts) > 0)
		return(above)
	else
		return(below)
}


pca.select.apts.area = function(data.prcomp, training.data=NULL, dims=c(1,2)) {
	plot.rotation(data.prcomp, training.data=training.data, dims=dims)
	loc.points = locator(type='l', col='red', lty=2)
	bools = pca.get.slope.bools(loc.points$x, loc.points$y, data.prcomp$rotation[,dims])
	cur.apts = pca.get.min.slope.list(data.prcomp, bools, row.names(data.prcomp$rotation))
	for (i in 2:(length(loc.points$x)-1)) {
		
	}

}

pca.select.samples = function(data.prcomp, training.data=NULL, dims=c(1,2)) {
	plot.projection(data.prcomp, training.data=training.data, dims=dims)
	loc.points = locator(n=2, type='n', col='red')
	bools = pca.get.slope.bools(loc.points$x, loc.points$y, data.prcomp$x[,dims])
	pca.get.min.slope.list(bools, row.names(data.prcomp$x))
}
get.projection.scores.apts = function(data.prcomp, training.data=NULL, dims=c(1,2), ...) {
	plot.projection(data.prcomp, training.data=training.data, dims=dims, ...)
	loc.points = locator(n=2, type='n', col='red')
	pca.get.slope.scores(loc.points$x, loc.points$y, data.prcomp$e[,dims])
}
get.rotation.scores.apts = function(data.prcomp, training.data=NULL, dims=c(1,2), ...) {
	plot.rotation(data.prcomp, training.data=training.data, dims=dims, ...)
	loc.points = locator(n=2, type='n', col='red')
	pca.get.slope.scores(loc.points$x, loc.points$y, data.prcomp$rotation[,dims])
}
pca.select.apts = function(data.prcomp, training.data=NULL, dims=c(1,2)) {
	plot.rotation(data.prcomp, training.data=training.data, dims=dims)
	loc.points = locator(n=2, type='n', col='red')
	bools = pca.get.slope.bools(loc.points$x, loc.points$y, data.prcomp$rotation[,dims])
	pca.get.min.slope.list(bools, row.names(data.prcomp$rotation))
}

opt.list = function(dataset, apt.list, reps=5, delta=5, start.size=30) {	
	start.size = min(start.size, length(apt.list))
	#dataset = center.scale.data(dataset, do.log=T, center=T, scale=F)
	apts = get.aptamers(dataset)
	apt.pool = apts
	cur.apts = match.seq.ids(apt.list[1:start.size], names(dataset))
	num.apts = length(cur.apts)
	layout(matrix(1:(reps*3), byrow=F, nrow=3))
	eigen.diff = c()
	for (i in 1:reps) {
		delta = (reps-i)/3+1
		weight.mask = as.numeric(apt.pool %in% cur.apts)
		pch = c(21, 19)[weight.mask+1] # 21 if False, 20 if True
		col = c('orange', 'blue')[weight.mask+1] # 21 if False, 20 if True
		cur.pca = supervised.peel(dataset[,c(apt.pool, "Response")], cur.apts, do.log=F, scale=T, do.plot=F)
		eigen.diff = c(eigen.diff, -diff(cur.pca$weighted$sdev[1:2]))
		screeplot.auc(cur.pca$unweighted, auc.train=dataset, main=sprintf("Round %i", i))
		plot.rotation(cur.pca$unweighted, col=col, pch=20, cex=.5)
		tab = cur.pca$unweighted$rotation[,1:2]
		load1 = cur.pca$unweighted$rotation[,1]
		load2 = cur.pca$unweighted$rotation[,2]
		pc.diff = as.matrix(abs(load1)-abs(load2))
		tab = cbind(load1, load2, pc.diff)


		tab.rem = tab[cur.apts,]
		#rem.apts = row.names(tab.rem[order(abs(tab.rem[,3]),decreasing=T),])[1:(reps-i)]
		tab.rem[,3] = apply(cur.pca$unweighted$rotation[cur.apts,2:3], 1, function (x) sum(x^2))
		#rem.apts = row.names(tab.rem[order(abs(tab.rem[,2]), decreasing=T),])[1:delta]
		rem.apts = row.names(tab.rem[order(tab.rem[,3], decreasing=T),])[1:delta]
		
		
		
		tab.add = tab[setdiff(apt.pool, cur.apts),]
		tab.add = tab.add[order(tab.add[,3], decreasing=T),]
		add.apts = row.names(tab.add)[1:delta]

		add.apts = row.names(tab.add[order(tab.add[,1], decreasing=F),])[1:length(rem.apts)]


		#plot.rotation(cur.pca$unweighted, col=col, pch=20, cex=.5)
		for (apt in add.apts)
			#text(cur.pca$unweighted$rotation[apt,1], cur.pca$unweighted$rotation[apt,2], apt, col='green')
			points(cur.pca$unweighted$rotation[apt,1], cur.pca$unweighted$rotation[apt,2], pch=4, col='green')
		for (apt in rem.apts) {
			#text(cur.pca$unweighted$rotation[apt,1], cur.pca$unweighted$rotation[apt,2], apt, col='red')
			points(cur.pca$unweighted$rotation[apt,1], cur.pca$unweighted$rotation[apt,2], pch=4, col='red')
		}
		plot.rotation(cur.pca$unweighted, col=col, pch=20, cex=.5, dims=c(2,3))
		for (apt in add.apts)
			#text(cur.pca$unweighted$rotation[apt,1], cur.pca$unweighted$rotation[apt,2], apt, col='green')
			points(cur.pca$unweighted$rotation[apt,2], cur.pca$unweighted$rotation[apt,3], pch=4, col='green')
		for (apt in rem.apts) {
			#text(cur.pca$unweighted$rotation[apt,1], cur.pca$unweighted$rotation[apt,2], apt, col='red')
			points(cur.pca$unweighted$rotation[apt,2], cur.pca$unweighted$rotation[apt,3], pch=4, col='red')
		}

		cur.apts = c(add.apts, setdiff(cur.apts, rem.apts))
		#apt.pool = setdiff(apt.pool, rem.apts)
		
		cat("Round ", i, "\n\tRemoved:\t", rem.apts, "\n\tAdded:\t", add.apts, "\n", sep=" ")

	}
	#x11()
	#plot(eigen.diff)
	eigen.diff
	cur.apts
}

pca.apt.table = function(pca.data, dims=1:5, num.apts=30) {
	apply(pca.data$rotation[,dims], 2, function(x) row.names(pca.data$rotation)[order(abs(x), decreasing=T)][1:num.apts])
}

vm.project = function(dataset, vec, norm=F, subtract.cal=F, smv.norm=F, sample.norm=F) {
	dat.apts = get.aptamers(dataset)
	matches = get.seq.ids.matches(names(vec), dat.apts)
	dataset.norm = apply.center.scale(dataset[,matches[,2]], dataset[matches[,2]], center=F, scale=F, do.log=T)
	vals = apply(dataset.norm, 1, function(sample) sum(sample * vec[matches[,1]]))
	if (smv.norm)
		vals = vals / norm(matrix(vec),"f")
	if (subtract.cal) {
		cal.proj = vm.project(dataset[dataset$ClassName %in% c("Calibrator", "cALIBRATOr"),], vec)
		vals = vals - median(cal.proj)
	}
	vals
}

plot.vm.project = function(dataset, vec, add=F, col=NULL, show.score=T, score.x=.1, score.y=.8, score.shift=1.25, adj=c(0,0), plot.fun=my.ecdf, do.plot=T, subtract.cal=T, ...) {
	if (is.null(col))
		col = as.numeric(add)+2

	vm.proj = vm.project(dataset, vec)

	score = median(vm.proj)
	if (do.plot)
		plot.fun(vm.proj, add=add, col=col, ...)
	if (any(dataset$ClassName %in% c("Calibrator", "cALIBRATOr"))) {
		if (subtract.cal) {
			cal.proj = vm.project(dataset[dataset$ClassName %in% c("Calibrator", "cALIBRATOr"),], vec)
			if (do.plot)
				plot.fun(cal.proj, lty=2, add=T, col=col, ...)
			score = score - median(cal.proj)
		}

	}


	if (show.score & do.plot) {
		if (add)
			adj[2] = adj[2] + score.shift*add
		my.text(score.x, score.y, sprintf("Score: %0.3f", score), col=col, cex=1.2, adj=adj)
	}
	invisible(return(score))
}


smv.fit = function(dataset, vec) {
	dat.apts = get.aptamers(dataset)
	dataset = log(dataset[,dat.apts])
	matches = get.seq.ids.matches(names(vec), dat.apts)
	orig.ratios = dataset[,matches[,2]] / vec[matches[,1]]
	coefs = trim(seq(min(orig.ratios)*1.2, max(orig.ratios)*.8, length.out=1000), type='right', percent=.2)
	print(range(coefs))
	coefs = apply(dataset[,matches[,2]], 1, function(sample) {
		mses = sapply(coefs,  function(coef) {
			mean(abs((sample / coef - vec[matches[,1]])^1))
		})	
#	plot(as.numeric(dataset[matches[,2]]),  vec[matches[,1]])
		min.coef = coefs[which( mses == min(mses))]
		min.coef
	})
	return(coefs)
	#par(mfrow=c(2,1))
	print(class(vec))
	plot(vec[matches[,1]], ylim=range(c(dataset[,matches[,2]] / min.coef, vec[matches[,1]])))
	print(dim(dataset[matches[,2]]))
	points(1:nrow(matches), dataset[,matches[,2]] / min.coef, col='red')

	#barplot(as.numeric(vec[matches[,1]] - dataset[,matches[,2]] / min.coef))
	dataset[matches[,2]] / min.coef

}



robust.pca.shrinkage = function(x, tau) {
	x.sign = sign(x)
	x.shrink = abs(x) - tau
	x.shrink[x.shrink < 0] = 0
	x.sign * x.shrink
}

robust.pca.svdThresh = function(x, tau) {
	x.svd = svd(x)
	x.svd$u %*% robust.pca.shrinkage(diag(x.svd$d), tau) %*% t(x.svd$v)
}

robust.pca = function(x, tolerance=1e-7, verbose=T, ...) {
	library(Matrix)
	x = as.matrix(x)
	S.cur = matrix(0, nrow=nrow(x), ncol=ncol(x))
	Y.cur = matrix(0, nrow=nrow(x), ncol=ncol(x))
	mu = nrow(x)*ncol(x)/4/norm(x, '1')

	lambda = nrow(x)^-.65

	cur.err = 100
	rounds = 1
	while (cur.err > tolerance*norm(x, 'F')) {
		L.next = robust.pca.svdThresh(x - S.cur + mu^-1*Y.cur, mu^-1)
		S.next = robust.pca.shrinkage(x - L.next + mu^-1*Y.cur, mu^-1*lambda)
		Y.next = Y.cur + mu*(x - L.next - S.next)

		L.cur = L.next
		S.cur = S.next
		Y.cur = Y.next
		cur.err = norm(x - L.cur - S.cur, 'F')
		print(sprintf("Round %i: %0.5f", rounds, cur.err))
		if (rounds > 100)
			break
		rounds = rounds + 1
	}
	dimnames(L.cur) = dimnames(S.cur)
	return(list(L=L.cur, S=S.cur))

}



