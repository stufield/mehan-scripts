#library(randomForest)
#library(e1071)
#library(kknn)
#library(pamr)

filter.rf.imp = function(rf.imp, count=nrow(rf.imp)) {
	as.matrix(rf.imp[order(rf.imp[,1], decreasing=T),][1:count])
}


cv.snippet = function(training, num.cv=10) {
	mses = lapply(1:10, function(i) { # for loop that returns a list(dictionary)
		withheld.set = seq(i, nrow(training), by=num.cv)
		lasso.fit = lasso(Response ~ ., training[-withheld.set,])
		lasso.pred = predict(lasso.fit, training[withheld.set,])
		# Some mse or performance calculation here, that gets returned 
	})
	unlist(mses) # returns a vector of mses. unlist concatenates the dictionary values
}



compare.nb = function(dataset_list, features, num.range=NULL, cv=T, labels=NULL, bw=F, lwd=2, ...) {
	if (is.null(num.range)) 
		num.range = 1:length(dataset_list)
	leg = c()
	aucs = c()
	if (bw) {
		col = rep(c("black", "black"), length(num.range))
		#lty = rep(c(1,2,3), length(num.range))
		lty = c(1, 3, 4, 5, 6, 1, 3, 1, 3, 1, 3)
		pch = c(rep(NA_integer_, 5), 1, 1, 2, 2, 3, 3)
		#pch = rep(c(rep(0,3), rep(3, 3), rep(2, 3), rep(4,3)), length(num.range))
	}
	else {
		col = seq(length(num.range))
		lty = as.integer(((1:length(num.range))-1)/8)+1
		pch = rep(NA_integer_, 11)
		#pch= rep(" ", length(num.range))
	}
	for (j in seq(length(num.range))) {
		i = num.range[j]
		#than this one true.roc(vote.data = class.cv(dataset_list[[i]][,features], dataset_list[[i]]$Response, model.type="nb", 10), pos.class=levels(dataset_list[[i]]$Response)[2], add=j-1, col=j, ...)
		if (cv)
			auc = true.roc(vote.data = class.cv(dataset_list[[i]][,features], dataset_list[[i]]$Response, model.type="nb", 10), pos.class=levels(dataset_list[[i]]$Response)[2], add=j-1, lty=lty[j], col=col[j], pch=pch[j], lwd=lwd, ...)
		else
			auc = true.roc(model = fitted.naiveBayes(dataset_list[[i]][,features], dataset_list[[i]]$Response, mad=T), test=dataset_list[[i]][,c(features, "Response")], add=j-1, col=col[j], lty=lty[j], pch=pch[j], lwd=lwd, ...)
			
		aucs = c(aucs, auc)
		leg = c(leg, names(dataset_list)[i])
	}

	if (!is.null(labels))
		leg = labels
	legend("bottomright", leg, col=col, lty=lty, pch=pch, lwd=lwd)
	aucs
}

compare.rf.nb = function(training, features, ...) {
	rf.model = randomForest(training[,features], training$Response)
	nb.cv.votes = class.cv(training[,features], training$Response, model.type='nb', 10)
	nb.model = fitted.naiveBayes(training[,features], training$Response)

	true.roc(model=rf.model, col='blue', ...)
	true.roc(vote.data=nb.cv.votes, add=T, col='red', pos.class=levels(training$Response)[2], ...)
	true.roc(model=nb.model, add=2, col='green', pos.class=levels(training$Response)[2], test=training[,c(features, "Response")], ...)
	legend("bottomright", c("Random Forest", "Naive Bayes (CV)", "Naive Bayes"), col=c("blue", "red", "green"), lty=1)
}

compare.rf.nb.scores = function(rf.model, nb.model, test, do.ident=F, min.prob=NULL, ...) {
	plot(create.vote.data(rf.model, test=test)$vote, create.vote.data(nb.model, test=test, min.prob=min.prob)$vote, ...)
	abline(h=.5, v=.5, lty=2, col=8)
	if (do.ident)
		identify(create.vote.data(rf.model)$vote, create.vote.data(nb.model, test=test)$vote, labels=row.names(test))
}

plot.nb.weird.single = function(model, data, sample, apt) {
	plot(-1000, xlim=log(range(data[,apt])), ylim=range(density(log(data[,apt]))$y), xlab="RFU", ylab="f(x)")
	points(density(log(data[data$Response==data$Response[1],apt])), col='blue', type='l')
	points(density(log(data[data$Response!=data$Response[1],apt])), col='red', type='l')
	x = log(seq(min(data[,apt]), max(data[,apt]), length.out=500))
	points(x, dnorm(x, mean=model$tables[[apt]][2,1], sd=model$tables[[apt]][2,2]), col='red', type='l', lty=2)
	points(x, dnorm(x, mean=model$tables[[apt]][1,1], sd=model$tables[[apt]][1,2]), col='blue', type='l', lty=2)
	# dnorms
	abline(v=log(data[sample,apt]))

}
check.nb.weird = function(model, data, apts=NULL, max.lr=10000) {
	data = log(data[,get.aptamers(data)])
	if (is.null(apts)) 
		apts = names(model$tables)
	sapply(row.names(data), function(sample) 
		sapply(apts, function(apt) {
			lr = max(min(max.lr, log( dnorm(data[sample,apt], mean=model$tables[[apt]][2,1], sd=model$tables[[apt]][2,2]) / dnorm(data[sample,apt], mean=model$tables[[apt]][1,1], sd=model$tables[[apt]][1,2]), base=2)), -max.lr)
			if (model$tables[[apt]][2,1] > model$tables[[apt]][1,1]) {# Disease Up
				if (data[sample,apt] > model$tables[[apt]][2,1])  # Sample Up
					return(lr < 0)
				if (data[sample,apt] < model$tables[[apt]][1,1])  # Sample Down
					return(lr > 0)
			}
			else { # Disease Down
				if (data[sample,apt] > model$tables[[apt]][1,1]) # Sample Up
					return(lr > 0)
				if (data[sample,apt] < model$tables[[apt]][2,1]) # Sample Down
					return(lr < 0)
			}
			return(F)
		})
	)
}


backwards.nb = function(data, cur.list=NULL, cut.prop=.10, ...) {
	if (is.null(cur.list))
		cur.list = names(data)[-dim(data)[2]]
	out.list = list()
	iter = 1
	while(length(cur.list) > 0) {
		cat(noquote(paste("Iteration ", iter, ": Building tree with ", length(cur.list), " features\n", sep="")))
		data.nb = naiveBayes(as.data.frame(data[,cur.list]), data[,"Response"], ...)
		out.list[[iter]] = data.rf
		num.to.keep = as.integer((1-cut.prop)*length(cur.list))
		o = order(data.rf$importance, decreasing=T)
		if (length(cur.list) == 1)  # Can't take range from 1 to 0, so need to break here
			break
		cur.list = names(data.rf$importance[o,])[1:num.to.keep]
		iter = iter+1
	}
	out.list
}

	
backwards.rf.pan  = function(data.list, cur.list=NULL, cut.prop=.10, min.size=1, do.nb=T, weights=rep(1, length(data.list)), ...) {
	if (is.null(cur.list))
		cur.list = names(data)[-dim(data)[2]]
	rf.list = list()
	nb.list = list()
	for (disease in names(data.list)) {
		rf.list[[disease]] = list()
		nb.list[[disease]] = list()
	}
	iter = 1
	while(length(cur.list) >= min.size) {
		cat(noquote(paste("Iteration ", iter, ": Building classifiers with ", length(cur.list), " features\n", sep="")))
		for (disease in names(data.list)) {
			cat(noquote(paste("\tRandom Forest on ", disease, "\n", sep="")))
			rf.list[[disease]][[length(cur.list)]] = randomForest(as.data.frame(data.list[[disease]][,cur.list]), data.list[[disease]][,"Response"], keep.forest=F, ...)
			if (do.nb) {
				cat(noquote(paste("\tNaive Bayes on ", disease, "\n", sep="")))
				nb.list[[disease]][[length(cur.list)]] = fitted.naiveBayes(as.data.frame(data.list[[disease]][,cur.list]), data.list[[disease]][,"Response"], mad=T)
			}
		}
		importance = matrix(0, nrow=length(cur.list), ncol=1)
		for (i in 1:length(weights)) {
			print(i)
			#print(rf.list[[i]])
			importance = as.matrix(as.matrix(importance) + as.matrix(rf.list[[i]][[length(cur.list)]]$importance*weights[i]))
		}


		num.to.keep = as.integer((1-cut.prop)*length(cur.list))
		o = order(importance, decreasing=T)
		if (length(cur.list) == 1)  # Can't take range from 1 to 0, so need to break here
			break
		
	
		cur.list = names(importance[o,])[1:num.to.keep]
		
		iter = iter+1
	}
	if(do.nb)
		return(list(rf.list, nb.list))
	else
		return(rf.list)
}

opt.rf = function(data, weights=c(.25,.5,1,2,4), mtrys=NULL, nodesizes=NULL, cur.list=NULL) {

	if (is.null(cur.list))
		cur.list = get.aptamers(data)
	if (is.null(mtrys))
		mtrys = unique(floor(seq(floor(sqrt(length(cur.list))/2), floor(sqrt(length(cur.list))*2), length.out=10)))
	if (is.null(nodesizes))
		nodesizes = unique(floor(seq(1, floor(nrow(data)/3), length.out=10)))
	auc.table = NULL
	par(mfrow=c(length(nodesizes), length(mtrys)), mar=c(3, 2, 2, 1))
	print(mtrys)
	print(nodesizes)
	print(weights)
	#mat = matrix(nrow=length(nodesizes), 
	for (nodesize in nodesizes) {
		for (mtry in mtrys) {
			count = 0
			for (weight in weights) {
				model=randomForest(data[,cur.list], data$Response, weight=c(1-weight, weight), nodesize=nodesize, mtry=mtry, ntree=1000)
				#auc = true.roc(model=model, add=count, main=sprintf("Nodesize: %i M-try: %i", nodesize, mtry), cex=.5, lwd=1, lty=1, auc.x=.4, auc.y=.4, )
				vote.data = create.vote.data(model)
				auc = colAUC(vote.data$vote, vote.data$class)
				ss = get.sens.spec(model, cutoff=.5)
				count = count + 1
				auc.table = rbind(auc.table, c(nodesize, mtry, weight, ss, auc))

			}
			#legend("bottomright", as.character(weights), col=1:length(weights), lwd=1)
		}
	}
	auc.table = as.data.frame(auc.table)
	names(auc.table) = c("nodesize", "mtry", "weight", "sens", "spec", "auc")
	auc.table
}

greedy.forward = function(data, cur.list, min.auc, max.size=20, max.classifiers=100, verbose=F, min.prob=NULL, ...) {
	cur.list = sort(cur.list)
	library(caTools)
	model = fitted.naiveBayes(data[,cur.list], data$Response, mad=T)
	auc.list = list()
	features.list = list()
	features.list[[2]] = list()
	cat("Building two feature classifiers\n")
	print(length(cur.list))
	features.table = data.frame(auc=NULL, index=NULL)
	for (i in 1:(length(cur.list)-1)) {
		for (j in (i+1):length(cur.list)) {
			vote.data = create.vote.data(model=model, test=data[,c(cur.list[i], cur.list[j], "Response")], min.prob=min.prob)
			roc.auc = colAUC(vote.data$vote, vote.data$class)
			if (roc.auc > min.auc) {
				features.list[[2]][[length(features.list[[2]])+1]] = list(c(cur.list[i], cur.list[j]), roc.auc)
				cat("Added", c(cur.list[i], cur.list[j], roc.auc), "\n")
				features.table = rbind(features.table, c(roc.auc, length(features.list[[2]])))
				#features.list[[2]][[length(features.list)+1]] = roc.auc
			}
		}
	}
	new.features.list = list()
	features.table = features.table[order(features.table[,1], decreasing=T),]
	for (i in 1:min(nrow(features.table), max.classifiers) )
		new.features.list[[i]] = features.list[[2]][[features.table[i,2]]]
	features.list[[2]] = new.features.list



	for (size in 3:max.size) {
		features.table = data.frame(auc=NULL, index=NULL)
		cat("Building ", size, " feature classifiers\n")
		features.list[[size]] = list()
		count = 1
		for (feature.set.roc in features.list[[size-1]]) {
			feature.set = feature.set.roc[[1]]
			cat(count, "/", length(features.list[[size-1]]), " ", paste(sapply(feature.set, remove.seq_id), compress=" "), "\n", sep="")
			remaining.apts = setdiff(cur.list, feature.set)
			for (feature in remaining.apts) {
				vote.data = create.vote.data(model=model, test=data[,c(feature.set, feature, "Response")], min.prob=min.prob)
				roc.auc = colAUC(vote.data$vote, vote.data$class)
				if (roc.auc > min.auc) {
					features.list[[size]][[length(features.list[[size]])+1]] = list(sort(c(feature.set, feature)), roc.auc)
					features.table = rbind(features.table, c(roc.auc, length(features.list[[size]])))
					cat("Added", c(sapply(c(feature), remove.seq_id), roc.auc), "\n")
				}
			}
			count = count + 1
		}
		new.features.list = list()
		if (nrow(features.table) < 1)
			break
		features.table = features.table[order(features.table[,1], decreasing=T),]
		features.table = features.table[!duplicated(features.table[,1]),]
		print(nrow(features.table))
		for (i in 1:min(nrow(features.table), max.classifiers) )
			new.features.list[[i]] = features.list[[size]][[features.table[i,2]]]
		features.list[[size]] = new.features.list
	}

	features.list

}





backwards.rf = function(data, cur.list=NULL, cut.prop=.10, feature.only=F, do.nb=F, min.size=1, verbose=F, keep.forest=F, ...) {
	if (is.null(cur.list))
		cur.list = names(data)[-dim(data)[2]]
	out.list = list()
	nb.list = list()
	iter = 1
	while(length(cur.list) >= min.size) {
		if (verbose)
			cat(noquote(paste("Iteration ", iter, ": Building tree with ", length(cur.list), " features\n", sep="")))
		data.rf = randomForest(as.data.frame(data[,cur.list]), data[,"Response"], keep.forest=keep.forest, ...)

		num.to.keep = as.integer((1-cut.prop)*length(cur.list))
		o = order(data.rf$importance, decreasing=T)
		
		if (feature.only)
			out.list[[length(cur.list)]] = names(data.rf$importance[o,])[num.to.keep+1:length(names(data.rf$importance[o,]))]
		else 
			if (length(cur.list) < 30)
				out.list[[length(cur.list)]] = data.rf

		if (do.nb) {
			if (verbose)
				cat(noquote(paste("\t Also performing Naive Bayes with ", length(cur.list), " features\n", sep="")))
			nb.list[[length(cur.list)]] = fitted.naiveBayes(as.data.frame(data[,cur.list]), data[,"Response"], mad=T)
		}
		if (length(cur.list) == 1)  # Can't take range from 1 to 0, so need to break here
			row.names(out.list[[1]]$importance) = cur.list
			#break
	
		cur.list = names(data.rf$importance[o,])[1:num.to.keep]
		
		iter = iter+1
	}
	if(do.nb)
		return(list(out.list, nb.list))
	else
		return(out.list)
}

list.recall = function(list.item) {
	confusion_stats(list.item$confusion)$stats[1]
}
list.specificity = function(list.item) {
	confusion_stats(list.item$confusion)$stats[2]
}
list.precision = function(list.item) {
	confusion_stats(list.item$confusion)$stats[3]
}
list.err.rate = function(list.item) {
	list.item$err.rate[length(list.item$err.rate)]
}
list.auc = function(x, test=NULL) { 
	vd =  create.vote.data(x, test=test)
	colAUC(vd$vote, vd$class)
}

naiveBayes.table = function(model) {
	init.table = as.data.frame(t(as.data.frame(model$tables)))
	means = init.table[seq(1,nrow(init.table), by=2),]
	sds = init.table[seq(2,nrow(init.table), by=2),]
	final.table = as.data.frame(cbind(means,sds))
	names(final.table) = c(sprintf("%s Mean", names(init.table)[1]), sprintf("%s Mean", names(init.table)[2]), sprintf("%s SD", names(init.table)[1]), sprintf("%s SD", names(init.table)[2]))
	final.table
}



plot.backwards.rocs = function(rf.list, sizes=2:17, rf.test=NULL, nb.list=NULL, nb.test=NULL, col='blue', mfrow=NULL, no.main=FALSE, ...) {

	if (length(rf.list) == 0) {
		plot(0,0, frame.plot=F, cex=0, axes=F, ylab="", xlab="")
		text(0,0,"No Significant Markers")
	}
	if (length(rf.list) == 2) 
		print("Warning: List has length 2. Did you forget that you included Naive Bayes?")

	else {
		rf.max.auc = 0
		rf.max.size = 0
		nb.max.auc = 0
		nb.max.size = 0
		if (is.null(mfrow))
			mfrow = get.row.col(length(sizes))
		par(mfrow=mfrow)
		for (num in sizes) {
			if (num > length(rf.list))
				break
			if (is.null(rf.list[[num]]))
				next
			if (class(rf.list[[num]]) == "naiveBayes") {
				classes = rf.list[[num]]$levels
				features = names(rf.list[[num]]$tables)
			}
			else if (class(rf.list[[num]]) == "randomForest") {
				classes = rf.list[[num]]$classes
				features = row.names(rf.list[[num]]$importance)
			}
			if (!no.main)
				main = sprintf("%i Markers", num)
			else
				main = ""
				
			rf.auc = true.roc(model=rf.list[[num]], main=main, test=rf.test[,c(features, "Response")], col=col, ...)
			if (!is.null(nb.list)) {
				nb.auc = true.roc(model=nb.list[[num]], test=nb.test[,c(features, "Response")], col='red', add=T, min.prob=.1, ...)
				if (nb.auc > nb.max.auc) {
					nb.max.auc = nb.auc
					nb.max.size = num
				}
			}

			if (rf.auc > rf.max.auc) {
				rf.max.auc = rf.auc
				rf.max.size = num
			}

		}
		#legend("bottomright", c("Random Forest", "Naive Bayes"), col=c("blue", "red"), lty=1)

		if (!is.null(nb.list)) 
			return(rf.max.size) # Removed NB maximum here
		else
			return(rf.max.size)
	}
}

backwards.rf.one.out = function(data, num.reps, cur.list=NULL, model.size=15, ...) {
	if (is.null(cur.list))
		cur.list = names(data)[-dim(data)[2]]
	ret.list = list()
	for (i in 1:num.reps) {
		if (length(cur.list) < 1)
			next
		back.rf = backwards.rf(data=data, cur.list=cur.list, do.nb=F, cut.prop=1/(model.size+2), ...)

		if (length(cur.list) < model.size) 
			model.size = length(cur.list)
		data.rf = back.rf[[model.size]]
		
		o = order(as.numeric(data.rf$importance[,1]), decreasing=T)
		print(sprintf("Removing %s", row.names(as.matrix(data.rf$importance[o,]))[1]))
		cur.list = cur.list[-grep(row.names(as.matrix(data.rf$importance[o,]))[1],  cur.list)]
		ret.list[[i]] = data.rf
	}
	ret.list
}

plot.backwards.rf.one.out = function(one.out.list) {
	aucs = c()
	for (i in 1:length(one.out.list)) 
		aucs = c(aucs, list.auc(one.out.list[[i]]))

	plot(aucs, type='l', col='blue')
}
get.backwards.rf.one.out = function(one.out.list, stop) {
	out.list = c()
	for (i in 1:stop)
		out.list = union(out.list, row.names(one.out.list[[i]]$importance))
	out.list
}
		
rf.table = function(model, remove.seq.ids=T, apt.data=NULL) {
	tab = as.data.frame(sorted.importance(model))
	names(tab) = "Gini Importance"
	if (!is.null(apt.data))
		tab = cbind(apt.data[rn(tab),c("GeneName", "Target", "GeneId", "SwissProt")], tab)
	if (remove.seq.ids)
		rownames(tab) = sapply(rownames(tab), remove.seq_id)
	tab
}
nb.table = function(model, remove.seq.ids=T) {
	control.mean.sd = t(as.data.frame(lapply(model$tables, function(x) x[1,])))
	disease.mean.sd = t(as.data.frame(lapply(model$tables, function(x) x[2,])))
	tab = as.data.frame(cbind(control.mean.sd, disease.mean.sd))
	names(tab) = c(paste(row.names(model$tables[[1]])[1], c("Mean", "SD")), paste(row.names(model$tables[[1]])[2], c("Mean", "SD")))
	if (remove.seq.ids)
		rownames(tab) = sapply(rownames(tab), remove.seq_id)
	tab
}


plot.backwards.stats = function(rf.list, func, label, add=F, col='blue', ylim=NULL, lty=lty, xlim=xlim, ...) {
	recall = c()
	lengths = c() 
	for (rf in rf.list) { 
		if (is.null(rf))
			next
		if (class(rf) == "randomForest")
			lengths = c(lengths, length(rf$importance))
		else {
			size = length(rf$tables)
			if (size == 1)
				next
			lengths = c(lengths, size)
		}
		recall = c(recall, func(rf, ...))

	}
	if (add) {
		points(lengths, recall, type='l', col=col, lty=lty)
	}
	else
		plot(lengths, recall, type='l', main=paste("RF", label, "vs. Num Biomarkers"), xlab="Number of Biomarkers in RF", ylab=paste("RF", label), col=col, ylim=ylim, lty=lty, xlim=xlim)
}

get.list.vals = function(rf.list, func) {
	vals = c()
	for (i in 1:length(rf.list)) { 
		vals = c(vals, func(rf.list[[i]]))
	}
	vals
}


plot.backwards.stats.all.both = function(both.list, ...) {
	rf.list = both.list[[1]]
	nb.list = both.list[[2]]
	par(mfrow = c(2,2))
	
	# Recall
	all.vals = c(get.list.vals(rf.list, list.recall), get.list.vals(nb.list, list.recall))
	plot.backwards.stats(rf.list, list.recall, "Recall", ylim=c(min(all.vals, na.rm=T), max(all.vals, na.rm=T)), ...)
	plot.backwards.stats(nb.list, list.recall, "Recall", add=T);
	abline(v=2, lty=2)
	legend("bottomright", c("RF", "NB"), col=c("blue", "red"), lty=1, inset=c(0,.0365))
	
	# Specificity 
	all.vals = c(get.list.vals(rf.list, list.specificity), get.list.vals(nb.list, list.specificity))
	plot.backwards.stats(rf.list, list.specificity, "Specificity", ylim=c(min(all.vals, na.rm=T), max(all.vals, na.rm=T)), ...)
	plot.backwards.stats(nb.list, list.specificity, "Specificity", add=T);
	abline(v=2, lty=2)
	legend("bottomright", c("RF", "NB"), col=c("blue", "red"), lty=1, inset=c(0,.0365))
	
	# Precision
	all.vals = c(get.list.vals(rf.list, list.precision), get.list.vals(nb.list, list.precision))
	plot.backwards.stats(rf.list, list.precision, "Precision", ylim=c(min(all.vals, na.rm=T), max(all.vals, na.rm=T)), ...)
	plot.backwards.stats(nb.list, list.precision, "Precision", add=T);
	abline(v=2, lty=2)
	legend("bottomright", c("RF", "NB"), col=c("blue", "red"), lty=1, inset=c(0,.0365))
	
	# RF Err
	plot.backwards.stats(rf.list, list.err.rate, "Error Rate");
}
plot.backwards.stats.all = function(rf.list, ...) {
	if (length(rf.list) == 2)
		plot.backwards.stats.all.both(rf.list)
	else {
		par(mfrow = c(2,2))
		plot.backwards.stats(rf.list, list.recall, "Recall", ...);
		plot.backwards.stats(rf.list, list.specificity, "Specificity", ...);
		plot.backwards.stats(rf.list, list.precision, "Precision", ...);
		plot.backwards.stats(rf.list, list.err.rate, "Error Rate", ...);
	}
}

my.tune.knn = function(training, response, k.range, num.cv=10, add=F, ...) {
	aucs = data.frame(k=NULL, auc=NULL)
	data = cbind(training, response)
	names(data) = c(names(training), "Response")
	for (k in k.range) {
		all.vote.data = my.kknn.cv(training, response, k, num.cv)
		auc = true.roc(model=NULL, test=NULL, vote.data=all.vote.data, col=k, auc=F, pos.class="stage1-3", add=add)
		aucs = rbind(aucs, data.frame(k=k, auc=auc))
		add=T
	}
	legend("bottomright", as.character(k.range), col=k.range, lwd=1)
	aucs
}


kknn.vote.data.n1 = function(data, k=NULL, kernel='gaussian') {
	vote.data = data.frame(vote=NULL, class=NULL)
	for (i in 1:nrow(data)) 
		vote.data = rbind(vote.data, kknn.vote.data(knn.data[-i,], k, knn.data[i,], kernel=kernel))
	vote.data
}


my.kknn.cv = function(training, response, k, num.cv=10, ...) {
	data = cbind(training, response)
	names(data) = c(names(training), "Response")
	all.vote.data = data.frame(vote=NULL, class=NULL)
	for (cv in 1:num.cv) {
		cv.group = seq(nrow(data))[seq(nrow(data))%%num.cv+1!=cv]
		training.knn = kknn(Response~., data[cv.group,], data[-cv.group,], k=k, ...)
		pred=data.frame(pred.class=training.knn$fitted.values, prob=training.knn$prob[,2])
		classes = levels(data$Response)
		vote.data = data.frame(vote=pred$prob, class=data$Response[-cv.group])
		all.vote.data = rbind(all.vote.data, vote.data)
	}
	all.vote.data
}

kknn.vote.data = function(data, k, test=NULL, ...) {
	if (is.null(test))
		test = data
	training.knn = kknn(Response~., data, test, k=k, ...)
	vote.data = data.frame(vote=training.knn$prob[,2], class=test$Response)
	vote.data
}



boost.cv = function(training.data, num.cv=10, apt.subset=NULL, ...) {
	if (!("Response" %in% names(training.data))) 
		stop("Missing Response in training.data")
	if (!is.null(apt.subset))
		training.data = training.data[,c(apt.subset, "Response")]
	rows = nrow(training.data)
	groups = ((1:rows) %% 10) + 1
	all.vote.data = NULL
	for (i in 1:num.cv) {
		train.group = which(groups!=i)
		test.group = which(groups==i)
		training.model = gbm.wrapper(training.data[train.group,], ...)
		vote.data = create.vote.data(training.model, training.data[test.group,])
		all.vote.data = rbind(all.vote.data, vote.data)
	}
	all.vote.data

}

	

class.cv = function(training, response, model.type, num.cv, min.prob=NULL, ...) {
	data = cbind(training, Response=response)
	names(data) = c(names(training), "Response")
	conf.table = as.table(matrix(0, ncol=2, nrow=2, dimnames=list(levels(response), levels(response))))
	if (ncol(training) == 1)
		return(list(confusion=conf.table, importance=0))
	all.vote.data = data.frame(vote=NULL, class=NULL)
	for (i in 1:num.cv) {
		cv.group = seq(length(response))[seq(length(response))%%num.cv+1!=i]
		if (model.type=="svm")
			training.model = svm(Response~., data=data, subset=cv.group,  probability=T)
		else if (model.type=="nb") {
			training.model = fitted.naiveBayes(training[cv.group,], response[cv.group], mad=T)
			#training.model = fitted.naiveBayes(training[cv.group,], response[cv.group])
		}
		else if (model.type=="rf")
			training.model = randomForest(Response~., data=data, subset=cv.group)
		else if (model.type=="kknn")
			training.model = kknn(Response~., data[cv.group,], data[-cv.group,], ...)
		else
			stop("Invalid model type, must be svm, rb, or nb")
		vote.data = create.vote.data(training.model, data[-cv.group,], min.prob=min.prob)
		all.vote.data = rbind(all.vote.data, vote.data)
	}
	all.vote.data[order(all.vote.data[,1], decreasing=T),]
}


nb.cv.old = function(x, y, num.cv) {
	conf.table = as.table(matrix(0, ncol=2, nrow=2, dimnames=list(c('ben_nod', 'stage1-3'), c('ben_nod', 'stage1-3'))))
	if (dim(x)[2] == 1)
		return(list(confusion=conf.table, importance=0))
	for (i in 1:num.cv) {
		cv.group = seq(length(y))[seq(length(y))%%num.cv+1!=i]
		x.nb = naiveBayes(x, y, subset=cv.group)
		x.nb.pred = predict(x.nb, x[-cv.group,])
		conf.table = conf.table + table(y[-cv.group], x.nb.pred)
	}
	print(dim(x))
	list(confusion=conf.table, importance=x[1,]) # This is a hack to create a vector of correct length for plotting
}


fitted.naiveBayes <- function(x, y, do.log=T, trim.perc=0, mad=F) {
	if (do.log)
		x = log(x)

    call <- match.call()
    Yname <- deparse(substitute(y))
    x <- as.data.frame(x)

    ## estimation-function
    est <- function(var)
        if (is.numeric(var)) {
            cbind(tapply(var, y, coefs.wrapper, trim.perc=trim.perc, mad=mad, coef=1),
                  tapply(var, y, coefs.wrapper, trim.perc=trim.perc, mad=mad, coef=2))
        } else {
            tab <- table(y, var)
            (tab + laplace) / (rowSums(tab) + laplace * nlevels(var))
        }

    ## create tables
    apriori <- table(y)
	apriori[1] = length(y)/2
	apriori[2] = length(y)/2
    tables <- lapply(x, est)

    ## fix dimname names
    for (i in 1:length(tables))
        names(dimnames(tables[[i]])) <- c(Yname, colnames(x)[i])
    names(dimnames(apriori)) <- Yname

    structure(list(apriori = apriori,
                   tables = tables,
                   levels = levels(y),
                   call   = call
                   ),

              class = "naiveBayes"
              )
}
predict.naiveBayes.bad <- function(object,
                               newdata,
                               type = c("class", "raw"),
                               threshold = 0.000001,
							   min.prob = NULL,
                               ...) {
	print('yup')
    type <- match.arg(type)
    newdata <- as.data.frame(newdata)
	if (length(object$tables) ==1 )
		stop("Naive Bayes must have more than 1 feature")
	#print(object$tables)
	#print(paste("test", names(newdata)))
    attribs <- which(names(object$tables) %in% names(newdata))
	print(attribs)
    isnumeric <- sapply(newdata, is.numeric)
    newdata <- data.matrix(newdata)
    L <- sapply(1:nrow(newdata), function(i) {
        ndata <- newdata[i,]
		print(ndata)
        L <- log(object$apriori) +
            apply(log(sapply(attribs, function(v) {
				print(paste("v", v))
                nd <- ndata[v]
                if(is.na(nd)) {
                    rep(1, length(object$apriori))
					print("bad")
				}
                else {
                    prob <- if (isnumeric[v]) {
                        msd <- object$tables[[v]]
                        msd[,2][msd[,2]==0] <- threshold
                        dnorm(nd, msd[,1], msd[,2])
                    } else
                        object$tables[[v]][,nd]
                    prob[prob == 0] <- threshold
					print(prob)
					if (!is.null(min.prob)) {
						prob[prob < min.prob] = min.prob
						prob[prob > 1-min.prob ] = 1-min.prob
					}
                    prob# + min.prob
                }
            })), 1, sum)
        if (type == "class")
            L
        else {
            L <- exp(L)
            L / sum(L)
        }
    })
    if (type == "class")
        factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
    else
        t(L)
}
predict.naiveBayes <- function(object,
                               newdata,
                               type = c("class", "raw"),
                               threshold = 0.000001,
							   min.prob = NULL,
                               ...) {
	#if (!is.null(min.prob))
		#print(sprintf( "thresholding %f", min.prob))
    type <- match.arg(type)
    newdata <- as.data.frame(newdata)
	if (length(object$tables) ==1 )
		stop("Naive Bayes must have more than 1 feature")
	#print(object$tables)
	#print(paste("test", names(newdata)))
    attribs <- intersect(names(object$tables), names(newdata))
	if (length(attribs) == 0) {
		stop("No intersection in features")
	}
    isnumeric <- sapply(newdata, is.numeric)
    newdata <- data.matrix(newdata)
    L <- sapply(1:nrow(newdata), function(i) {
        #ndata <- newdata[i,]
		#print(ndata)
        L <- log(object$apriori) +
            apply(log(sapply(attribs, function(v) {
				#print(paste("v", v))
                nd <- newdata[i,v]
                if(is.na(nd)) {
                    rep(1, length(object$apriori))
					print("bad")
				}
                else {
                    prob <- if (isnumeric[v]) {
                        msd <- object$tables[[v]]
                        msd[,2][msd[,2]==0] <- threshold
                        dnorm(nd, msd[,1], msd[,2])
                    } else
                        object$tables[[v]][,nd]
                    prob[prob == 0] <- threshold
					if (!is.null(min.prob)) {
						prob[prob < min.prob] = min.prob
						prob[prob > 1-min.prob ] = 1-min.prob
					}
                    prob# + min.prob
                }
            })), 1, sum)
        if (type == "class")
            L
        else {
            L <- exp(L)
            L / sum(L)
        }
    })
    if (type == "class")
        factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
    else
        t(L)
}

coefs.wrapper = function(vec, trim.perc, mad, coef) {
	get.norm.coefs(vec, trim.perc, mad)[coef]
}


fitted.naiveBayes.old = function(training, response, do.log=T, plot=F, trim.perc=.01, mad=F) {
	if (do.log)
		training = log(training)
	data = cbind(training, response)
	names(data) = c(names(training), "Response")
	class1 = training[response==levels(response)[1],]
	class2 = training[response==levels(response)[2],]
	out = list()
	out[["apriori"]] = table(all.ben_smk.stages$Response)
	out[["tables"]] = list()
	if (plot) {
		par(mfrow=get.row.col(ncol(training)), mar=c(2,2,2,2))
		print(get.row.col(ncol(training)))
	}
	for (feature in names(training)) {
		coefs1 = get.norm.coefs(class1[,feature], trim.perc, mad, plot, add=F, feature=feature, class.label=levels(response)[1], col='green')
		#coefs1 = get.norm.coefs(class1[,feature], trim.perc, !mad, plot, col='green', add=T, feature=feature, class.label=levels(response)[1])
		coefs2 = get.norm.coefs(class2[,feature], trim.perc, mad, plot, add=T, feature=feature, class.label=levels(response)[2], col='red')
		#coefs2 = get.norm.coefs(class2[,feature], trim.perc, !mad, plot, col='green', add=T, feature=feature, class.label=levels(response)[2])
		out[["tables"]][[feature]] = matrix( c(coefs1, coefs2), byrow=T, nrow=2)
		row.names(out[["tables"]][[feature]]) = levels(response)
	}
	out[["levels"]] = levels(response)
	class(out) = "naiveBayes"
	out
}



get.norm.coefs = function(data, trim.perc=0, mad=F, plot=F, add=T, xlab="log(RFU)", col='red', feature="", class.label="", trim.seq_id=F, main=paste(feature, class.label), lwd=1, ...) {
	x = trim(data, trim.perc)
	y = rank(x, tie="max")/length(x)
	if (mad) 
		params = c( median(x), IQR(x)/1.349)
	else {
		nls.data = data.frame(x=x,y=y)
		data.nls = nls(y ~ pnorm(x, mean=m, sd=s), data=nls.data, start=list(m=mean(x), s=sd(x)), control=nls.control(maxiter=2000, minFactor=1/1024, warnOnly=T))
		params = coef(data.nls)
	}
	if (plot) {
		if (add) 
			points(x,y, xlab="", ylab="", col=col, type='s', lwd=lwd, ...)
		else {
			if (trim.seq_id)
				feature = remove.seq_id(feature)
			plot(x,y, main=main, xlab=xlab, ylab="F(x)", col=col, type='s', lwd=lwd, ...)
		}

		points(seq(min(x), max(x), length.out=100), pnorm(seq(min(x), max(x), length.out=100), mean=params[1], sd=params[2]), type='l', col=8, cex=.8, lty=2, lwd=lwd)
	}
	params 	
}

my.norm.mle = function(data, mad=T, do.log=F) {
	#x = trim(data, trim.perc)
	if (do.log)
		x=log(data)
	else
		x=data
	y = rank(x, tie="max")/length(x)
	if (mad) 
		params = c( median(x), IQR(x)/1.349)
	else {
		nls.data = data.frame(x=x,y=y)
		data.nls = nls(y ~ pnorm(x, mean=m, sd=s), data=nls.data, start=list(m=mean(x), s=sd(x)), control=nls.control(maxiter=2000, minFactor=1/1024, warnOnly=T))
		params = coef(data.nls)
	}
	mean((y-pnorm(x,params[1], params[2]))^2)
}


apt.diffs = function(data, apts=NULL, force=F, do.log=F) {
	if (is.null(apts))
		apts = get.aptamers(data)
	if (length(apts) > 50 & !force)
		stop("Too many aptamers. Use force to use more than 50")
	apt.data = data[,apts]
	num.apts = length(apts)
	if (do.log)
		apt.data = log(apt.data)
	dists = apply(apt.data, 1, function(data.row) {
		unlist(sapply(1:(num.apts-1), function(i) {
			sapply((i+1):num.apts, function(j) {
				data.row[i] - data.row[j]
			})
		}))
	})
	ret.data = as.data.frame(t(as.data.frame(dists)))
	names(ret.data) = unlist(sapply(1:(num.apts-1), function(i) {
			sapply((i+1):num.apts, function(j) {
				paste(apts[i], apts[j], sep=".")
			})
		}))
	meta = setdiff(get.meta(data), "Response")
	if (length(meta) > 0)
		ret.data = cbind(data[,meta], ret.data)
	if ("Response" %in% get.meta(data))
		ret.data$Response = data$Response
	ret.data
}

get.apt.pairs = function(data, apts, within=F) {
	if (within)
		intersect(sapply(apts, function(apt1) sapply(apts, function(apt2) paste(apt1, apt2, sep="_"))), get.aptamers(data))
					
	else
		unique(unlist(sapply(apts, function(apt) grep(apt, get.aptamers(data), value=T))))
}

apt.diffs.pairs = function(data, apts.down, apts.up, force=F, do.log=F) {
	apts.down = match.seq.ids(apts.down, names(data))
	apts.up = match.seq.ids(apts.up, names(data))
	apt.data = data[,c(apts.down, apts.up)]
	if (do.log)
		apt.data = log(apt.data)
	dists = apply(apt.data, 1, function(data.row) {
		unlist(sapply(apts.down, function(apt.down) {
			sapply(apts.up, function(apt.up) {
				data.row[apt.up] - data.row[apt.down]
			})
		}))
	})
	ret.data = as.data.frame(t(as.data.frame(dists)))
	names(ret.data) = unlist(sapply(apts.down, function(apt.down) {
			sapply(apts.up, function(apt.up) {
				paste(apt.up, apt.down, sep="_")
			})
		}))
	meta = setdiff(get.meta(data), "Response")
	if (length(meta) > 0)
		ret.data = cbind(data[,meta], ret.data)
	if ("Response" %in% get.meta(data))
		ret.data$Response = data$Response
	ret.data
}

