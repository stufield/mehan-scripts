# This is called by true.roc
# Given a vote.data structure and a positive predicted class, returns the x,y coordinates of a ROC plot 

# Creates a ROC curve using the ordering of the predictions. Called true roc because I was doing it differently before
# It always needs a positive class label (from create.training.data) and it can be called in two ways:
# 1) With a model and possibly test data. If no test data is provided the training performance from the model is printed. Note this will not work for models like Naive Bayes and KNN which need to be applied to the training or test data
# 		Ex: true.roc("stage1-3", model=nsclc.rf)
# 		Ex: true.roc("stage1-3", model=nsclc.knn, test=nsclc.test)
# 2) Provide vote.data data object, which has all information required to create a ROC plot
# 		Ex: true.roc("state1-3", vote.data=nsclc.vote)

# Other Arguments:
# auc: Should the AUC be printed on the plot. If t
# add: Dictates whether new plot is drawn or ROC is added to current plot. Should be passed as an integer, so that if this is the second roc added, add=2. This will properly align the AUCs so they don't get written over the top of each other
# auc.x: Float in [0,1] that places the AUC on the x-axis
# auc.y: Float in [0,1] that places the AUC on the y-axis
# col: color for curve and AUC. Didn't pass it through as ... because I need to put it different places
# adj: used to align AUC, see help for 'text' function
# auc.shift: The vertical shift between AUC values when multiple ROCs are plotted
true.roc <- function(model=NULL, test=NULL, vote.data=NULL, pos.class=NULL, add=F, 
							auc=T, auc.x=.5, auc.y=.5, col=NULL, adj=c(0,0), auc.shift=1.25, 
							conf.int=F, switch=F, lty=NULL, lwd=2, pch = NULL, file=NULL, 
							min.prob=NULL, cutoff=.5, skip.nb.log=F, do.grid=T, boxes=T, box.alpha=.35, ...) {

	if (is.null(vote.data))
		vote.data = create.vote.data(model, test, switch=switch, min.prob=min.prob, skip.nb.log=skip.nb.log)
	else {
		if (is.null(pos.class)) 
			pos.class = vote.data[nrow(vote.data),2]
	}
	if (is.null(col))
		col = as.numeric(add)+1
	if (is.null(lty))
		lty = as.integer((as.numeric(add))/8)+1
	vote.data.sorted = vote.data[order(vote.data[,1], decreasing=T),]
	if (is.null(pos.class))
		pos.class = get.pos.class(model)
	plot.labels=F
	line.col = col
	roc.data = roc.xy(vote.data.sorted, pos.class)
	num.samples = nrow(vote.data)
	if (!is.null(file))
		png(file=file, width=6, height=6, res=300, units='in')
	if (!add) {
		plot(-10,-10, xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,1), ylim=c(0,1), yaxs='i', xaxs='i', lwd=lwd, ...)
		if (do.grid) {
			abline(v=seq(.2,.8,by=.2), lty=2, col=8)
			abline(h=seq(.2,.8,by=.2), lty=2, col=8)
			abline(0,1, col=8, lty=2, lwd=1)
		}
	}
	if (conf.int) {
		full.vote.data = full.vote.table(vote.data=vote.data, conf.int=T)
		if (boxes) {
			dist = abs(full.vote.data$DiseaseProb - cutoff)
			w = which(dist == min(dist))[1]
			apply(full.vote.data[w,c('Spec.Upper', 'Sens.Lower', 'Spec.Lower', 'Sens.Upper')], 1, function(x) {
			#apply(full.vote.data[seq(1, num.samples, length.out=min(num.samples, 40)),c('Spec.Upper', 'Sens.Lower', 'Spec.Lower', 'Sens.Upper')], 1, function(x) {
				rect(1-x[1], x[2], 1-x[3], x[4], col=color.alpha(col, box.alpha), border=NA)
			})
			#line.col='black'
		}

		else {
			apply(full.vote.data[seq(1, num.samples, length.out=min(num.samples, 20)),c('Sensitivity', 'Specificity', 'Spec.Upper', 'Sens.Lower', 'Spec.Lower', 'Sens.Upper')], 1, function(x) {
				segments(1-x[3], x[1], x1=1-x[5], col=col, lty=2) # Spec
				segments(1-x[2], x[4], y1=x[6], col=col, lty=2) # Sens
				points(1-x[2], x[1], col=col, cex=.5, pch=19)
			})
		}

	}
	if (conf.int & boxes)
		points(roc.data$x[1:(length(roc.data$x))], roc.data$y, type='S', col='black',  lwd=lwd*2, lty=lty, ...)
	points(roc.data$x[1:(length(roc.data$x))], roc.data$y, type='S', col=line.col,  lwd=lwd, lty=lty, ...)
		#plot(stepfun(roc.data$x[1:(length(roc.data$x)-1)], roc.data$y), xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,1), col=col, yaxs='i', xaxs='i', lwd=2, do.points=F, lty=lty,  ...)
		#plot(roc.data$x,roc.data$y, type='l', xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,1), col=col, yaxs='i', xaxs='i', lwd=2, ...)
	if (!is.null(pch))
		points(roc.data$x[1:(length(roc.data$x))], roc.data$y, cex = .75, pch=pch)
	if (plot.labels)
		text(roc.data$x, roc.data$y, as.character(vote.data.sorted[,1]), adj=c(0,2))
	#roc.auc = auc2(roc.data$x, roc.data$y)

	if (is.numeric(cutoff)) {
		sens.spec = get.sens.spec(vote.data=vote.data, cutoff=cutoff)
		points(1-sens.spec[2], sens.spec[1], cex=1.2, bg='white', col=line.col, pch=23)
	}
	#roc.auc = colAUC(vote.data.sorted$vote, vote.data.sorted$class)
	auc.data = delong.auc.var(vote.data)
	#roc.auc2 = wilcox.test(vote.data[vote.data$class==levels(vote.data$class)[2], 'vote'], vote.data[vote.data$class==levels(vote.data$class)[1], 'vote'])$statistic
	if (auc) {
		if (add)
			adj[2] = adj[2] + auc.shift*add
		text(auc.x, auc.y, col=col, labels=paste("AUC:", sprintf("%0.2f (%0.2f,%0.2f)", auc.data$auc, auc.data$lower.limit, auc.data$upper.limit)), cex=1.1, adj=adj)
	}

	if (!is.null(file))
		dev.off()
	auc.data$auc
}





	
# Vote data structure is a two column data frame containing the model prediction score and the true class label for the sample. This is called by true.roc but is often called by itself as well. See comments for true.roc for explanation of the model and test. 
create.vote.data = function(model, test=NULL, switch=F, sorted=F, correct.test=F, skip.nb.log=F, ...) {

	if (length(class(model)) > 1)
		model.class = class(model)[1] # Fix for svm
	else
		model.class = class(model)
	if (!is.null(test))  {
		if (names(test)[length(names(test))] != "Response")
			stop("Must include Response with test in create.vote.data")

	if  (model.class == "glm" ) {
		if (is.null(test)) 
			stop("Must provide test set for glm")
		common.names = intersect(attributes(model$terms)$term.labels, names(test))
		if (length(common.names) == 0)
			stop("No common features in model and test.data")
		pred.set = test[,common.names]
		if (!skip.nb.log) {
			pred.set = log(test[,common.names])
		}
		pred = predict(model, pred.set, type="response")
		vote.data = data.frame(vote=pred, class=test$Response)
		row.names(vote.data) = row.names(test)
	}
	else if  (model.class == "naiveBayes" ) {
		if (is.null(test)) 
			stop("Must provide test set for naive bayes")
		common.names = intersect(names(model$tables), names(test))
		if (length(common.names) == 0)
			stop("No common features in model and test.data")
		pred.set = test[,common.names]
		if (!skip.nb.log)
			pred.set = log(test[,common.names])
		pred = predict(model, pred.set, type="raw", ...)
		vote.data = data.frame(vote=pred[,2], class=test$Response)
		row.names(vote.data) = row.names(test)
	}
	else if  (model.class == "randomForest") {
		if (is.null(model$votes))
			stop("No votes in random forest. Did you forget to make Response a factor")
		if (is.null(test)) {
			vote.data = data.frame(vote=model$votes[,2], class=model$y)
		}
		else {
			if (!all(levels(test$Response) == model$classes)) {
				class.table = table(test$Response)
				if (all(test$Response[1:class.table[test$Response[1]]] == test$Response[1])) {
					warning(sprintf("Renaming class labels to match model. Controls: %s has become %s. Cases: %s has become %s.", test$Response[1], model$classes[1], test$Response[nrow(test)], model$classes[2]))
					test$Response = as.character(test$Response)
					test$Response[test$Response == test$Response[1]] = model$classes[1]
					test$Response[test$Response != test$Response[1]] = model$classes[2]
				}
			}
			common.names = intersect(row.names(model$importance), names(test))
			if (length(common.names) == 0)
				stop("No common features in model and test.data")

			old = intersect(row.names(test), row.names(model$votes))
			old.vote = data.frame(vote=model$votes[old,2], class=test[old,]$Response)
			vote.data = old.vote


			new = setdiff(row.names(test), row.names(model$votes))
			if (length(new) > 0) {
				new.pred = predict(model, test[new,common.names], type="vote", ...)
				new.vote = data.frame(vote=new.pred[,2], class=test[new,]$Response)
				vote.data = rbind(old.vote, new.vote)
			}
		}
	}
	else if  (model.class == "svm") {
		if (is.null(test)) {
			stop("Must provide test set for svm")
		}
		else {
			pred = predict(model, test[,-ncol(test)], probability=T, ...)
			pred2 = attr(pred, "probabilities")
			vote.data = data.frame(vote=pred2[,2], class=test$Response)
		}
	}
	else if  (model.class == "gbm") {
		if (is.null(test)) {
			vote.data = data.frame(vote=model$fit, class=model$classes, row.names = names(model$classes))
		}
		else {
			common.names = intersect(model$var.names, names(test))

			old = intersect(row.names(test), names(model$classes))
			row.names.bool = names(model$classes) %in% old
			old.vote = data.frame(vote=model$fit[row.names.bool], class=test[old,]$Response)
			vote.data = old.vote


			new = setdiff(row.names(test), names(model$classes))
			if (length(new) > 0) {
				new.pred = predict(model, test[new,common.names], n.trees=model$n.trees, ...)
				new.vote = data.frame(vote=new.pred, class=test[new,]$Response)
				vote.data = rbind(old.vote, new.vote)
			}
		}
	}
	else if (model.class == "kknn") {
		if (is.null(test)) {
			vote.data = data.frame(vote=model$prob[,2], class=model$Response, row.names = row.names(model$prob))
		}
		else {
			test.knn = my.kknn(Response ~ ., eval(parse(text=model$train)), test[,c(attr(model$terms, "term.labels"), "Response")], k=model$k, kernel=model$kernel, distance=model$distance)
			vote.data = data.frame(vote=test.knn$prob[,2], class=test.knn$Response, row.names = row.names(test.knn$prob))

		}

	}
	else 
		stop(paste("invalid class to create.vote.data: ", model.class))
	vote.data = vote.data[!is.nan(vote.data$vote),]
	if (switch) {
		vote.data$vote = -vote.data$vote
	}
	if (sorted)
		vote.data = vote.data[order(vote.data$vote, decreasing=T),]
	vote.data
}















############################################
################## Out-Dated ###############
############################################

# Not used
plot.roc = function(roc.data, add=F, plot.labels=F, ...) {
	if (is.null(roc.data)) {
		print("No data provided to plot.roc")
		if (!add)
			plot(-1, -1, type='l', xlab="1-Specificity", ylab="Sensitivity", ...)
		return(NULL)
	}
	if (add)
		points(1-roc.data$tnr, roc.data$tpr, type='l', ...)
	else
		plot(1-roc.data$tnr, roc.data$tpr, type='l', xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,1), ...)
	abline(0,1, col='black', lty=2, lwd=2)
	if (plot.labels)
		text(1-roc.data$tnr, roc.data$tpr, as.character(roc.data$cutoff), adj=c(0,2))
}
# Not used
plot.pr = function(roc.data, add=F, plot.labels=F, ...) {
	if (is.null(roc.data)) {
		print("No data provided to plot.pr")
		if (!add)
			plot(-1, -1, type='l', xlab="1-Specificity", ylab="Sensitivity", ...)
		return(NULL)
	}
	if (add)
		points(roc.data$tpr, roc.data$ppv, type='l', ...)
	else
		plot(roc.data$tpr, roc.data$ppv, type='l', xlab="Recall", ylab="Precision", xlim=c(0,1), ...)
	if (plot.labels)
		text(roc.data$tpr, roc.data$ppv, as.character(roc.data$cutoff), adj=c(1,2))
}

# Not used
plot.both = function(roc.data, plot.labels=F, ...) {
	par(mfrow=c(2,1))
	plot.roc(roc.data, main="ROC Plot", col='blue', plot.labels=plot.labels, ...)
	plot.pr(roc.data, main="PR Plot", col='red', plot.labels=plot.labels, ...)
}

# Not Used
auc = function(roc.data) {

	x = c(1, 1-roc.data$tnr, 0)
	y = c(1, roc.data$tpr, 0)
	total = 0
	for (i in seq(2, length(x))) 
		total = total + trapez(x[(i-1):i], y[(i-1):i])
	total
}

# Not Used 
add.roc.stats = function(vote.data, pos.class="stage1-3", neg.class="ben_smk") {
	o = order(vote.data[,1], decreasing=T)
	vote.data = vote.data[o,]
	num.pos = sum(as.numeric(vote.data$class==pos.class))
	num.neg = sum(as.numeric(vote.data$class!=pos.class))
	ret.table = data.frame(Discr=NULL, Class=NULL, Sensitivity=NULL, CI_Low=NULL, CI_High=NULL, Specificity=NULL, CI_Low=NULL, CI_High=NULL, Precision=NULL, CI_Low=NULL, CI_High=NULL)
	for (row in 1:(nrow(vote.data)-1)) {
		calls = c(rep(pos.class, row), rep(neg.class, nrow(vote.data)-row))
		confusion.mat = table(vote.data$class, calls)
		confusion.stats = confusion.stats(confusion.mat)
		ret.table = rbind(ret.table, data.frame(Discr=vote.data[row,"vote"], Class=vote.data[row, "class"], Sensitivity=confusion.stats$stats[1], CI_Low=confusion.stats$intervals$Sensitivity[1], CI_High=confusion.stats$intervals$Sensitivity[2], Specificity=confusion.stats$stats[2], CI_Low=confusion.stats$intervals$Specificity[1], CI_High=confusion.stats$intervals$Specificity[2], Precision=confusion.stats$stats[3], CI_Low=confusion.stats$intervals$Precision[1], CI_High=confusion.stats$intervals$Precision[2]))
	}
	ret.table[, 3:ncol(ret.table)] = apply(ret.table[,3:ncol(ret.table)], 2, function(x) {as.integer(.5+x*100)})
	ret.table
}


# Single use, won't need
plot.feature.pair = function(training.data, features, vote.data, vote.data2=NULL, option=1) {
	par(mfrow=c(length(features), length(features)))
	for (i in seq(length(features))) {
		feature1 = features[i]
		for (j in seq(length(features))) {
			feature2 = features[j]

			if (feature1 == feature2) {
				plot(ecdf(training.data[training.data$Response==levels(training.data$Response)[1],feature1]), verticals=T, do.points=F, col.hor='green', col.vert='green', main=feature1)
				plot(ecdf(training.data[training.data$Response==levels(training.data$Response)[2],feature1]), verticals=T, do.points=F, col.hor='red', col.vert='red', add=T)
				#plot(-1,-1,axes=F, xlim=c(0,1), ylim=c(0,1), xlab=" ", ylab=" ")
				#text(.5,.5,feature1, cex=1.5)
			}
			if (i < j) {
				if (option==1)
					plot(training.data[,feature1], training.data[,feature2], pch=c(19, 17)[as.numeric(training.data$Response)], col=c('green', 'red')[as.integer(vote.data[,1]>.5)+1], xlab="RFU", ylab="RFU", main=" ")
				if (option==2)
					plot(training.data[,feature1], training.data[,feature2], pch=c(19, 17)[as.numeric(training.data$Response)],  col=rainbow(1000, start=.4, end=)[as.integer(nsclc10.rf.class_data[,1]*1000)], xlab="RFU", ylab="RFU", main=" ")
				if (option==3)
					plot(training.data[,feature1], training.data[,feature2], pch=c(4, 20)[as.numeric(training.data$Response)],  col=heat.colors(1000)[as.integer(nsclc10.rf.class_data[,1]*1000)], xlab="RFU", ylab="RFU", main=" ")
				if (option==4)
					plot(training.data[,feature1], training.data[,feature2], pch=c(19, 17)[as.numeric(training.data$Response)],  col=topo.colors(1000)[as.integer(nsclc10.rf.class_data[,1]*1000)], xlab="RFU", ylab="RFU", main=" ")

			}
			if (j < i) {
				if (option==1)
					plot(training.data[,feature2], training.data[,feature1], pch=c(19, 17)[as.numeric(training.data$Response)], col=c('green', 'red')[as.integer(vote.data2[,1]>.5)+1], xlab="RFU", ylab="RFU", main=" ")
				if (option==2)
					plot(training.data[,feature2], training.data[,feature1], pch=c(19, 17)[as.numeric(training.data$Response)],  col=rainbow(1000, start=.4, end=)[as.integer(nsclc10.rf.class_data[,1]*1000)], xlab="RFU", ylab="RFU", main=" ")
				if (option==3)
					plot(training.data[,feature2], training.data[,feature1], pch=c(4, 20)[as.numeric(training.data$Response)],  col=heat.colors(1000)[as.integer(nsclc10.rf.class_data[,1]*1000)], xlab="RFU", ylab="RFU", main=" ")
				if (option==4)
					plot(training.data[,feature2], training.data[,feature1], pch=c(19, 17)[as.numeric(training.data$Response)],  col=topo.colors(1000)[as.integer(nsclc10.rf.class_data[,1]*1000)], xlab="RFU", ylab="RFU", main=" ")

			}
		}
		
	}
}




# not used
nb.roc.wrapper = function(nb.data, training, response, pos.class, do.plot=T, ...) {
	pred = predict(nb.data, training, type="raw")
	vote.data = data.frame(vote=pred[,2], class=response)
	vote.data = vote.data[!is.nan(vote.data$vote),]
	roc.data = eval.roc.pr(vote.data, pos.class, ...)
	if (do.plot)
		plot.both(roc.data, T)
	list(vote.data=vote.data, roc.data=roc.data)
}

# not used
rf.roc.wrapper = function(rf.data, pos.class, do.plot=T, ...) {
	vote.data = create.vote.data(rf.data)
	roc.data = eval.roc.pr(vote.data, pos.class, ...)
	if (do.plot)
		plot.both(roc.data, T)
	list(vote.data=vote.data, roc.data=roc.data)
}


