# forest helper functions 
#
# Revision History
#---------------------
# $Id: $: $Date: $ 
#
 
sorted.importance = function(rf.data, remove.seq.ids=F) {
	if (nrow(rf.data$importance) == 1) {
		print(c("imp", row.names(rf.data$importance)))
		return(rf.data$importance)
	}
	table = as.matrix(rf.data$importance[order(rf.data$importance[,1], decreasing=T),])
	if (remove.seq.ids) 
		row.names(table) = sapply(row.names(table), remove.seq_id)
	colnames(table) = "Gini Importance"
	table


}

traverse.tree = function(model, tree, sample) {
	sample = sample[row.names(model$importance)]
	status = 1
	cur.row = 1
	while (tree[cur.row,5] == 1) {
		if (sample[tree[cur.row,3]] < tree[cur.row,4])
			cur.row = tree[cur.row,1]
		else
			cur.row = tree[cur.row,2]
	}
	return(list(node=cur.row, pred=model$classes[tree[cur.row,6]]))
}

print_trees = function(rf.data, trees, path) {
	for (i in trees) 
		write.csv(getTree(rf.data, i, labelVar=T), file=paste(path, "/tree", i, ".csv", sep=""))
}

reforest = function(model, data, ...) {
	randomForest(data[,row.names(model$importance)], data$Response, ntree=model$ntree, mtry=model$mtry, ...)
}

iqr.dataset = function(data, apts=NULL, range.length=50, p1=.2, p2=.8, ...) {
	if (is.null(apts))
		apts = get.aptamers(data)
	ranges = lapply(apts, function(apt) seq(quantile(data[,apt], p=p1), quantile(data[,apt], p=p2), length.out=range.length))
	iqr.mat = t(apply(all.apt.combos.recur(length(apts), range.length), 1, function(x) sapply(1:length(x), function(i) ranges[[i]][x[i]])))
	iqr.data = as.data.frame(iqr.mat, row.names=as.character(1:nrow(iqr.mat)))
	names(iqr.data) = apts
	iqr.data
}

rf.dec.boundary = function(data, model, apts) {
	stop("Bad Idea")
	if (length(apts) != 2)	
		stop("Must be two apts for RF boundary")
	total = 0
	model.apts = row.names(model$importance)
	step.size=.1
	#totals =  sapply(model.apts, function(apt) length(seq(min(data[,apt]), max(data([,apt])), by=step.size)))
	#xrange = range(data[,apts[1]])
	#yrange = range(data[,apts[2]])

}


test.site.diff = function(rf.data,  ...) {
	imp.data = data.frame(matrix(NA, nrow=dim(rf.data)[2]-1, 4))

	# Without site data
	print("Building Forest")
	rf.without = randomForest(Response ~ .-StudyName, data=rf.data,...)
	print("Adding imp data")
	blank.study = matrix(0, nrow=1, ncol=1, dimnames=list("StudyName"))
	imp = rbind(blank.study, rf.without$importance)
	imp.data[,1] = imp
	imp.data[,2] = rank(imp, ties.method="min")


	# With site Data
	print("Building Forest")
	rf.with = randomForest(Response ~ ., data=rf.data, ...)
	print("Adding imp data")
	imp = rf.with$importance
	imp.data[,3] = imp
	imp.data[,4] = rank(imp, ties.method="min")
	
	names(imp.data) = c("Without Sites Imp", "Without Sites Rank", "With Sites Imp", "With Sites Rank")
	row.names(imp.data) = names(rf.data)[1:(dim(rf.data)[2]-1)]
	imp.data

}

proxMDS = function(rf.data, data, k=2, classes=NULL, pch=NULL, col=NULL, iso=F, scores=rf.data$votes[,2], dims=1:k, ...) {
	if (iso)
		rf.mds <- isoMDS(1 - rf.data$proximity, k = k)
	else
		rf.mds <- myMDS(1 - rf.data$proximity, eig = TRUE, k = k)
		#rf.mds <- stats:::cmdscale(1 - rf.data$proximity, eig = TRUE, k = k)
	colnames(rf.mds$points) <- paste("Dim", 1:k)
	

	if (is.null(col))
		col = map.color(scores)
	
	if (is.null(classes))
		classes = data$Response
	pch=class.pch(classes)
	plot(as.data.frame(rf.mds$points)[,dims], col=col, pch=pch, ...)
	invisible(rf.mds$points)
}

create.full.proximity.matrix = function(data, model) {
	#final.mat = matrix(0, nrow=nrow(data), ncol=nrow(data)) 
	final.mat = foreach (i = 1:model$ntree, .combine="+") %dopar% {
		#final.mat = final.mat + create.proximity.matrix(data, model, i)
		cat('.')
		create.proximity.matrix(data, model, i)
	}

	#print("hi")
	#for (i in 1:length(final.mat.list)) {
		#cat('.')
		#print(dim(final.mat))
		#print(dim(final.mat.list[[i]]))
		#final.mat = final.mat + final.mat.list[[i]]
	#}
	#print("hi")


	for (i in 1:nrow(final.mat))
		final.mat[i,i] = 1
	
	ntree.mat = matrix(0, nrow=nrow(data), ncol=nrow(data)) + diag(rep(1, nrow(data)))
	for (n in 1:model$ntree) {
		inbag = which(model$inbag[,n]==1)
		outbag = setdiff(1:nrow(data), inbag)
		for (ib in inbag)  {
			for (ob in outbag) {
				ntree.mat[ib,ob] = ntree.mat[ib,ob] + 1
				ntree.mat[ob,ib] = ntree.mat[ob,ib] + 1
			}
		}
	}
	cat("\n")

	ntree.mat[ntree.mat==0] = 1

	final.mat / ntree.mat
}
create.proximity.matrix = function(data, model, tree.num) {
	if(is.null(model$inbag)) 
		stop("RF model must be run with keep.inbag=T")
	inbag = which(model$inbag[,tree.num]==1)
	corner = create.proximity.matrix.corner(data[inbag,], data[-inbag,], model, tree.num)
	top.part = cbind(diag(rep(1, nrow(data)-length(inbag))), corner)
	bottom.part = cbind(t(corner), diag(rep(1, length(inbag))))
	mat = rbind(top.part, bottom.part)
	mat2 = mat
	mapping = c(setdiff(1:nrow(data), inbag), inbag)
	mat2[mapping, mapping] = mat
	mat2

}

myMDS.project = function(model, train.data, test.prox) {
	test.data = myMDS.apply.center(model$proximity, test.prox)
	test.proj = test.data %*% myMDS(1-model$proximity)$e.vectors

	invisible(test.proj)

}

myMDS.predict = function(model, train.data, test.data, do.plot=T) {
	par(mfrow=c(1,2))
	if (!all(rn(model$votes) == rn(train.data)))
		stop("Mismatched model and training.data")
	train.mds = cmdscale(1-model$proximity)[,1]
	# Calculate proximities between test group and forest
	test.prox = test.proximity(model, train.data, test.data)
	# Project test proximities into MDS space
	test.proj = myMDS.project(model, train.data, test.prox)


	if (do.plot) {
		proxMDS(model, data=train.data, main="MDS Plot")
		points(test.proj, col=1, cex=.8, pch=19)

		true.roc(model=model, col='blue', main="Performance Comparison")
		if ("Response" %in% names(test.data)) 
			true.roc(model=model, test=test.data, col='blue', add=T, lty=2)
		switch = cor(model$votes[,2], cmdscale(1-model$proximity)[,1]) < 0
		if (switch)
			true.roc(vote.data=data.frame(vote=-train.mds, class=train.data$Response), col='red', add=2)
		else
			true.roc(vote.data=data.frame(vote=train.mds, class=train.data$Response), col='red', add=2)

		if ("Response" %in% names(test.data)) {
			if (switch)
				true.roc(vote.data=data.frame(vote=-test.proj[,1], class=test.data$Response), add=3, col='red', lty=2)
			else
				true.roc(vote.data=data.frame(vote=test.proj[,1], class=test.data$Response), add=3, col='red', lty=2)
			legend('bottomright', c("RF - Train", "RF - Test", "MDS - Train", "MDS - Test"), col=c(4,4,2,2), lty=c(1,2,1,2))
		}
		else
			legend('bottomright', c("RF - Train", "MDS - Train"), col=c(4,2), lty=c(1,1))

	}
	invisible(list(test.prox=test.prox, test.proj=test.proj))
	

}
myMDS.cv = function(x, y, num.cv=10, do.plot=F, ...) {
	full.data = cbind(x,Response=y)
	full.rf = randomForest(x,y,proximity=T,...)
	cv.vote.data = NULL
	rf.cv.vd = NULL
	#cv.vote.data.list = lapply(1:num.cv, function(i) {
	for (i in 1:num.cv) {
	#cv.vote.data = foreach (i = 1:num.cv, .combine='rbind') %do% {
		cv.indices = seq(i, nrow(x), by=num.cv)
		rf = randomForest(x[-cv.indices,], y[-cv.indices], proximity=T, keep.inbag=T, ...)
		mds = myMDS.predict(rf, full.data[-cv.indices,], full.data[cv.indices,], do.plot)
		vd = data.frame(vote = mds$test.proj[,1], class=full.data$Response[cv.indices])
		if (cor(rf$votes[,2], cmdscale(1-rf$proximity)[,1]) < 0) 
			vd[,1] = -vd[,1]
		vd
		cv.vote.data = rbind(cv.vote.data, vd)
		rf.cv.vd = rbind(rf.cv.vd, create.vote.data(model=rf, test=full.data[cv.indices,]))
	}
	#for (i in 1:length(cv.vote.data.list))
		#cv.vote.data = rbind(cv.vote.data, cv.vote.data.list[[i]])

	proxMDS(full.rf, full.data)
	true.roc(vote.data=rf.cv.vd, col='blue')
	true.roc(vote.data=cv.vote.data, add=T, col='red')
	legend('bottomright', c("RF CV", "MDS CV"), col=c('blue', 'red'), lwd=2)
	invisible(cv.vote.data)
}


balance.rf = function(x, y, num.splits, do.plot=T) {
	if (missing(num.splits))
		stop("Must provide num.splits")
	control.label = y[1]
	models = lapply(1:num.splits, function(i) {
		x.controls = x[y==control.label,]
		y.controls = y[y==control.label]
		x.cases = x[y!=control.label,]
		y.cases = y[y!=control.label]
		which.controls = seq(i, as.integer(nrow(x.controls)/num.splits)*num.splits, num.splits)
		which.cases = seq(1, nrow(x.cases))

		print(nrow(x.controls[which.controls,]))
		print(nrow(x.cases[which.cases,]))

		temp.x = rbind(x.controls[which.controls,], x.cases[which.cases,])
		temp.y = factor(c(y.controls[which.controls], y.cases[which.cases]))

		rf = randomForest(temp.x, temp.y)
		if (do.plot)
			true.roc(model=rf, col=col.string[i], add=i-1)
		rf
	})
	return(models)
	final.rf = combine(models[[1]], models[[2]])
	for (i in 3:num.splits)
		final.rf = combine(final.rf, models[[i]])
	if (do.plot)
		true.roc(model=final.rf, add=num.splits, col=col.string[num.splits+1])
	final.rf


}


test.proximity = function(model, train.data, test.data) {
	#final.mat = matrix(0, nrow=nrow(test.data), ncol=nrow(train.data))
	cat(sprintf("Calculating between group proximities for %i trees\n", model$ntree))
	final.mat = foreach (i = 1:model$ntree, .combine="+") %dopar% {
		cat('.')
		test.proximity.single(train.data, test.data, model, i)
	}
	cat('\n')
	final.mat / model$ntree
}
test.proximity.single = function(train.data, test.data, model, tree.num) {
	inbag = which(model$inbag[,tree.num]==1)
	corner = create.proximity.matrix.corner(train.data[inbag,], test.data, model, tree.num)
	mat = matrix(0, nrow=nrow(test.data), ncol=nrow(train.data))
	mat[,inbag] = corner
	mat
}

create.proximity.matrix.corner = function(in.bag.samples, oob.samples, model, tree.num) {
	ib.node2sample = get.node.list(in.bag.samples, model, tree.num)
	oob.sample2node = reverse.list(get.node.list(oob.samples, model, tree.num))
	mat = matrix(0, ncol=nrow(in.bag.samples), nrow=nrow(oob.samples))
	for (i in 1:nrow(oob.samples)) {
		oob.node = oob.sample2node[[i]]
		#print(sprintf("oob sample %i in node %i", i, oob.node))
		#print(sprintf("ib samples in node %i: %s", oob.node, paste(ib.node2sample[[oob.node]], collapse=" ")))
		#sprintf("
		mat[i, ib.node2sample[[oob.node]] ] = rep(1, length(ib.node2sample[[oob.node]]))
	}
	mat

}

get.node.list = function(data, model, tree.num) {
	tree = getTree(model, k=tree.num)
	node.list = vector("list", nrow(tree))
	for (i in 1:nrow(data)) 
		node.list[[traverse.tree(model, tree, data[i,])$node]] = c(node.list[[traverse.tree(model, tree, data[i,])$node]], i)
	node.list
}
get.sym.center = function(mat) {
	mean1 = apply(mat, 2, mean)
	mat2 = scale(mat, scale=F)
	mean2 = apply(mat2, 1, mean)
	mean.mat = matrix(0, nrow=nrow(mat), ncol=ncol(mat))
	for (i in 1:nrow(mat))
		for (j in 1:nrow(mat))
			mean.mat[i,j] = mean1[j] + mean2[i]
	mean.mat
}
myMDS.getx2 = function (d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE) {
	x = as.matrix(d^2)
	storage.mode(x) = 'double'
	x = t(scale(t(scale(x, scale=F)), scale=F))
	-x/2
}
myMDS.apply.center = function (prox.matrix.orig, prox.matrix.new, k=2) {
	d.orig = 1-prox.matrix.orig
	d.new = 1-prox.matrix.new
	x.orig = as.matrix(d.orig^2)
	x.new = as.matrix(d.new^2)
	storage.mode(x.orig) = 'double'
	storage.mode(x.new) = 'double'
	col.means = apply(x.orig, 2, mean)
	ret.matrix = t(t(x.new) - col.means)
	ret.matrix = t(scale(t(ret.matrix), scale=F))
	-ret.matrix/2
}
myMDS = function (d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE) 
{
    if (any(is.na(d))) 
        stop("NA values not allowed in 'd'")
    if (is.null(n <- attr(d, "Size"))) {
        if (add) 
            d <- as.matrix(d)
        x <- as.matrix(d^2)
        if ((n <- nrow(x)) != ncol(x)) 
            stop("distances must be result of 'dist' or a square matrix")
        rn <- rownames(x)
    }
    else {
        x <- matrix(0, n, n)
        if (add) 
            d0 <- x
        x[row(x) > col(x)] <- d^2
        x <- x + t(x)
        if (add) {
            d0[row(x) > col(x)] <- d
            d <- d0 + t(d0)
        }
        rn <- attr(d, "Labels")
    }
    if ((k <- as.integer(k)) > n - 1 || k < 1) 
        stop("'k' must be in {1, 2, ..  n - 1}")
    storage.mode(x) <- "double"
    .C(stats:::R_dblcen, x, as.integer(n), DUP = FALSE)
    if (add) {
        i2 <- n + (i <- 1L:n)
        Z <- matrix(0, 2L * n, 2L * n)
        Z[cbind(i2, i)] <- -1
        Z[i, i2] <- -x
        Z[i2, i2] <- .C(R_dblcen, x = 2 * d, as.integer(n))$x
        e <- eigen(Z, symmetric = FALSE, only.values = TRUE)$values
        add.c <- max(Re(e))
        x <- matrix(double(n * n), n, n)
        non.diag <- row(d) != col(d)
        x[non.diag] <- (d[non.diag] + add.c)^2
    }
    e <- eigen(-x/2, symmetric = TRUE)
    ev <- e$values[1L:k]
    if (any(ev < 0)) 
        warning(gettextf("some of the first %d eigenvalues are < 0", 
            k), domain = NA)
    points <- e$vectors[, 1L:k, drop = FALSE] %*% diag(ev, k)
    #points <- e$vectors[, 1L:k, drop = FALSE] %*% diag(sqrt(ev), k)
    dimnames(points) <- list(rn, NULL)
    if (eig || x.ret || add) {
        evalus <- e$values[-n]
        list(points = points, eig = if (eig) ev, x = if (x.ret) x, 
            ac = if (add) add.c else 0, GOF = sum(ev)/c(sum(abs(evalus)), 
                sum(evalus[evalus > 0])), e.vectors=e$vectors[, 1L:k, drop = FALSE])
    }
    else list(points=points, e.values=ev, e.vectors=e$vectors[, 1L:k, drop = FALSE])
}



myMDS.getx = function (d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE) 
{
    if (any(is.na(d))) 
        stop("NA values not allowed in 'd'")
    if (is.null(n <- attr(d, "Size"))) {
        if (add) 
            d <- as.matrix(d)
        x <- as.matrix(d^2)
        if ((n <- nrow(x)) != ncol(x)) 
            stop("distances must be result of 'dist' or a square matrix")
        rn <- rownames(x)
    }
    else {
        x <- matrix(0, n, n)
        if (add) 
            d0 <- x
        x[row(x) > col(x)] <- d^2
        x <- x + t(x)
        if (add) {
            d0[row(x) > col(x)] <- d
            d <- d0 + t(d0)
        }
        rn <- attr(d, "Labels")
    }
    if ((k <- as.integer(k)) > n - 1 || k < 1) 
        stop("'k' must be in {1, 2, ..  n - 1}")
    storage.mode(x) <- "double"
    .C(stats:::R_dblcen, x, as.integer(n), DUP = FALSE)
	return(-x/2)
}

all.apt.combos.recur = function(num.apts, num.values, cur.vec=rep(1, num.apts)) {
	if (all(cur.vec == rep(num.values, num.apts)))
		return(cur.vec)
	new.vec = cur.vec
	new.vec[num.apts] = cur.vec[num.apts] + 1
	for (i in seq(num.apts, 1, -1)) {
		if (new.vec[i] <= num.values) # successfully incremented
			break
		new.vec[i] = 1
		new.vec[i-1] = new.vec[i-1] + 1
	}
	rbind(cur.vec, all.apt.combos.recur(num.apts, num.values, new.vec))
}
		

#all.apt.combos = function(num.apts, num.values) {
	#mat = NULL
	#depth = 0
	#for (i in 1:num.values) {
		#


		
