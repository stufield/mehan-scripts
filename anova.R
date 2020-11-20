
create.lm.table = function(data, response, as.factor=F) {
	if (as.factor)
		data[,response] = as.numeric(factor(data[,response]))
	else if (is.factor(data[,response])) 
		data[,response] = as.numeric(as.character(data[,response]))
	data.anova = as.data.frame(t(apply(as.data.frame(t(get.aptamers(data))), 2, function(x) get.lm.data(data, response, x))))
	row.names(data.anova) = get.aptamers(data)
	names(data.anova) = c("Intercept", "Slope", "p-value")
	data.anova$qvalue = qvalue(data.anova[,3])$qvalues
	cors = sapply(get.aptamers(data), function(x) cor(data[,response], data[,x]))
	names(data.anova[4]) = "q-value"
	data.anova$Correlation = cors
	data.anova = data.anova[order(data.anova[,3]),]
	data.anova
}

create.anova.table = function(data, response) {
	data.anova = as.data.frame(t(as.data.frame(lapply(get.aptamers(data), function(apt) {
		form = as.formula(sprintf("%s ~ %s", apt, response))
		unlist(anova(lm(form, data=data)))
	}))))
	row.names(data.anova) = get.aptamers(data)
	data.anova$qvalue = qvalue(data.anova[,"Pr(>F)1"])$qvalues
	data.anova = data.anova[, c("Pr(>F)1", "qvalue")]
		
	data.hsd = as.data.frame(t(as.data.frame(lapply(get.aptamers(data), function(apt) {
		form = as.formula(sprintf("%s ~ %s", apt, response))
		hsd.data = as.data.frame(TukeyHSD(aov(form, data=data))[[response]])
		adj.p = hsd.data[['p adj']]
		adj.p = as.numeric(adj.p)
		names(adj.p) = row.names(hsd.data)
		adj.p
	}))))
	row.names(data.hsd) = get.aptamers(data)
	final.data = cbind(data.anova, data.hsd)	
	final.data = final.data[order(final.data[,"Pr(>F)1"]),]
	names(final.data) = c("ANOVA pv", "ANOVA qv", paste(names(final.data)[3:ncol(final.data)], "pv"))
	final.data
}

fit.lm.model = function(data, response, filename=F, ...) {
	table = get.sig.table(create.lm.table(data, response), ...)
	apts = row.names(table)
	print('blah')
	if (length(apts) == 0)
		return(NULL)

	num.apts = length(apts)
	errs = c()
	min.err = 10000000
	min.apts = apts
	if (is.factor(data[,response])) {
		data[,response] = as.numeric(as.character(data[,response]))
	}
	for (size in seq(num.apts, 1, by=-1)) {
		lm.summary = summary(lm(as.formula(sprintf("%s ~ .",response)), data=data[,c(response, apts)]))
		table = lm.summary[[4]]
		table = table[order(table[,4], decreasing=T),]
		apts = get.aptamers.list(row.names(table[2:nrow(table),]))
		err = lm.summary[[6]]
		errs = c(errs, err)
		if (lm.summary[[6]] < min.err) {
			min.err = err
			min.apts = apts
		}
	}

	if (filename)
		pdf(width=7, height=7, file=sprintf("plots/%s_err.pdf", response))
	plot(seq(num.apts, 1, by=-1), errs, main=response, xlab="Model Size", ylab="Residual Standard Error", type='l')
	points(length(min.apts), min.err, pch=4)
	if (filename)
		dev.off()

	min.apts
}

get.lm.data = function(data, response, apt) {
	data[[response]] = as.numeric(as.character(data[[response]]))
	data.lm = lm(as.formula(sprintf("%s ~ %s", response, apt)), data=data)
	pv = unlist(summary(data.lm))[["coefficients8"]]
	#print(data.lm$coef)
	#data.frame(Intercept=data.lm$coef[1], Slope=data.lm$coef[2], pv)
	c(data.lm$coef[1], data.lm$coef[2], pv)
}

