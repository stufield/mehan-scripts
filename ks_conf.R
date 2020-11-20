


ks.conf.single = function(training, do.t=F) {
	out = data.frame(matrix(0, ncol=reps, nrow=length(get.aptamers(training))))
	cat('.')
	training$Response = sample(training$Response, nrow(training), replace=F)
	calc.ks(training, fdr=F, do.t=do.t)[,1]
}





ks.conf.int.cdfs = function(training, num.samples=200, no.shuffle=F, do.t=F, num.cores=7) {
	fake.training = training
	reps.per.core = num.samples/num.cores
	if (reps.per.core > 0) {
		#out.list = as.data.frame(collect(lapply(1:num.cores, function(x) mcparallel(ks.conf.sub(training, reps.per.core, do.t=do.t), mc.set.seed=T))))
		out.final = foreach(i=1:num.samples, .combine=cbind) %dopar% ks.conf.sub(training, 1, do.t=do.t)
	}
	out.final.sorted = out.final[,order(apply(out.final, 2, median))]
	cushion = max(as.integer(num.samples*.025), 1)
	conf.int = t(out.final.sorted[,c(num.samples-cushion, cushion)])
	#lower = apply(out.final, 1, function(x) quantile(p=.025))
	return(conf.int)
	ret.dat = as.data.frame(rbind(upper, lower))
	cat("\n")
	list(ret.dat, out.final)
}
ks.conf.int.mc = function(training, num.samples=200, no.shuffle=F, do.t=F, num.cores=7) {
	fake.training = training
	reps.per.core = num.samples/num.cores
	if (reps.per.core > 0) {
		#out.list = as.data.frame(collect(lapply(1:num.cores, function(x) mcparallel(ks.conf.sub(training, reps.per.core, do.t=do.t), mc.set.seed=T))))
		out.final = foreach(i=1:num.samples, .combine=cbind) %dopar% ks.conf.sub(training, 1, do.t=do.t)
	}
	conf.int = apply(out.final, 1, quantile, p=c(1-.025, .025))
	#lower = apply(out.final, 1, function(x) quantile(p=.025))
	return(conf.int)
	ret.dat = as.data.frame(rbind(upper, lower))
	cat("\n")
	list(ret.dat, out.final)
}


