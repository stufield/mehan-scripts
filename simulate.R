



simulate.dataset = function(nrows, ncols, means=NULL, sds=NULL, noise=0) {
	if (is.null(means))
		means = rep(1000, ncols)
	if (is.null(sds))
		sds = rep(0, ncols)
	data = as.data.frame(sapply(1:ncols, function(i) rnorm(nrows, mean=means[i], sd=((sds[i]+noise)*means[i]))))
	names(data) = paste(sprintf("%04i", 1:ncols), "00", sep="-")
	data
}


add.disease.effect = function(sim.data, rows, cols, percent) {
	sim.data[rows,cols] = apply(sim.data[rows,cols], 2, function(x) x + x*percent)
	sim.data
}
	
