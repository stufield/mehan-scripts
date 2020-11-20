



get.standard.curve.data = function(clone.screen.data, apt, ...) {
	apt.clone.data = clone.screen.data[grep(apt2seqid(apt), clone.screen.data$APTAMERID),]
	lapply(1:nrow(apt.clone.data), function(i) {
		apt.clone.row = apt.clone.data[i,]
		concentration = as.numeric(strsplit(as.character(apt.clone.row$STDCURVE_CONCENTRATION), ',')[[1]])
		rfu = as.numeric(strsplit(as.character(apt.clone.row$STDCURVE_RFU), ',')[[1]])
		plasma_dil = as.numeric(strsplit(as.character(apt.clone.row$PLASMA_DILUTIONFACTOR), ',')[[1]])
		plasma_rfu = as.numeric(strsplit(as.character(apt.clone.row$PLASMA_RFU), ',')[[1]])
		apt.plateau = 10^as.numeric(as.character(apt.clone.row$PLATEAU))
		apt.baseline = 10^as.numeric(as.character(apt.clone.row$BASELINE))
		kd = as.numeric(as.character(apt.clone.row$KD))
		log.concentration = log(concentration, base=10)
		#print(concentration)
		#for (i in 1:length(log.concentration)) {
			#if (!is.finite(log.concentration[i])) {
				##decay = concentration[length(log.concentration)-1]/concentration[length(log.concentration)]
				#print(decay)
				#log.concentration[i] = log(concentration[length(log.concentration)]/decay, base=10)
			#}
		#}
		log.concentration[concentration == 0] = min(log.concentration, na.rm=TRUE)

		data = data.frame(RFU=log(rfu, base=10), conc=concentration);
		
		#start = list(baseline=min(rfu), slope=1, inflec=10^-10, plateau=max(rfu))
		start = list(baseline=min(rfu), inflec=10^-10, plateau=max(rfu))
		#print(start)

		plot(log(data$conc, base=10), data$RFU, xlab="Concentration", ylab="RFU", ...)

		slope = 1
		apt.nls = nls( RFU ~  log((plateau-baseline)*( conc^slope / (conc^slope + inflec^slope) ) + baseline, base=10),  data=data, start=start)
		coefs = coefficients(apt.nls)
		if (!"slope" %in% names(coefs)) {
			coefs = c(coefs, slope)
			names(coefs)[length(coefs)] = 'slope'
		}
		else
			slope = coefs['slope']
		conc.range =  10^seq(min(log(data$conc+10^-30, base=10)), max(log(data$conc, base=10)), length.out=5000)
		points( log(conc.range, base=10), log((coefs['plateau']-coefs['baseline'])*( conc.range^slope / (conc.range^slope + coefs['inflec']^slope) ) + coefs['baseline'], base=10), col='red', type='l')
		
		abline(v=log(coefs['inflec'], base=10), col='blue', lty=3)

		#print(coefs)
		coefs	

	})
}

rfu2conc = function(rfu, coefs, do.log=T) {
	plateau = coefs['plateau']
	baseline = coefs['baseline']
	slope = coefs['slope']
	kd = coefs['inflec']
	Y = (rfu-baseline)/(plateau-baseline)
	ret = (Y/(1-Y))^(1/slope) * kd
	ifelse(do.log, log(ret, base=10), ret)
}

conc2rfu = function(conc,coefs) {
	plateau = coefs['plateau']
	baseline = coefs['baseline']
	slope = coefs['slope']
	kd = coefs['inflec']
	(plateau-baseline)*( conc^slope / (conc^slope + kd^slope) ) + baseline
}

mw.to.M = function(pg.ml, kDa) {
	pg.ml / 10^12 / kDa 
}




