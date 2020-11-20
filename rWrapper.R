


cbind.rownames = function(data) {
	data2 = as.data.frame(cbind(Aptamers=row.names(data), as.data.frame(data)))
	return(data2)
}



na.to.none = function(adat, ...) {
	for (meta in get.meta(adat)) {
		if (class(adat[,meta]) == "factor")
			adat[,meta] = as.character(adat[,meta])
		if (class(adat[,meta]) == "character") {
			adat[,meta][is.na(adat[,meta])] = "None"
			adat[,meta] = factor(adat[,meta])
		}
	}
	adat
}


