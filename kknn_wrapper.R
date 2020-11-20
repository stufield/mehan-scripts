


# My wrapper for kknn that allows use with the Mike in a Box
# See help for kknn for more details
my.kknn = function(form, train, test, k=15, distance=2, kernel='triangular', ...) {
	ret.knn = kknn(form, train, test, k=min(k, nrow(train)), distance=distance, kernel=kernel, ...)
	ret.knn$Response = test$Response
	ret.knn$train = deparse(substitute(train))
	row.names(ret.knn$prob) = row.names(test)
	ret.knn$k = k
	#ret.knn$distance=distance
	ret.knn$kernel = kernel
	ret.knn
}
