
cv.snippet = function(X, y, num.cv=10) {
    mses = lapply(1:num.cv, function(i) { # for loop that returns a list (dictionary) object
        withheld.set = seq(i, nrow(X), by=num.cv)
        lasso.fit = lars(X[-withheld.set,], y[-withheld.set], type="lasso", use.Gram=FALSE)
        lasso.pred = predict(lasso.fit, X[withheld.set,])
        resid = lasso.pred- y[withheld.set];
        sse = sum(resid*resid);
        # Some mse or performance calculation here, that gets returned
    })
    unlist(mses) # returns a vector of mses. unlist concatenates the dictionary values
}
