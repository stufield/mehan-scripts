#library(klaR)
#library(e1071)






total.auc.old = function(training, test1, test2, test3) {

	training_data.rf = randomForest(Response ~ . -StudyName, data=training, ntree=10)
	par(mfrow=c(2,2))
	get.auc(training_data.rf, test1) + get.auc(training_data.rf, test2) + get.auc(training_data.rf, test3)
}
total.auc = function(formula, training, test) {

	training_data.rf = randomForest(as.formula(formula), data=training, ntree=10)
	get.auc(training_data.rf, test)*3 
}

nb.wrapper = function(dataset, features=NULL) {
	if (is.null(features))
		features = 2:(ncol(dataset)-1)
	naiveBayes(dataset[,features], dataset$Response)
}
	
get.auc = function(model, testing_data, cols, type) {
	if (type=="rf") {
		testing_data.pred = predict(model, testing_data[,cols], type="vote")
		pos.class = model$classes[2]
	}
	else if (type == "nb") {
		pos.class = model$levels[2]
		testing_data.pred = predict(model, testing_data[,cols], type="raw")
		#testing_data.pred = predict(model, testing_data[,cols])$posterior
	}
	else {
		print("Invalid type to get.auc")
		return(NULL)
	}
	vote.data = data.frame(vote=testing_data.pred[,2], class=testing_data$Response)
	vote.data = vote.data[!is.nan(vote.data$vote),]
	roc.data = eval.roc.pr(vote.data, pos.class=pos.class)
	auc(roc.data)
}


calc_post = function(model, row, type='raw') {
	features = names(row)

	feature1 = features[1]
	feature2 = features[2]

	prior1 = model$apriori[1]
	prior2 = model$apriori[2]


	smk.mean1 = model$tables[[feature1]][1,1]
	smk.sd1 = model$tables[[feature1]][1,2]
	stage.mean1 = model$tables[[feature1]][2,1]
	stage.sd1 = model$tables[[feature1]][2,2]

	smk.mean2 = model$tables[[feature2]][1,1]
	stage.mean2 = model$tables[[feature2]][2,1]
	smk.sd2 = model$tables[[feature2]][1,2]
	stage.sd2 = model$tables[[feature2]][2,2]

#	print(paste("feature1", feature1))
#	print(paste("smk.mean1", smk.mean1))
#	print(paste("smk.sd1", smk.sd1))
#	print(paste("stage.mean1", stage.mean1))
#	print(paste("stage.sd1", stage.sd1))
#	print(paste("feature2", feature2))
#	print(paste("smk.mean2", smk.mean2))
#	print(paste("smk.sd2", smk.sd2))
#	print(paste("stage.mean2", stage.mean2))
#	print(paste("stage.sd2", stage.sd2))
#
	val1 = row[,feature1]
	val2 = row[,feature2]
	#print(paste("vals", c(val1, val2) ) )
	post1 = prior1 * dnorm(val1, smk.mean1, smk.sd1) * dnorm(val2, smk.mean2, smk.sd2)
	post2 = prior2 * dnorm(val1, stage.mean1, stage.sd1) * dnorm(val2, stage.mean2, stage.sd2)
	p1 = post1/(post1+post2)
	p2 = post2/(post1+post2)
	#print(c(p1,p2))
	a = data.frame(log=log(p2/p1), val1=val1, val2=val2)
	print(a)



}
