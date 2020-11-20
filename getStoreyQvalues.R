#
#  Wee function to compute Qvalues and write out to matlab file
#
require(qvalue)
require(R.matlab)
getStoreyQvalues <- function(matfileBaseName) {
  mdata = readMat(paste(matfileBaseName,".mat", sep=""))
  qobj = qvalue(as.numeric(mdata$p))
  qwrite(qobj,filename=paste(matfileBaseName,"StoreyQvalues.txt", sep="_"))
  writeMat(paste(matfileBaseName,"pqvalues.mat", sep="_"), q=qobj$qvalues, p=qobj$pvalues,pi0=qobj$pi0, lambda=qobj$lambda)
}

