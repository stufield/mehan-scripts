"CoxglmnetStabilityselection" <- function(dframe,t,delta,nbootstrap=100,nsteps=20,lambda=3,alpha=0.8, Pw=0,plotme=FALSE) 
{
  # Stability selection in the spirit of Meinshausen&Buhlman using glmnet for Cox model regression
  # Input: 
  #     dframe -  nxp data frame (asumes data standardized)
  #     t      -  nx1 numeric vector of event times
  #     delta  -  nx1 numeric vector for censoring indicator (delta=0 -> censored)
  #
  # result is a list containing: 
  #     $freq          - nsteps x p selection frequency matrix 
  #     $lambda        - vector of penalty values for L1 penalty
  #     $selectionProb - vector containing the probability that each feature is selected 
  #                      during the first nsteps steps of the Lasso path when half of the 
  #                      samples are used and the features are reweigthed by a random 
  #                      weight uniformaly sampled in [alpha,1]. This probability is 
  #                      estimated by nbootstrap bootstrap samples
  require(glmnet)
  require(survival)
 	dimZ <- dim(dframe) 
	n <- dimZ[1]
	p <- dimZ[2]
	halfsize <- as.integer(n/2)
	freq <- as.data.frame(matrix(0,nsteps,p));
  targetNames <-gsub(" ", ".", colnames(dframe))
  colnames(freq) = targetNames
  y =Surv(t, delta)
  survival = cbind(time=y[,1],status=y[,2])
  Z = as.matrix(dframe)
   
  # Generate a **fixed** lambda sequence using the whole data
  gfit = glmnet(Z,survival, family="cox", nlambda=nsteps, alpha=1, standardize=FALSE, lambda.min.ratio=.1)
  lambda = gfit$lambda;
  
	for (i in seq(nbootstrap)) {
	  
  		# Ramdomly split the sample in two sets
    		perm <- sample(dimZ[1])
    		i1 <- perm[1:halfsize]
    		i2 <- perm[(halfsize+1):n]
        
    		# Randomly reweight each variable
    		perturb = rep(1,length=p);
    		draw = runif(p,0,1)
    		perturb[draw<Pw] = alpha;
    		
    		# Randomly reweight (perturb) each variable with bernoulli weights
    		Zs <- t(t(Z)*perturb);
    		
      # Generate a penalized fit
    		#gfit = glmnet(Z[i1,],survival, family="cox", lambda=lambda, alpha=1, standardize=FALSE) # fixed lambda
    		gfit = glmnet(Zs[i1,],survival, family="cox", nlambda=nsteps, alpha=1, standardize=FALSE)
    		#lambda1 = gfit$lambda  # Debuggin
    
  		 # Find the variable with non-zero coeffs
     	 	 tmp = as.matrix(gfit$beta)  # p x Nlambda matrix
     	   sel.index = (1:dim(Zs)[2])[abs(tmp[,dim(tmp)[2]])>0]
     		 #targetNames[sel.index]
    		 #freq[i,targetNames[sel.index]] = freq[i,targetNames[sel.index]] + 1;  # Update proteins selected along the way
  	   
        # Increment the counters in each row (lambda value) of the freq matrix where the jth feature appeared in the model
         for (j in 1:length(sel.index)) {
  	      lamidx = (1:dim(tmp)[2])[abs(tmp[sel.index[j],])>0];
  	      freq[lamidx,sel.index[j]] = freq[lamidx,sel.index[j]]+1;
  	     }
        
    		#gfit = glmnet(Z[i2,],survival, family="cox", lambda=lambda, alpha=1, standardize=FALSE)  #Fixed lambda
    		gfit = glmnet(Zs[i2,],survival, family="cox", nlambda=nsteps, alpha=1, standardize=FALSE)
    		#lambda2 = gfit$lambda
    		
    	  # Find the variable with non-zero coeffs
    	  tmp = as.matrix(gfit$beta)  # p x Nlambda matrix
    	  sel.index = (1:dim(Z)[2])[abs(tmp[,dim(tmp)[2]])>0]
        
    		for (j in 1:length(sel.index)) {
    		  lamidx = (1:dim(tmp)[2])[abs(tmp[sel.index[j],])>0];
    		  freq[lamidx,sel.index[j]] = freq[lamidx,sel.index[j]]+1;
    		}
    		
	}
		
	# normalize frequence in [0,1]
	freq <- freq/(2*nbootstrap)
	
	if (plotme) {
		matplot(lambda, freq,type='l',xlab="Lambda",ylab="Frequency")
	}
	
	# the final stability score is the maximum frequency over the steps
	#result <- apply(freq,2,max)
 
  result <- list(freq=freq, lambda=lambda, selectionProb=apply(freq,2,max), numSamples=2*nbootstrap)
  
}

#
#
#
"CoxglmStabilityPlot" <- function(object, targetNames=colnames(object$freq), labelsize = 0.6, pLabelThresh=0.3, ymax=1, normalize=FALSE,...)
{
  
freq   = object$freq;
if ( normalize ) {
  freq = freq / apply(freq,1,sum);
}

lambda = object$lambda;
remove = apply(freq, 2, function(x) all(x==0)) 
if (all(remove)) 
  stop("all frequencies are zero for all values of lambda in this object")
label <- "lambda1/lambda_max"

labwidth <- ifelse(labelsize > 0, max(strwidth(colnames(freq[,!remove]), "inches", labelsize)), 0)
margins <- par("mai")
par(mai = c(margins[1:3], max(margins[4], labwidth * 1.4)), xaxs="i", yaxs="i")
if ( normalize ) {
  matplot((lambda)/max(lambda), freq[,!remove,  drop = FALSE], type = "l", axes=F,
          ylab="Prob [Select]",xlab = label, col = rainbow(sum(!remove)), 
          xlim = range(lambda/max(lambda)), ylim=c(0,ceiling(10*ymax)/10+1e-5), main="Normalized Cox glm Regression Stability Paths",...)
  
}
else {
  matplot((lambda)/max(lambda), freq[,!remove,  drop = FALSE], type = "l", axes=F,
          ylab="Prob [Select]",xlab = label, col = rainbow(sum(!remove)), 
          xlim = rev(range(lambda/max(lambda))), ylim=c(0,ceiling(10*ymax)/10+1e-5), main="Cox glm Regression Stability Paths",...)
}
axis(1,at=seq(0,1,0.2));
axis(2,at=seq(0,1,0.1))
#grid(NA,ny=ceiling(ymax*10))
abline(h=seq(0.1,floor(ymax*10)/10,0.1), lty=2, col=8)

if (labelsize > 0 && !is.null(colnames(freq))) {
  # Identify features which selection probabilities that exceed labeling threshold
  take <- which(!remove )
  selectNames = character(0); 
  selectIndex = vector(length=length(take))
  for (i in 1:length(take)) {
    j <- take[i]
    if (max(freq[,j]) >= pLabelThresh) {
      selectNames = paste(selectNames, colnames(freq)[j]);
      selectIndex[i] = j;
      if ( normalize ) { 
        axis(4, at = freq[1,j], labels = targetNames[j], 
             las = 1, cex.axis = labelsize, col.axis = rainbow(sum(!remove))[i], 
             lty = (i - 1)%%5 + 1, col = rainbow(sum(!remove))[i])
        
      }
      else {
        axis(4, at = freq[nrow(freq),j], labels = targetNames[j], 
             las = 1, cex.axis = labelsize, col.axis = rainbow(sum(!remove))[i], 
             lty = (i - 1)%%5 + 1, col = rainbow(sum(!remove))[i])
      }
    }
  }
}
par(mai = margins) 


# return the selected features in decreasing order of selection probability
if( length(which(selectIndex>0))) {
ord = order(freq[nrow(freq),selectIndex], decreasing=T)
s = strsplit(selectNames," ")[[1]];
s = s[2:length(s)][ord];
topFreq = freq[,s];
}
else topFreq=NA  # No features satified labeling threshold so return NA

return(topFreq)
}

