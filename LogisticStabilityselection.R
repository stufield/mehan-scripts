"LogisticStabilitySelectionTest1" <- function()
{ 
  #X = mvrnorm(100, mu=c(10,20,30), diag(c(1,1,1)))
  #stabSel = LogisticStabilitySelection(X, class, nbootstrap=1000, nsteps=50, maxCoeff=200)
  #LogisticStabilityPlot(stabSel, labelsize=.8, ymax=1, pLabelThresh=0.3)
  
  test1 = loadMatlabData("stabSelTest1-LDA.mat");
  stabSel = LogisticStabilitySelection(test1$Zframe, test1$class, nbootstrap=50, nsteps=50, maxCoeff=200, fam="gaussian")
  LogisticStabilityPlot(stabSel, labelsize=.8, pLabelThresh=0.3)
}
 
"LogisticStabilitySelection" <- function(dFrame,class,alpha=0.85,nbootstrap=100,nsteps=20, maxCoeff=25,plotme=TRUE, lamMax=-1, fam="binomial")
{
	# Stability selection in the spirit of Meinshausen&Buhlman using glmnet for logistic Reg
  # Input: 
	#     dFrame -  nxp data frame (asumes data standardized)
  #     y      -  nx1 numeric class label
	# result is a list containing: 
  #     $freq          - nsteps x p selection frequency matrix sw
  #     $lambda        - vector of penalty values for L1 penalty
  #     $selectionProb - vector containing the probability that each feature is selected 
  #                      during the first nsteps steps of the Lasso path when half of the 
  #                      samples are used and the features are reweigthed by a random 
  #                      weight uniformaly sampled in [alpha,1]. This probability is 
  #                      estimated by nbootstrap bootstrap samples
	require(glmnet)
 	n <- dim(dFrame)[1]
	p <- dim(dFrame)[2]
	halfsize <- as.integer(n/2)
  targetNames = colnames(dFrame) 
	numSel = 0;
  
  # First fit all of the data with the requested number of steps to generate a sequence of 
  # lambda values to use for subsequent fits
  gfit = glmnet(as.matrix(dFrame),class, family=fam, alpha=1, nlambda=nsteps, standardize=TRUE)  # alpha=1 --> lasso
  if ( lamMax == -1) {
    lam=gfit$lambda;  
  }
  else{
   lam=seq(lamMax,0.0001,-1/nsteps)
  }  

  nsteps = length(lam);
  freq <- as.data.frame(matrix(0,nsteps,p));
  colnames(freq) = colnames(dFrame)

  lambdaMax = vector(); 
  X = as.matrix(dFrame);
	for (i in seq(nbootstrap)) {
				
  		# Ramdomly splits the sample in two setsC
  		perm <- sample(n)
  		i1 <- perm[1:halfsize]
  		i2 <- perm[(halfsize+1):n]
      
      # Random perturbation to the predictors (multiply each column by unif[alpha,1] r.v.)
      perturb = runif(p,alpha,1)
    
      gfit = glmnet(X[i1,],class[i1], family=fam, alpha=1, pmax = maxCoeff, lambda=lam,standardize=TRUE, penalty.factor=perturb)  # alpha=1 --> lasso
      #gfit = glmnet(X[i1,],class[i1], family="binomial", alpha=1, pmax = maxCoeff, nlambda=nsteps,standardize=TRUE, penalty.factor=perturb)  # alpha=1 --> lasso
      lambdaMax = c(lambdaMax, gfit$lambda[1]);
    
      # Find the the variable that were selected at some point on the path
      tmp = as.matrix(gfit$beta)     
      sel.index = (1:length(targetNames))[abs(tmp[,dim(tmp)[2]])>0]
            #targetNames[sel.index]
     for (j in 1:length(sel.index)) {
         lidx = (1:nsteps)[abs(tmp[sel.index[j],])>0];
         freq[lidx,targetNames[sel.index[j]]] = freq[lidx,targetNames[sel.index[j]]]+1;
      }
      numSel = numSel+1;

      gfit = glmnet(X[i2,],class[i2], family=fam, pmax=maxCoeff, alpha=1, lambda=lam, standardize=TRUE, penalty.factor=perturb)  # alpha=1 --> lasso
      #gfit = glmnet(X[i2,],class[i2], family="binomial", pmax=maxCoeff, alpha=1, nlambda=nsteps, standardize=TRUE, penalty.factor=perturb)  # alpha=1 --> lasso
      lambdaMax = c(lambdaMax, gfit$lambda[1]);

      # Find the the variables that were selected at some point on the path
      # and increment their count at all lambda values where the are included
      # in the model
      tmp = as.matrix(gfit$beta)     
      sel.index = (1:length(targetNames))[abs(tmp[,dim(tmp)[2]])>0]
      for (j in 1:length(sel.index)) {
         lidx = (1:nsteps)[abs(tmp[sel.index[j],])>0];
         freq[lidx,targetNames[sel.index[j]]] = freq[lidx,targetNames[sel.index[j]]]+1;
      }      
      numSel = numSel+1;

 		}
		
	# normalize frequence in [0,1]
  #browser();
      
	freq <- freq/(2*nbootstrap)
  print(numSel)
  numSel = numSel/(2*nbootstrap)
  print(numSel)
  
	if (plotme) {
		matplot(freq,type='l',xlab="Lambda Step",ylab="Pr[select]")    
	}
	
	# the final stability score is the maximum frequency over the steps
	#result <- apply(freq,2,max)
 print(max(lambdaMax))
result <- list(freq=freq, lambda=lam, selectionProb=apply(freq,2,max), numSamples=2*nbootstrap, lambdaMax = lambdaMax)
}

#
#
#
"LogisticStabilityPlot" <- function(object, targetNames=colnames(object$freq), labelsize = 0.6, pLabelThresh=0.3, ymax=1, normalize=FALSE,...)
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
            xlim = range(lambda/max(lambda)), ylim=c(0,ceiling(10*ymax)/10+1e-5), main="Normalized Logistic Regression Stability Paths",...)
      
    }
    else {
      matplot((lambda)/max(lambda), freq[,!remove,  drop = FALSE], type = "l", axes=F,
          ylab="Prob [Select]",xlab = label, col = rainbow(sum(!remove)), 
          xlim = rev(range(lambda/max(lambda))), ylim=c(0,ceiling(10*ymax)/10+1e-5), main="Logistic Regression Stability Paths",...)
    }
    axis(1,at=seq(0,1,0.2));
    axis(2,at=seq(0,1,0.1))
    #grid(NA,ny=ceiling(ymax*10))
    abline(h=seq(0.1,floor(ymax*10)/10,0.1), lty=2, col=8)

    if (labelsize > 0 && !is.null(colnames(freq))) {
        # Identify features which selection probabilities that exceed labeling threshold
        take <- which(!remove )
        selectNames = ""; 
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
    ord = order(freq[nrow(freq),selectIndex], decreasing=T)
    s = strsplit(selectNames," ")[[1]];
    s = s[2:length(s)][ord]
    subFrame = freq[,s];
    colnames(subFrame) = targetNames[selectIndex][ord];
    return(subFrame)
}
