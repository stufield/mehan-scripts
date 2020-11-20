"LinearglmnetStabilitySelectionTest1" <- function(n=100, p=1000, sigma=.1)
{  
  n1 = floor(n/2)
  n2 = n-n1-1;
  X = mvrnorm(n, mu=rep(0,p), diag(rep(1,p)))
  y = X[,1]
  y =     X[,3]*.666 + X[,33]*.33 + 0.5*rnorm(n,0,sigma);
  y = y + X[,11]*.33 + X[,517]*.666 + 0.5*rnorm(n,0,sigma);
  stabSel = LinearStabilitySelection(data.frame(X), y, nbootstrap=100)
  LinearStabilityPlot(stabSel, labelsize=.8, ymax=1, pLabelThresh=0.3)
  
}
 
"LinearglmnetStabilitySelection" <- function(dFrame,y,alpha=0.85,Pw=0,nbootstrap=100,nsteps=-1, plotme=TRUE, lamMax=-1, logScale=FALSE)
{
	# Stability selection in the spirit of Meinshausen&Buhlman using glmnet for Linear Reg
  # Input: 
	#     dFrame -  nxp data frame (asumes data standardized)
  #     y      -  nx1 response vector
  #     alpha  -  Weight scale factor (weights set to alpha on success in Bernoulli trial)
  #     Pw     -  bernoulli probability for variable reweighting (default = 0 for no-reweighting)
  #
  # result is a list containing: 
  #     $freq          - nsteps x p selection frequency matrix sw
  #     $lambda        - vector of penalty values for L1 penalty
  #     $selectionProb - vector containing the probability that each feature is selected 
  #                      during the first nsteps steps of the Lasso path when half of the 
  #                      samples are used and the features are reweigthed by a random 
  #                      weight uniformaly sampled in [alpha,1]. This probability is 
  #                      estimated by nbootstrap bootstrap samples
  #
	require(glmnet)
 	n <- dim(dFrame)[1]
	p <- dim(dFrame)[2]
	halfsize <- as.integer(n/2)
  targetNames = colnames(dFrame) 
	numSel = 0;
  
  # First fit all of the data with the requested number of steps to generate a sequence of 
  # lambda values to use for subsequent fits
  if ( nsteps > 0) {
     glmnet.fit = glmnet(as.matrix(dFrame),y, family="gaussian",standardize=F, nlambda=nsteps)
  }
  else {
    glmnet.fit = glmnet(as.matrix(dFrame),y, family="gaussian",standardize=F)
  }
  lam=glmnet.fit$lambda;    # Fix lambda steps for all subsequent regularizations

  nsteps = length(lam);
  freq <- as.data.frame(matrix(0,nsteps,p));
  colnames(freq) = colnames(dFrame)
  
  lambdaMat = matrix(0,nbootstrap,nsteps); # Debuggin
  X = as.matrix(dFrame);
	for (i in seq(nbootstrap)) {
				
  		# Ramdomly splits the sample in two sets
  		perm <- sample(n)
  		i1 <- perm[1:halfsize]
  		i2 <- perm[(halfsize+1):n]
      		
      # Generate random perturbation: weights given by outcome of bernoulli trial with sucess prob Pw. On success weight=alpha, else weight=1 
  		perturb = rep(1,length=p);
  		draw = runif(p,0,1)
  		perturb[draw<Pw] = alpha;
  		
  		# Randomly reweight (perturb) each variable with bernoulli weights
  		Xs <- as.matrix(t(t(X)*perturb));
  		
      glmnet.fit = glmnet(Xs[i1,],y[i1], lambda = lam, standardize=F, family="gaussian")
      #lambdaMat[i,] = glmnet.fit$lambda / glmnet.fit$lambda[1]; # Debuggin'...Store the regulaization path so we can check that they are the same...
    
      # Find the the variable that were selected at some point on the path as those with non-zero coefficients
  		tmp = t(as.matrix(glmnet.fit$beta))    # nsteps x p
  		sel.index = (1:dim(Xs)[2])[apply(tmp,2,sum)>0]  # Column indices for selected variables
  		for (j in 1:length(sel.index)) {
  		  lidx = (1:nsteps)[abs(tmp[1:nsteps,sel.index[j]])>0];
        print(lidx)
        print(sel.index[j])
  		  if (length(lidx) > 0) {
  		    freq[lidx,sel.index[j]] = freq[lidx,sel.index[j]]+1;
  		  }
  		}
      numSel = numSel+1;
      glmnet.fit = glmnet(Xs[i2,],y[i2], lambda = lam, standardize=F, family="gaussian")

      # Find the the variables that were selected at some point on the path
      # and increment their count at all lambda values where they are included
      # in the model
  		tmp = t(as.matrix(glmnet.fit$beta))    # nsteps x p
  		sel.index = (1:dim(Xs)[2])[apply(tmp,2,sum)>0]  # Column indices for selected variables
  		for (j in 1:length(sel.index)) {
  		  lidx = (1:nsteps)[abs(tmp[1:nsteps,sel.index[j]])>0];
        if (length(lidx) > 0) {
          freq[lidx,sel.index[j]] = freq[lidx,sel.index[j]]+1;
        }
  		  
  		}
  		numSel = numSel+1;


 		}
		
  # Compute Frequency
	freq <- freq/(2*nbootstrap)
  print(numSel)
  numSel = numSel/(2*nbootstrap)
  print(numSel)
  
	if (plotme) {
		matplot(freq,type='l',xlab="Lambda Step",ylab="Pr[select]")    
	}
	

 result <- list(freq=freq, lambda=lam, selectionProb=apply(freq,2,max), numSamples=2*nbootstrap, nsteps = nsteps, alpha=alpha, Pw = Pw)
	
}

#
#
#
"LinearglmnetStabilityPlot" <- function(object, targetNames=colnames(object$freq), labelsize = 0.6, pLabelThresh=0.5, ymax=1, normalize=FALSE,...)
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
            xlim = range(lambda/max(lambda)), ylim=c(0,ceiling(10*ymax)/10+1e-5), main="Normalized Linear Regression Stability Paths",...)
      
    }
    else {
      matplot((lambda)/max(lambda), freq[,!remove,  drop = FALSE], type = "l", axes=F,
          ylab="Prob [Select]",xlab = label, col = rainbow(sum(!remove)), 
          xlim = rev(range(lambda/max(lambda))), ylim=c(0,ceiling(10*ymax)/10+1e-5), main="Linear Regression Stability Paths",...)
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
    ord = order(freq[nrow(freq),selectIndex], decreasing=T)
    s = strsplit(selectNames," ")[[1]];
    s = s[2:length(s)][ord]
    return(freq[,s])
}
