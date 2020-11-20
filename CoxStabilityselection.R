"CoxStabilityselection" <- function(Z,t,delta,nbootstrap=100,nsteps=15,lambda=2,plotme=FALSE, alpha=0.8, Pw=0)
{
  # Stability selection in the spirit of Meinshausen&Buhlman using penalized package for Cox model regression
  # Input: 
  #     dFrame -  nxp data frame (asumes data standardized)
  #     t      -  nx1 numeric vector of event times
  #     delta  -  nx1 numeric vector for censoring indicator (delta=0 -> censored)
  #     lambda -  minimum value of lambda
  #     alpha  -  Weight scale factor (weights set to alpha on success in Bernoulli trial)
  #     Pw     -  bernoulli probability for variable reweighting (default = 0 for no-reweighting)
  #
  # result is a list containing: 
  #     $freq          - nsteps x p selection frequency matrix 
  #     $lambda        - typical vector of penalty values for L1 penalty
  #     $selectionProb - vector containing the probability that each feature is selected 
  #                      during the first nsteps steps of the Lasso path when half of the 
  #                      samples are used and the features are reweigthed by a random 
  #                      weight uniformly sampled in [alpha,1]. This probability is 
  #                      estimated by nbootstrap bootstrap samples
  #
  #
  require(penalized)
  require(survival)

  dimZ <- dim(Z) 
	n <- dimZ[1]
	p <- dimZ[2]
	halfsize <- as.integer(n/2)
	freq <- as.data.frame(matrix(0,nsteps,p));
  targetNames <-gsub(" ", ".", colnames(Z))
  colnames(freq) = targetNames
  check = 0;
	for (i in seq(nbootstrap)) {
				
  		# Ramdomly split the sample in two disjoint sets
  		perm <- sample(dimZ[1])
  		i1 <- perm[1:halfsize]
  		i2 <- perm[(halfsize+1):n]
      
  		# Generate random perturbation: weights given by outcome of bernoulli trial with sucess prob Pw. On success weight=alpha, else weight=1 
   		perturb = rep(1,length=p);
  		draw = runif(p,0,1)
      perturb[draw<Pw] = alpha;
      
  		# Randomly reweight (perturb) each variable with bernoulli weights
  	  Zs <- as.data.frame(t(t(Z)*perturb));
  		#browser();
  		
  		# run the penalized lasso fit for each sub-sample and record which variables are selected at each lambda
      pen <- penalized(Surv(t[i1], delta[i1]), penalized = Zs[i1,],data=Zs[i1,],steps=nsteps, lambda2=0, lambda1=lambda,standardize=F)
  		for (j in 2:length(pen)) {
  		  theNames = names(coefficients(pen[[j]]));
  		  if (length(unique(theNames)) < length(theNames) ) {
  		    browser();
  		  }
         freq[j,theNames] = freq[j,theNames]+1; 
   		}
      check = check+1;

      #rm(pen)
      #gc(verbose=TRUE)
      
      pen <- penalized(Surv(t[i2], delta[i2]), penalized = Zs[i2,],data=Zs[i2,],steps=nsteps, lambda2=0, lambda1=lambda,standardize=F)
     	for (j in 2:length(pen)) {
         theNames = names(coefficients(pen[[j]]));
         if (length(unique(theNames)) < length(theNames) ) {
           browser();
         }
         #print(theNames)
         freq[j,theNames] = freq[j,theNames]+1; 
  		}
      check = check+1;

      #Print lambda values for the path  ( pathes for each subsample are different )
      #lam=unlist(lapply(pen,penalty))[seq(from=1,to=2*length(pen),by=2)];
      #print(lam)
    
      #browser()
 		}
		
	# normalize frequence in [0,1]
	freq <- freq/(2*nbootstrap)
  print(check)
  check = check/(2*nbootstrap)
  print(check)

  lam  = unlist(lapply(pen,penalty))[seq(from=1,to=2*length(pen),by=2)];
  
	if (plotme) {
		matplot(lam, freq,type='l',xlab="Lambda Step",ylab="Frequency")    
	}
	
	# the final stability score is the maximum frequency over the steps
 result <- list(freq=freq, lambda=lam, selectionProb=apply(freq,2,max), numSamples=2*nbootstrap, nsteps = nsteps, alpha=alpha, Pw = Pw)
  
}

#
#
#
"CoxStabilityPlot" <- function(object, targetNames=colnames(object$freq), labelsize = 0.6, pLabelThresh=0.3, ymax=1,xrng=c(0,1),...)
{
    freq   = object$freq;
    lambda = object$lambda;
    remove = apply(freq, 2, function(x) all(x==0)) 
    if (all(remove)) 
        stop("all frequencies are zero for all values of lambda in this object")
    label <- "lambda1/lambda_max"
 
    labwidth <- ifelse(labelsize > 0, max(strwidth(colnames(freq[,!remove]), "inches", labelsize)), 0)
    margins <- par("mai")
    par(mai = c(margins[1:3], max(margins[4], labwidth * 1.4)))
    
    ratio = rev(lambda)/max(lambda)
    pts = ratio >= xrng[1] & ratio <= xrng[2];
    matplot(ratio[pts], freq[pts,!remove,  drop = FALSE], type = "l", 
    ylab = "Prob[ Select ] ", xlab = label, col = rainbow(sum(!remove)), 
    xlim = rev(range(lambda/max(lambda))), ylim=c(0,ymax), main=sprintf("Penalized Cox Model Stability Paths\n (%d permutations, alpha=%.2f, Pw=%.2f)",object$numSamples, object$alpha, object$Pw),...)
    grid()
    if (labelsize > 0 && !is.null(colnames(freq))) {
        # Identify features which selection probabilities that exceed labeling threshold
        take <- which(!remove )
        selectNames = ""; 
        selectIndex = vector(length=length(take))
        fidx = nrow(freq);
        
        for (i in 1:length(take)) {
            j <- take[i]
            if (max(freq[,j]) >= pLabelThresh) {
              selectNames = paste(selectNames, colnames(freq)[j]);
              selectIndex[i] = j;
              axis(4, at = freq[fidx,j], labels = targetNames[j], 
                  las = 1, cex.axis = labelsize, col.axis = rainbow(sum(!remove))[i], 
                  lty = (i - 1)%%5 + 1, col = rainbow(sum(!remove))[i])
            }
       }
     }
     par(mai = margins)
     s = strsplit(selectNames," ")[[1]];
     subFrame = freq[,s[2:length(s)]];
     colnames(subFrame) = targetNames[selectIndex];
    return(subFrame)
}
