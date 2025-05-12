


#'@param tollam2  1- tollam2 is the largested value for alpha

sparseSmoothFitCV <- function(dat, n.k, stepSize=0.1, lambda = NULL, nlam = 100, lam2 = NULL, nlam2 = 10, maxInt = 500,
                              epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                              accelrt = TRUE, nfolds = 5, mc.cores,hugeCont =100000,tollam2=0.01,
                              method = "ordinarySSP", w1=NULL, w2=NULL){
  
  dat$ID <- as.numeric(as.character(dat$ID))
  dat = dat[with(dat, order(ID, Position)),]
  
  
  meth <- dat$Meth_Counts
  total <- dat$Total_Counts
  unmeth <-total-meth
  
  numCovs = ncol(dat)-4
  
  colnames(dat)[-c(1:4)] = paste0("X", 1:numCovs)
  
  myp = (numCovs+1)*n.k
  lambda.min.ratio = ifelse(nrow(dat)<myp,0.01,0.0001)
 # if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
  #---------------------------------
  # The sequence of lambda2 : ulam2
  #--------------------------------
  if (is.null(lam2)) {
    lambda_max <- 1- tollam2
    ulam2 <- seq(0, lambda_max, length.out = nlam2)
    
  } else { # user provided lambda values
    user_lambda2 = TRUE
    if (any(lam2 < 0)) stop("lambdas should be non-negative")
    ulam2 = as.double((sort(lam2)))
    nlam2 = as.integer(length(lam2))
  }

  initOut = extractMatsBD(dat,n.k=n.k) # basisMat0 basisMat1 and sparse and smooth penalty matrix
  
  # Step 1: decompose of H------------
  # Obtain the decomposition for the penalty matrix H 
  # Different H (Linv) corresponds to different penalty functions
  
 
  if(method == "ordinarySSP"){
    Linv = lapply(as.list(ulam2), function(x){
      getSeqLinv(x, sparOmega=initOut$sparOmega, 
                 smoOmega1=initOut$smoOmega1)
    })
    # Step 2: designmat in tilda------------  
    bigDesignAll <- extractDesignMatAllalphs(numCovs,initOut$basisMat1,as.matrix(dat[,-c(1:4)]),initOut$basisMat0,n.k, nlam2, Linv) # store in a list
  }
  if(method == "adapSSP"){
    Linv = lapply(as.list(ulam2), function(x){
      getSeqLinv(x, sparOmega=initOut$sparOmega, 
                 smoOmega1=initOut$smoOmega1, method = method, w1=w1, w2=w2)
    })
    bigDesignAll <- 
    extractDesignMatAllalphsAdp( numCovs=numCovs, basisMat1=initOut$basisMat1,
                                 snpdat=as.matrix(dat[,-c(1:4)]), basisMat0=initOut$basisMat0,
                                 nk=n.k,nlam2=nlam2,Linv=Linv)
    
  }
  # Step 3: lambdaMax and matrix lamGrid------------   this is the same for both ordinarySSP and adapSSP
  
 if(!is.null(lambda)) { # user provided lambda values
    user_lambda = TRUE
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
 }else{
   lamGrid <- 
   vapply(seq_along(ulam2), function(ii){
     lambda_max <- lambdaMaxTilda(y=meth, x=total, bigDesign=bigDesignAll[[ii]],  numCovs=numCovs, basisMat0=initOut$basisMat0, n.k = n.k)
     exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
             length.out = nlam))
   }, rep(0.0,nlam))
  }

  lamGrid[1,] <-  lamGrid[1,] + hugeCont

 # >   all.equal(lamGrid, lamGridRaw)
 # [1] "Mean relative difference: 0.0001396023"
  #--------------
  # Step 1, CV fold
  #----------------

  nsamp = nrow(dat)
  temp = nsamp-(floor(nsamp/nfolds)*nfolds)
  
  if(temp==0){
    folds <- c(rep(1:nfolds,   floor(nsamp/nfolds)))
  }
  if(temp>0){
    folds <- c(rep(1:nfolds,   floor(nsamp/nfolds)), 1: (nsamp-(floor(nsamp/nfolds)*nfolds)) )
  }

  #folds <- c(rep(1:nfolds,   floor(nsamp/nfolds)), 1: (nsamp-(floor(nsamp/nfolds)*nfolds)) )
  
  foldIndex <- 
  lapply(1:nfolds, function(x){
    which(folds==x)
  })
  
  
  # leave off the first and last points from any fold, and always include them in the training sets in each iteration of CV
  
  exm <- c(min(dat$Position), max(dat$Position))
  exmid <- which(dat$Position %in% exm)
  
  #foldIndex <- lapply(foldIndex, function(x){ x[-which(x %in%exmid)]}) -- modify Sept 7
  foldIndex <- lapply(foldIndex, function(x){ 
    if(any(x %in%exmid)){
      x[-which(x %in% exmid)]
    }else{
      x
    }
  })

  
  AllOut = parallel::mclapply(seq(nfolds), function(ijk){
 
    testID <- foldIndex[[ijk]]
    designTrain <- lapply(bigDesignAll, function(xx){xx[-testID,]})
    designTest <-lapply(bigDesignAll, function(xx){xx[testID,]})
    
   # Step 1: perform the CV----
    trainFit =sparseSmoothGridRaw_useRPathCpp(meth=meth[-testID], total=total[-testID], unmeth=unmeth[-testID], n.k=n.k, nlam2=nlam2, lamGrid=lamGrid, 
                                              theta=initTheta, stepSize=stepSize, shrinkScale=shrinkScale,
                                               bigDesignList = designTrain , numCovs=numCovs,
                                                maxInt = maxInt,  epsilon = epsilon,  truncation = truncation)
   
    testPred = sparseSmoothPredPassGrid(trainFit =trainFit,meth =meth[testID], total=total[testID], unmeth=unmeth[testID], 
                                       bigDesign=designTest , nlam=nlam, nlam2=nlam2,n.k=n.k, numCovs=numCovs,truncation=truncation, nsamp=length(testID))

    
   list(testPred, trainFit)

  }, mc.cores=mc.cores)

  
  if(nlam2==1){
    testAlldev =  lapply(seq(nfolds), function(i){
      AllOut[[i]][[1]][1,]
    }) 
    
  }else{
  
  testAlldev =  lapply(seq(nfolds), function(i){
    AllOut[[i]][[1]][1,,]
  }) 
}
  

  # What to export from the cv function
  # best lambda2, best lambda1
  
  #1. calculate the cross validation average/SD matrix of testPred
  
  testPredMean = Reduce("+", testAlldev) / length(testAlldev)
  
  if(nlam2>1){
    testPredSD = apply(simplify2array(testAlldev), 1:2, sd )
  }else{
    testPredSD =apply(simplify2array(testAlldev), 1, sd )
  }
  
  
  
  #which(testPredMean == min(testPredMean), arr.ind = TRUE)
  
  if(nlam2 >1){
  # overall best
  bestInd = which(testPredMean == min(testPredMean), arr.ind = TRUE)
  
  bestLambda1 = lamGrid[bestInd]
  bestLambda2 = ulam2[bestInd[2]]
  
  # best lambda for different alpha values
  bestIndAllLam2 = apply(testPredMean, 2, which.min )
  bestLambda1vec = lamGrid[cbind(bestIndAllLam2,seq(nlam2))]
  
  #------- implement the lambda.1se as well --- May 7, 2021
  
  bestIndAllLam2_1SE=
  vapply(1:ncol(testPredMean), function(x){
    min(which(testPredMean[,x]<testPredMean[bestIndAllLam2[x],x] + testPredSD[bestIndAllLam2[x],x]/sqrt(nfolds)))
  }, 1.0)
  
  bestLambda1vec_1SE = lamGrid[cbind( bestIndAllLam2_1SE,seq(nlam2))]


  bestInd_1SE = bestIndAllLam2_1SE[bestInd[2]]
  bestLambda1_1SE = lamGrid[bestInd_1SE, bestInd[2]]
  
  }
  if(nlam2 == 1){
    bestInd = which(testPredMean == min(testPredMean), arr.ind = TRUE)
    bestLambda1 = lamGrid[bestInd,1]
    bestLambda2 = ulam2
    
    
    #------- implement the lambda.1se as well --- May 7, 2021
    
    
    # For different alpha values
    
    
    bestInd_1SE=  min(which(testPredMean<testPredMean[bestInd] + testPredSD[bestInd]/sqrt(nfolds)))
    bestLambda1_1SE = lamGrid[bestInd_1SE,1]

    
  }

  
 # if(bestmethod == "1se")
    if(nlam2 >1){
      
      bestFitAll <- parallel::mclapply(seq(ulam2), function(i){
        sparseSmoothPathRawOneLam1(meth=meth, total=total,unmeth=unmeth, initTheta = AllOut[[1]][[2]][[i]][,bestIndAllLam2_1SE[i]],
                                   intStepSize= stepSize, ulam= bestLambda1vec_1SE[i], bigDesign = bigDesignAll[[i]],
                                   n.k = n.k, numCovs=numCovs, maxInt=maxInt, epsilon = epsilon, shrinkScale=shrinkScale,truncation=truncation,
                                   Linv=Linv[[i]], datPosition = dat$Position, basisMat1 = initOut$basisMat1, basisMat0=initOut$basisMat0,
                                   method=method)
                                  
      }, mc.cores=mc.cores)
      
      
      bestFit <- bestFitAll[[bestInd[2]]]
      
      bestFitAllmin <- parallel::mclapply(seq(ulam2), function(i){
       
        sparseSmoothPathRawOneLam1(meth=meth, total=total,unmeth=unmeth, initTheta = AllOut[[1]][[2]][[i]][,bestIndAllLam2[i]],
                                   intStepSize= stepSize, ulam= bestLambda1vec[i], bigDesign = bigDesignAll[[i]],
                                   n.k = n.k, numCovs=numCovs, maxInt=maxInt, epsilon = epsilon, shrinkScale=shrinkScale,truncation=truncation,
                                   Linv=Linv[[i]], datPosition = dat$Position, basisMat1 = initOut$basisMat1, basisMat0=initOut$basisMat0,
                                   method=method)
      }, mc.cores=mc.cores)
      
      bestFitmin <- bestFitAllmin[[bestInd[2]]]
    }
    if(nlam2 == 1){
      bestFit <- 
      sparseSmoothPathRawOneLam1(meth=meth, total=total,unmeth=unmeth, 
                                 initTheta = AllOut[[1]][[2]][[1]][,bestInd_1SE[1]],
                                 intStepSize= stepSize, 
                                 ulam= bestLambda1_1SE, 
                                 bigDesign = bigDesignAll[[1]],
                                 n.k = n.k, numCovs=numCovs, maxInt=maxInt, 
                                 epsilon = epsilon, shrinkScale=shrinkScale,
                                 truncation=truncation,
                                 Linv=Linv[[1]], datPosition = dat$Position, basisMat1 = initOut$basisMat1, basisMat0=initOut$basisMat0,
                                 method=method)
      bestFitmin <- 
      sparseSmoothPathRawOneLam1(meth=meth, total=total,unmeth=unmeth, 
                                 initTheta = AllOut[[1]][[2]][[1]][,bestInd[1]],
                                 intStepSize= stepSize, 
                                 ulam= bestLambda1, 
                                 bigDesign = bigDesignAll[[1]],
                                 n.k = n.k, numCovs=numCovs, maxInt=maxInt, 
                                 epsilon = epsilon, shrinkScale=shrinkScale,
                                 truncation=truncation,
                                 Linv=Linv[[1]], datPosition = dat$Position, basisMat1 = initOut$basisMat1, basisMat0=initOut$basisMat0,
                                 method=method)
      
    }
    
  if(nlam2==1){
  return(out = list(bestFit = bestFit,testPredMean=testPredMean,testPredSD=testPredSD,
                    ulam2=ulam2, lamGrid=lamGrid, bestLambda1=bestLambda1, bestLambda2=bestLambda2,
                    bestInd=bestInd, testAlldev=testAlldev, bestLambda1_1SE=bestLambda1_1SE,
                    bestFitmin = bestFitmin,  sparOmega = initOut$sparOmega,smoOmega1 = initOut$smoOmega1))
  }else{
    return(out = list(bestFit = bestFit,testPredMean=testPredMean,testPredSD=testPredSD,
                      ulam2=ulam2, lamGrid=lamGrid, bestLambda1=bestLambda1, bestLambda2=bestLambda2,
                      bestInd=bestInd, bestFitAll=bestFitAll, testAlldev=testAlldev,
                      bestLambda1_1SE=bestLambda1_1SE,bestFitmin = bestFitmin, bestFitAllmin=bestFitAllmin,
                      sparOmega = initOut$sparOmega,smoOmega1 = initOut$smoOmega1)) 
    }
}


#plotSScv <- function(cvOut){
#  testLossOut[4,,,]
#}
  
#library(colortools)
#library(gplots)
#library(colorspace)
plotCV <- function(cvout, hugeCont =100000, mycols = c('#a6cee3','#1f78b4','#b2df8a',
                                                       '#33a02c','#fb9a99','#e31a1c','#fdbf6f',
                                                       '#ff7f00','#cab2d6','#6a3d9a',
                                                       "#E16A86" ,"#D37A4A", "#B88A00")){
  
  lamGridPlot = cvout$lamGrid
  lamGridPlot[1,] = lamGridPlot[1,]-hugeCont

  
  #mycols = qualitative_hcl(n = length(cvout$ulam2))
 # cols = sequential("steelblue", plot = FALSE)
  
  #mycols <- cols[1:length(cvout$ulam2)]
  plotPred = cvout$testPredMean
  
  if(length(cvout$ulam2)>1){
  
    
  #plotCI(x=log(lamGridPlot[,1]),  y=plotPred[,1], uiw= cvout$testPredSD[,1],
  #       xlim = log(c(min(lamGridPlot), max(lamGridPlot))),
  #       ylim = c(min( plotPred)-max(cvout$testPredSD), max( plotPred)+max(cvout$testPredSD)),
  #       pch = 19,xlab = expression(log(lambda)), 
  #       ylab = "mean of CV prediction error", col = mycols[1], cex = 0.5)
  plot(log(lamGridPlot[,1]),  plotPred[,1], xlim = log(c(min(lamGridPlot), max(lamGridPlot))),
       ylim = c(min( plotPred), max( plotPred)), pch = 19,xlab = expression(log(lambda)), 
       ylab = "mean of CV prediction error", col = mycols[1], cex = 0.5)
  if(length(cvout$ulam2)>1){
  for(i in 2:length(cvout$ulam2)){
    points(log(lamGridPlot[,i]),  plotPred[,i], col = mycols[i], pch = 19, cex = 0.5)
  }
  }
  }
  
  if(length(cvout$ulam2)==1){
    
    plot(log(lamGridPlot[,1]),  plotPred, xlim = log(c(min(lamGridPlot), max(lamGridPlot))),
         ylim = c(min( plotPred), max( plotPred)), pch = 19,xlab = expression(log(lambda)), 
         ylab = "mean of CV prediction error", col = mycols[1], cex = 0.5)
  }
  legend('topleft', legend = paste( "alpha =", cvout$ulam2), bty = "n", cex = 0.7, fill = mycols)
  
 
}
  
#plotCV(cvout)


