sparseSmoothGridRaw_useRPathCpp <- function(meth, total, unmeth, n.k, nlam2, lamGrid, theta, stepSize, shrinkScale,
                                            bigDesignList,  numCovs,
                                            maxInt = 10^5,  epsilon = 1E-20,truncation = TRUE){
  
  AllOut <- parallel::mclapply(seq(nlam2), function(i){
    sparseSmoothPathCpp(ulam=lamGrid[,i], theta, stepSize, meth, unmeth, total,
                                             bigDesignList[[i]], n.k, maxInt, epsilon, shrinkScale, numCovs, truncation)
  }, mc.cores = 1)
  return(AllOut)
}



# Given the one lambda1 

sparseSmoothPathRawOneLam1 <- function(meth, total, unmeth, initTheta, intStepSize, ulam, bigDesign, 
                                       n.k, numCovs, maxInt, epsilon, shrinkScale, truncation, Linv, datPosition, basisMat1, basisMat0,
                                       method=method){
  
  thetaMat <- fitProxGradCppClean1NewArmaThetaOnly(initTheta, intStepSize, ulam[1], meth, total,unmeth, bigDesign, n.k, numCovs,
                                       maxInt, epsilon, shrinkScale, truncation);
  
  
  thetaSep <- getSeparateThetaCpp(thetaMat, n.k, numCovs)
  
  zeroCovsBool<-unlist(lapply( thetaSep[-1], function(x){all(x==0)}))
  
  #fit1$thetaEstSep
  #checkall[,1] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[1], Hp, L, Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
  
  if(method == "ordinarySSP"){
    thetaMatOriginalSep <-(lapply( thetaSep, function(x){Linv%*%x}))
  }
  if(method =="adapSSP"){
    thetaMatOriginalSep <- thetaTildaToOriginal(numCovs, Linv, thetaSep)
  }
  

  thetaMatOriginal <- unlist( thetaMatOriginalSep)
  uniqPos = unique(datPosition)
  
  uni_rows <- match(uniqPos, datPosition)
  
  #penalBetas <- basisMat1[uni_rows,]%*% thetaMatOriginal
  
  #thetaMatOriginalSep <- getSeparateThetaCpp(thetaMatOriginal, n.k, numCovs)
  
  penalBetas <-  lapply(thetaMatOriginalSep[-1][!zeroCovsBool], function(x){basisMat1[uni_rows,]%*% x})
  
  penalBetas0 <-basisMat0[uni_rows,]%*%thetaMatOriginalSep[[1]]
  
  return(out = list(thetaMat=Matrix::Matrix(thetaMat, sparse=TRUE), ulam= ulam, 
                    thetaMatOriginal=Matrix::Matrix(thetaMatOriginal, sparse=TRUE),
                    zeroCovsBool=zeroCovsBool, 
                    penalBetas = penalBetas, uniqPos = uniqPos, penalBetas0=penalBetas0))
  
}

