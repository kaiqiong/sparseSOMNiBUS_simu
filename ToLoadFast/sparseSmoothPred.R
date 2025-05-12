

sparseSmoothPredPassGrid <- function(trainFit,meth, total,  unmeth, bigDesign, nlam, nlam2,n.k,numCovs,truncation, nsamp){
  if(nlam2>1){
    lossvals <- 
      vapply(seq(nlam2), function(i){
        vapply(seq(nlam), function(j){
          
        binomObjectCppLossOnlyPassMat(theta=trainFit[[i]][,j], bigDesign[[i]], meth, total, unmeth, truncation=truncation, nsamp)%>% unlist()
          
          
         #mypi =   1/(1+exp(-bigDesign[[i]]%*%trainFit[[i]][,j]))
          
        # mean(-2*(meth*log(mypi)+unmeth*log(1-mypi)))
         
         #mypi_as <- asin(2*mypi-1)
         #obsmypi <- asin(2*meth/total-1)
         
#sqrt(mean((mypi_as-obsmypi)^2))
         
         # c(testOut$neg2loglik/nrow(testDat), testOut$sqrt_of_sum_of_dif/sqrt(nrow(testDat)))
        }, FUN.VALUE = rep(0.0, 2))
      }, FUN.VALUE = matrix(1, nrow = 2, ncol =  nlam[[1]]))
  }
  
  if(nlam2 ==1 ){
    
    lossvals= vapply(seq(nlam), function(j){
      
      binomObjectCppLossOnlyPassMat(theta=trainFit[[1]][,j], bigDesign[[1]], meth, total, unmeth, truncation=truncation, nsamp)%>% unlist()
      
    }, FUN.VALUE = rep(0.0, 2))
  }
  
  return(lossvals)
} 

#' @title Prediction of a sparse-smoothness fit
#' @param fit an object from \code{sparseSmoothFit} or \code{fitProxGrad}
#' @param dat test dataset for prediction
#' @param n.k number of knots for all covariates (including intercept);
#' curretnly, we assume the same n.k for all covariates
#' @param lengthUniqueDataID number of samples in the data
#' @param numCovs number of covariates
#' @param lambda1 penalization parameter for the L2 norm
#' @param lambda2 penalization parameter for the weight between two penalities
#' @return This function return a vector including objects:
#' \itemize{
#' \item \code{testLoss} var of alpha
#' \item \code{binomLoss} var of alpha0
#' \item \code{penTerms} var of alpha_p, p = 1, 2, P
#' }
#' @author Kaiqiong Zhao
#' @noRd
# one theta update from the prox(theta_old - t gradient(theta_old))
#' We can directly supply the basisMat0, designMat1 for the testDataset
sparseSmoothPredPassTest <- function(trainFit,testDat, basisMat0test, designMat1test, n.k, numCovs, truncation){
  
  
  meth = testDat$Meth_Counts
  total = testDat$Total_Counts
  unmeth = total - meth
  
  nlam2 = length(trainFit$thetaOut)
  
  
  nlam = lapply(trainFit$thetaOut, ncol)
  
  
  if(nlam2>1){
    lossvals <- 
      vapply(seq(nlam2), function(i){
        vapply(seq(nlam[[i]]), function(j){
          
          testOut <-binomObjectCppLossOnlyVec(theta=trainFit$thetaOutOri[[i]][,j], basisMat0=basisMat0test,
                                              numCovs=numCovs,designMat1=designMat1test, truncation=truncation, meth=meth, unmeth=unmeth,total=total, nk = n.k 
          )
          
          c(testOut$neg2loglik/nrow(testDat), testOut$sqrt_of_sum_of_dif/sqrt(nrow(testDat)))
        }, FUN.VALUE = rep(0.0, 2))
      }, FUN.VALUE = matrix(1, nrow = 2, ncol =  nlam[[1]]))
  }
  
  if(nlam2 ==1 ){
    
    lossvals= vapply(seq(nlam[[1]]), function(j){
      
      testOut <-binomObjectCppLossOnlyVec(theta=trainFit$thetaOutOri[[1]][,j], basisMat0=basisMat0test,
                                          numCovs=numCovs,designMat1=designMat1test, truncation=truncation,meth=meth, unmeth=unmeth,total=total,nk=n.k
      )
      c(testOut$neg2loglik/nrow(testDat), testOut$sqrt_of_sum_of_dif/sqrt(nrow(testDat)))
    }, FUN.VALUE = rep(0.0, 2))
  }
  
  return(lossvals)
  
}


#'The hidden requirement for the function 'sparseSmoothPredOneSetting' is that the order of components
#'in thetaEst should correspond its order in the fitDat as well as the testDat
#' i.e the order of columns matter
#' the orders of rows of fitDat does not matter. 
sparseSmoothPredOneSetting <- function(thetaEst, fitDat, testDat, n.k, numCovs, truncation){
  initOut = extractMats(fitDat,n.k=n.k)
  
  basisMat0 = initOut$basisMat0
  basisMat1 = initOut$basisMat1
  
  trainDatPos = fitDat$Position
  
  
  #sparseSmoothPred <- function(trainFit, trainDatPos, testDat, basisMat0, basisMat1, n.k, numCovs, truncation){
  
  rowsID <- match(testDat$Position, trainDatPos)
  
  
  if (any(is.na(rowsID))){message(paste0("Some positions in the test dataset are not present in the train fit;
                                         predictions at those postions are not available"))}
  # calculate the design matrix for the test dataset
  basisMat0 <- basisMat0[rowsID,]
  basisMat1 <- basisMat1[rowsID,]
  designMat1 <- extractDesignMat1(numCovs, basisMat1, testDat)
  
  # calculate criterion
  # 1, loss function value / number of observations
  # 2, mean prediction errors
  
  meth = testDat$Meth_Counts
  total = testDat$Total_Counts
  unmeth = total - meth
  
  
  testOut <-binomObjectCppLossOnlyVecExpOutcome(theta=thetaEst, basisMat0=basisMat0,
                                                numCovs=numCovs,designMat1=designMat1, truncation=truncation,
                                                meth=meth, unmeth=unmeth,total=total,nk=n.k)
  
  trueOut <- asin(2* (meth/total)-1)
  preOut <- asin(2* (testOut$pi_ij_est)-1)
  
  
  c(deviance = testOut$neg2loglik/nrow(testDat), 
                     RMSE = sqrt(mean((trueOut-preOut)^2)),
                     corRaw = cor(meth/total, testOut$pi_ij_est),
                     corTran = cor(trueOut, preOut))
  #list(predError=c(deviance = testOut$neg2loglik/nrow(testDat), 
  #                 RMSE = sqrt(mean((trueOut-preOut)^2)),
  #                 corRaw = cor(meth/total, testOut$pi_ij_est),
  #                 corTran = cor(trueOut, preOut)),
  #     pred=testOut$pi_ij_est, 
  #     meth = meth, total=total)
  
}




