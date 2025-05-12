

source("../../ToLoad/sparseSmoothPred.R")
library(Rcpp)
sourceCpp("../../ToLoad/sparseOmegaCr.cpp")
sourceCpp("../../ToLoad/fitProxGradCpp.cpp")
source("../../ToLoad/utils.R")
source("../../ToLoad/Simu.R")
load("../../ToLoad/5betaShapes.RData")

#load("../ToLoad/5betaShapes.RData")
#load("../../ToLoad/4Step_1_curve.RData")


#load("../../ToLoad/5step_funsShapes.RData")
my.samp=50
n.snp=1000
n.sig.snp=5
#--------------
M = 100
n.k = 10
numCovs = n.snp
truncation = TRUE


lam2_vec <- c(seq(0, 0.9, 0.1), 0.95, 0.99)
nlam2 = length(lam2_vec)

#-------------------
# Estimation results
#--------------------

PosAll <- matrix(0, nrow = length(pos), ncol = M)

BetaAll <- array(NA, c(length(pos), M, 100), dimnames = list(NULL, NULL, paste0("beta", 1:100)))

beta_all_alp0 <- BetaAll
beta0best <- beta0alp0 <- matrix(NA, nrow = length(pos),ncol = M)

#-------------------
# Prediction results
#--------------------
trueDev <- rep(NA, M)

pred <- matrix(NA, nrow = M, ncol = 4)
#colnames(pred) <- c("deviance", "rmse_t")

pred_alp0 <- pred
pred_all_alp <- array(NA, c(M, 4, nlam2), dimnames= list(NULL, c("deviance", "tran_rmse", "corRaw", "corTrans"),
                                                         paste0("alp", lam2_vec)))
trueDev_c <- rep(NA, M)

bestLams <- matrix(NA, nrow = M, ncol = 2)
colnames(bestLams) <- c("lam1", "lam2")
#------------------------
# Selection results
#-------------------------


SelectResBest <- matrix(NA, nrow= 100, ncol = M)

SlectResBestAlp0 <- SelectResBest


for( i in 1:M){
  
  try({
    load(paste0("Res", i, ".RData"))
    source("../../ToLoad/Simu.R")
    library(Rcpp)
    sourceCpp("../../ToLoad/sparseOmegaCr.cpp")
    sourceCpp("../../ToLoad/fitProxGradCpp.cpp")
    #------------------
    # estimation
    #--------------------
    PosAll[,i] = cvout$bestFit$uniqPos
    
    selCovs =  which(!cvout$bestFit$zeroCovsBool)
    
    for(jj in selCovs[selCovs<=n.snp]){
      BetaAll[,i,jj] = cvout$bestFit$penalBetas[[match(jj, selCovs)]]
    }
    BetaAll[,i,which(cvout$bestFit$zeroCovsBool)] <- 0
    
    
    selCovs =  which(!cvout$bestFitAll[[1]]$zeroCovsBool)
    for(jj in selCovs[selCovs<=n.snp]){
      beta_all_alp0[,i,jj] = cvout$bestFitAll[[1]]$penalBetas[[match(jj, selCovs)]]
    }
    beta_all_alp0[,i,which(cvout$bestFitAll[[1]]$zeroCovsBool)] <- 0
    
    
    
    
    
    # estimating beta.0
    
    set.seed(32314242+i)
    
    dat = sparseSimu(my.samp, n.snp, n.sig.snp, beta.0, BETAs)$dat
    
    
    dat$ID <- as.numeric(as.character(dat$ID))
    dat = dat[with(dat, order(ID, Position)),]
    initOut = extractMats(dat,n.k=n.k)
  
    basisMat0 = initOut$basisMat0
    uniqPos = unique(dat$Position)
    
    #all.equal(uniqPos, cvout$bestFit$uniqPos)
    
    uni_rows <- match(uniqPos, dat$Position)
    
    #penalBetas <- basisMat1[uni_rows,]%*% thetaMatOriginal
    
    thetaMatOriginalSep <- getSeparateThetaCpp(cvout$bestFit$thetaMatOriginal[,1], n.k, numCovs)

    beta0best[,i]<-basisMat0[uni_rows,]%*%thetaMatOriginalSep[[1]]
    
    thetaMatOriginalSep <- getSeparateThetaCpp(cvout$bestFitAll[[1]]$thetaMatOriginal[,1], n.k, numCovs)
    
    beta0alp0[,i]<-basisMat0[uni_rows,]%*%thetaMatOriginalSep[[1]]
    
    
 
    #-----------
    # selection
    #--------------
    
    SelectResBest[,i]  <- cvout$bestFit$zeroCovsBool
    
    SlectResBestAlp0[,i] <- cvout$bestFitAll[[1]]$zeroCovsBool
    
    bestLams[i,] <- c(cvout$bestLambda1, cvout$bestLambda2)
    
    #-------------------------
    # Out-of-sample Prediction
    #---------------------------
    
    set.seed(902314242+i)
    out = sparseSimu(my.samp, n.snp, n.sig.snp, beta.0, BETAs)
    testDat = out$dat
    truePi = out$truePi
    
    trueDev_c[i]  = -2*mean(testDat$Meth_Counts*log(truePi)+ (testDat$Total_Counts-testDat$Meth_Counts)*log(1-truePi))
    
    
    
    pred_alp0 [i,] = sparseSmoothPredOneSetting(thetaEst=cvout$bestFitAll[[1]]$thetaMatOriginal[,1], 
                                                fitDat=dat, 
                                                testDat=testDat, 
                                                n.k=n.k, numCovs=numCovs, truncation=truncation)
    pred[i,] = sparseSmoothPredOneSetting(thetaEst=cvout$bestFit$thetaMatOriginal[,1], 
                                          fitDat=dat, 
                                          testDat=testDat, 
                                          n.k=n.k, numCovs=numCovs, truncation=truncation)
    
    
    
    pred_all_alp[i,,] <- t(vapply(seq(nlam2), function(ii){
      sparseSmoothPredOneSetting(thetaEst=cvout$bestFitAll[[ii]]$thetaMatOriginal[,1], 
                                 fitDat=dat, 
                                 testDat=testDat, 
                                 n.k=n.k, numCovs=numCovs, truncation=truncation)
    }, rep(0, 4)))
  
    #trueDev_c[i] <- trueDev
  })
}



save.image(file="ResAll.RData")



# before calculating IMSE, we need to properly order the array BetaAll


#---------------
# for all alphas
#--------------

IMSE <- matrix(NA, nrow = M, ncol = n.snp)


BetaAllposOrd = BetaAll  # The rownames of BetaAllposOrd corresponds exactly to the the saved object pos

for(mm in 1:M){
  
  BetaAllposOrd[,mm,] <- BetaAll[match(pos,PosAll[,mm]),mm,]
  
}

for(bb in 1:n.sig.snp){
  
  see = (BetaAllposOrd[  ,,bb]-BETAs[,bb])^2 # substract truth for each columns
  
  IMSE[,bb]= apply(see, 2, sum)
}


for(bb in (n.sig.snp+1):n.snp){
  
  see = (BetaAllposOrd[,,bb]-0)^2 # substract truth for each columns
  
  IMSE[,bb]= apply(see, 2, sum)
}


# for the intercept

beta0bestOrd <- beta0best


for(mm in 1:M){
  beta0bestOrd[,mm] <- beta0best[match(pos,PosAll[,mm]),mm]
  
}
IMSEbeta0 <- rep(NA, M)
for( i in 1:ncol(beta0bestOrd)){
  IMSEbeta0[i] <- sum((beta0bestOrd[,i]-beta.0)^2)
}



save.image(file="ResAll.RData")



#---------------
# for alpha 0 only
#--------------

IMSE0 <- matrix(NA, nrow = M, ncol = n.snp)


BetaAllposOrd = beta_all_alp0  # The rownames of BetaAllposOrd corresponds exactly to the the saved object pos

for(mm in 1:M){
  
  BetaAllposOrd[,mm,] <- beta_all_alp0[match(pos,PosAll[,mm]),mm,]
  
}

for(bb in 1:n.sig.snp){
  
  see = (BetaAllposOrd[  ,,bb]-BETAs[,bb])^2 # substract truth for each columns
  
  IMSE0[,bb]= apply(see, 2, sum)
}


for(bb in (n.sig.snp+1):n.snp){
  
  see = (BetaAllposOrd[,,bb]-0)^2 # substract truth for each columns
  
  IMSE0[,bb]= apply(see, 2, sum)
}


# for the intercept


beta0bestOrd <- beta0alp0

for(mm in 1:M){
  beta0bestOrd[,mm] <- beta0alp0[match(pos,PosAll[,mm]),mm]
  
}
IMSEalp0 <- rep(NA, M)
for( i in 1:ncol(beta0bestOrd)){
  IMSEalp0[i] <- sum((beta0bestOrd[,i]-beta.0)^2)
}



save.image(file="ResAll.RData")



#---------
# correct some mistakes



beta0bestOrd <- beta0alp0

for(mm in 1:M){
  beta0bestOrd[,mm] <- beta0alp0[match(pos,PosAll[,mm]),mm]
  
}
IMSEalp0 <- rep(NA, M)
for( i in 1:ncol(beta0bestOrd)){
  IMSEalp0[i] <- sum((beta0bestOrd[,i]-beta.0)^2)
}


#---------------------------------------------
# calculate integrated squared bias & itegrated variance
#------------------------------------

ibias2 <- rep(NA, n.snp)

BetaAllposOrd = BetaAll  # The rownames of BetaAllposOrd corresponds exactly to the the saved object pos

for(mm in 1:M){
  BetaAllposOrd[,mm,] <- BetaAll[match(pos,PosAll[,mm]),mm,]
}

for(bb in 1:n.sig.snp){
  ibias2[bb] <- sum((apply(BetaAllposOrd[  ,,bb], 1, mean, na.rm = T) - BETAs[,bb])^2)
}
for(bb in (n.sig.snp+1):n.snp){
  ibias2[bb] <- sum((apply(BetaAllposOrd[  ,,bb], 1, mean, na.rm = T))^2)
}

# for the intercept
beta0bestOrd <- beta0best
for(mm in 1:M){
  beta0bestOrd[,mm] <- beta0best[match(pos,PosAll[,mm]),mm]
  
}
ibias2b0 <- sum((apply(beta0bestOrd, 1, mean, na.rm = T) - beta.0)^2)


#ivarMax <- matrix(NA, nrow = M, ncol = n.snp)

ivar <- rep(NA, n.snp)
for(bb in 1:n.snp){
  see <-  (BetaAllposOrd[,,bb]- apply(BetaAllposOrd[  ,,bb], 1, mean, na.rm = T))^2
  ivar[bb] <- sum(apply(see, 1, mean, na.rm = T))
}

# for the intercept
ivarb0 <- sum(apply((beta0bestOrd-apply(beta0bestOrd, 1, mean, na.rm = T))^2, 1, mean, na.rm = T))


all.equal(ivarb0+ibias2b0, mean(IMSEbeta0, na.rm = T))
all.equal(ivar+ibias2, apply(IMSE, 2, mean, na.rm = T))




#--------------
# Merge results to a cleaner format
#--------------

est_summary <- data.frame(matrix(NA, nrow = n.snp+1, ncol = 3))
colnames(est_summary) <- c("ibias2", "ivar", "imse")
rownames(est_summary) <- paste0("beta", 0:n.snp)

est_summary$ibias2 <- c(ibias2b0, ibias2)

est_summary$ivar <- c(ivarb0, ivar)

est_summary$imse <- c( mean(IMSEbeta0, na.rm = T), apply(IMSE, 2, mean, na.rm = T))

all.equal(est_summary[,3], est_summary[,1]+est_summary[,2])
#--------------------
# for ssp(alpha = 0)
#-------------------

ibias20 <- rep(NA, n.snp)

BetaAllposOrd = beta_all_alp0  # The rownames of BetaAllposOrd corresponds exactly to the the saved object pos

for(mm in 1:M){
  BetaAllposOrd[,mm,] <- beta_all_alp0[match(pos,PosAll[,mm]),mm,]
}

for(bb in 1:n.sig.snp){
  ibias20[bb] <- sum((apply(BetaAllposOrd[  ,,bb], 1, mean, na.rm = T) - BETAs[,bb])^2)
}
for(bb in (n.sig.snp+1):n.snp){
  ibias20[bb] <- sum((apply(BetaAllposOrd[  ,,bb], 1, mean, na.rm = T))^2)
}


#ivarMax <- matrix(NA, nrow = M, ncol = n.snp)

ivar0 <- rep(NA, n.snp)
for(bb in 1:n.snp){
  see <-  (BetaAllposOrd[,,bb]- apply(BetaAllposOrd[  ,,bb], 1, mean, na.rm = T))^2
  ivar0[bb] <- sum(apply(see, 1, mean, na.rm = T))
}


#------------------
# for the intercept
#---------
beta0bestOrd <- beta0alp0
for(mm in 1:M){
  beta0bestOrd[,mm] <- beta0alp0[match(pos,PosAll[,mm]),mm]
  
}
ibias2b0alp0 <- sum((apply(beta0bestOrd, 1, mean, na.rm = T) - beta.0)^2)
# for the intercept
ivarb0alp0 <- sum(apply((beta0bestOrd-apply(beta0bestOrd, 1, mean, na.rm = T))^2, 1, mean, na.rm = T))


all.equal(ivarb0alp0+ibias2b0alp0, mean(IMSEalp0, na.rm = T))


all.equal(ivar0+ibias20, apply(IMSE0, 2, mean, na.rm = T))




#-------------


est_summary_alp0 <- est_summary

est_summary_alp0$ibias2 <- c(ibias2b0alp0, ibias20)

est_summary_alp0$ivar <- c(ivarb0alp0, ivar0)

est_summary_alp0$imse <- c( mean(IMSEalp0, na.rm = T), apply(IMSE0, 2, mean, na.rm = T))

all.equal(est_summary_alp0[,3], est_summary_alp0[,1]+est_summary_alp0[,2])






save.image(file="ResAll.RData")



