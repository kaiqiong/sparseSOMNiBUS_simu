



#load("../../ToLoad/5SmoothBetas.RData")

load("../ToLoad/5betaShapes.RData")
#load("../../ToLoad/4Step_1_curve.RData")


#load("../../ToLoad/5step_funsShapes.RData")
#----------
my.samp=50
n.snp=1000
n.sig.snp=5
#--------------
M = 50
n.k = 10
numCovs = n.snp



#-------------------
# Estimation results
#--------------------

PosAll <- matrix(NA, nrow = length(pos), ncol = M)

BetaAll <- array(NA, c(length(pos), M, n.snp), dimnames = list(NULL, NULL, paste0("beta", 1:n.snp)))
beta0best <- matrix(NA, nrow = length(pos),ncol = M)
#beta_all_alp0 <- BetaAll
#-------------------
# Prediction results
#--------------------

pred <- matrix(NA, nrow = M, ncol = 4)
colnames(pred) <- c("deviance", "tran_rmse", "corRaw", "corTrans")
trueDev_c <- rep(NA, M)
#pred_alp0 <- pred


bestLams <- matrix(NA, nrow = M, ncol = 2)
colnames(bestLams) <- c("lam1", "lam2")
#------------------------
# Selection results
#-------------------------


SelectResBest <- matrix(NA, nrow= n.snp, ncol = M)

#SlectResBestAlp0 <- SelectResBest


for( i in 1:M){
  
  try({
  load(paste0("Res", i, ".RData"))
  
  # estimation
  PosAll[,i] = cvout$bestFit$uniqPos
 selCovs =  which(!cvout$bestFit$zeroCovsBool)

 for(jj in selCovs[selCovs<=n.snp]){
   
     BetaAll[,i,jj] = cvout$bestFit$penalBetas[[match(jj, selCovs)]]
 }
 BetaAll[,i,which(cvout$bestFit$zeroCovsBool)] <- 0
 
 beta0best[,i] <- cvout$bestFit$penalBetas0
 # prediction
# pred_alp0 [i,] <- ind_pred_alp0
 pred[i,] <- ind_pred
 trueDev_c[i] <- trueDev
 # selection
 
 SelectResBest[,i]  <- cvout$bestFit$zeroCovsBool
 
 #SlectResBestAlp0[,i] <- cvout$bestFitAll[[1]]$zeroCovsBool
 
 bestLams[i,] <- c(cvout$bestLambda1, cvout$bestLambda2)
 
 
  })
}

save.image(file="ResAll.RData")



# before calculating IMSE, we need to properly order the array BetaAll



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



#IMSEbeta0 <- apply((beta0bestOrd-beta.0)^2, 2, sum)



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


save.image(file="ResAll.RData")


save.image(file="ResAll.RData")
