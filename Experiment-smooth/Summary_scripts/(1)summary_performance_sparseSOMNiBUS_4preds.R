



load("../../ToLoad/5betaShapes.RData")

#load("../ToLoad/5betaShapes.RData")
#load("../../ToLoad/4Step_1_curve.RData")


#load("../../ToLoad/5step_funsShapes.RData")
my.samp=50
n.snp=100
n.sig.snp=5
#--------------
M = 100
n.k = 10
numCovs = n.snp

lam2_vec <- c(seq(0, 0.9, 0.1), 0.95, 0.99)
nlam2 = length(lam2_vec)
#-------------------
# Estimation results
#--------------------

PosAll <- matrix(0, nrow = length(pos), ncol = M)

BetaAll <- array(0, c(length(pos), M, n.snp), dimnames = list(NULL, NULL, paste0("beta", 1:n.snp)))

beta_all_alp0 <- BetaAll
#-------------------
# Prediction results
#--------------------

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
  
  # estimation
  
  PosAll[,i] = cvout$bestFit$uniqPos
  
 selCovs =  which(!cvout$bestFit$zeroCovsBool)
  
 
 
 for(jj in selCovs[selCovs<=n.snp]){
   
     BetaAll[,i,jj] = cvout$bestFit$penalBetas[[match(jj, selCovs)]]
 }
 
 selCovs =  which(!cvout$bestFitAll[[1]]$zeroCovsBool)
 for(jj in selCovs[selCovs<=n.snp]){
   beta_all_alp0[,i,jj] = cvout$bestFitAll[[1]]$penalBetas[[match(jj, selCovs)]]
 }
 
 # prediction
 pred_alp0 [i,] <- ind_pred_alp0
 pred[i,] <- ind_pred
 
 # selection
 
 SelectResBest[,i]  <- cvout$bestFit$zeroCovsBool
 
 SlectResBestAlp0[,i] <- cvout$bestFitAll[[1]]$zeroCovsBool
 
 bestLams[i,] <- c(cvout$bestLambda1, cvout$bestLambda2)
 
 trueDev_c[i] <- trueDev
 
 pred_all_alp[i,,] <- t(ind_pred_all_alp)
 
 
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

save.image(file="ResAll.RData")

