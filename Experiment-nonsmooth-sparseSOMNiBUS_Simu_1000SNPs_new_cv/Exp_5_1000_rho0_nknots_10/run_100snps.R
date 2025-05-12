

sss <- commandArgs(TRUE)


sss <- as.numeric(sss)
#library(sparseSOMNiBUS)
source("../../ToLoad/sparseSmoothFitGrid.R")
source("../../ToLoad/sparseSmoothPred.R")
source("../../ToLoad/utils.R")
source("../../ToLoad/sparseSmoothFitCV.R")

source("../../ToLoad/Simu.R")

#source("../../ToLoad/sparseSmoothFitCV.R")
#setwd("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("../../ToLoad/sparseOmegaCr.cpp")
sourceCpp("../../ToLoad/fitProxGradCpp.cpp")


library(magrittr)
library(mvtnorm)
library(doParallel)
library(foreach)


#-------------------
# Step1 : Load a real data set 
#---------------------

load("../../ToLoad/5betaShapes.RData")
#beta.1 <- 3* beta.1


#----------
my.samp=50
n.snp=1000
n.sig.snp=5
#--------------
M = 50
n.k = 10
numCovs = n.snp


lambda = NULL
nlam = 100

lam2_vec <- seq(0, 0.9, 0.1)

gridIDs <- expand.grid(1:10, 1:M); colnames(gridIDs) <- c("lam2", "simu")

mm <- gridIDs[sss,2]
lam2 <- lam2_vec[gridIDs[sss,1]]



#lam2 = 0

nlam2 = length(lam2)

set.seed(1213)
initTheta <- rnorm(n.k*(numCovs+1))
stepSize=0.1
shrinkScale=0.5

maxInt = 10^5
epsilon = 1E-6
accelrt = FALSE

truncation = TRUE
mc.cores = 1

nfolds = 5


#SelectRes <- matrix(NA, nrow = M, ncol = n.snp)

#setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/sparseSOMNiBUS_Simu_no3times/Exp_1_5_rho0")

rho = 0
Sigma <- matrix(rep(rho, n.snp^2), nrow = n.snp, ncol = n.snp)
diag(Sigma) <- 1


time0 = Sys.time()
#for ( mm in 1:M){
set.seed(32314242+mm)

dat = sparseSimu(my.samp, n.snp, n.sig.snp, beta.0, BETAs)



cvout = sparseSmoothFitCV(dat, n.k, stepSize, lambda = lambda, nlam = nlam, lam2 = lam2, nlam2=nlam2, maxInt = maxInt,
                          epsilon, printDetail = TRUE, initTheta, shrinkScale=0.5,
                          accelrt = FALSE, nfolds = 5, mc.cores=mc.cores)

save(cvout, file = paste0("Res", mm, "_alph", lam2, ".RData"))
set.seed(902314242+mm)

testDat = sparseSimu(my.samp, n.snp, n.sig.snp, beta.0, BETAs)


ind_pred = sparseSmoothPredOneSetting(thetaEst=cvout$bestFit$thetaMat[,1], 
                                      fitDat=dat, 
                                      testDat=testDat, 
                                      n.k=n.k, numCovs=numCovs, truncation=TRUE)


# Evaluate the test set prediction




save(cvout,ind_pred, file = paste0("Res", mm, "_alph", lam2, ".RData"))


print(mm)
print(Sys.time()-time0)


