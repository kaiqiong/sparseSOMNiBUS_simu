




sss <- 14
#library(sparseSOMNiBUS)
#source("../../ToLoad/sparseSmoothFitPath.R")
source("../../ToLoadFast/sparseSmoothGrid.R")
source("../../ToLoadFast/sparseSmoothPred.R")
source("../../ToLoadFast/utils.R")
source("../../ToLoadFast/sparseSmoothFitCV.R")





source("../../ToLoad/Simu.R")

#source("../../ToLoad/sparseSmoothFitCV.R")
#setwd("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("../../ToLoadFast/fitProxGradCpp.cpp")
sourceCpp("../../ToLoad/sparseOmegaCr.cpp")



library(magrittr)
library(mvtnorm)



#-------------------
# Step1 : Load a real data set 
#---------------------

load("../../ToLoad/5betaShapes.RData")
#beta.1 <- 3* beta.1


#----------
my.samp=50
n.snp=100
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

rho = 0.3

block_n_snp <- 20
Sigma <- matrix(rep(rho,block_n_snp^2), nrow = block_n_snp, ncol = block_n_snp)
diag(Sigma) <- 1

library(Matrix)

temp <- vector("list", n.snp/block_n_snp)
Sigma <- lapply(1:length(temp), function(i){Sigma})

Sigma <- bdiag(Sigma)


#for ( mm in 1:M){
set.seed(32314242+mm)

dat = sparseSimuCor(my.samp, n.snp, n.sig.snp, beta.0, BETAs, Sigma=Sigma)


time0 = Sys.time()
cvout = sparseSmoothFitCV(dat$dat, n.k, stepSize, lambda = lambda, nlam = nlam, lam2 = lam2, nlam2=nlam2, maxInt = maxInt,
                          epsilon, printDetail = TRUE, initTheta, shrinkScale=0.5,
                          accelrt = FALSE, nfolds = 5, mc.cores=mc.cores)
print(Sys.time()-time0)

print(paste0("Sample Size", my.samp ))