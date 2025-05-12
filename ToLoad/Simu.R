library(magrittr)
library(mvtnorm)

BSMethSim_bbinom <-function(n, posit, theta.0, beta, phi, random.eff = F, mu.e=0,
                            sigma.ee=1,p0 = 0.003, p1 = 0.9, X, Z,binom.link="logit"){
  if( !is.matrix((Z)) ) message ("covariate Z is not a matrix")
  #  if( !is.matrix(beta) ) message ("the covariate effect parameter beta is not a matrix")
  
  if( !(nrow(X)==nrow(Z) & nrow(X) == n) ) message("Both X and Z should have n rows")
  if( !(ncol(X)==length(theta.0) & ncol(X) ==nrow (beta) & ncol(X)==length(posit) )) message ("The columns of X should be the same as length of beta theta.0 and posit; They all equals to the number of CpGs")
  if( ncol(beta)!= ncol(Z)) message("beta and Z should have the same dimentions")
  
  # the random effect term
  if(random.eff == T){
    my.e <- rnorm(n, mean=mu.e, sd = sqrt(sigma.ee))
  }else{
    my.e <- rep(mu.e, n)
  }
  
  
  
  my.theta <- t(sapply(1:n, function(i){
    theta.0 + rowSums(sapply(1:ncol(Z), function(j){Z[i,j] * beta[,j]})) + my.e[i]
  }))
  
  
  # Transform my.theta to my.pi for each (i, j)
  
  my.pi <- t(sapply(1:nrow(my.theta), function(i){
    #exp(my.theta[i,])/(1+exp(my.theta[i,]))
    binomial(link=binom.link)$linkinv(my.theta[i,])
  }))
  #  Generate S-ij based on the my.pi and my.Z
  #---------------------------------------------------#
  my.S <- my.pi
  for ( i in 1:nrow(my.S)){
    for ( j in 1:ncol(my.S)){
      #my.S[i,j] <- rbinom (1, size = X[i,j], prob= my.pi[i,j])
      my.S[i,j] <-VGAM::rbetabinom(1, size = X[i,j], prob = my.pi[i,j], rho = (phi[j]-1)/(X[i,j]-1) )
    }
  }
  #---------------------------------------------------#
  # Generate Y-ij based on the S-ij and the error rate (1-p1) and p0
  #---------------------------------------------------#
  my.Y <- my.S
  for ( i in 1:nrow(my.Y)){
    for ( j in 1:ncol(my.Y)){
      my.Y[i,j] <- sum(rbinom(my.S[i,j], size =1, prob=p1)) +
        sum(rbinom(X[i,j]-my.S[i,j], size = 1, prob=p0))
    }
  }
  out = list(S = my.S, Y = my.Y, theta = my.theta, pi = my.pi)
}

sparseSimu <- function(my.samp=50, n.snp=5, n.sig.snp=1, beta.0, BETAs, phi = rep(1,length(pos)), 
                       add_read_depth =0){
  
  # betas_non_zeros <- matrix(NA, nrow = length(beta.0), ncol = n.sig.snp)
  betas_non_zeros <- BETAs
  #for (i in 1:ncol(betas_non_zeros)){
  #  betas_non_zeros[,i] <- beta.1
  # }
  
  beta_all <- cbind(beta.0, betas_non_zeros)
  
  Z <- data.frame(matrix(NA, nrow= my.samp, ncol = n.snp))
  
  
  mafs <- runif(n.snp,0.1, 0.5)
  for ( i in 1:ncol(Z)){
    Z[,i] <- sample(c(0,1,2), size = my.samp, prob = c((1-mafs[i])^2, 2*mafs[i]*(1-mafs[i]),  mafs[i]^2), replace = T)
  }
  
  
  Z <-as.matrix(Z);rownames(Z)<- NULL
  
  samp.Z <- Z
  
  
  #---------------------------
  # Step 3.2: Read depth matrix X
  #---------------------------
  # Build a read-depth matrix which sort of preserve the dependence structure in read-depths
  
  my.X <- matrix(sample(0:1, my.samp*length(pos), replace = T), nrow = my.samp, ncol = length(pos))
  
  ff = smooth.spline(pos, apply(totalMat, 1, median), nknots = 10)
  
  spacial_shape <- round(predict(ff, pos)$y)
  for ( i in 1:my.samp){
    my.X[i,] <- my.X[i,] + spacial_shape 
    #+ round(10*beta.0) + round(beta.1 * Z[i,1] + beta.2 * Z[i, 2])
  } 
  
  my.X <-  (my.X + add_read_depth)
  
  
  #---------------------------
  # Step 3.2: Simulate the methylated count matrix
  #---------------------------
  
  sim.dat<-BSMethSim_bbinom(n= my.samp, posit = pos, theta.0 =beta_all[,1], beta= as.matrix(beta_all[,-1]), phi=phi, 
                            X = my.X, Z =as.matrix(Z[, 1:n.sig.snp], ncol = n.sig.snp),p0 = 0, p1 = 1,random.eff = F)
  
  
  
  plot(pos, sim.dat$pi[1,], ylim = c(0,1))
  
  for( i in 1:my.samp){
    points(pos, sim.dat$pi[i,], pch = 19, cex =0.5, col = i)
  }
  
  
  
  #-------------------------------------
  # Generate loss function; proximal gradient descent; and etc.
  
  
  #-- Step 1: generate the matrix of Omega1 and Omega2.
  
  
  # n.k: K number of knots ---- 
  # we choose the same k and same basis function B for all predictors 
  # so, the same mat_omega2 for all p = 1, 2, ... P
  length(pos)
  
  # because curretnly, I didn't add smoothness penalty for the intercept, I will use n.k = 5 for the intercept
  n.k = 10 # for the rest of Zs (non-intercept predictors)
  
  
  
  
  
  #--- Organize the data before EM-smooth ---#
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  
  dat.use <- data.frame(Meth_Counts=as.vector(t(Y)), 
                        Total_Counts=as.vector(t(X)), 
                        Position = rep(pos, samp.size),
                        ID = rep(1:samp.size, each=my.p))
  
  truePi <- as.vector(t(sim.dat$pi))
  covs_use <- colnames(samp.Z)
  for( j in 1:length(covs_use)){
    dat.use <- data.frame(dat.use, rep(samp.Z[,j], each = my.p))
  }
  colnames(dat.use)[-c(1:4)] <- covs_use
  
  dat.use <- dat.use[dat.use$Total_Counts>0,]
  
  #my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))
  
  #Z <- dat.use[,-c(1:4)]
  
  
  # pre-set parameters
  dat = dat.use
  
  return(list(dat=dat, truePi=truePi))
}


sparseSimuCor <- function(my.samp=50, n.snp=5, n.sig.snp=1, beta.0, BETAs, phi = rep(1,length(pos)), 
                          add_read_depth =0, Sigma){
  
  betas_non_zeros <- BETAs
  
  beta_all <- cbind(beta.0, betas_non_zeros)
  
  Z <- data.frame(matrix(1, nrow= my.samp, ncol = n.snp))
  
  # Generate genotypes with correlation 
  
  normRand <- rmvnorm(my.samp, sigma=as.matrix(Sigma) )
  
  #cor(normRand)
  
  mafs <- runif(n.snp,0.1, 0.5)
  
  cut1 = qnorm(p= (1-mafs)^2)  # The 2 cutoff
  cut2 = qnorm(p = mafs^2, lower.tail = FALSE)
  
  for(i in 1:my.samp){
    Z[i, normRand[i,]<cut1]=0
    Z[i, normRand[i,]>cut2]=2
  }
  
  Z <-as.matrix(Z);rownames(Z)<- NULL
  
  samp.Z <- Z
  
  
  #---------------------------
  # Step 3.2: Read depth matrix X
  #---------------------------
  # Build a read-depth matrix which sort of preserve the dependence structure in read-depths
  
  my.X <- matrix(sample(0:1, my.samp*length(pos), replace = T), nrow = my.samp, ncol = length(pos))
  
  ff = smooth.spline(pos, apply(totalMat, 1, median), nknots = 10)
  
  spacial_shape <- round(predict(ff, pos)$y)
  for ( i in 1:my.samp){
    my.X[i,] <- my.X[i,] + spacial_shape 
    #+ round(10*beta.0) + round(beta.1 * Z[i,1] + beta.2 * Z[i, 2])
  } 
  
  my.X <-  (my.X + add_read_depth)
  
  
  #---------------------------
  # Step 3.2: Simulate the methylated count matrix
  #---------------------------
  
  sim.dat<-BSMethSim_bbinom(n= my.samp, posit = pos, theta.0 =beta_all[,1], beta= as.matrix(beta_all[,-1]), phi=phi, 
                            X = my.X, Z =as.matrix(Z[, 1:n.sig.snp], ncol = n.sig.snp),p0 = 0, p1 = 1,random.eff = F)
  
  
  
  plot(pos, sim.dat$pi[1,], ylim = c(0,1))
  
  for( i in 1:my.samp){
    points(pos, sim.dat$pi[i,], pch = 19, cex =0.5, col = i)
  }
  
  
  
  #-------------------------------------
  # Generate loss function; proximal gradient descent; and etc.
  
  
  #-- Step 1: generate the matrix of Omega1 and Omega2.
  
  
  # n.k: K number of knots ---- 
  # we choose the same k and same basis function B for all predictors 
  # so, the same mat_omega2 for all p = 1, 2, ... P
  #length(pos)
  
  # because curretnly, I didn't add smoothness penalty for the intercept, I will use n.k = 5 for the intercept
  #n.k = 10 # for the rest of Zs (non-intercept predictors)
  
  
  #--- Organize the data before EM-smooth ---#
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  
  dat.use <- data.frame(Meth_Counts=as.vector(t(Y)), 
                        Total_Counts=as.vector(t(X)), 
                        Position = rep(pos, samp.size),
                        ID = rep(1:samp.size, each=my.p))
  
  truePi <- as.vector(t(sim.dat$pi))
  
  covs_use <- colnames(samp.Z)
  for( j in 1:length(covs_use)){
    dat.use <- data.frame(dat.use, rep(samp.Z[,j], each = my.p))
  }
  colnames(dat.use)[-c(1:4)] <- covs_use
  
  dat.use <- dat.use[dat.use$Total_Counts>0,]
  
  #my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))
  
  #Z <- dat.use[,-c(1:4)]
  
  
  # pre-set parameters
  dat = dat.use
  
  return(list(dat=dat, truePi=truePi))
}