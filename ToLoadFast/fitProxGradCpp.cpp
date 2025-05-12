#include <RcppArmadillo.h>
#include "/usr/local/opt/libomp/include/omp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double getPrec(double x) {
  return nextafter(x, std::numeric_limits<double>::infinity()) - x;}


NumericVector myseq(int &first, int &last) {
  NumericVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}

// [[Rcpp::export]]
List getSeparateThetaCpp(const NumericVector& theta,
                          const int& nk, 
                          const int& numCovs){
  int size=numCovs+1;
  List thetaSep(size);
  for (int i=0; i < size; ++i){
    int lower = (i*nk);
    int upper = (i+1)*nk-1;
    thetaSep[i] = theta[myseq(lower, upper)];
  }
  return(thetaSep);
}


// [[Rcpp::export]]
arma::vec meanFromThetaArma(const arma::vec& theta,
                            const arma::mat& bigDesign, 
                            const bool& truncation){
  
  arma::vec lp_ij = bigDesign * theta;
  lp_ij = 1/(1+trunc_exp(-lp_ij));
  //NumericVector pi_ij1 = NumericVector(lp_ij.begin(), lp_ij.end());
  double eps=getPrec(1) * 10;
  if(truncation==TRUE){
    lp_ij.elem(arma::find(lp_ij > (1 - eps))).fill(1 - eps);
    lp_ij.elem(arma::find(lp_ij < eps)).fill(eps);
  }
  return(lp_ij);
}










// I removed the PathRaw and GridRaw, which are now inside the folder 'src_new_big_func'


// I will further modify the series of functions fitProxGradCppClean1New; meanFromTheta; newGrad; thetaUpdate and oneUpdate
// to avoid transforming from NumericVector to arma::vec or vice versas
// to avoid making copies of objects




// [[Rcpp::export]]
arma::vec newGradArma(const arma::vec& meth,
                      const arma::vec& total,
                      const arma::vec& pi_ij,
                      const arma::mat& bigDesign){
  arma::vec gPi_ij = -2 *(meth - total % pi_ij);
  arma::vec grad = trans(bigDesign) * gPi_ij;
  return grad;
}



// [[Rcpp::export]]
arma::colvec proximalOperatorCppArma(const double& t, 
                                  const double& lambda1, 
                                  const arma::vec& u_p,
                                  const int& nk){
  
  double res1 = std::inner_product(u_p.begin(), u_p.end(), u_p.begin(), 0.0);
  // double res=sqrt(dot(temp1, u_p));
  res1=sqrt(res1);
  res1 =1-(t*lambda1)/res1;
  if(res1>=0){
    arma::colvec out =res1 *u_p;
    return out;
  }else{
    arma::colvec out(nk);
    out.fill(0);
    return out;
  }
}

//'@title thetaUpdateCpp
//'@description One proximal gradient descent update given the stepSize, current theta, and current gradient

// [[Rcpp::export]]

arma::vec thetaUpdateCppNewArma(const double& stepSize,
                           const arma::vec& theta, 
                           const arma::vec& gBinomLoss, 
                           const int& nk,
                           const int& numCovs, 
                           const double& lambda1){
  
  arma::vec theta_l =  theta - stepSize * gBinomLoss;
  
  for (unsigned int i=1; i < numCovs+1; ++i){ // only iterate across the covariates with penalization
    theta_l(span(i*nk, (i+1)*nk-1)) = proximalOperatorCppArma(stepSize, lambda1, theta_l(span(i*nk, (i+1)*nk-1)), nk);
  }
  return(theta_l);
  
}

// [[Rcpp::export]]
//'@title oneUpdateCpp
//'@description One proximal gradient descent update given the stepSize



List oneUpdateCppNewArma(const arma::vec& theta,
                         const double& current_lossval,
                         const arma::vec& gBinomLossNum,
                         double stepSize,    // Not adding the pointer for stepSize variable is important, other wise, stepSize will be changed after running the line searching
                         const double& lambda1,
                         const arma::vec& meth,
                         const arma::vec& unmeth,
                         const arma::vec& total,
                         const arma::mat bigDesign,
                         const int& nk,
                         const int& numCovs,
                         const double& shrinkScale,
                         const bool& truncation){
  
  arma::vec thetaNew = thetaUpdateCppNewArma(stepSize,theta, gBinomLossNum,
                                   nk, numCovs, lambda1);
  arma::vec  Gttheta = (theta - thetaNew)/stepSize;
  
  // Step1: estimate pi_ij
  arma::vec pi_ij =  meanFromThetaArma(thetaNew,bigDesign, truncation);
  // Step2: loss function
  double newloss = (-2)*arma::sum(meth % arma::log(pi_ij) + unmeth % arma::log(1-pi_ij));
  
  
  double innerdot1 = std::inner_product(gBinomLossNum.begin(), gBinomLossNum.end(), Gttheta.begin(), 0.0);
  double innerdot2 = std::inner_product(Gttheta.begin(), Gttheta.end(), Gttheta.begin(), 0.0);
  double linearSupport =  current_lossval - stepSize * innerdot1 + stepSize/2*innerdot2;
  
  
  bool   shrinkCondi = newloss >linearSupport;
  while(shrinkCondi == TRUE){
    stepSize=stepSize*shrinkScale;
    thetaNew = thetaUpdateCppNewArma(stepSize,theta, gBinomLossNum,
                                nk, numCovs, lambda1);
    Gttheta = (theta - thetaNew)/stepSize;


    pi_ij =  meanFromThetaArma(thetaNew,bigDesign, truncation);
    newloss = (-2)*arma::sum(meth % arma::log(pi_ij) + unmeth % arma::log(1-pi_ij));
    
    //--------------------//
    
    double innerdot1 = std::inner_product(gBinomLossNum.begin(), gBinomLossNum.end(), Gttheta.begin(), 0.0);
    double innerdot2 = std::inner_product(Gttheta.begin(), Gttheta.end(), Gttheta.begin(), 0.0);
    linearSupport =  current_lossval - stepSize * innerdot1 + stepSize/2*innerdot2;
    
    shrinkCondi = newloss >linearSupport;
  }
  
  
  List output=List::create(Named("theta_l_proximal")=thetaNew, 
                           Named("stepSize")=stepSize,
                           Named("pi_ij_new") = pi_ij,
                           Named("current_lossval")=newloss);
  return(output);
  
}



// [[Rcpp::export]]
double twoPenaltiesCppArmaCleaner(const arma::vec& theta,
                           const double& lambda1, 
                           const int& numCovs, 
                           const int& nk){
  
  NumericVector squaredIndPen(numCovs);
  arma::vec temp(nk);
  for(int i=0; i<numCovs; ++i){
    temp =  theta(span((i+1)*nk, (i+2)*nk-1));
    squaredIndPen[i] = std::inner_product(temp.begin(), temp.end(), temp.begin(), 0.0);
  }
  return( lambda1 * sum(sqrt(squaredIndPen )));
}







// [[Rcpp::export]]
//'@title fitProxGradCpp
//'@description use proximial gradient descent with backtracking line search to minimize
//'a penalized negative binomial likelihood with a group LASSO penalty. This function is
//'for the tilda version of our objective function. So matrix Hp is not needed
//'@name fitProxGradCpp
//'@param theta the initial value for theta tilda
//'@param intStepSize the initial step size used
//'@param lambda1 the penalty parameter lambda in our paper
//'@param dat data frame with two columns named as "Meth_Counts" and "Total_Counts"
//'@param basisMat0 design matrix for the intercept. Its row equals to the row of dat
//'@param nk number of knots used for the covariate --- currently it is the same for all
arma::vec fitProxGradCppClean1NewArmaThetaOnly(const arma::vec& Inittheta, 
                                 double& intStepSize,
                                 const double& lambda1,
                                 const arma::vec& meth,
                                 const arma::vec& total,
                                 const arma::vec& unmeth,
                                 const arma::mat& bigDesign, 
                                 const int& nk, 
                                 const int& numCovs, 
                                 const double& maxInt,
                                 const double& epsilon,
                                 const double& shrinkScale,
                                 const bool& truncation){
  //Initialization
  int iter=0;
  int nsamp=bigDesign.n_rows;
  
  arma::vec theta=Inittheta;
  // Step1: estimate pi_ij
  arma::vec pi_ij =  meanFromThetaArma(Inittheta,bigDesign, truncation);
  // Step2: loss function
  double current_lossval = (-2)*arma::sum(meth % arma::log(pi_ij) + unmeth %  arma::log(1-pi_ij));
  // Step 3: calculate gradients
  arma::vec gBinomLossNum=newGradArma(meth,total,pi_ij,bigDesign);
  
  
  
  // Proximal gradient descent update
  List out=oneUpdateCppNewArma(theta, current_lossval,gBinomLossNum,intStepSize, lambda1, meth, unmeth, total,
                               bigDesign, nk, numCovs,shrinkScale,truncation);
  
  arma::vec theta_new = out["theta_l_proximal"];
  double lossnew = out["current_lossval"];
  arma::vec pi_ij_new = out["pi_ij_new"];
  
  
  // calculating stopping rules   -- criterion 1: theta between two iterations are very close
  arma::vec diff=theta_new-theta;
  double innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double tol = pow(innerdot, 0.5);
  
  // criterion 2: objective function(binomloss + penalty) are too close
  
  double penTerms = twoPenaltiesCppArmaCleaner(theta_new, lambda1, numCovs, nk);
  double lossSumOld = lossnew + penTerms;
  double lossSumNew = lossSumOld;
  // 
  double tol2 = 100; // the second tolerance is the lossSum/observation are very close
  
  
  while(((tol > epsilon) & (iter < maxInt-1)) & (tol2 >epsilon)){
    
    lossSumOld = lossSumNew;
    theta = theta_new;
    current_lossval = lossnew;
    
    
    // Step 3: calculate gradients
    gBinomLossNum=newGradArma(meth,total,pi_ij_new,bigDesign);
    out = oneUpdateCppNewArma(theta, current_lossval,gBinomLossNum,intStepSize, lambda1, meth, unmeth, total,
                              bigDesign, nk, numCovs,shrinkScale,truncation);
    arma::vec temp1 = out["theta_l_proximal"];
    theta_new = temp1;
    lossnew = out["current_lossval"];
    arma::vec temp2 = out["pi_ij_new"];  
    pi_ij_new =temp2;
    
    
    // check condition
    diff=theta_new-theta;
    innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    tol = pow(innerdot, 0.5);
    
    // Add another tolerance checking for the values of objective functions  
    penTerms = twoPenaltiesCppArmaCleaner(theta_new, lambda1, numCovs, nk);
    lossSumNew = lossnew + penTerms;
    
    tol2 = (lossSumOld-lossSumNew)/nsamp;
    iter = iter+1;
    
    //   Rcout << "tol1 : " << tol << "\n";
    //   Rcout << "tol2 : " << tol2 << "\n";
    //   Rcout << "tol2 : " << iter << "\n";
  }
  
 
  return(theta_new);
  
}






// [[Rcpp::export]]
arma::mat sparseSmoothPathCpp(const arma::vec& ulam,
                                              const arma::vec& intTheta, 
                                              double& intStepSize,
                                              const arma::vec& meth,
                                              const arma::vec& unmeth,
                                              const arma::vec& total,
                                              const arma::mat& bigDesign, 
                                              const int& nk, 
                                              const double& maxInt,
                                              const double& epsilon,
                                              const double& shrinkScale,
                                              const int& numCovs, 
                                              const bool& truncation){
  int myp = (numCovs+1)*nk;
  int nlam = ulam.n_elem;
  
  arma::mat thetaMat(myp, nlam);
  thetaMat.col(0)  = fitProxGradCppClean1NewArmaThetaOnly(intTheta, intStepSize, ulam[0], meth, total,unmeth, bigDesign, nk, numCovs,
                                                          maxInt, epsilon, shrinkScale, truncation);
  for (unsigned int i=1; i < nlam; ++i){
    
    thetaMat.col(i)  = fitProxGradCppClean1NewArmaThetaOnly(thetaMat.col(i-1), intStepSize, ulam[i], meth, total,unmeth, bigDesign, nk, numCovs,
                 maxInt, epsilon, shrinkScale, truncation);
  
  }
  return(thetaMat);
}






// [[Rcpp::export]]
arma::mat extractDesignMat1Cppmat_userep(const int&numCovs,
                                         const arma::mat& basisMat1,
                                         const arma::mat& snpdat,
                                         const int&nk){
  
  arma::mat X = repmat(basisMat1, 1, numCovs);
  arma::uvec indices(nk);
  for(unsigned int i=0; i<numCovs; i++){
    
    for(unsigned int j=0; j < nk;j++){
      indices(j) = i*nk+j;
    }
    X.each_col(indices) %= snpdat.col(i);
  }
  return X;
}








//'@description faster way to do the sweep for each snp
//'@param
// [[Rcpp::export]]
arma::mat extractDesignMat1Cppmat_userep_all(const int& numCovs,
                                             const arma::mat& basisMat1,
                                             const arma::mat& snpdat,
                                             const arma::mat& basisMat0,
                                             const int&nk){
  
  arma::mat X = repmat(basisMat1, 1, numCovs+1);
  X.cols(0, nk-1) = basisMat0;
  arma::uvec indices(nk);
  for(unsigned int i=0; i<numCovs; i++){
    for(unsigned int j=0; j < nk;j++){
      indices(j) = (i+1)*nk+j;
    }
    X.each_col(indices) %= snpdat.col(i);
  }
  return X;
}


//'@description faster way to do the sweep for each snp
//'@param
// [[Rcpp::export]]
List extractDesignMatAllalphs(const int& numCovs,
                              const arma::mat& basisMat1,
                              const arma::mat& snpdat,
                              const arma::mat& basisMat0,
                              const int&nk,
                              const int&nlam2,
                              const List&Linv){
  List out(nlam2);
  
  for(unsigned int k=0; k<nlam2;k++){
  arma::mat temp = Linv[k];
    
  arma::mat temp2 = basisMat1 * temp;
  
  arma::mat X = repmat(temp2, 1, numCovs+1);
  X.cols(0, nk-1) = basisMat0 * temp;
  
  arma::uvec indices(nk);
  for(unsigned int i=0; i<numCovs; i++){
    for(unsigned int j=0; j < nk;j++){
      indices(j) = (i+1)*nk+j;
    }
    X.each_col(indices) %= snpdat.col(i);
  }
  out[k] = X;
  }
  
  return out;
}


//'@description faster way to do the sweep for each snp
//'@param
// [[Rcpp::export]]
List extractDesignMatAllalphsAdp(const int& numCovs,
                              const arma::mat& basisMat1,
                              const arma::mat& snpdat,
                              const arma::mat& basisMat0,
                              const int&nk,
                              const int&nlam2,
                              const List&Linv){
  List out(nlam2);
  int nrows = basisMat0.n_rows;
  int ncols = nk * (numCovs+1);
  
  arma::uvec indices(nk);
  arma::mat X(nrows, ncols);
  
  for(unsigned int k=0; k<nlam2;k++){
    List Linvind = Linv[k]; // individual Linv for each snp (the first one is for the intercept) -- numCovs+1
    
    arma::mat Linvi =  Linvind[0];
    X.cols(0, nk-1) = basisMat0 * Linvi;
    

    for(unsigned int i=0; i<numCovs; i++){
      
      arma::mat Linvi1 =  Linvind[i+1];
      
      for(unsigned int j=0; j < nk;j++){
        indices(j) = (i+1)*nk+j;
      }
      
      X.cols((i+1)*nk, (i+2)*nk-1) = basisMat1 *  Linvi1 ;
      X.each_col(indices) %= snpdat.col(i);
    }
    out[k] = X;
  }
  
  return out;
}                            


// [[Rcpp::export]]
List binomObjectCppLossOnlyPassMat(const arma::vec& theta,
                                   const arma::mat& bigDesign, 
                                   const arma::vec& meth,
                                   const arma::vec& total,
                                   const arma::vec& unmeth,
                                   const bool& truncation,
                                   const int& nsamp){
  arma::vec pi_ij =  meanFromThetaArma(theta,bigDesign, truncation);
  double average_dev = (-2)*arma::mean(meth % arma::log(pi_ij) + unmeth %  arma::log(1-pi_ij));
  // RMSE in the arcsin sclae
  
  arma::vec obs_pi=meth/total;
  arma::vec obs_pi_as = asin(2*obs_pi-1);
  arma::vec pred_pi_as = asin(2*pi_ij-1);
  
  
  arma::vec diff = pred_pi_as-obs_pi_as ;
  double innerdot1 = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double rmse = sqrt(innerdot1/nsamp);
  // double corraw = arma::cor(obs_pi, pi_ij);
  //  double cortran = arma::cor(obs_pi_as,pred_pi_as);
  
  List output=List::create(Named("average_dev")=average_dev, 
                           Named("rmse")=rmse);
  //   Named("corraw") = corraw,
  //   Named("cortran")=cortran);
  return(output);
  
}

//'@title fitProxGradGridRaw.cpp
//'@description given a grid of lambda and a sequence of alpha fit the sequence for each alpha
//'@param lamGrid A grid values of lambda1, each column corresponds to a different ulam2 
//'@param ulam2 a sequence of alpha/lambda2
//'@param Linv a list for the matrix Linv; its length equal to length(ulam2); inverse of the cholesky decomp L; H = t(L)%*%L
//'@param designMat1 design matrix for the original version


// [[Rcpp::export]]


arma::cube sparseSmoothGridRawCppArmaList(const arma::mat& lamGrid,
                                          const arma::vec& theta, 
                                          double& intStepSize,
                                          const arma::vec& meth,
                                          const arma::vec& total,
                                          const arma::vec& unmeth,
                                          const List& bigDesignList, 
                                          const int& nk, 
                                          const double& maxInt,
                                          const double& epsilon,
                                          const double& shrinkScale,
                                          const int& numCovs,
                                          const bool& truncation){
  
  int nlam2 = lamGrid.n_cols;
  int myp = (numCovs+1)*nk;
  int nlam = lamGrid.n_rows;
  
  //omp_set_num_threads(nthr);
  
  // Initialization:
  
  arma::cube thetaMat(myp, nlam, nlam2);
  
  // #pragma omp parallel for
  for (int ii=0; ii < nlam2; ++ii){ // loop over different ulam2, different columns of lamGrid
    
    
    thetaMat.subcube(0,0,ii, myp-1,nlam-1,ii) = sparseSmoothPathCpp(lamGrid.col(ii),theta, intStepSize,meth, unmeth, total, bigDesignList[ii],
                     nk, maxInt,epsilon, shrinkScale, numCovs, truncation);
  }
  return( thetaMat);
  
}


// [[Rcpp::export]]
List thetaTildaToOriginal(const int& numCovs,
                          const List&Linv,
                          const List thetaSep){
  List out(numCovs+1);
  for(int ii=0; ii< numCovs+1; ++ii){
    arma::mat linvnow = Linv[ii];
    arma::vec temp = thetaSep[ii];
    out[ii] = linvnow * temp;
    }
  return(out);
  }

