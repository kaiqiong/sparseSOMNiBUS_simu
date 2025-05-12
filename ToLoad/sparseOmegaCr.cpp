#include <RcppArmadillo.h>
#include "/usr/local/opt/libomp/include/omp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat sparseOmegaCr (const arma::vec& myh,
                        const int& K,
                        const arma::mat& matF
                          ){
  // initialize matrices A11 A12 A22 to 0
  arma::sp_mat A11(K, K) ;  
  //A11.fill(0) ;
  
  arma::sp_mat A12(K, K) ;  
  //A12.fill(0) ;
  
  arma::sp_mat A22(K, K) ; 
  //A22.fill(0) ;
  // compute A11 A12 and A22
 // double test = (4/315)*pow(myh(0), 5); // no report
 // double test2 = (pow(myh(0), 5) *4)/315; // with report
 // double test3 = (4/315)*pow(myh(0), 5);
  for ( int i = 0; i<K; i++){
    if( i == 0){
      A11(i,i) = myh(i)/3;
      A12(i,i) = -pow(myh(i), 3)/45;
      A22(i,i) = (pow(myh(i), 5)* 4 ) /315 ;
    }
    if (i == K-1){
      A11(i,i) = myh(i-1)/3;
      A12(i,i) = -pow(myh(i-1), 3)/45;
      A22(i,i) = (pow(myh(i-1), 5)* 4 ) /315;
    }
    if((i >0) & (i < (K-1))){
      A11(i,i) = myh(i)/3 + myh(i-1)/3;
      A12(i,i) = -pow(myh(i), 3)/45-pow(myh(i-1), 3)/45;
      A22(i,i) = (pow(myh(i), 5)* 4 ) /315 + (pow(myh(i-1), 5)* 4 ) /315;
    }
    if(i > 0){
      A11(i-1, i) = myh(i-1)/6;
      A11(i, i-1) = A11(i-1, i);
      A12(i-1, i) = (pow(myh(i-1), 3) * (-7))/360;
      A12(i, i-1) = A12(i-1, i) ;
      A22(i-1, i) = (pow(myh(i-1), 5) * 31)/15120;
      A22(i, i-1) = A22(i-1, i) ;
    }
   
  }
  arma::mat A12F = A12 * matF;
  arma::mat FtA12t = arma::trans(A12F);
  arma::mat Omega1 = A11 + FtA12t + A12F + arma::trans(matF) * A22* matF;
  return(Omega1);
}
