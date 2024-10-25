#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace std;
using namespace arma; 


//[[Rcpp::export]]
double assigFun_g(double x, double y, double gamma = 2){
  if(x < 0.000001 && x > -0.000001){
    return 1; 
  }
  if(x < 1.000001 && x > 0.999999){
    return 0; 
  }
  else{
    double t1 = y * pow((y / x), gamma); 
    double t2 = (1 - y) * pow((1 - y) / (1 - x), gamma); 
    double p = t1 / (t1 + t2); 
    return p; 
  }
}

//[[Rcpp::export]]
double target_alloc(double muA, double sigmaA, 
                    double muB, double sigmaB, 
                    Rcpp::String target = "Neyman", 
                    double TB = 4){
  //Rcpp::Rcout<<"mu1, mu0, sigma1, sigma0: "<<muA<<' '<<muB<<" "<<sigmaA * sigmaA <<" " <<sigmaB*sigmaB<<std::endl;
  double rho = 0; 
  if(sigmaA<1e-5)sigmaA=0.1;
  if(sigmaB<1e-5)sigmaB=0.1;
  if(target == "Neyman"){
    rho = sigmaA /(sigmaA + sigmaB);
  }else if(target == "RSIHR"){
    if(muB<=0)muB=1e-5;
    if(muA<=0)muA=1e-5;
    rho = sigmaA * sqrt(muB) / (sigmaA * sqrt(muB) + sigmaB * sqrt(muA)); 
  }else if(target == "BandBis"){
    rho = arma::normcdf((muA - muB) / TB); 
  }else if(target == "ZhangRosenberger"){
    double Rstar = sigmaA * sqrt(muB) / (sigmaB * sqrt(muA));
    double s; 
    if((muA < muB && Rstar > 1) || (muA > muB && Rstar < 1)){
      s = 1; 
    }else{
      s = 0; 
    }
    if(s == 1){
      rho = sigmaA * sqrt(muB) / (sigmaA * sqrt(muB) + sigmaB * sqrt(muA)); 
    }else{
      rho = 0.5; 
    }
  }else if(target == "New"){
    double rho_Neyman = sigmaA / (sigmaA + sigmaB);
    double rho_TB = (TB - muB) / (muA - muB);
    if( rho_Neyman*muA+(1-rho_Neyman)*muB <= TB){
      rho = rho_Neyman;
    }else if(rho_TB > 0){
      rho = rho_TB;
    }else{
      rho = rho_Neyman;
    }
  }
  //threshold
  if(rho>0.9)rho=0.9;
  if(rho<0.1)rho=0.1;
  return rho; 
}


arma::uvec Csample(int n, int num, bool replace, arma::vec proba) {
  Rcpp::NumericVector a(n);
  for(int i =0;i<n;i++){
    a(i) = i + 1;
  }
  Rcpp::NumericVector ret = Rcpp::RcppArmadillo::sample(a, num, replace, proba);
  return Rcpp::as<arma::uvec>(ret) - 1;
}

//[[Rcpp::export]]
arma::vec n0Rand(int n0){
  arma::uvec ind = Csample(2 * n0, n0, FALSE, linspace<arma::vec>(1, 1, 2 * n0)); 
  arma::vec Tvec(2 * n0, arma::fill::zeros); Tvec(ind).fill(1); 
  return Tvec; 
}

//[[Rcpp::export]]
arma::vec DBCD(arma::mat Ymat, double gamma = 2, int n0 = 10, 
               Rcpp::String target = "BandBis", double TB = 4){
  int n = Ymat.n_rows; 
  arma::vec Tvec(n); Tvec.fill(arma::datum::nan); 
  Tvec.subvec(0, 2 * n0 - 1) = n0Rand(n0); 
  
  for(int i = 2 * n0; i < n; i++){
    
    arma::vec Tvec_t = Tvec.subvec(0, i - 1); 
    
    //double N1 = sum(Tvec_t); 
    //int N2 = i - N1; 
    double x = arma::mean(Tvec_t); 
    
    arma::uvec ind1 = arma::find(Tvec_t > 0.99); 
    arma::uvec ind0 = arma::find(Tvec_t < 0.01); 
    
    arma::vec Y1 = Ymat.col(1).rows(0, i - 1), 
      Y0 = Ymat.col(0).rows(0, i - 1); 
    
    double Y1bar = mean(Y1(ind1)), 
      Y0bar = mean(Y0(ind0)); 
    
    double sigma1 = arma::stddev(Y1(ind1)),
      sigma0 = arma::stddev(Y0(ind0));
    double rhohat = target_alloc(Y1bar, sigma1, 
                                 Y0bar, sigma0, target, TB);
    
    double g = assigFun_g(x, rhohat, gamma); 
    
    arma::vec Temp = arma::randu<arma::vec>(1); 
    Tvec(i) = sum(Temp > 1 - g); 
  }
  
  return Tvec; 
}

//[[Rcpp::export]]
arma::vec CARA(List model_output, double gamma = 2, int n0 = 10, 
               std::string target = "Neyman", double TB = 1) {
  
  arma::mat Ymat = as<arma::mat>(model_output["Ymat"]);
  arma::vec strata = as<arma::vec>(model_output["strata"]);
  
  int n = Ymat.n_rows;
  arma::vec An(n);
  
  arma::vec strata_set = arma::unique(strata);
  
  for(unsigned int s = 0; s < strata_set.n_elem; ++s) {
    arma::uvec strata_ind = arma::find(strata == strata_set(s));
    arma::mat Ymat_strata = Ymat.rows(strata_ind);
    
    An.elem(strata_ind) = DBCD(Ymat_strata, gamma, n0, target, TB);
  }
  
  return An;
}

// [[Rcpp::export]]
arma::vec CADBCD(List model_output, double gamma = 2, int n0 = 30, 
                 Rcpp::String target = "Neyman", double TB = 30) {
  arma::mat Ymat = as<arma::mat>(model_output["Ymat"]);
  arma::vec strata = as<arma::vec>(model_output["strata"]);
  
  int n = Ymat.n_rows;
  arma::vec An(n, fill::zeros);
  arma::vec pi(max(strata), fill::zeros);
  
  // Initialize the first 2*n0 elements of An with a random sample of 0's and 1's
  arma::vec initial_sample = n0Rand(n0);
  for (int i = 0; i < 2 * n0; ++i) {
    An[i] = initial_sample[i];
  }
  
  // Main loop
  for (int i = 2 * n0; i < n; ++i) {
    double rho = 0.0;
    
    for (int s = 1; s <= max(strata); ++s) {
      arma::uvec pre_seq = arma::regspace<uvec>(0, i - 1);
      arma::uvec ind = find(strata(pre_seq) == s);
      arma::uvec ind1 = find(An(pre_seq) == 1 && strata(pre_seq) == s);
      arma::uvec ind0 = find(An(pre_seq) == 0 && strata(pre_seq) == s);
      
      arma::vec Y1 = Ymat.col(1);
      arma::vec Y0 = Ymat.col(0);
      double Y1bar = 0;
      double Y0bar = 0;
      
      double sigma1 = 0;
      double sigma0 = 0;
      if(ind1.n_elem==0 || ind0.n_elem==0){
        ;
      }else{
        Y1bar = mean(Y1(ind1));
        Y0bar = mean(Y0(ind0));
        
        sigma1 = arma::stddev(Y1(ind1));
        sigma0 = arma::stddev(Y0(ind0));
      }
      
      
      // Call the target_alloc function (you need to define it properly)
      pi[s - 1] = target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB);
      
      rho = rho + ind.n_elem*1.0/i * pi[s - 1];
      //printf("%f\n\n",sum(strata(pre_seq) == s)/i);
      
    }
    
    uint N1 = arma::sum(An.subvec(0, i - 1));
    uint j = strata[i] - 1; // Adjust indexing for C++ (0-based)
    
    double prob = pi[j] * pow( (rho / N1 * (i - 1)), gamma) /
      (pi[j] * pow( (rho / N1 * (i - 1)), gamma) +
        (1 - pi[j]) * pow( (1 - rho) / (1 - N1*1.0 / (i - 1)), gamma));
    
 //printf("%f\n%f\n\n",pi[j],prob);
  prob = std::min(std::max(prob, 0.1), 0.9);  // Ensure prob is between 0.1 and 0.9
  
  arma::vec Temp = arma::randu<arma::vec>(1); 
  An(i) = sum(Temp > 1 - prob);
  }
  
  return An;
}



//[[Rcpp::export]]
arma::vec CRand(arma::mat Ymat, double delta = 0.5)
{
  int n = Ymat.n_rows; 
  arma::vec Tvec(n); Tvec.fill(arma::datum::nan); 
  for(int i = 0; i < n; i++){
    arma::vec Temp = arma::randu<arma::vec>(1); 
    Tvec(i) = sum(Temp > 1 - delta); 
  }
  
  return Tvec; 
}

