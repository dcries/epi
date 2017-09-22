#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) { 
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma); 
}

arma::cube rinvwish(int n, int v, arma::mat S){
  //RNGScope scope;
  int p = S.n_rows;
  arma::mat L = chol(inv_sympd(S), "lower");
  arma::cube sims(p, p, n, arma::fill::zeros);
  for(int j = 0; j < n; j++){
    arma::mat A(p,p, arma::fill::zeros);
    for(int i = 0; i < p; i++){
      int df = v - (i + 1) + 1; //zero-indexing
      A(i,i) = sqrt(R::rchisq(df)); 
    }
    for(int row = 1; row < p; row++){
      for(int col = 0; col < row; col++){
        A(row, col) = R::rnorm(0,1);
      }
    }
    arma::mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
    sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}


const double log2pi = std::log(2.0 * M_PI);



//// [[Rcpp::export]]
double dmvnrm_arma(arma::vec x,
                   arma::rowvec mean,
                   arma::mat sigma,
                   bool logd = false) {
  //int n = x.n_rows;
  
  int xdim = x.size();
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  //for (int i=0; i < n; i++) {
  ////std::cout << x.row(i) << "\n" << mean;
  
  arma::vec z = rooti * arma::trans( x.t() - mean) ;
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;
  //}
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

double logl_b(arma::mat y, arma::vec tstar, arma::vec beta, arma::mat Sigma, 
              arma::vec priormean, arma::mat priorcov){
  int n;
  double ll = 0.0;
  n = y.n_rows;
  arma::mat means(n,7);
  means.col(0) = beta[3]-beta[0]/(1.0+exp(-beta[1]*(tstar-beta[2])));
  means.col(1) = beta[7]-beta[4]/(1.0+exp(-beta[5]*(tstar-beta[6])));
  means.col(2) = beta[11]-beta[8]/(1.0+exp(-beta[9]*(tstar-beta[10])));
  means.col(3) =  beta[17]-beta[14]/(1.0+exp(-beta[15]*(tstar-beta[16])));
  means.col(4) = beta[18] + beta[19]*tstar;
  means.col(5) = beta[12] + beta[13]*tstar;
  means.col(6) = beta[20] + beta[21]*tstar;
  for(int i=0;i<n;i++){
    ll += dmvnrm_arma(y.row(i).t(),means.row(i),Sigma,true);
  }

  ll += dmvnrm_arma(beta,priormean.t(),priorcov,true);
  return ll;
}

arma::mat calc_mean(arma::vec tstar, arma::vec beta){
  int n;
  n = tstar.size();
  arma::mat means(n,7);
  means.col(0) = beta[3]-beta[0]/(1.0+exp(-beta[1]*(tstar-beta[2])));
  means.col(1) = beta[7]-beta[4]/(1.0+exp(-beta[5]*(tstar-beta[6])));
  means.col(2) = beta[11]-beta[8]/(1.0+exp(-beta[9]*(tstar-beta[10])));
  means.col(3) =  beta[17]-beta[14]/(1.0+exp(-beta[15]*(tstar-beta[16])));
  means.col(4) = beta[18] + beta[19]*tstar;
  means.col(5) = beta[12] + beta[13]*tstar;
  means.col(6) = beta[20] + beta[21]*tstar;
  return means;
}

// [[Rcpp::export]]
List mcmc_epi(arma::mat y, arma::mat tstar, List start, List prior, int nsim,int burn){
  
  arma::vec currentbeta         = as<arma::vec>(start["currentbeta"]);
  arma::mat currentSigma        = as<arma::mat>(start["currentSigma"]);
  arma::mat propcov             = as<arma::mat>(start["propcov"]);
  
  arma::vec bm            = as<arma::vec>(prior["bm"]);
  arma::mat bcov            = as<arma::mat>(prior["bcov"]);
  arma::mat D               = as<arma::mat>(prior["D"]);
  double d            = as<double>(prior["d"]);
  
  int n = y.n_rows;
  int p = bm.size();
  int k = y.n_cols;
  int ntstar = tstar.n_rows;
  
  //storage
  arma::mat beta(nsim,p);
  arma::cube Sigma(k,k,nsim);
  
  //arma::ivec index(nsim);
  // arma::ivec index2(ntstar);
  // arma::vec probs(ntstar);
  // 
  // for(int i=0;i<ntstar;i++){
  //   index2[i] = i;
  //   probs[i] = 1.0/ntstar;
  // }
  Rcpp::IntegerVector index;
  Rcpp::IntegerVector frame = Rcpp::Range(1, ntstar);
  Rcpp::NumericVector wts = Rcpp::rep(1.0/ntstar,ntstar);//Rcpp::runif(n, 0.0, 1.0);
  index = Rcpp::RcppArmadillo::sample(frame, nsim, false, wts / Rcpp::sum(wts));
  //std::cout << "1\n";
  
  arma::vec propb(p);
  double logr;
  double logrprop;
  double logacceptprob;
  arma::mat currentmeans;
  
  for(int i=0;i<nsim;i++){
    //update tstar
    propb = (mvrnormArma(1,currentbeta,0.2618*propcov)).row(0).t();
    
    
    logrprop = logl_b(y,tstar.row(index[i]).t(),propb,currentSigma,bm,bcov); 
    //std::cout << "2\n";
    
    logr = logl_b(y,tstar.row(index[i]).t(),currentbeta,currentSigma,bm,bcov);
    //std::cout << "3\n";
    
    logacceptprob = logrprop - logr;
    if(log(R::runif(0,1))<logacceptprob){
      currentbeta = propb;
    }
    if((i<burn) && (i>20) && (i%20==0)){
      propcov = cov(beta.rows(0,i-1));
    }
    //std::cout << "4\n";
    
    currentmeans = calc_mean(tstar.row(index[i]).t(),currentbeta);
    //std::cout << "5\n";
    
    currentSigma = rinvwish(1,n+d,D+(y-currentmeans).t()*(y-currentmeans));
    //std::cout << "6\n";
    
    beta.row(i) = currentbeta.t();
    //std::cout << "7\n";
    
    Sigma.slice(i) = currentSigma;
    
    if(i % 100==0){
      std::cout << "i= " << i << "\n";
    }
  }
  
  return List::create(
    Named("beta") = beta,
    Named("Sigma") = Sigma,
    Named("propcov") = propcov);
} 

