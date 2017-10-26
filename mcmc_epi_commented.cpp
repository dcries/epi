#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Description: random number generator for K dimensional multivariate normal
// n: number of random generations
// mu: mean vector of length K
// sigma: covariance matrix of dim KxK
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) { 
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma); 
}

// Description: random number generator for dirichlet
// a: vector of values for dirichlet parameter a
arma::vec rdirich(arma::vec a){
  int n = a.size();
  double total;
  arma::vec out(n);
  arma::vec gam(n);
  for(int i=0;i<n;i++){
    gam[i] = R::rgamma(a[i],1);
  }
  total = arma::accu(gam);
  out = gam/total;
  return out;
}

// Description: random number generator for inverse wishart
// n: number of random generations
// v: degrees of freedom
// S: scale matrix
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



// Description: density of multivariate normal
// x: vector of values you want density for
// mean: mean vector for observations
// sigma: covariance matrix for all observations
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

// Description: similar to which in R, finds indicies for which x==val
// x: vector of integers
// val: value you want to know for which indicies in x it equals
arma::uvec which(arma::ivec x, double val) {
  arma::uvec ids = find(x == val); // Find indices
  return ids;
}

// Description: subsets data x according to which rows in index==val
// x: data to be subsetted
// index: value determinnig which ``group" each row of x is in
// val: indicates which ``group" you want the subset for
arma::mat subset(arma::mat x, arma::ivec index, int val){
  arma::uvec ids = find(index == val);
  return x.rows(ids);
}

// Description: subsets data x according to which rows in index==val
// x: data to be subsetted
// index: value determinnig which ``group" each value of x is in
// val: indicates which ``group" you want the subset for
arma::vec subset(arma::vec x, arma::ivec index, int val){
  arma::uvec ids = find(index == val);
  return x(ids);
}


// Description: log likelihood + prior for gamma coefficients for SUR regression
// y: matrix of data
// tstar: usual minutes in MVPA
// beta: current values of gamma in regression
// Sigma: covariance matrix for y
// lambda: intercept values for y for all groups
// zeta: zeta values, indicating which group observations are in
// pi: proportion for each group of the mixture
// priormean: prior mean for gamma coefficients
// priorcov: prior cov matrix for gamma coefficients
double logl_b(arma::mat y, arma::vec tstar, arma::vec beta, arma::cube Sigma, 
              arma::mat lambda, arma::ivec zeta, arma::vec pi, arma::vec priormean, arma::mat priorcov){
  int n;
  double ll = 0.0;
  n = y.n_rows;
  arma::mat means(n,7);

  for(int i=0;i<n;i++){
    means(i,0) = lambda(0,zeta[i])-beta[0]/(1.0+exp(-beta[1]*(tstar[i]-beta[2])));
    means(i,1) = lambda(1,zeta[i])-beta[3]/(1.0+exp(-beta[4]*(tstar[i]-beta[5])));
    means(i,2) = lambda(2,zeta[i])-beta[6]/(1.0+exp(-beta[7]*(tstar[i]-beta[8])));
    means(i,3) =  lambda(3,zeta[i])-beta[9]/(1.0+exp(-beta[10]*(tstar[i]-beta[11])));
    means(i,4) =  lambda(4,zeta[i]) + beta[12]*tstar[i];
    means(i,5) =  lambda(5,zeta[i]) + beta[13]*tstar[i];
    means(i,6) =  lambda(6,zeta[i]) + beta[14]*tstar[i];
    
    ll += log(pi[zeta[i]]) + dmvnrm_arma(y.row(i).t(),means.row(i),Sigma.slice(zeta[i]),true);
  }
  
  ll += dmvnrm_arma(beta,priormean.t(),priorcov,true);
  return ll;
}

//Description: calculates mean vectors for each individual
// tstar: values of usual MVPA
// beta: current gamma regression coefficients
// lambda: current intercept terms
arma::mat calc_mean(arma::vec tstar, arma::vec beta, arma::vec lambda){
  int n;
  n = tstar.size();
  arma::mat means(n,7);
  means.col(0) = lambda[0]-beta[0]/(1.0+exp(-beta[1]*(tstar-beta[2])));
  means.col(1) = lambda[1]-beta[3]/(1.0+exp(-beta[4]*(tstar-beta[5])));
  means.col(2) = lambda[2]-beta[6]/(1.0+exp(-beta[7]*(tstar-beta[8])));
  means.col(3) =  lambda[3]-beta[9]/(1.0+exp(-beta[10]*(tstar-beta[11])));
  means.col(4) =  lambda[4] + beta[12]*tstar;
  means.col(5) =  lambda[5] + beta[13]*tstar;
  means.col(6) =  lambda[6] + beta[14]*tstar;
  return means;
}

//Description: calculates mean vectors for each individual minues their respective intercept
// tstar: values of usual MVPA
// beta: current gamma regression coefficients
arma::mat calc_mean_noint(arma::vec tstar, arma::vec beta){
  int n;
  n = tstar.size();
  arma::mat means(n,7);
  means.col(0) = -beta[0]/(1.0+exp(-beta[1]*(tstar-beta[2])));
  means.col(1) = -beta[3]/(1.0+exp(-beta[4]*(tstar-beta[5])));
  means.col(2) = -beta[6]/(1.0+exp(-beta[7]*(tstar-beta[8])));
  means.col(3) =  -beta[9]/(1.0+exp(-beta[10]*(tstar-beta[11])));
  means.col(4) =   beta[12]*tstar;
  means.col(5) =   beta[13]*tstar;
  means.col(6) =   beta[14]*tstar;
  return means;
}

//Description: calculates specific mean vectors for each individual
// tstar: values of usual MVPA
// beta: current gamma regression coefficients
// lambda: matrix of all lambda intercept values for all groups of mixture
// zeta: vector of current zeta values indicating which group each individual belongs to
arma::mat calc_specific_mean(arma::vec tstar, arma::vec beta, arma::mat lambda,arma::ivec zeta){
  int n;
  n = tstar.size();
  arma::mat means(n,7);
  for(int i=0;i<n;i++){
    means(i,0) = lambda(0,zeta[i])-beta[0]/(1.0+exp(-beta[1]*(tstar[i]-beta[2])));
    means(i,1) = lambda(1,zeta[i])-beta[3]/(1.0+exp(-beta[4]*(tstar[i]-beta[5])));
    means(i,2) = lambda(2,zeta[i])-beta[6]/(1.0+exp(-beta[7]*(tstar[i]-beta[8])));
    means(i,3) =  lambda(3,zeta[i])-beta[9]/(1.0+exp(-beta[10]*(tstar[i]-beta[11])));
    means(i,4) =  lambda(4,zeta[i]) + beta[12]*tstar[i];
    means(i,5) =  lambda(5,zeta[i]) + beta[13]*tstar[i];
    means(i,6) =  lambda(6,zeta[i]) + beta[14]*tstar[i];
  }
  
  return means;
}

// Description: function to sample zeta, the group indicator latent variables
// y: matrix of data, nxk matrix, n individuals, k  length vector of obs
// pi: vector of probabilities for each group of mixture
// meanmat: mean values for each individual for each group of mixture
// Sigma: covariance matrix for y for each group of mixture
arma::ivec sample_zeta(arma::mat y, arma::vec pi, arma::cube meanmat, arma::cube Sigma){
  int k = pi.size();
  int n = meanmat.slice(0).n_rows;
  arma::vec probs(k);
  arma::ivec out(n);
  arma::ivec temp(1);
  Rcpp::IntegerVector frame = Rcpp::Range(0,k-1);
  
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      probs[j] = log(pi[j]) + dmvnrm_arma(y.row(i).t(),meanmat.slice(j).row(i),Sigma.slice(j),true);
    }
    probs = exp(probs-log(sum(exp(probs))));
    temp = Rcpp::RcppArmadillo::sample(frame, 1, true, probs);
    out[i] = temp[0];
  }
  return out;
}

// Description: calculates probabilities that individual will be in each of the K groups of the mixture
// y: matrix of data, nxk matrix, n individuals, k  length vector of obs
// pi: vector of probabilities for each group of mixture
// meanmat: mean values for each individual for each group of mixture
// Sigma: covariance matrix for y for each group of mixture
arma::mat calc_pmat(arma::mat y, arma::vec pi, arma::cube meanmat, arma::cube Sigma){
  int k = pi.size();
  int n = meanmat.slice(0).n_rows;
  arma::vec probs(k);
  arma::mat out(n,k);
  
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      probs[j] = log(pi[j]) + dmvnrm_arma(y.row(i).t(),meanmat.slice(j).row(i),Sigma.slice(j),true);
    }
    out.row(i) = exp(probs-log(sum(exp(probs)))).t();
    
  }
  return out;
}

// Description: samples the vector p giving proportions to each group of the mixture
// zeta: vector indicating which group of mixture each individual belongs to
// a: prior for values of p
arma::vec sample_pi(arma::ivec zeta, arma::vec a){
  int k = a.size();
  int n = zeta.size();
  arma::ivec counts(k);
  arma::vec out;
  
  for(int i=0;i<k;i++){
    counts[i] = 0;
  }
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      if(zeta[i]==j){
        counts[j] += 1;
      }
    }
  }
  //std::cout << counts << "\n";
  out = rdirich(a+counts);
  //std::cout << out << "\n";
  return out;
}

// Description: samples covariance matrix for y for a given group
// y: matrix of nxk data
// meanmat: matrix of mean values corresponding to y
// D: prior scale matrix
// d: prior degrees of freedom
// n: nrow(y)
arma::mat sample_Sigma(arma::mat y, arma::mat meanmat, arma::mat D, double d, int n){
  arma::mat out = rinvwish(1,n+d,D+(y-meanmat).t()*(y-meanmat));
  return out;
}

//Description: samples lambda, the intercept terms for all regression models
// y: matrix of nxk data
// meannoint: matrix of mean values of y without intercept
// Sigma: covariance of y
// priormean: vector of prior means for lambda
// priorcov: prior covariance matrix for lambda
arma::vec sample_lambda(arma::mat y, arma::mat meannoint, arma::mat Sigma, arma::vec priormean, arma::mat priorcov){
  int n = y.n_rows;
  arma::vec out;
  arma::vec ybar = (mean(y-meannoint,0)).t();
  arma::mat postcov = arma::inv(arma::inv(priorcov)+n*arma::inv(Sigma));
  arma::vec postmean = postcov*(arma::inv(priorcov)*priormean + n*arma::inv(Sigma)*ybar);
  out = mvrnormArma(1,postmean,postcov).row(0).t();
  
  return out;
}

// Description: samples gamma from the SUR mixture regression model using random walk metropolis
// y: nxk data
// tstar: vector of lengh n corresponding to covariates, or usual min of MVPA
// beta: current values of gamma
// lambda: matrix of lambda values, or intercepts for each group of mixture
// zeta: vector indicating which group each individual belongs to
// Sigma: covariance matricies for each of the groups
// pi: proportions for each of the groups
// propcov: proposal covariance matrix
// bm: prior means for gamma
// bcov: prior covariance for gamma
arma::vec sample_beta(arma::mat y, arma::vec tstar, arma::vec beta, arma::mat lambda,
                      arma::ivec zeta, arma::cube Sigma,
                      arma::vec pi, arma::mat propcov, arma::vec bm, arma::mat bcov,double mult){
  int p = beta.size();
  arma::vec propb(p);
  double logr = 0.0;
  double logrprop = 0.0;
  double logacceptprob;
  //int k = pi.size();
  
  
  propb = (mvrnormArma(1,beta,mult*propcov)).row(0).t();
  //tstar.row(index[i]).t()
  
  logrprop = logl_b(y,tstar,propb,Sigma,lambda,zeta,pi,bm,bcov); 
  //std::cout << "2\n";
  
  logr = logl_b(y,tstar,beta,Sigma,lambda,zeta,pi,bm,bcov);
  //std::cout << "3\n";
  
  
  
  logacceptprob = logrprop - logr;
  if(log(R::runif(0,1))<logacceptprob){
    beta = propb;
  }
  
  return beta;
}

// Description: calculates full log likelihood for y
// y: observed data
// zeta: vector indicating which group each individual belongs to
// specificmeanmat: matrix of specific means for given group
// Sigma: covariance matricies for all groups
double calc_full_ll(arma::mat y, arma::ivec zeta, arma::mat specificmeanmat, arma::cube Sigma){
  int n = y.n_rows;
  double ll=0.0;
  for(int i=0;i<n;i++){
    ll += dmvnrm_arma(y.row(i).t(),specificmeanmat.row(i),Sigma.slice(zeta[i]),true);
  }
  return ll;
}

// Description: calculates the mean value of Sigma for a specific group
// sds: values of standard deviations for each group
// covmat: off diagonal values for each group
arma::cube calc_meanSigma(arma::cube sds, arma::cube covmat){
  int K = sds.n_slices;
  int p = sds.n_cols;
  arma::cube out(p,p,K);
  arma::mat sdmat;
  int count;
  for(int i=0;i<K;i++){
    sdmat = pow(sds.slice(i),2.0);
    out.slice(i).diag() = mean(sdmat,0);
    count = 0;
    for(int j=0;j<(p-1);j++){
      for(int k=(j+1);k<p;k++){
        out(j,k,i) = out(k,j,i) = mean(covmat.slice(i).col(count));
        count++;
      }
    }
  }
  return out;
}

// Description: calculates mean values of lambda for each group
// lambda: values of lambda for each group
arma::mat calc_meanlambda(arma::cube lambda){
  int p = lambda.n_cols;
  int K = lambda.n_slices;
  arma::mat out(p,K);
  
  for(int i=0;i<K;i++){
    for(int j=0;j<p;j++){
      out(j,i) = mean(lambda.slice(i).col(j));
    }
  }
  return out;
}

// [[Rcpp::export]]
List mcmc_epi_mixture(arma::mat y, arma::mat tstar, List start, List prior, 
                      int K, int nsim,int burn, int thin=1, double mult=0.1){
  
  arma::vec currentbeta           = as<arma::vec>(start["currentbeta"]);
  arma::ivec currentzeta          = as<arma::ivec>(start["currentzeta"]);
  arma::vec currentpi             = as<arma::vec>(start["currentpi"]);
  arma::mat currentlambda         = as<arma::mat>(start["currentlambda"]);
  arma::mat Sigmadiag             = as<arma::mat>(start["Sigmadiag"]);
  arma::mat propcov               = as<arma::mat>(start["propcov"]);
  
  arma::vec bm                    = as<arma::vec>(prior["bm"]);
  arma::mat bcov                  = as<arma::mat>(prior["bcov"]);
  arma::vec lm                    = as<arma::vec>(prior["lm"]);
  arma::mat lcov                  = as<arma::mat>(prior["lcov"]);
  arma::mat D                     = as<arma::mat>(prior["D"]);
  double d                        = as<double>(prior["d"]);
  arma::vec a                     = as<arma::vec>(prior["a"]);

  
  int n = y.n_rows;
  int p = bm.size();
  int k = y.n_cols;
  int ntstar = tstar.n_rows;
  int keep = (nsim-burn)/thin;
  
  arma::cube currentSigma = arma::zeros(k,k,K);
  for(int i=0;i<K;i++){
    currentSigma.slice(i).diag() = Sigmadiag.col(i);
  }
  
  //storage
  arma::mat beta(keep,p);
  arma::mat betaburn(burn,p);
  arma::cube Sigma(k,k,keep);
  arma::cube sds(keep,k,K);
  arma::mat pi(keep,K);
  arma::cube lambda(keep,k,K);
  arma::imat zeta(keep,n);
  arma::cube cormat(keep,0.5*k*(k-1),K);
  arma::cube covmat(keep,0.5*k*(k-1),K);
  arma::cube pmat(n,K,keep);
  
  //for dic
  arma::vec full_ll(keep);
  double penalty;
  double dic;
  arma::ivec medianzeta(n);
  arma::cube meanSigma(k,k,K);
  arma::mat meanlambda;
  arma::mat finalmeanmat;
  

  Rcpp::IntegerVector index;
  Rcpp::IntegerVector frame = Rcpp::Range(0, ntstar-1);
  Rcpp::NumericVector wts = Rcpp::rep(1.0/ntstar,ntstar);
  index = Rcpp::RcppArmadillo::sample(frame, nsim, true, wts / Rcpp::sum(wts));

  arma::cube currentmeans(n,k,K);
  arma::mat currentspecificmean;
  arma::mat currentmeansnoint;
  int count;
  int ind = 0;
  
  for(int j=0;j<K;j++){
    currentmeans.slice(j) = calc_mean(tstar.row(index[0]).t(),currentbeta,currentlambda.col(j));
  }
  
  for(int i=0;i<nsim;i++){
    currentmeansnoint = calc_mean_noint(tstar.row(index[i]).t(),currentbeta);

    currentzeta = sample_zeta(y, currentpi, currentmeans, currentSigma);

    currentpi = sample_pi(currentzeta,a);

    currentspecificmean = calc_specific_mean(tstar.row(index[i]).t(),currentbeta,currentlambda,currentzeta);

    for(int j=0;j<K;j++){
      currentSigma.slice(j) = sample_Sigma(subset(y,currentzeta,j),subset(currentspecificmean,currentzeta,j),D,d,(which(currentzeta,j)).size());
    }
    
    for(int j=0;j<K;j++){
      if((which(currentzeta,j)).size() > 0){
        currentlambda.col(j) = sample_lambda(subset(y,currentzeta,j),subset(currentmeansnoint,currentzeta,j),currentSigma.slice(j),lm,lcov);
      }
    }
    
    
    currentbeta = sample_beta(y,tstar.row(index[i]).t(),currentbeta,currentlambda,
                              currentzeta,currentSigma,currentpi,propcov,bm,bcov,mult);
    
    for(int j=0;j<K;j++){
      currentmeans.slice(j) = calc_mean(tstar.row(index[i]).t(),currentbeta,currentlambda.col(j));
    }
    
    if((i<burn) && (i>20) && (i%20==0)){
      propcov = cov(betaburn.rows(0,i-1));
    }
    
    currentspecificmean = calc_specific_mean(tstar.row(index[i]).t(),currentbeta,currentlambda,currentzeta);
    
    
    //storage 
    if(i < burn){
      betaburn.row(i) = currentbeta.t();
    }

    if((i >= burn) && (i%thin==0)){
      pmat.slice(ind) = calc_pmat(y, currentpi, currentmeans, currentSigma);
      full_ll[ind] = calc_full_ll(y,currentzeta,currentspecificmean,currentSigma);
      beta.row(ind) = currentbeta.t();
      zeta.row(ind) = currentzeta.t();
      pi.row(ind) = currentpi.t();
      
      for(int j=0;j<K;j++){
        lambda.slice(j).row(ind) = currentlambda.col(j).t();
        sds.slice(j).row(ind) = sqrt(currentSigma.slice(j).diag().t());
        count = 0;
        for(int l=0;l<(k-1);l++){
          for(int ll=(l+1);ll<k;ll++){
            cormat(ind,count,j) = currentSigma(l,ll,j)/sqrt(currentSigma(l,l,j)*currentSigma(ll,ll,j));
            covmat(ind,count,j) = currentSigma(l,ll,j);
            count++;
          }
        }
      }
      
      ind++;
    }

    if(i % 1000==0){
      std::cout << "i= " << i << "\n";
    }
  }
  
  // to calculate dic penalty term and DIC
  meanlambda = calc_meanlambda(lambda);
  finalmeanmat = calc_specific_mean(mean(tstar,0).t(),mean(beta,0).t(),meanlambda,median(zeta,0).t());
  meanSigma = calc_meanSigma(sds,covmat);
  penalty = calc_full_ll(y,median(zeta,0).t(),finalmeanmat,meanSigma);
  
  dic = -4*mean(full_ll) + 2*penalty;
  
  return List::create(
    Named("beta") = beta,
    Named("lambda") = lambda,
    Named("pi") = pi,
    Named("zeta") = zeta,
    Named("sds") = sds,
    Named("cormat") = cormat,
    Named("dic") = dic,
    Named("meanll") = mean(full_ll),
    Named("pentalty") = penalty,
    Named("pmat") = pmat,
    Named("propcov") = propcov);
}

