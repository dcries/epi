#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) { 
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma); 
}

// [[Rcpp::export]]
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


arma::uvec which(arma::ivec x, double val) {
  
  arma::uvec ids = find(x == val); // Find indices
  //x.elem(ids).fill(replace);       // Assign value to condition
  return ids;
}

arma::mat subset(arma::mat x, arma::ivec index, int val){
  arma::uvec ids = find(index == val);
  return x.rows(ids);
}

arma::vec subset(arma::vec x, arma::ivec index, int val){
  arma::uvec ids = find(index == val);
  return x(ids);
}



double logl_b(arma::mat y, arma::vec tstar, arma::vec beta, arma::cube Sigma, 
              arma::mat lambda, arma::ivec zeta, arma::vec pi, arma::vec priormean, arma::mat priorcov){
  int n;
  double ll = 0.0;
  n = y.n_rows;
  arma::mat means(n,7);
  // means.col(0) = lambda[0]-beta[0]/(1.0+exp(-beta[1]*(tstar-beta[2])));
  // means.col(1) = lambda[1]-beta[3]/(1.0+exp(-beta[4]*(tstar-beta[5])));
  // means.col(2) = lambda[2]-beta[6]/(1.0+exp(-beta[7]*(tstar-beta[8])));
  // means.col(3) =  lambda[3]-beta[9]/(1.0+exp(-beta[10]*(tstar-beta[11])));
  // means.col(4) =  lambda[4] + beta[12]*tstar;
  // means.col(5) =  lambda[5] + beta[13]*tstar;
  // means.col(6) =  lambda[6] + beta[14]*tstar;
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
    probs = probs/sum(probs);
    temp = Rcpp::RcppArmadillo::sample(frame, 1, true, probs);
    out[i] = temp[0];
  }
  return out;
}

arma::vec sample_pi(arma::ivec zeta, arma::vec a){
  int k = a.size();
  int n = zeta.size();
  arma::ivec counts(k);
  arma::vec out;
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      if(zeta[i]==j){
        counts[j] += 1;
      }
    }
  }
  out = rdirich(a+counts);
  return out;
}

arma::mat sample_Sigma(arma::mat y, arma::mat meanmat, arma::mat D, double d, int n){
  arma::mat out = rinvwish(1,n+d,D+(y-meanmat).t()*(y-meanmat));
  return out;
}

arma::vec sample_lambda(arma::mat y, arma::mat meannoint, arma::mat Sigma, arma::vec priormean, arma::mat priorcov){
  int n = y.n_rows;
  arma::vec out;
  arma::vec ybar = (mean(y-meannoint,0)).t();
  arma::mat postcov = arma::inv(arma::inv(priorcov)+n*arma::inv(Sigma));
  arma::vec postmean = postcov*(arma::inv(priorcov)*priormean + n*arma::inv(Sigma)*ybar);
  out = mvrnormArma(1,postmean,postcov).row(0).t();

  return out;
}

arma::vec sample_beta(arma::mat y, arma::vec tstar, arma::vec beta, arma::mat lambda,
                      arma::ivec zeta, arma::cube Sigma,
                      arma::vec pi, arma::mat propcov, arma::vec bm, arma::mat bcov){
  int p = beta.size();
  arma::vec propb(p);
  double logr = 0.0;
  double logrprop = 0.0;
  double logacceptprob;
  //int k = pi.size();
  
  
  propb = (mvrnormArma(1,beta,0.15*propcov)).row(0).t();
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


// [[Rcpp::export]]
List mcmc_epi_mixture(arma::mat y, arma::mat tstar, List start, List prior, int K, int nsim,int burn){
  
  std::cout << "1\n";
  arma::vec currentbeta         = as<arma::vec>(start["currentbeta"]);
  arma::ivec currentzeta         = as<arma::ivec>(start["currentzeta"]);
  arma::vec currentpi         = as<arma::vec>(start["currentpi"]);
  arma::mat currentlambda             = as<arma::mat>(start["currentlambda"]);
  std::cout << "2\n";
  
  arma::vec Sigmadiag        = as<arma::vec>(start["Sigmadiag"]);
  arma::mat propcov             = as<arma::mat>(start["propcov"]);
  
  std::cout << "3\n";
  
  arma::vec bm            = as<arma::vec>(prior["bm"]);
  arma::mat bcov            = as<arma::mat>(prior["bcov"]);
  arma::vec lm            = as<arma::vec>(prior["lm"]);
  arma::mat lcov            = as<arma::mat>(prior["lcov"]);
  arma::mat D               = as<arma::mat>(prior["D"]);
  double d            = as<double>(prior["d"]);
  arma::vec a            = as<arma::vec>(prior["a"]);
  std::cout << "4\n";
  
  
  int n = y.n_rows;
  int p = bm.size();
  int k = y.n_cols;
  int ntstar = tstar.n_rows;
  
  arma::cube currentSigma = arma::zeros(k,k,K);
  for(int i=0;i<K;i++){
    currentSigma.slice(i).diag() = Sigmadiag;
  }
  std::cout << "5\n";
  
  
  //storage
  arma::mat beta(nsim,p);
  arma::cube Sigma(k,k,nsim);
  arma::mat pi(nsim,K);
  arma::cube lambda(nsim,k,K);
  arma::imat zeta(nsim,n);
  std::cout << "6\n";
  
  //arma::ivec index(nsim);
  // arma::ivec index2(ntstar);
  // arma::vec probs(ntstar);
  // 
  // for(int i=0;i<ntstar;i++){
  //   index2[i] = i;
  //   probs[i] = 1.0/ntstar;
  // }
  Rcpp::IntegerVector index;
  Rcpp::IntegerVector frame = Rcpp::Range(0, ntstar-1);
  Rcpp::NumericVector wts = Rcpp::rep(1.0/ntstar,ntstar);//Rcpp::runif(n, 0.0, 1.0);
  index = Rcpp::RcppArmadillo::sample(frame, nsim, true, wts / Rcpp::sum(wts));
  //std::cout << "1\n";
  std::cout << "7\n";
  
  arma::cube currentmeans(n,k,K);
  arma::mat currentspecificmean;
  arma::mat currentmeansnoint;// = calc_mean_noint(tstar.row(index[0]).t(),currentbeta);
  std::cout << "8\n";
  
  for(int j=0;j<K;j++){
    currentmeans.slice(j) = calc_mean(tstar.row(index[0]).t(),currentbeta,currentlambda.col(j));
  }
  
  for(int i=0;i<nsim;i++){
    currentmeansnoint = calc_mean_noint(tstar.row(index[i]).t(),currentbeta);
    std::cout << "9\n";
    
    //update tstar
    currentzeta = sample_zeta(y, currentpi, currentmeans, currentSigma);
    std::cout << "10\n";
    
    currentpi = sample_pi(currentzeta,a);
    std::cout << "11\n";
    
    currentspecificmean = calc_specific_mean(tstar.row(index[i]).t(),currentbeta,currentlambda,currentzeta);
    std::cout << "12\n";
    
    for(int j=0;j<K;j++){
      currentSigma.slice(j) = sample_Sigma(subset(y,currentzeta,j),subset(currentspecificmean,currentzeta,j),D,d,(which(currentzeta,j)).size());
    }
    std::cout << "13\n";
    
    for(int j=0;j<K;j++){
      currentlambda.col(j) = sample_lambda(subset(y,currentzeta,j),subset(currentmeansnoint,currentzeta,j),currentSigma.slice(j),lm,lcov);
    }
    std::cout << "14\n";
    
    //std::cout << "4\n";
    
    for(int j=0;j<K;j++){
      currentmeans.slice(j) = calc_mean(tstar.row(index[i]).t(),currentbeta,currentlambda.col(j));
    }
    std::cout << "15\n";
    
    
    currentbeta = sample_beta(y,tstar.row(index[i]).t(),currentbeta,currentlambda,
                              currentzeta,currentSigma,currentpi,propcov,bm,bcov);
      
    //std::cout << "5\n";
    if((i<burn) && (i>20) && (i%20==0)){
      propcov = cov(beta.rows(0,i-1));
    }
    //std::cout << "6\n";
    
    beta.row(i) = currentbeta.t();
    zeta.row(i) = currentzeta.t();
    //std::cout << "7\n";
    std::cout << "16\n";
    
    Sigma.slice(i) = currentSigma.slice(0);
    pi.row(i) = currentpi.t();
    
    for(int j=0;j<K;j++){
      lambda.slice(j).row(i) = currentlambda.col(j).t();
    }
    std::cout << "17\n";
    
    if(i % 1000==0){
      std::cout << "i= " << i << "\n";
    }
  }
  
  return List::create(
    Named("beta") = beta,
    Named("Sigma") = Sigma,
    Named("lambda") = lambda,
    Named("pi") = pi,
    Named("zeta") = zeta,
    Named("propcov") = propcov);
} 

