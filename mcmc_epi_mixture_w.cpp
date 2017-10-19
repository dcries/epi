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
              arma::mat lambda, arma::ivec zeta, arma::vec pi, arma::vec priormean, arma::mat priorcov,arma::vec weights){
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
    
    ll += weights[i]*( dmvnrm_arma(y.row(i).t(),means.row(i),Sigma.slice(zeta[i]),true));
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

arma::ivec sample_zeta(arma::mat y, arma::vec pi, arma::cube meanmat, arma::cube Sigma,arma::vec weights){
  int k = pi.size();
  int n = meanmat.slice(0).n_rows;
  arma::vec probs(k);

  arma::ivec out(n);
  arma::ivec temp(1);
  Rcpp::IntegerVector frame = Rcpp::Range(0,k-1);
  
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      // if((i==1 )| (i==15) | (i==600) | (i== 2000)){
      //   std::cout << "for i=" << i << " log pi=" << log(pi[j]) << "\n";
      //   std::cout << "dmvnorm = " << dmvnrm_arma(y.row(i).t(),meanmat.slice(j).row(i),Sigma.slice(j),true) << "\n";
      // }
      probs[j] = weights[i]*(log(pi[j]) + dmvnrm_arma(y.row(i).t(),meanmat.slice(j).row(i),Sigma.slice(j),true));
    }
    probs = exp(probs-log(sum(exp(probs))));
    temp = Rcpp::RcppArmadillo::sample(frame, 1, true, probs);
    out[i] = temp[0];
  }
  return out;
}

arma::mat calc_pmat(arma::mat y, arma::vec pi, arma::cube meanmat, arma::cube Sigma){
  int k = pi.size();
  int n = meanmat.slice(0).n_rows;
  arma::vec probs(k);
  arma::mat out(n,k);
  
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      // if((i==1 )| (i==15) | (i==600) | (i== 2000)){
      //   std::cout << "for i=" << i << " log pi=" << log(pi[j]) << "\n";
      //   std::cout << "dmvnorm = " << dmvnrm_arma(y.row(i).t(),meanmat.slice(j).row(i),Sigma.slice(j),true) << "\n";
      // }
      probs[j] = log(pi[j]) + dmvnrm_arma(y.row(i).t(),meanmat.slice(j).row(i),Sigma.slice(j),true);
    }
    out.row(i) = exp(probs-log(sum(exp(probs)))).t();
    
  }
  return out;
}

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

arma::mat sample_Sigma(arma::mat y, arma::mat meanmat, arma::mat D, double d, int n,arma::vec weights){
  int k = y.n_cols;
  arma::mat total=arma::zeros(k,k);

  for(int i=0;i<n;i++){
      total += weights[i]*((y.row(i)-meanmat.row(i)).t()*(y.row(i)-meanmat.row(i)));

  }
  arma::mat out = rinvwish(1,n+d,D+total);

  return out;
}

arma::vec sample_lambda(arma::mat y, arma::mat meannoint, arma::mat Sigma, arma::vec priormean, arma::mat priorcov,arma::vec weights){
  int n = y.n_rows;
  int k = y.n_cols;
  arma::vec out;
  
  arma::vec total=arma::zeros(k);
  arma::vec total2=arma::zeros(k);
  
  for(int i=0;i<n;i++){
    total += (weights[i]*arma::inv(Sigma))*(y.row(i)-meannoint.row(i)).t();
    total2 += (arma::inv(Sigma))*(y.row(i)-meannoint.row(i)).t();
    
  }
  //std::cout << "total weights = " << total << "\n" << "total2 = " << total2 << "\n";
  //arma::vec ybar = (mean(y-meannoint,0)).t();
  arma::mat postcov = arma::inv(arma::inv(priorcov)+n*arma::inv(Sigma));
  arma::vec postmean = postcov*(arma::inv(priorcov)*priormean + total2);
  arma::vec postmean2 = postcov*(arma::inv(priorcov)*priormean + total2);
  
  //std:: cout << "postmean = " << postmean << "\n" << "postmean2 = " << postmean2 << "\n";
  out = mvrnormArma(1,postmean,postcov).row(0).t();
  
  return out;
}

arma::vec sample_beta(arma::mat y, arma::vec tstar, arma::vec beta, arma::mat lambda,
                      arma::ivec zeta, arma::cube Sigma,
                      arma::vec pi, arma::mat propcov, arma::vec bm, arma::mat bcov,double mult,arma::vec weights){
  int p = beta.size();
  arma::vec propb(p);
  double logr = 0.0;
  double logrprop = 0.0;
  double logacceptprob;
  //int k = pi.size();
  
  
  propb = (mvrnormArma(1,beta,mult*propcov)).row(0).t();
  //tstar.row(index[i]).t()
  
  logrprop = logl_b(y,tstar,propb,Sigma,lambda,zeta,pi,bm,bcov,weights); 
  //std::cout << "2\n";
  
  logr = logl_b(y,tstar,beta,Sigma,lambda,zeta,pi,bm,bcov,weights);
  //std::cout << "3\n";
  
  
  
  logacceptprob = logrprop - logr;
  if(log(R::runif(0,1))<logacceptprob){
    beta = propb;
  }
  
  return beta;
}

double calc_full_ll(arma::mat y, arma::ivec zeta, arma::mat specificmeanmat, arma::cube Sigma){
  int n = y.n_rows;
  double ll=0.0;
  for(int i=0;i<n;i++){
    ll += dmvnrm_arma(y.row(i).t(),specificmeanmat.row(i),Sigma.slice(zeta[i]),true);
  }
  return ll;
}

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
List mcmc_epi_mixture_w(arma::mat y, arma::mat tstar, List start, List prior, arma::vec weights,
                      int K, int nsim,int burn, int thin=1, double mult=0.1){
  
  //std::cout << "1\n";
  arma::vec currentbeta         = as<arma::vec>(start["currentbeta"]);
  arma::ivec currentzeta         = as<arma::ivec>(start["currentzeta"]);
  arma::vec currentpi         = as<arma::vec>(start["currentpi"]);
  arma::mat currentlambda             = as<arma::mat>(start["currentlambda"]);
  //std::cout << "2\n";
  
  arma::mat Sigmadiag        = as<arma::mat>(start["Sigmadiag"]);
  arma::mat propcov             = as<arma::mat>(start["propcov"]);
  
  //std::cout << "3\n";
  
  arma::vec bm            = as<arma::vec>(prior["bm"]);
  arma::mat bcov            = as<arma::mat>(prior["bcov"]);
  arma::vec lm            = as<arma::vec>(prior["lm"]);
  arma::mat lcov            = as<arma::mat>(prior["lcov"]);
  arma::mat D               = as<arma::mat>(prior["D"]);
  double d            = as<double>(prior["d"]);
  arma::vec a            = as<arma::vec>(prior["a"]);
  //std::cout << "4\n";
  
  
  int n = y.n_rows;
  int p = bm.size();
  int k = y.n_cols;
  int ntstar = tstar.n_rows;
  int keep = (nsim-burn)/thin;
  
  arma::cube currentSigma = arma::zeros(k,k,K);
  for(int i=0;i<K;i++){
    currentSigma.slice(i).diag() = Sigmadiag.col(i);
  }
  //std::cout << "5\n";
  
  
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
  
  //std::cout << "6\n";
  
  //for dic
  arma::vec full_ll(keep);
  double penalty;
  double dic;
  arma::ivec medianzeta(n);
  arma::cube meanSigma(k,k,K);
  arma::mat meanlambda;
  arma::mat finalmeanmat;
  
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
  //std::cout << "7\n";
  
  arma::cube currentmeans(n,k,K);
  arma::mat currentspecificmean;
  arma::mat currentmeansnoint;// = calc_mean_noint(tstar.row(index[0]).t(),currentbeta);
  int count;
  int ind = 0;
  
  for(int j=0;j<K;j++){
    currentmeans.slice(j) = calc_mean(tstar.row(index[0]).t(),currentbeta,currentlambda.col(j));
  }
  
  for(int i=0;i<nsim;i++){
    currentmeansnoint = calc_mean_noint(tstar.row(index[i]).t(),currentbeta);
    //std::cout << "9\n";
    
    //update tstar
    currentzeta = sample_zeta(y, currentpi, currentmeans, currentSigma,weights);
    //std::cout << "10\n";
    
    currentpi = sample_pi(currentzeta,a);
    //std::cout << "11\n";
    
    currentspecificmean = calc_specific_mean(tstar.row(index[i]).t(),currentbeta,currentlambda,currentzeta);
    //std::cout << "12\n";
    
    for(int j=0;j<K;j++){
      currentSigma.slice(j) = sample_Sigma(subset(y,currentzeta,j),subset(currentspecificmean,currentzeta,j),D,d,(which(currentzeta,j)).size(),subset(weights,currentzeta,j));
    }
    
    //std::cout << "13\n";
    
    for(int j=0;j<K;j++){
      if((which(currentzeta,j)).size() > 0){
        currentlambda.col(j) = sample_lambda(subset(y,currentzeta,j),subset(currentmeansnoint,currentzeta,j),currentSigma.slice(j),lm,lcov,subset(weights,currentzeta,j));
      }
    }
    
    //std::cout << "14\n";
    
    
    
    currentbeta = sample_beta(y,tstar.row(index[i]).t(),currentbeta,currentlambda,
                              currentzeta,currentSigma,currentpi,propcov,bm,bcov,mult,weights);
    //std::cout << "15\n";
    
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
    //std::cout << "8\n";
    
    if((i >= burn) && (i%thin==0)){
      pmat.slice(ind) = calc_pmat(y, currentpi, currentmeans, currentSigma);
      full_ll[ind] = calc_full_ll(y,currentzeta,currentspecificmean,currentSigma);
      beta.row(ind) = currentbeta.t();
      zeta.row(ind) = currentzeta.t();
      //Sigma.slice(i) = currentSigma.slice(0);
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
    //std::cout << "17\n";
    
    if(i % 1000==0){
      std::cout << "i= " << i << "\n";
    }
  }
  
  meanlambda = calc_meanlambda(lambda);
  finalmeanmat = calc_specific_mean(mean(tstar,0).t(),mean(beta,0).t(),meanlambda,median(zeta,0).t());
  meanSigma = calc_meanSigma(sds,covmat);
  penalty = calc_full_ll(y,median(zeta,0).t(),finalmeanmat,meanSigma);
  
  dic = -4*mean(full_ll) + 2*penalty;
  
  return List::create(
    Named("beta") = beta,
    //Named("Sigma") = Sigma,
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

