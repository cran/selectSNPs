#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector EScore(NumericVector f){
  int n = f.size();
  NumericVector entropy(n);
  
  for (int i=0; i<n;  ++i){
    double p = f[i];
    if(p==0){
      p=0.000001;
    }
    double q = 1-p;
    entropy[i]=(-1)*(p*log2(p)+q*log2(q));
  }
  return(entropy);
}


// [[Rcpp::export]]
double UScore(NumericVector x){
  int n = x.size();
  
  double sumx2 = 0;
  for (int i=0; i<(n-1); ++i){
    sumx2 += pow((x[(i+1)]-x[i]),2.0);
  }
  double spacing = (x[(n-1)]-x[0])/(n-1);
  double sumz2 = (n-1)*pow(spacing,2.0);
    
  double uscore = sqrt(sumz2/sumx2);
  return(uscore);
}
  
// [[Rcpp::export]]
double localScore(NumericVector f, NumericVector x, 
                  double w1 = 0.5, double w2 = 0.5, 
                  double t1 = 1, double t2 = 1){
  NumericVector ent = EScore(f);
  int n = ent.size();
  double mu = 0;
  for (int i=0; i<n; ++i){
    mu += ent[i]/n;
  }
  
  double u = UScore(x); 
  
  w1 = w1/(w1+w2);
  w2 = w2/(w1+w2);
  double local = (w1*pow(mu,t1)) + (w2*pow(u,t2));
  return local;
}


