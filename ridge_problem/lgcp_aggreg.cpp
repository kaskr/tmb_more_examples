#include <TMB.hpp>
//#include "../atomic_sort.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(N);
  PARAMETER_VECTOR(mu);
  PARAMETER(logsd);
  PARAMETER_ARRAY(x); // ngroup x nobs
  int ngroup=x.rows();
  int nobs=x.cols();

  array<Type> lam=x;
  lam=exp(x);
  Type ans=0;
  for(int i=0;i<nobs;i++){
    ans-=dnorm(vector<Type>(x.col(i)),mu,exp(logsd),true).sum();
    ans-=dpois(N(i),lam.col(i).sum(),true);
    //ans-=dpois(N(i),max(vector<Type>(lam.col(i))),true);
    //ans-=dpois(N(i),atomic::max(vector<Type>(lam.col(i))),true);
  }
  
  return ans;

}

