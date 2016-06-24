#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(N);
  PARAMETER_VECTOR(mu);
  PARAMETER(logsd);
  PARAMETER_VECTOR(x);
  vector<Type> lambda = exp(x);
  
  Type ans = 0;

  ans -= dnorm(x, mu, exp(logsd), true).sum();
  ans -= dpois(N, lambda.sum(), true);

  return ans;
}

