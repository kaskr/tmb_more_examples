#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs);
  PARAMETER(phi);
  PARAMETER(logSigma);
  PARAMETER_VECTOR(u);
  using namespace density;
  Type sigma= exp(logSigma);
  Type f = 0;
  f += SCALE(AR1(phi), sigma)(u);
  vector<Type> p = invlogit(u);
  f -= dbinom(obs, Type(1), p, true).sum();
  return f;
}
