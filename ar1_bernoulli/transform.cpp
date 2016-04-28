#include <TMB.hpp>

template<class Type> 
Type mydbinom(Type k, Type size, Type prob, int give_log=0)
{
	Type logres = k*log(prob)+(size-k)*log(1-prob);
	if(!give_log) return exp(logres);
	else return logres;
}


template<class Type>
struct gauss_bernoulli {
  typedef Type scalartype;
  Type mu, sigma, x;
  Type operator()(Type u){
    Type p = invlogit(u);
    Type ans = 
      dnorm(u, mu, sigma, false) * 
      mydbinom(x, Type(1), p, false);
    return ans;
  }
};


/*
 F(x) := integrate(f, -inf, x)
 F0   := integrate(f, -inf, inf)

 CDF of f
 G(x) := F(x) / F0

 Solve G(x)=y for som 0<y<1.
 H(x) := G(x)-y

 Newton:
 H(x0) + H'(x0)(x-x0) = 0
 x = -H(x0)/H'(x0) + x0

 ====>

 H(x) = G(x)-y = F(x)/F0 - y
 H'(x) = G'(x) = f(x)/F0

 =======

 G(G^-1(y)) = y
 G'(G^-1(y)) * (G^-1(y))' = 1
 (G^-1(y))' = 1 / G'(G^-1(y))
            = 1 / (f(x)/F0)
	    = F0/f(x)
*/
template<class distribution>
struct quantile_function {
  typedef typename distribution::scalartype Type;
  distribution f; // 1D un-normalized density
  Type deriv;
  Type operator()(Type y){
    int niter = 10;
    Type inf = 20;
    using romberg::integrate;
    Type F0 = integrate(f, -inf, inf, 9 );
    Type x = 0; // Initial guess
    for(int i=0; i<niter; i++){
      Type Fx = integrate(f, -inf, x, 9 );
      Type dx = -( Fx/F0 - y) / ( f(x) / F0 );
      x = x + dx;
    }
    deriv = F0 / f(x);
    return x;
  }
};

template<class Type>
vector<Type> my_transform(vector<Type> parms){
  Type u=parms[0], mu = parms[1], sigma = parms[2], x = parms[3], f=0;
  vector<Type> ans(2);
  gauss_bernoulli<Type> gb = {mu, sigma, x};
  quantile_function<gauss_bernoulli<Type> > qf = {gb};
  f -= dnorm(u, Type(0.0), Type(1.0), true);
  u = pnorm(u, Type(0.0), Type(1.0));
  u = qf(u);
  f -= log(qf.deriv);
  ans(0) = u; ans(1) = f;
  return ans;
}
REGISTER_ATOMIC(my_transform)
template<class Type>
void my_transform(Type &u, Type mu, Type sigma, Type x, Type &f, Type uaux, Type tiny){
  vector<Type> parms(4);
  parms << u, mu, sigma, x;
  vector<Type> ans = my_transform(parms);
  f -= dnorm(uaux, ans(0), tiny, true); // Softlink 'copy'
  f += ans(1);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(do_transform);
  DATA_VECTOR(obs);
  PARAMETER(phi);
  PARAMETER(logSigma);
  PARAMETER_VECTOR(u);     // Actual states
  PARAMETER_VECTOR(uaux);  // Transformed states
  DATA_SCALAR(tiny);
  using namespace density;
  Type sigma = exp(logSigma);
  int n = obs.size();
  Type f = 0;

  for(int i=0; i<n; i++){
    switch (do_transform) {
    case 0: // No transform
      f -= dnorm(uaux(i), Type(0), Type(1), true);
      break;
    case 1: // Forward recursive transform
      if(i==0) {
	my_transform(uaux(i), Type(0), sigma, obs(i), f, u(i), tiny);
      } else {
	my_transform(uaux(i), phi*u(i-1), sigma * sqrt(1.0 - phi*phi),
		     obs(i), f, u(i), tiny);
      }
      break;
    case 2: // Full conditionals
      if(i==0) {
	my_transform(uaux(i),
		     phi*u(i+1),
		     sigma * sqrt(1.0 - phi*phi),
		     obs(i), f, u(i), tiny);
      }
      else if
	(i==n-1) {
	my_transform(uaux(i),
		     phi*u(i-1),
		     sigma * sqrt(1.0 - phi*phi),
		     obs(i), f, u(i), tiny);
      }
      else {
	my_transform(uaux(i),
		     phi / (Type(1)+phi*phi) * (u(i-1) + u(i+1)),
		     sigma * sqrt((1.0 - phi*phi) / (1.0 + phi*phi)),
		     obs(i), f, u(i), tiny);
      }
      break;
    default:
      error("Invalid do_transform");
      break;
    }
  }

  f += SCALE(AR1(phi), sigma)(u);
  vector<Type> p = invlogit(u);
  f -= dbinom(obs, Type(1), p, true).sum();
  return f;
}
