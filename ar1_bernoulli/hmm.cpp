// Inference in AR1-bernulli using HMM filter
#include <TMB.hpp>

/* Matrix vector multiply */
template<class Type>
vector<Type> multiply(matrix<Type> x, vector<Type> y){
  return atomic::matmul(x, matrix<Type>(y.matrix())).vec();
}

/* yobs  */
template<class Type>
void update(vector<Type> &px, vector<Type> xm, matrix<Type> P, Type &nll, int yobs){
  // Update joint distribution (J) of state and measurement, and
  // extract 'yobs'-slice:
  vector<Type> Ppx = multiply(P, px);
  vector<Type> Ppx1 = Ppx * xm;
  vector<Type> Ppx0 = Ppx * (Type(1) - xm);
  // Joint: cbind( Ppx * xm , Ppx * (1-xm) )
  // Integrated slice is yobs-likelihood:
  Type Ly = (yobs == 1 ? Ppx1.sum() : Ppx0.sum() );
  nll -= log(Ly);
  // Normalized slice is updated px:
  px = (yobs == 1 ? Ppx1 : Ppx0 );
  px = px / Ly;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Grid of the state space
  DATA_VECTOR(grid);
  DATA_IVECTOR(obs);
  PARAMETER(phi);
  PARAMETER(logSigma);
  Type sigma= exp(logSigma);
  // Construct discrete generator
  // Note: We formulate the dynamics in 'probability-domain'
  //       I.e.: grid is a discrete unit interval
  int n = grid.size() - 1;
  matrix<Type> P(n,n);
  Type sigma_increm = sigma * sqrt(1.0 - phi*phi);
  vector<Type> xm = Type(.5) * (grid.head(n) + grid.tail(n));
  vector<Type> xm_logit = logit(xm);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      P(i,j) = 
	pnorm(logit(grid[i+1]), phi * logit(xm[j]), sigma_increm) -
	pnorm(logit(grid[i])  , phi * logit(xm[j]), sigma_increm);
  // Normalize
  vector<Type> cs = P.colwise().sum();
  for(int i=0; i<n; i++) P.col(i) /= cs[i];
  // Likelihood
  vector<Type> px(n);
  for(int i=0; i<n; i++)
    px(i) =
      pnorm(logit(grid[i+1]), Type(0), sigma) -
      pnorm(logit(grid[i])  , Type(0), sigma);
  px = px / px.sum();
  Type nll = 0;
  for(int k=0; k<obs.size(); k++){
    update(px, xm, P, nll, obs(k));
  }
  return nll;
}
