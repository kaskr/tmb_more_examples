## AR1-bernoulli

obs(t) | p ~ Binomial(1, p(t))

logit(p(t)) ~ AR1(phi, sigma)

Compare accuracy of different methods:

1. HMM filter ('exact')
2. Laplace approximation
3. Transformed Laplace approximation
