rm(list=ls());gc()

library(pomp)
library(ggplot2)

simulate(t0=0, times=1:20,
         params=c(r=1.2,K=200,sigma=0.1,N_0=50),
         rinit=function (N_0, ...) {
           c(N=N_0)
         },
         rprocess=discrete_time(
           function (N, r, K, sigma, ...) {
             eps <- rnorm(n=1,mean=0,sd=sigma)
             c(N=r*N*exp(1-N/K+eps))
           },
           delta.t=1
         )
) -> sim1

# spy(sim1)
# plot(sim1)

ggplot(data=as.data.frame(sim1),aes(x=time,y=N))+
  geom_line() +
  theme_bw()

Csnippet("
  pop = rpois(b*N);  
  ") -> rmeas

Csnippet("
  N = N_0;
  ") -> rinit

Csnippet("
  double eps = rnorm(0,sigma);
  N = r*N*exp(1-N/K+eps);
") -> rickstepC

Csnippet("
  double dW = rnorm(0,sqrt(dt));
  N += r*N*(1-N/K)*dt+sigma*N*dW;
") -> vpstepC

Csnippet("
  lik = dpois(pop,b*N,give_log);
") -> dmeas

data("parus")
