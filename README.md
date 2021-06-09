# fitexamples

## software packages

* [Nimble](https://r-nimble.org) has [HTML docs](https://r-nimble.org/html_manual/cha-welcome-nimble.html) here.
* [pomp](https://kingaa.github.io/pomp/index.html) is here.
* [mcstate](https://mrc-ide.github.io/mcstate/) does SMC and PMCMC for state space models.
* [Stan](https://mc-stan.org) is here.
* [Greta](https://greta-stats.org) is here, fits graphical models.
* [BayesianTools](https://github.com/florianhartig/BayesianTools) just a collection of MCMC and SMC algorithms.
* [stemr](https://github.com/fintzij/stemr): Baysian Inference for Stochastic Epidemic Models via the Linear Noise Approximation
* [LibBi](http://libbi.org) state space modeling and inference, with R package [rbi](https://CRAN.R-project.org/package=rbi).
* [Torsten](https://github.com/metrumresearchgroup/Torsten) provides helper functions to fit PK/PD models in Stan

## examples, etc.

* Nimble, Stan, BUGS code for discrete time stochastic epi models: https://github.com/wzmli/hybridx
* Stan ecology examples: https://stanecology.github.io
* Intro to PMCMC for R via LibBi: https://github.com/akira-endo/Intro-PMCMC 
* Estimating HIV incidence via Stan: https://github.com/frbrz25/Thesis_Codes
* TSIR for CHIK/ZIKV in Stan: https://github.com/jriou/comp_chik_zika
* More TSIR in Stan: https://github.com/jriou/bora
* List of likelihood-free examples, mostly with code, including some real interesting ones (ERGMs): https://github.com/dennisprangle/LFexamples
* chaotic voles: https://github.com/mfasiolo/volesModel uses the [synlik](https://CRAN.R-project.org/package=synlik) package to fit data with synthetic likelihood MCMC and accompanies this paper https://arxiv.org/abs/1511.02644
* pomp and LiBBi (via rbi) for simple flu data: https://github.com/sbfnk/inference_pomp_rbi 

## new research, miscellanea

* discontinuous HMC (for discrete parameters): https://github.com/aki-nishimura/discontinuous-hmc
* general profiling for DDEs: https://github.com/cran/gpDDE 
* profile likelihood for DDEs: https://github.com/VulpeculaZ/DDE
* inference for stochastic delays: https://github.com/cbskust/DDE_BD
* pseudo-marginal inference in Julia via LNA/CLE: https://github.com/davidwarne/Warne2019_GuideToPseudoMarginal
* fits SDEs with time-varying parameters: https://github.com/TheoMichelot/smoothSDE

## specific papers
* (Oliver Maclaren) Likelihood based methods for misspecified epidemic models: https://arxiv.org/abs/2103.16058
* (Simon Wood) fit COVID-19 infection trajectories with GAMs https://onlinelibrary.wiley.com/doi/full/10.1111/biom.13462
* (Simon Wood) fit COVID-19 ODE models with spline beta(t) function https://www.medrxiv.org/content/10.1101/2021.02.03.21251112v2
