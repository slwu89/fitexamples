## (c) Simon N Wood
## occupancy likelihood corrected
## Version allowing the generation time to be manipulated while keeping the
## disease progression times fixed. 


lfun <- function(theta,link) {
## applies link to theta...
  if (is.matrix(link)) {
    a1 <- link[,2];a2 <- link[,3];link <- link[,1]
  } else {
    a1 <- rep(0,length(theta));a2 <- rep(1,length(theta))
  }
  ii <- link==1
  if (length(ii)) theta[ii] <- log(theta[ii])
  ii <- link==2
  if (length(ii)) theta[ii] <- binomial()$linkfun((theta[ii]-a1[ii])/(a2[ii]-a1[ii]))
  theta
}

ilink <- function(eta,link) {
## apply link inverse and compute dtheta/deta - eta is the unconstrained parameter
  if (is.matrix(link)) {
    a1 <- link[,2];a2 <- link[,3];link <- round(link[,1])
  } else {
    a1 <- rep(0,length(eta));a2 <- rep(1,length(eta))
  }
  theta <- eta
  dtheta <- eta*0+1
  ii <- link==1
  if (length(ii)) {
    dtheta[ii] <- theta[ii] <- exp(eta[ii])
  }  
  ii <- link==2
  if (length(ii)) {
    theta[ii] <- binomial()$linkinv(eta[ii])*(a2[ii]-a1[ii]) + a1[ii]
    dtheta[ii] <- binomial()$mu.eta(eta[ii])*(a2[ii]-a1[ii])
  }  
  list(theta=theta,dtheta=dtheta)
}

ll <- function(theta,gamma,S0,iaux,daux,dat,Sb=NULL,link=rep(0,length(theta)),plot=FALSE,
               step=1,nbk=rep(2,5),getR=FALSE,lab="",xlab="day") {
## log likelihood for the rep41 covid SEIR+hospital model.  nbk is a vector of overdispersion parameters.
## rep41 uses neg bin kappa = 2 for all - this is very hard to justify for deaths.
  nb <- iaux[1];np <- iaux[2];nc <- iaux[3];ns <- iaux[4]
  nx <- ns*nc*(1+np);
  nyo <- 14
  ny <- nyo*(1+np)
  x <- rep(0,nx)
  for (i in 1:nc) x[(i-1)*ns+1] <- S0[i]
  ii <- 1:length(gamma)+(nc-1)^2
  daux[ii] <- gamma
  n <- max(dat$day)
  th <- ilink(theta,link)
  thetaw <- theta ## original scale
  theta <- th$theta ## parameters on dynamic model specification scale
  if (any(!is.finite(theta))) return(list(l=NA,dl=NA))
  vm <- .C("RK4",yout = as.double(rep(0,ny*n)),xout=as.double(rep(0,nx*n)),x=as.double(x),
          theta=as.double(theta),h=as.double(1/step),t0=as.double(0),n=as.integer(n),nx=as.integer(nx),ny=as.integer(ny),
          ostep =as.integer(step),work=as.double(rep(0,5*nx)),iaux=as.integer(iaux),daux=as.double(daux))
  yout <- matrix(vm$yout,ny,n)

  ## hospital deaths...
  
  k <- nbk[4]
  y <- dat$death3;
  ynb <- matrix(NA,max(dat$day),7);munb <- matrix(NA,max(dat$day),9)
  ii <- !is.na(y); y<- y[ii];ii <- dat$day[ii]
  mu <- yout[4,ii]
  ynb[ii,4] <- y;munb[ii,4] <- mu
  n.death <- sum(mu)
  ld <- sum(k*log(k/(k+mu)) + y*log(mu/(mu+k)) + lgamma(k+y) - lgamma(k) - lgamma(y+1))
  dld.dmu <- y/mu - (y+k)/(mu+k)
  dld <- drop(yout[1:np*nyo+4,ii]%*%dld.dmu)
  if (plot&&length(y)) {
    ylim=c(range(y,mu));xlim <- range(ii)
    maj <- 2.2;mas <- .7
    pch=19;pex=.4
    plot(ii,y,ylab="",col="grey",ylim=ylim,xlim=xlim,xlab="",pch=pch,cex=pex);
    mtext("deaths",2,maj,cex=mas);mtext(xlab,1,maj,cex=mas)
    lines(ii,mu)
    text(220,ylim[2]*.86,lab)
  }
  
  ## care home deaths...

  k <- nbk[5]
  y <- dat$death2-dat$death3;ii <- !is.na(y); y<- y[ii];ii <- dat$day[ii]
  mu <- yout[5,ii];n.death <- n.death+ sum(mu)
  ynb[ii,5] <- y;munb[ii,5] <- mu
  lg <- sum(k*log(k/(k+mu)) + y*log(mu/(mu+k)) + lgamma(k+y) - lgamma(k) - lgamma(y+1))
  dlg.dmu <- y/mu - (y+k)/(mu+k)
  dlg <- drop(yout[1:np*nyo+5,ii]%*%dlg.dmu)
  if (plot&&length(y)) {
    points(ii,y,col=2,pch=pch,cex=pex);lines(ii,mu,col=4)
  }
  
  ## occupancy...

  k <- 100000

  ## total occupancy... 
  y0 <- dat$phe_patients;i0 <- dat$day
  ii <- !is.na(y0);y <- sum(y0[ii]);i0 <- i0[ii]
  mu <- sum(yout[2,i0])

  ## log likelihood total occupancy....
  lot <- sum(k*log(k/(k+mu)) + y*log(mu/(mu+k)) + lgamma(k+y) - lgamma(k) - lgamma(y+1))
  dlot.dmu <- y/mu - (y+k)/(mu+k)
  dlot <- rowSums(yout[1:np*nyo+2,i0]) * dlot.dmu 
  
  ## get change in occupancy over the day - difference between in and outflow
  k <- nbk[2] ## overdispersion parameter
  y <- diff(dat$phe_patients);ii <- !is.na(y); y <- y[ii];ii <- dat$day[ii]
  mu1 <- yout[11,ii] # general ward inflow
  mu2 <- yout[12,ii] # general ward outflow
  wi <- rep(1:ceiling(length(y)/7),each=7)[1:length(y)]
  y <- tapply(y,wi,sum) ## aggregate differences by week to avoid unmodelled weekly discharge cycle
  mu1 <- tapply(mu1,wi,sum);mu2 <- tapply(mu2,wi,sum)
  
  iw <- ii[1]+1:length(y)*7 - 4 ## output to allow over-dispersion update
  ynb[iw,6] <- y;munb[iw,6] <- mu1;munb[iw,7] <- mu2; 
  alpha <- y-mu1+mu2;sig <- mu1+mu2

  lo <- sum(-alpha^2/(2*k*sig) - log(k*sig)/2) ## log likelihood occupancy change
  dlo.dmu1 <- alpha/(k*sig) + alpha^2/(2*k*sig^2) - .5/sig
  dlo.dmu2 <- -alpha/(k*sig) + alpha^2/(2*k*sig^2) -.5/sig
  dlo <- t(apply(yout[1:np*nyo+11,ii],1,function(y,wi) tapply(y,wi,sum),wi=wi))%*%dlo.dmu1 +
	 t(apply(yout[1:np*nyo+12,ii],1,function(y,wi) tapply(y,wi,sum),wi=wi))%*%dlo.dmu2 
  lo <- lo + lot;
  dlo <- drop(dlo) + dlot
  
  ## occupancy itself for plotting and return..
  y <- dat$phe_patients;ii <- !is.na(y); y<- y[ii];ii <- dat$day[ii]
  mu <- yout[2,ii]
  ynb[ii,2] <- y;munb[ii,2] <- mu

  if (plot&&length(y)) { 
    plot(ii,y,col="grey",ylim=range(c(y,mu)),xlim=xlim,ylab="",xlab="",pch=pch,cex=pex)
    mtext("occupancy",2,maj,cex=mas);mtext(xlab,1,maj,cex=mas)
    lines(ii,mu)
  }
  
  ## ICU occupancy...
  k <- 20000

  y0 <- dat$phe_occupied_mv_beds;i0 <- dat$day
  ii <- !is.na(y0);y <- sum(y0[ii]);i0 <- i0[ii]
  mu <- sum(yout[3,i0])
  ## log likelihood total icu occupancy
  lut <- sum(k*log(k/(k+mu)) + y*log(mu/(mu+k)) + lgamma(k+y) - lgamma(k) - lgamma(y+1))
  dlut.dmu <- y/mu - (y+k)/(mu+k)
  dlut <- rowSums(yout[1:np*nyo+3,i0]) * dlut.dmu

  ## get likelihood for icu occupancy changes...
  k <- nbk[3] ## overdispersion parameter
  y <- diff(dat$phe_occupied_mv_beds);ii <- !is.na(y); y<- y[ii];ii <- dat$day[ii]
  mu1 <- yout[13,ii] # ICU inflow
  mu2 <- yout[14,ii] # ICU outflow
  wi <- rep(1:ceiling(length(y)/7),each=7)[1:length(y)]
  y <- tapply(y,wi,sum) ## aggregate differences by week to avoid unmodelled weekly discharge cycle
  mu1 <- tapply(mu1,wi,sum);mu2 <- tapply(mu2,wi,sum)
  
  iw <- ii[1]+1:length(y)*7 - 4
  ynb[iw,7] <- y;munb[iw,8] <- mu1;munb[iw,9] <- mu2;

  alpha <- y-mu1+mu2;sig <- mu1+mu2
  lu <- sum(-alpha^2/(2*k*sig) - log(k*sig)/2) ## log likelihood - icu occupancy change
  dlu.dmu1 <- alpha/(k*sig) + alpha^2/(2*k*sig^2) - .5/sig
  dlu.dmu2 <- -alpha/(k*sig) + alpha^2/(2*k*sig^2) -.5/sig
  dlu <- t(apply(yout[1:np*nyo+13,ii],1,function(y,wi) tapply(y,wi,sum),wi=wi))%*%dlu.dmu1 +
	 t(apply(yout[1:np*nyo+14,ii],1,function(y,wi) tapply(y,wi,sum),wi=wi))%*%dlu.dmu2 
  #dlu <- drop(yout[1:np*nyo+13,ii]%*%dlu.dmu1)+drop(yout[1:np*nyo+14,ii]%*%dlu.dmu2)
  lu <- lu + lut
  dlu <- dlu + dlut
  # ICU occupancy itself for plotting and return
  y <- dat$phe_occupied_mv_beds;ii <- !is.na(y); y<- y[ii];ii <- dat$day[ii]
  mu <- yout[3,ii]
  ynb[ii,3] <- y;munb[ii,3] <- mu
  if (plot&&length(y)) {
    points(ii,y,col=2,pch=pch,cex=pex)#xlab="day",ylab="ICU",col="grey",ylim=range(c(y,mu)));
    lines(ii,mu,col=4)
  }
  
  ## hospitalization...

  y <- dat$phe_admissions; ii <- !is.na(y); y<- y[ii]
  ii <- dat$day[ii] ## cols of yout corresponding to observed data
  k <- nbk[1] ## theta parameter
  mu <- yout[1,ii]
  ynb[ii,1] <- y;munb[ii,1] <- mu
  lh <- sum(k*log(k/(k+mu)) + y*log(mu/(mu+k)) + lgamma(k+y) - lgamma(k) - lgamma(y+1))
  dlh.dmu <- y/mu - (y+k)/(mu+k)
  dlh <- drop(yout[1:np*nyo+1,ii]%*%dlh.dmu)
  #if (plot) par(mfrow=c(2,3));
  if (plot&&length(y)) {
    plot(ii,y,col="grey",ylim=range(c(y,mu)),xlim=xlim,ylab="",xlab="",pch=pch,cex=pex);
    mtext("admissions",2,maj,cex=mas);mtext(xlab,1,maj,cex=mas)
    lines(ii,mu)
  }
 
  ## serology samples...
  
  N <- sum(S0[4:13])
  y <- dat$n_positive;ii <- !is.na(y); y<- y[ii];
  n <- dat$total_samples[ii]
  ii <- dat$day[ii]
  mu <- yout[6,ii]/N
  ls <- sum(dbinom(y,n,mu,log=TRUE))
  dls.dmu <- y/mu - (n-y)/(1-mu)
  dls <- drop(yout[1:np*nyo+6,ii]%*%dls.dmu)/N
  if (plot&&length(y)) {
    plot(ii,y/n,ylab="",ylim=range(c(y/n,mu)),col="grey",xlim=xlim,xlab="",pch=pch,cex=pex)
    mtext("P(sero)",2,maj,cex=mas);mtext(xlab,1,maj,cex=mas)
    lines(ii,mu,col=1)
  }
  
  ## REACT PCR samples...

  N <- sum(S0[2:18]) 
  y <- dat$react_pos;ii <- !is.na(y); y<- y[ii];
  n <- dat$react_samples[ii]
  ii <- dat$day[ii]
  mu <- yout[7,ii]/N
  lr <- sum(dbinom(y,n,mu,log=TRUE))
  dlr.dmu <- y/mu - (n-y)/(1-mu)
  dlr <- drop(yout[1:np*nyo+7,ii]%*%dlr.dmu)/N
  if (plot&&length(y)) {
    plot(ii,y/n,ylab="",col="grey",xlim=xlim,xlab="",pch=pch,cex=pex);lines(ii,mu,col=1)
    mtext("P(pcr)",2,maj,cex=mas);mtext(xlab,1,maj,cex=mas)
  }
  
  l=lh+lo+lu+ld+lg+lr+ls ## total log likelihood
  dl=dlh+dlo+dlu+dld+dlg+dlr+dls ## corresponding derivative vector

  if (!is.null(Sb)) { ## add smoothing penalty
    Sb <- drop(Sb %*% theta[1:nb])
    l <- l - sum(theta[1:nb]*Sb)/2
    dl[1:nb] <- dl[1:nb] - Sb
  }
  R <- yout[9,]; R[R<0] <- NA; Rd <- NULL
  if (plot||getR) { ## get the Diekmann version of R...
    gamma_A <- gamma[2];gamma_C <- gamma[3];p_c <- gamma[4]
    Rd <- R
    beta <- yout[10,]
    C <- matrix(daux[1:(nc-1)^2],nc-1,nc-1)[-(nc-1),-(nc-1)]
    n <- max(dat$day)
    xout <- matrix(vm$xout,nx,n)[1:(nc*ns),]
    S <- xout[0:(nc-3)*ns+1,]
    for (i in 1:n) {
      C0 <- S[,i]*C
      K <- beta[i]*rbind(cbind((1-p_c)*C0/gamma_A,(1-p_c)*C0/gamma_C),cbind(p_c*C0/gamma_A,p_c*C0/gamma_C))
      Rd[i] <- abs(eigen(K,only.values = TRUE)$values[1])
    }
  } else xout <- NULL
  if (ncol(link)>3) { ## prior on dynamic params too
     ii <- (nb+1):np
     sig.p <- link[ii,5]
     l <- l - sum((thetaw[ii]-link[ii,4])^2/(2*sig.p^2)) 
     dprior <- c(rep(0,nb),-(thetaw[ii]-link[ii,4])/sig.p^2)
  } else dprior <- 0
  list(l=l,dl=dl*th$dtheta+dprior,R=R,Rd=Rd,fi=yout[8,],xout=xout,n.death=n.death,munb=munb,ynb=ynb)
} ## ll


thfix <- function(theta,fixi,fixv) {
## if fixi!=NULL then it indexes parameters that have been dropped from the optimization of
## theta, while fixv gives their fixed values. This routine inserts the fixed values into
## theta in the correct place...
  if (!is.null(fixi)) {
    if (length(fixi)!=length(fixv)) stop("fixv wrong length")
    th1 <- rep(0,length(theta)+length(fixi))
    ii <- rep(FALSE,length(th1))
    ii[fixi] <- TRUE
    th1[ii] <- fixv
    th1[!ii] <- theta
    theta <- th1
  }
  theta
} ## thfix

nll <- function(theta,gamma,S0,iaux,daux,dat,Sb,link,step=1,nbk=rep(2,5),fixi=NULL,fixv=NULL) {
## if fixi!=NULL then it indexes parameters that have been dropped from the optimization of
## theta, while fixv gives their fixed values.
 iaux[7] <- 0; ## no derivatives needed
 theta <- thfix(theta,fixi,fixv)
 -ll(theta,gamma,S0,iaux,daux,dat,Sb=Sb,link=link,plot=FALSE,step=step,nbk=nbk)$l
}

nll0 <- function(theta,gamma,S0,iaux,daux,dat,Sb,link,step=1,nbk=rep(2,5),fixi=NULL,fixv=NULL) {
## if fixi!=NULL then it indexes parameters that have been dropped from the optimization of
## theta, while fixv gives their fixed values.
 iaux[7] <- 0; ## no derivatives needed
 theta <- thfix(theta,fixi,fixv)
 ll(theta,gamma,S0,iaux,daux,dat,Sb=Sb,link=link,plot=FALSE,step=step,nbk=nbk)
}

nllg <- function(theta,gamma,S0,iaux,daux,dat,Sb,link,step=1,nbk=rep(2,5),fixi=NULL,fixv=NULL) {
## grad function corresponding to nll
 theta <- thfix(theta,fixi,fixv)
 g <- -ll(theta,gamma,S0,iaux,daux,dat,Sb=Sb,link=link,plot=FALSE,step=step,nbk=nbk)$dl
 if (!is.null(fixi)) g <- g[-fixi]
 g
}

hessian <- function(theta,gamma,S0,iaux,daux,dat,Sb,link,step=1,nbk=rep(2,5),fixi=NULL,fixv=NULL) {
  ## compute Hessian by finite differencing 
  nt <- length(theta)
  eps <- 1e-7
  g0 <- nllg(theta,gamma,S0,iaux,daux,dat,Sb,link,step,nbk,fixi,fixv)
  H <- matrix(0,nt,nt)
  for (i in 1:nt) {
    theta1 <- theta;theta1[i] <- theta[i] + eps
    g1 <- nllg(theta1,gamma,S0,iaux,daux,dat,Sb,link,step,nbk,fixi,fixv)
    H[,i] <- (g1-g0)/eps;
  }
  .5*(t(H)+H)
} ## hessian

nll1 <-  function(theta,gamma,S0,iaux,daux,dat,Sb,link,step=1,nbk=rep(2,5),fixi=NULL,fixv=NULL) {
## unused nlm version
 theta <- thfix(theta,fixi,fixv)
 l0 <- ll(theta,gamma,S0,iaux,daux,dat,Sb=Sb,link=link,plot=FALSE,step=step,nbk=nbk)
 nll <- -l0$l
 attr(nll,"gradient") <- -l0$dl
 nll
}

ginv <- function(V,warn=FALSE) { ## generalized inverse
  if (!all(is.finite(V))) return(NULL)
  ev <- eigen(V,symmetric=TRUE)
  d <- ev$values;ii <- d>0
  d[ii] <-1/d[ii];d[!ii] <- 0
  if (warn&&sum(!ii)) warning("indefinite V\n")
  ev$vectors%*%(d*t(ev$vectors))
}

tf <- function(x,a1=.00,a2=.1) { ## transform as applied to b(t)
  ii <- x<0
  x[ii] <- exp(x[ii])
  x[!ii] <- exp(-x[!ii])
  x[ii] <- a1 + (a2-a1)*x[ii]/(1+x[ii])
  x[!ii] <- a1 + (a2-a1)*1/(1+x[!ii])
  x
}

itf <- function(y,a1=.00,a2=.1) { ## transform as applied to b(t)
  y <- (y-a1)/(a2-a1)
  log(y/(1-y)) 
}


setup.model <- function(location="london",Xtype=2,beta.off = .25,rep41C=FALSE) {
## sets up the information required for a model fit, including fixed parameters,
## data, initial parameter values and limits.
 
  source("params.r") ## some fixed parameters

  region <- c("east_of_england","london","midlands","north_east_and_yorkshire",
  "north_west","south_east","south_west")

  if (!location %in% region) stop("location wrong")

  source("data.r",local=TRUE) ## obtain the fit data and some other parameters

  load(paste("C_",location,".rda",sep=""))
  ## or use m from sircovid package - problem is that it only matches what is documented
  ## up to age 70 (parameters change but results unchanged)... 
  if (rep41C) C <- ch$m[1:18,1:18] 
  
  n <- max(rtr$day) ## number of output steps
  step <- 1 ## steps between output
  h <- 1/step ## step length
  hx <- h/2 ## step at which basis evaluation required
  t <- h*0:(2*step*n)/2 ## when basis evaluation required 
  nb <- 80 ## number of b(t) coeffs
  if (Xtype==1) {## adaptive
    t0 <- 45 ## constant b(t) before here
    ii <- t >= t0 ## set up adaptive smooth from t0
    smb <- smoothCon(s(t,k=nb,bs="ad"),data=data.frame(t=t[ii]))[[1]]
    ni <- sum(!ii) ## need to pad model matrix to fill in constant prior to t0
    X0 <- matrix(smb$X[1,],ni,ncol(smb$X),byrow=TRUE)
    Xb <- rbind(X0,smb$X) ## b(t) basis matrix
    Sb <- smb$S ## penalty matrix
    attr(Sb,"rank") <- smb$rank
  } else if (Xtype==2) { ## force b(t) to be constant at start
    ne <- 5
    smb <- smoothCon(s(t,k=nb+ne,bs="bs"),data=data.frame(t=t))[[1]]
    T <- diag(nb+ne)[,-(1:ne)];T[1:ne,1] <- 1
    Xb <- smb$X%*%T ## b(t) basis matrix  
    Sb <- t(T)%*%smb$S[[1]]%*%T ## penalty matrix
  } else { ## rep 41
    nb <- length(bk)
    Xb <- matrix(0,length(t),nb)
    Xb[,1] <- as.numeric(t<=bk[1])
    for (i in 1:nb) {
      yk <- bk*0;yk[i] <- 1
      Xb[,i] <- approx(bk,yk,t,rule=2)$y
    }
    Sb <- NULL
  }


  b0 <- c(rep(-2.4,5),-2.404542,-2.468239,-2.676757,-3.068759,-3.510927,
  -3.833448,-3.971020,-3.999736,rep(-4,10),-3.999977,-3.987959,-3.920592,
  -3.752347,-3.506797,-3.277703,-3.147313,-3.104178,rep(-3.5,14))
  t0 <- 0:44*7.5+7
  b0 <- exp(approx(t0,b0,t,rule=2)$y)
  a1 <- 0; a2 <- 0.1
  b0[b0>a2] <- a1 + (a2-a1)*.95
  b0[b0<a1] <- a1 + (a2-a1)*.05
  b0 <- itf(b0,a1,a2)
  beta <- coef(lm(b0~Xb-1)) ## initial b(t) parameters
  beta <- beta - beta.off

  ## set up param vector
  eps <- 0.1; m_chw <- 4e-6; m_chr <- 4e-6
  t0 <- 20;p_g <- 0.6;p_h_max <- 0.45
  mu_ic <- .9; mu_d <- .9; p_ic_max <- .3
  p_icd_max <- .6;p_wd_max <- .3;p_hd_max <- .37
  gamma_G <- .4;
  gamma_IC_Wr <- .0641;gamma_Wr <- .16;gamma_Hr<- 0.0935

  theta <- c(beta,eps,m_chw,m_chr,t0,p_h_max,p_g,mu_ic,p_ic_max,p_icd_max,mu_d,p_wd_max,p_hd_max,gamma_G,
             gamma_IC_Wr,gamma_Wr,gamma_Hr)
  names(theta)[(nb+1):length(theta)] <- c("eps","m_chw","m_chr","t0", "p_h_max","p_g","mu_ic","p_ic_max",
                             "p_icd_max","mu_d","p_wd_max","p_hd_max","gamma_G","gamma_IC_Wr","gamma_Wr","gamma_Hr")

  link <- c(rep(0,nb),1,1,1,0,2,2,2,2,2,2,2,2,1,1,1,1) ## link types: 0 - I(), 1 - log, 2 - logit
  a1 <-   c(rep(0,nb),0,0,0,0,0,0,.66,.14,.57,.66,.25,.36,0,0,0,0) ## lower logistic limit
  a2 <-   c(rep(0,nb),0,0,0,0,1,1, 1 ,.36,.77,  1,.46,.56,0,0,0,0)  ## upper logistic limit
  ### prior mean on working scale...
  mup <- c(rep(0,nb),-.7,-12.7,-12.7,30,1.1,1.1,-1,-.2,0,-1,0,0,-1,-2.75,-1.8,-2.37)
  sdp <- c(rep(NA,nb),2,2,2,10,2,2,2,2,2,2,1,1,2,2,2,2) ## prior s.d. on working scale
  link <- cbind(link,a1,a2,mup,sdp)
  ## fixed param vector
  p_c <- 0.6 # .6 ## Knock et al
  inoc_sig <- 7.5 # arbitrary
  p_star <- 0.25 ## Knock et al
  f_start <- 92 ## improvement start 
  f_end <- 152  ## improvement end
  gamma_U <- .5 ## Table S7
  ## gamma_E, gamma_C and gamma_ph should keep mean time from infection to hospital at
  ## 12.3 days to match cited literature 
  #gamma_E <- .44 ## Table S2 or use ch$gamma_E based on early lauer paper
  #gamma_E <- 1/2.9 ## based on McAloon meta-analysis!
  gamma_E <- 1/1.5 ## mean 3 days to infective
  #gamma_C <- .25 ## Table S2 - based on equating the median and mean of an exponential!
  #gamma_C <- 1/6.5 ## Based on cited Docherty paper *mean* for males and females (see supp)
  gamma_C <- 1/4 ## mean 4 days infective => generation time 6.2 (brought down by asymptomatics) 
  gamma_ph <- 1/5.3 ## 5.3 is the length of the post I_C pre hospital phase (2.8 without correction)
  gamma_A <- ch$gamma_asympt
  gamma_IC_Wr <- .0641 ## Table S2
  gamma_IC_Wd <- .14 ## Table S2
  gamma_IC_d <- .17 ## Table S2
  gamma_Wr <- .16  ## Table S2
  gamma_Wd <- .12  ## Table S2
  gamma_Hr <- .0935 ## Table S2
  gamma_Hd <- .19 ## Table S2
  gamma_S <- 1/13 ## Table S7
  p_S <- .85 ## Table S7

  gamma0 <- c(gamma_E,gamma_A,gamma_C,p_c,inoc_sig,gamma_ph,p_star,
           f_start,f_end,ch$gamma_ICU_D,gamma_U,gamma_IC_Wr,gamma_IC_Wd,gamma_IC_d,
	   gamma_Wr,gamma_Wd,gamma_Hr,gamma_Hd,ch$gamma_PCR_pre,
           ch$gamma_PCR_pos,gamma_S,p_S)

  names(gamma0) <- c("gamma_E","gamma_A","gamma_C","p_c","inoc_sig","gamma_ph","p_star",
  "f_start","f_end","gamma_IC_pre","gamma_U","gamma_IC_Wr","gamma_IC_Wd","gamma_IC_d",
  "gamma_Wr","gamma_Wd","gamma_Hr","gamma_Hd","gamma_P","gamma_Po","gamma_S","p_S")

  gamma <- c(gamma0,psi.h,psi.ic,psi.icd,psi.wd,psi.hd) ## 19 vector
  ng <- length(gamma)
  np <- length(theta)
  nc <- 19 ## classes
  ns <- 37 ## states per class
  nk <- prod(dim(Xb))
  deriv <- 1
  iaux <- c(nb,np,nc,ns,ng,nk,deriv,-1)
  daux <- c(C,gamma,t(Xb),hx,rep(0,nb+2*(nc+nc*np)))
  S0=ch$N_tot
  thetaw <- lfun(theta,link=link);fixi <- fixv <- NULL
  arg <- list(location=location,gamma=gamma,S0=S0,thetaw=thetaw,iaux=iaux,daux=daux,dat=rtr,Sb=Sb,link=link,step=step)
  arg
} ## setup.model

nlnb <- function(theta,y,mu) {
## negative log lik for NB
  ii <- is.na(y)|is.na(mu)
  y <- y[!ii];mu <- mu[!ii]
  -sum(theta*log(theta/(theta+mu)) + y*log(mu/(mu+theta)) + lgamma(theta+y) - lgamma(theta))
}


fit.model <- function(arg) {
## fit model using Fellner Schall sp updates...
  if (is.null(arg$Sb)) { m <- 0; Sb <- NULL} else
  if (is.list(arg$Sb)) {
    lambda <- rep(10,length(arg$Sb))
    Sb <- arg$Sb[[1]]*lambda[1]
    m <- length(arg$Sb)
    if (m>1) for (i in 2:m) Sb <- Sb +  arg$Sb[[i]]*lambda[i]
    if (m==1) arg$Sb <- arg$Sb[[1]]
  } else {
    m <- 1
    Sb <- arg$Sb * 10 
  }
  thetaw <- if (is.null(arg$fixi)) arg$thetaw else arg$thetaw[-arg$fixi]
  er <- optim(thetaw,nll,gr=nllg,method="BFGS",control=list(maxit=400),gamma=arg$gamma,S0=arg$S0,
            iaux=arg$iaux,daux=arg$daux,dat=arg$dat,Sb=Sb,link=arg$link,step=arg$step,
	    nbk=rep(2,5),fixi=arg$fixi,fixv=arg$fixv)
  cat(arg$location," ")
  
  er$par -> thetaw -> thetaw0

  if (is.null(arg$Sb)) { m <- 0; Sb <- NULL} else
  if (is.list(arg$Sb)) {
    lambda <- rep(10,length(arg$Sb))
    Sb <- arg$Sb[[1]]*lambda[1]
    m <- length(arg$Sb)
    if (m>1) for (i in 2:m) Sb <- Sb +  arg$Sb[[i]]*lambda[i]
  } else {
    m <- 1; lambda <- 10
    Sb <- arg$Sb * lambda
  }
  nrep <- if (m>0) 20 else 1 ## based on prelim run

  nbk <- c(2,2,2,2000,50) ## admissions, ward and icu occupancy, hospital deaths, care home deaths
  k1 <- c(40,100,100,2000,2000) ## upper limits on overdispersion param
  for (i in 1:nrep) {
    er0 <- er
    er <- optim(thetaw,nll,gr=nllg,method="BFGS",control=list(maxit=400),gamma=arg$gamma,S0=arg$S0,
            iaux=arg$iaux,daux=arg$daux,dat=arg$dat,Sb=Sb,link=arg$link,step=arg$step,
	    nbk=nbk,fixi=arg$fixi,fixv=arg$fixv,hessian=FALSE)
    ## get Hessian without Sb first - check +ve semidef	    
    H <- hessian(er$par,arg$gamma,arg$S0,arg$iaux,arg$daux,arg$dat,Sb=NULL,link=arg$link,
                 step=arg$step,nbk=nbk,arg$fixi,arg$fixv)
    eh <- eigen(H)
    if (any(eh$values<0)) { ## find nearest +ve semi definite
      eh$values[eh$values<0] <- 0
      H <- eh$vectors %*% (eh$values*t(eh$vectors))
    }
    Hp <- H ## penalized version
    if (!is.null(Sb)) {
      nb <- ncol(Sb)
      Hp[1:nb,1:nb] <- Hp[1:nb,1:nb] + Sb
    }
    
    ## update nb theta...
    l0 <- nll0(thetaw,gamma=arg$gamma,S0=arg$S0,iaux=arg$iaux,daux=arg$daux,
            dat=arg$dat,Sb=Sb,link=arg$link,step=arg$step)
    for (i in c(1,5)) { ## optimize log lik
      nbk[i] <- optimize(nlnb,interval=c(2,k1[i]),y=l0$ynb[,i],mu=l0$munb[,i])$minimum 
    }
    for (i in 2:3) { ## update occupancy rates overdispersion
      y <- l0$ynb[,i+4];ii <- !is.na(y);y <- y[ii]
      mu1 <- l0$munb[ii,i*2+2]
      mu2 <- l0$munb[ii,i*2+3]
      k.mle <- mean((y-mu1+mu2)^2/(mu1+mu2)) 
      k.mle <- if (k.mle < 1) 1 else { if (k.mle>k1[i]) k1[i] else k.mle}
      nbk[i] <- k.mle
    }
    er$par -> thetaw
    thetaw <- thfix(er$par,arg$fixi,arg$fixv)
    Vb <- ginv(Hp,warn=TRUE)
    if (is.null(Vb)) break
    nb <- arg$iaux[1]
    b <- thetaw[1:nb]
    if (m>1) { 
      Si <- ginv(Sb)
      if (is.null(Si)) { Vb <- NULL; break}
      for (i in 1:m) {
        fac <- drop((sum(Si*arg$Sb[[i]]) - sum(Vb[1:nb,1:nb]*arg$Sb[[i]]))/t(b)%*%arg$Sb[[i]]%*%b)
	if (fac<0.2) fac <- .2
	if (fac>5) fac <- 5 
        lambda[i] <- lambda[i]*fac
      }
     
      Sb <- arg$Sb[[1]]*lambda[1]
      for (i in 2:m) Sb <- Sb +  arg$Sb[[i]]*lambda[i]
    } else if (m==1) {
      fac <- drop(((nb-2)/lambda - sum(Vb[1:nb,1:nb]*arg$Sb))/t(b)%*%arg$Sb%*%b)
      if (fac<0.2) fac <- .2
      if (fac>5) fac <- 5
      lambda <- lambda*fac
      Sb <- lambda*arg$Sb
    }
    if (m>0) cat(arg$region,": ",lambda,"\n")
  }
  if (is.null(Vb)) {
    er <- er0
    er$par -> thetaw
    Vb <- ginv(er$hessian)
  }
  list(thetaw=thetaw,theta=ilink(thetaw,arg$link)$theta,Vb=Vb,H=Hp,lambda=lambda,nbk=nbk)
} ## fit.model


check.model <- function(arg,plot=TRUE,deriv=FALSE) {
## checks if initial params gtive something remotely reasonable
  arg$iaux[7] <- as.numeric(deriv) 
  l0 <- ll(arg$thetaw,arg$gamma,arg$S0,arg$iaux,arg$daux,arg$dat,Sb=NULL,link=arg$link,plot=plot,step=arg$step)
}

check.derivs <- function(thetaw,arg,eps = 1e-6) {
## checks derivatives of penalized log likelihood. Note currently not checking with smoothing penalty.
  l0 <- ll(thetaw,arg$gamma,arg$S0,arg$iaux,arg$daux,arg$dat,Sb=NULL,link=arg$link,plot=FALSE,step=arg$step)
  fd <- l0$dl*0
  arg$iaux[7] <- 0
  for (i in 1:length(fd)) {
    theta1 <- thetaw;theta1[i] <- theta1[i] + eps
    l1 <- ll(theta1,arg$gamma,arg$S0,arg$iaux,arg$daux,arg$dat,Sb=NULL,link=arg$link,plot=FALSE,step=arg$step)
    fd[i] <- (l1$l-l0$l)/eps
    cat(".")
  }
  cat("\n")
  plot(l0$dl,fd);abline(0,1)
  list(dl=l0$dl,fd=fd,l=l0$l)
}

dRfi <- function(arg) {
## lazy evaluation of derivatives of R and incidence trajectories w.r.t. parameters, assuming thetaw and Vb added to arg
  eps = 1e-5
  thetaw <- arg$thetaw
  l0 <- ll(thetaw,arg$gamma,arg$S0,arg$iaux,arg$daux,arg$dat,Sb=NULL,link=arg$link,plot=FALSE,step=arg$step,getR=TRUE)
  dR <- df <- matrix(0,length(l0$Rd),length(thetaw))
  for (i in 1:length(thetaw)) {
    theta1 <- thetaw;theta1[i] <- theta1[i] + eps
    l1 <- ll(theta1,arg$gamma,arg$S0,arg$iaux,arg$daux,arg$dat,Sb=NULL,link=arg$link,plot=FALSE,step=arg$step,getR=TRUE)
    dR[,i] <- (log(l1$Rd)-log(l0$Rd))/eps
    df[,i] <- (log(l1$fi)-log(l0$fi))/eps
  }
  list(dR=dR,df=df,se.r=rowSums((dR%*%arg$Vb) * dR)^.5,se.f=rowSums((df%*%arg$Vb) * df)^.5,fi=l0$fi,Rd=l0$Rd)
} ## dRfi

plot.fit <- function(fit,arg) {
## quick plot of fit to data + R and incidence.
  l0 <- ll(fit$thetaw,arg$gamma,arg$S0,arg$iaux,arg$daux,arg$dat,Sb=NULL,link=arg$link,plot=TRUE,step=arg$step)
  X11()
  par(mfrow=c(2,2))
  plot(log(l0$Rd),xlim=c(50,150));abline(0,0);abline(v=84)
  plot(log(l0$Rd));abline(1,0);abline(v=84)
  plot(l0$Rd,xlim=c(50,150));abline(1,0);abline(v=84)

  plot(l0$fi,xlim=c(50,150));abline(v=84)
  
} ## plot.fit

###############################################################################
setwd("~sw283/research.notes/covid19/rep41/code/albus") ## Data & Code location
###############################################################################

library(mgcv);library(parallel)
dyn.load("albus-gl.so")

region <- c("east_of_england","london","midlands","north_east_and_yorkshire",
  "north_west","south_east","south_west")

rlab <- c("East","London","Midlands","N.E. & Yorkshire", "North West", "South East", "South West")

## setup regional fits...
arg <- list()
for (i in 1:7) arg[[i]] <-setup.model(region[i],Xtype=1,beta.off=.75,rep41C=FALSE)
names(arg) <- rlab
load("thetaw0-g.rda")
for (i in 1:7) {
  arg[[i]]$region <- rlab[i]
  arg[[i]]$thetaw <- thetaw0 ## initialize from old fit to midlands to speed up
}


## par(mfrow=c(2,3));check.model(arg[[1]])
## plot.fit(er,arg[[1]])
## system.time(fd <- check.derivs(arg[[1]]$thetaw,arg[[1]]));fd$dl[-(1:80)];fd$fd[-(1:80)]

system.time(er <- mclapply(arg,fit.model,mc.cores=7))


for (i in 1:7) { ## add fit results to arg list
  arg[[i]]$thetaw <- er[[i]]$thetaw
  arg[[i]]$Vb <- er[[i]]$Vb
}
## get derivatives of R and fi (incidence) w.r.t. parameters
system.time(dd <- mclapply(arg,dRfi,mc.cores=7)) 

# for (i in 1:7) arg[[i]]$iaux[8] <- 83 ## freeze b(t) at b(83)

## check all fits and compute total incidence and average R -
## R is number of new infections per existing infection, so
## must be weighted by existing infectives.
ns <- arg[[1]]$iaux[4]
nc <- arg[[1]]$iaux[3]
ii <- rep((0:(nc-3))*ns,each=2)+4:5 ## index of all infectives outside care homes
n.deaths <- 0

ps <- FALSE
if (ps) postscript("fits.eps",horizontal=F,height=9,width=11)
par(mfrow=c(7,5),mar=c(3.5,4,0,0))# par(mfcol=c(5,7),mar=c(4,4,1,1))
xlab=""
for (k in 1:7) {
  arg[[k]]$iaux[7] <- 0 ## turn off derivatives
  if (k==7) xlab="day"
  l0 <- ll(er[[k]]$thetaw,arg[[k]]$gamma,arg[[k]]$S0,arg[[k]]$iaux,arg[[k]]$daux,arg[[k]]$dat,Sb=NULL,
        link=arg[[k]]$link,plot=TRUE,step=arg[[k]]$step,getR=TRUE,lab=arg[[k]]$region,xlab=xlab)
  n.deaths <- n.deaths + l0$n.death	
  arg[[k]]$iaux[7] <- 1
  Ik <- colSums(l0$xout[ii,])
  rep <- 1000;n <- length(l0$Rd)
  if (k==1) {
    R <- exp(log(dd[[1]]$Rd) + matrix(rnorm(n*rep),n,rep)*dd[[1]]$se.r)*Ik 
    f <- exp(log(dd[[1]]$fi) + matrix(rnorm(n*rep),n,rep)*dd[[1]]$se.f)
    Rd <- dd[[k]]$Rd*Ik 
    ftot <- dd[[k]]$fi
    Itot <- Ik
  } else {
    R <- R + exp(log(dd[[k]]$Rd) + matrix(rnorm(n*rep),n,rep)*dd[[k]]$se.r)*Ik
    f <- f + exp(log(dd[[k]]$fi) + matrix(rnorm(n*rep),n,rep)*dd[[k]]$se.f)
    Itot <- Itot + Ik
    Rd <- Rd+ dd[[k]]$Rd*Ik
    ftot <-ftot + dd[[k]]$fi
  }
  #text(200,.01,names(arg)[k])
  #X11()
}
if (ps) dev.off()

Rd <- Rd/Itot
Rav <- R/Itot 
Rci <- apply(Rav,1,quantile,prob=c(.025,.975))
fci <- apply(f,1,quantile,prob=c(.025,.975))

for (j in 1:2) {
if (j==1) {
  xlim <- c(50,150);xpos <- 110
  fn <- "fw-gc.eps"
} else {  
  xlim <- c(60,320);xpos <-150
  fn <- "w2-gc.eps"
}  
vt <- c(84,216,245,258,288,310)

if (ps) postscript(fn) else X11()
par(mfrow=c(2,2),mar=c(5,5,1,1))
t <- 1:length(Rd)
c1 <- 1.3
plot(ftot,xlim=xlim,type="l",ylim=c(0,max(fci[2,50:150])),xlab="day",ylab="incidence",cex.lab=c1)
lines(t,fci[1,],lty=2);lines(t,fci[2,],lty=2)
abline(v=vt)

plot(dd[[1]]$fi,xlim=xlim,ylim=c(0,70000),type="l",lwd=3,xlab="day",ylab="incidence",cex.lab=c1)
lines(t,exp(log(dd[[1]]$fi)+dd[[1]]$se.f*2),lty=2)
lines(t,exp(log(dd[[1]]$fi)-dd[[1]]$se.f*2),lty=2)
text(xpos,20000,arg[[1]]$region)
for (k in 2:7) {
  lines(dd[[k]]$fi,col=k,lwd=3)
  text(xpos,10000+k*10000,arg[[k]]$region,col=k)
  lines(t,exp(log(dd[[k]]$fi)+dd[[k]]$se.f*2),lty=2,col=k)
  lines(t,exp(log(dd[[k]]$fi)-dd[[k]]$se.f*2),lty=2,col=k)
}
abline(v=vt)
plot(Rd,xlim=xlim,ylim=c(0,4),type="l",lwd=1,xlab="day",ylab="R",cex.lab=c1)
lines(t,Rci[1,],lty=2);lines(t,Rci[2,],lty=2)

abline(1,0);abline(v=vt)


plot(log(dd[[1]]$Rd),xlim=xlim,ylim=c(-2,2),type="l",lwd=3,xlab="day",ylab="log(R)",cex.lab=c1)
lines(t,log(dd[[1]]$Rd)+dd[[1]]$se.r*2,lty=2)
lines(t,log(dd[[1]]$Rd)-dd[[1]]$se.r*2,lty=2)
text(xpos,.35,arg[[1]]$region)
for (k in 2:7) {
  lines(log(dd[[k]]$Rd),col=k,lwd=3)
  lines(t,log(dd[[k]]$Rd)+dd[[k]]$se.r*2,lty=2,col=k)
  lines(t,log(dd[[k]]$Rd)-dd[[k]]$se.r*2,lty=2,col=k)
  text(xpos,.1+k*.25,arg[[k]]$region,col=k)
}
abline(0,0);abline(v=vt)
if (ps) dev.off()
}
