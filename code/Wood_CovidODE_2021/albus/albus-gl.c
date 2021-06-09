/* (c) Simon N Wood (2021) Alternative Light But Useful Simulation of rep 41 model... 
   cd ~sw283/research.notes/covid19/rep41/code/albus
   R CMD SHLIB albus-gl.c
   version in which four of the gamma parameters are fit parameters
   also there is an extra stage between I_c and hospital to allow 
   the generation time to be manipulated without changing disease 
   progression average timings.
   This version provides flows needed to correct hospital likelihood 
   to something defensible.
   knarf dna lezah rof
*/
#include <math.h>
#include "albus.h"

void betaX(double *beta,double *b,double *theta,double t,double *baux,int nk,int nb,int t_freeze) {
/* the basis for the infection modifier - simple case: assume baux is a model matrix, with time step
   h appended at end. Integration starts from zero.    
   infection modifier is t(sum_i(beta[i]*theta[i])) for an exponential of logistic transform, t
*/
  int i,j,n,r;
  double h,a1=.0,a2=.1,d,dd,x;
  n = nk/nb;
  h = baux[nk]; // the model matrix timestep between rows
  if (t_freeze>0&&t>t_freeze) t = (double) t_freeze; // return b(t_freeze) for any t>t_freeze
  t = t/h;r = floor(t);if (t-r > .5) r++;
  *beta =0.0;
  for (i=0,j=r*nb;i<nb;i++,j++) {
    //b[i] = baux[r+i*n];
    b[i] = baux[j]; // assuming X^T supplied for faste access
    *beta += b[i]*theta[i]; 
  }
  if (0) { // exponential transform
    *beta = exp(*beta);
    for (i=0;i<nb;i++) {
      b[i] *= *beta; // derivative of beta w.r.t. ith coef 
    }  
  } else { // a1 + (a2-a1)*exp(beta)/(1+exp(beta)) transform
    if (*beta>0) {
      x = exp(- *beta);d = 1/(x+1);dd = x*d*d;
    } else {
      x = exp(*beta);d = x/(1+x);dd = d/(1+x);
    } // stable evaluation of exp(beta)/(1+exp(beta))
    *beta = a1 + (a2-a1)*d; // beta' - transformed beta
    d = dd*(a2-a1); // dbeta'/dbeta
    for (i=0;i<nb;i++) b[i] *= d; // dbeta'/dtheta
  }  
} // betaX  

void lambda(double *beta,double *lam,double *dlam,double *x,double *theta,double *C,
	    double t,int nb,int np,int nc, int ns,int nk,int deriv,
	    double *baux,double *work,int finf,int t_freeze) {
/* Evaluates the force of infection at t for each of the nc classes, 
   returning the result in lam. If deriv!=0, then return derivatives 
   of lamda w.r.t. first deriv (>=nb+4) params in dlam. C is the symmetric supplied 
   scaled mixing matrix and is nc-1 by nc-1 - its final element is ignored, 
   and its final row/col is c_CHW.    
   x[i*ns+j] is jth state variable in ith class
   x[(k+1)*ns*nc+i*ns+j] is deriv of the above w.r.t. param k.

   lam is an nc vector
   dlam is nc*np dlam[i+k*nc] is d lam[i]/d theta[k]
   work is nb+nc+nc*np    

   This is basically the interesting bit of the dynamics, the rest is tedious 
   follow on compartments. 

   NOTE: in fact lambda only depends on parameters beta, eps, m_chw, m_chr and t0,
         all other sensitivities will be zero and can be ignored. So typically 
         deriv set to nb+4, and will be reset to this if non-zero but <nb+4. 

   if finf!=0 then routine computes the total infections caused by asymptomatics in lam
   and the total infections caused by symptomatics in dlam. No derivatives are then computed.

*/
  double *bX,*I,z,z1,*Ci,eps,m_chw,m_chr,*dI,*dum;
  int i,j,j1,k,nsc,drop_ch=0;
  if (finf) {
    deriv=0;
    drop_ch = 1; // drop care homes from force of infections, so only returning community transmission
  }  
  if (deriv&&deriv<nb+4) deriv = nb+4; //
  nsc = ns * nc;
  eps = theta[nb];
  m_chw = theta[nb+1];
  m_chr = theta[nb+2];
  bX = work; work += nb;
  I = work;work += nc;
  betaX(beta,bX,theta,t,baux,nk,nb,t_freeze); // evaluate basis for infection modifier at t
  // add up the infectives in each class...
  if (finf) for (i=0;i<nc;i++) I[i] = x[i*ns + 3];
  else for (i=0;i<nc;i++) I[i] = x[i*ns + 3] + x[i*ns + 4];
  if (deriv) { // collect the sensitivities of the infectives
    dum = dI = work; work += nc * np;
    for (j=0;j<deriv;j++) { // param loop
      for (i=0;i<nc;i++,dum++) { // class loop
	k = (j+1)*nsc+i*ns;
	*dum = x[k+3] + x[k+4]; // d I_i/d theta_j in dI[j*nc+i]
      }
    }  
  }  
 
  for (i=0;i<nc-2;i++) { // age class loop
    for (z=0,j=0,Ci=C+(nc-1)*i;j<nc-1-drop_ch;j++,Ci++) z += *Ci * I[j];
    if (!drop_ch) z += C[(nc-1)*i + nc-3] * I[nc-1] * eps;
    lam[i] = z * *beta;
    if (finf) {
      lam[i] *= x[i*ns]; // multiply lam[i] by S_i - rate of asymptomatic caused infections
      for (z=0,j=0,Ci=C+(nc-1)*i;j<nc-1-drop_ch;j++,Ci++) z += *Ci * x[j*ns+4]; // 
      if (!drop_ch) z += C[(nc-1)*i + nc-3] * x[(nc-1)*ns+4] * eps;
      dlam[i] = z * *beta * x[i*ns]; //daily rate of infections caused by symptomatics 
    }  
    if (deriv) {
      // first the sensitivities from the I state dependencies...
      for (k=0;k<deriv;k++) { // parameter loop
        for (z1=0,j=k*nc,j1=j+nc-1,Ci=C+(nc-1)*i;j<j1;j++,Ci++) z1 += *Ci * dI[j];
	z1 += C[(nc-1)*i + nc-3] * dI[k*nc+nc-1] * eps;
	dlam[i+k*nc] = z1 * *beta; 
      }
      // now add the direct beta dependencies
      for (k=0;k<nb;k++) {
        dlam[i+k*nc] += bX[k]*z;
      } 	
      // and the epsilon (theta[nb]) dependency
      dlam[i+nb*nc] += *beta * C[(nc-1)*i + nc-3] * I[nc-1];//dI[k*nc+nc-1];
    } // deriv 
  }  
  // now care home workers...
  i = nc-2; // CHW row/col
  for (z=0,j=0,Ci=C+(nc-1)*i;j<nc-2;j++,Ci++) z += *Ci * I[j];
  lam[nc-2] = *beta * z + m_chw * (I[nc-2]+I[nc-1]);
  if (finf) {
     lam[nc-2] *= x[(nc-2)*ns]; // multiply lam[i] by S_i - rate of symptomatic caused infections
     for (z=0,j=0,Ci=C+(nc-1)*i;j<nc-2;j++,Ci++) z += *Ci * x[j*ns+4]; // 
     dlam[nc-2] = (z * *beta + m_chw * (x[(nc-2)*ns+4]+x[(nc-1)*ns+4]))*x[(nc-2)*ns]; //daily rate of infections caused by symptomatics 
  }
  if (deriv) {
    // first the sensitivities from the I state dependencies...
    for (k=0;k<deriv;k++) { // parameter loop
      for (z1=0,j=0,Ci=C+(nc-1)*i;j<nc-2;j++,Ci++) z1 += *Ci * dI[k*nc+j];
      dlam[i+k*nc] = z1 * *beta + m_chw*(dI[k*nc+nc-2]+dI[k*nc+nc-1]); 
    }
    // now add the direct beta dependencies
    for (k=0;k<nb;k++) {
        dlam[i+k*nc] += bX[k]*z;
    }
    // and the m_chw (theta[nb+1]) dependency
    dlam[i+nc*(nb+1)] += (I[nc-2]+I[nc-1]);
    
  } // deriv 
  // and care home residents
  i = nc-3; // 80+ row/col
  for (z=0,j=0,Ci=C+(nc-1)*i;j<nc-2;j++,Ci++) z += *Ci * I[j];
  lam[nc-1] = eps * *beta * z + m_chw * I[nc-2] + m_chr * I[nc-1];
  if (finf) {
     lam[nc-1] *= x[(nc-1)*ns]; // multiply lam[i] by S_i - rate of symptomatic caused infections
     for (z=0,j=0,Ci=C+(nc-1)*i;j<nc-2;j++,Ci++) z += *Ci * x[j*ns+4]; // 
     dlam[nc-1] = (eps* z * *beta + m_chw * x[(nc-2)*ns+4]+ m_chr*x[(nc-1)*ns+4])*x[(nc-1)*ns]; //daily rate of infections caused by symptomatics 
  }
  if (deriv) {
    // first the sensitivities from the I state dependencies...
    for (k=0;k<deriv;k++) { // parameter loop
      for (z1=0,j=0,Ci=C+(nc-1)*i;j<nc-2;j++,Ci++) z1 += *Ci * dI[k*nc+j];
      dlam[nc-1+k*nc] = z1 * *beta*eps + m_chw * dI[k*nc+nc-2] + m_chr * dI[k*nc+nc-1];
    }
    // now add the direct beta dependencies
    for (k=0;k<nb;k++) {
        dlam[nc-1+k*nc] += bX[k]*z*eps;
    }
    // direct eps (theta[nb])
    dlam[nc-1+nc*nb] += *beta*z;
    // direct m_chw (theta[nb+1])
    dlam[nc-1+nc*(nb+1)] += I[nc-2];
    // direct m_chr (theta[nb+2])
    dlam[nc-1+nc*(nb+2)] += I[nc-1];
  }  
} // lambda 

double egduf(double *dmu, double mu,double t,double t0,double t1) {
/* This is the care improvement term, multiplying mortality, basically. It drops linearly 
   from 1 on April 1 to mu on June 1. The justification given is a trial of Dexamethasone
   which gave about a 15% RR. The trial completed in July.  
*/
  double h,dt;
  if (t<t0) {
    h=1.0;*dmu=0.0;
  } else if (t>t1) {
    h=mu;*dmu=1.0;
  } else {
    dt = t1-t0;
    h = 1 + (1-mu)*(t0-t)/dt;
    *dmu=(t-t0)/dt;
  }    
  return(h);
} // egduf  

double p_s(double t) {
/* p_star is problematic - it is set as constant at .25, meaning that 25% of people are 
   assumed arriving with a confirmed diagnosis even in March. The cited source of the 
   figure...
https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
  ... does not appear to contain it. 
  this function simply ramps it up slowly.  
*/
  if (t<90) return(0.0);
  if (t>200) return(0.25);
  else return(0.25*(t-90)/110);
	      
}  

void f(double *g,double *x, double t,double *theta,int *iaux,double *daux) {
/* evaluate gradient of system 
   [nb,np,nc,ns,ng,nk,deriv,t_freeze] = iaux; 
   [C,gamma,knots,work] = daux;
     - C is (nc-1)*(nc-1) matrix/vector 
     - gamma is ng vector
     - baux is nk+1 vector of basis eval info
     - work is nb+2*(nc+nc*np) vector (includes lam, dlam)
 
   x[i*ns+j] is jth state variable in ith class
   x[(k+1)*ns*nc+i*ns+j] is deriv of the above w.r.t. param k.
*/
  int nb,np,nc,ns,ng,nk,deriv,deriv1,ins,nsc,k1,k1nsc,i,k,t_freeze,k1nscins;
  double *C,*work,*lam,*dlam,*gamma,gamma_E,gamma_A,gamma_C,gamma_G,p_c,dinoc,beta,
    *baux,inoc_sig,a1,a2,t0,*psi_h,*psi_ic,*psi_icd,*psi_wd,*psi_hd,p_h_max,p_g,dmu_ic,mu_ic,f_ic,dmu_d,mu_d,f_icd,
    p_star,f_start,f_end,gamma_IC_pre,gamma_U,gamma_IC_Wr,gamma_IC_Wd,gamma_IC_d,gamma_Wd,gamma_Wr,gamma_Hr,gamma_Hd,
    gamma_P,gamma_Po,gamma_S,p_S,gamma_ph,
    p_ic_max,p_icd_max,p_wd_max,p_hd_max,
    c_a,c_b,c_0,c_1,c_2,c_3,c_4,c_5;
  nb=iaux[0]; // number of beta params controlling b(t)
  np=iaux[1]; // length of theta - free param vector
  nc=iaux[2]; // number of classes
  ns=iaux[3]; // number of states per class
  ng=iaux[4]; // number of fixed parameters
  nk= iaux[5]; // number of knots
  deriv=iaux[6]; // compute derivatives?
  t_freeze=iaux[7]; // when,  if ever, to freeze b(t)
  nsc = ns*nc;
  C = daux; daux += (nc-1)*(nc-1);
  gamma = daux; daux += ng;
  baux = daux; daux += nk+1;
  lam = daux; daux += nc;
  dlam = daux; daux += nc*np;
  work = daux;
  // unpack fixed coeffs
  gamma_E = gamma[0];gamma_A = gamma[1];gamma_C=gamma[2];p_c=gamma[3];
  inoc_sig=gamma[4];gamma_ph=gamma[5];//gamma_G=gamma[5];
  p_star = p_s(t);// gamma[6];
  f_start=gamma[7];f_end = gamma[8];
  gamma_IC_pre= gamma[9];gamma_U= gamma[10];//gamma_IC_Wr = gamma[11];
  gamma_IC_Wd = gamma[12];
  gamma_IC_d = gamma[13];//gamma_Wr = gamma[14];
  gamma_Wd = gamma[15];//gamma_Hr=gamma[16];
  gamma_Hd = gamma[17];gamma_P = gamma[18];gamma_Po=gamma[19];gamma_S = gamma[20];p_S=gamma[21];
  k = 22;
  psi_h = gamma + k;psi_ic = gamma + k + nc;psi_icd = gamma + k + 2*nc;psi_wd=gamma + k + 3*nc;
  psi_hd = gamma + k + 4*nc;
  a1=1/(sqrt(4*asin(1.0))*inoc_sig);a2=1/(2*inoc_sig*inoc_sig); // inoc params
  // Obtain the force of infection vector and (if deriv!=0) sensitivities...
  if (deriv) deriv1 = nb+4; else deriv1=0;
  lambda(&beta,lam,dlam,x,theta,C,t,nb,np,nc,ns,nk,deriv1,baux,work,0,t_freeze);
  // unpack free parameters needed here...
  t0 = theta[nb+3];p_h_max = theta[nb+4];p_g = theta[nb+5];
  mu_ic=theta[nb+6];p_ic_max=theta[nb+7];p_icd_max=theta[nb+8];
  mu_d=theta[nb+9];p_wd_max=theta[nb+10];p_hd_max=theta[nb+11];
  gamma_G=theta[nb+12]; // care home progression to death - cited reference did not contain given rate, and timing wrong.
  gamma_IC_Wr=theta[nb+13];gamma_Wr=theta[nb+14];gamma_Hr=theta[nb+15]; // recovery rates in hospital - most clinical latitude
  /* The main SEIR dynamics... */
  for (i=0;i<nc;i++) { // class loop
    ins = i * ns;
    // The main SEIR dynamics...
    g[ins] = -lam[i]*x[ins]; // susceptible
    g[ins+1] = lam[i]*x[ins] - gamma_E*x[ins+1]; // Exposed 1
    g[ins+2] = gamma_E*(x[ins+1]-x[ins+2]); // Exposed 2
    g[ins+3] = (1-p_c)*gamma_E*x[ins+2] - gamma_A*x[ins+3]; // I_A infectious asymptomatic
    g[ins+4] = p_c*gamma_E*x[ins+2] - gamma_C*x[ins+4]; // I_C infectious symptomatic
    // pre-hospital state to allow separation of mean time to hospital and mean time
    // infectious... 
    g[ins+36] = gamma_C*x[ins+4] - gamma_ph*x[ins+36]; 
    if (i>2&&i<13) { // initialization of 1 into each class (10) from 15-65 I_A
      dinoc = a1*exp(-a2*(t-t0)*(t-t0));g[ins+3] += dinoc;
    }
    // The SEIR sensitivities...
    if (deriv) {
      for (k=0;k<nb+4;k++) { // parameter loop
	k1=k+1;k1nsc=k1*nsc;k1nscins=k1nsc+ins;
	// lambda involved...
	g[k1nscins] = -dlam[k*nc+i]*x[ins] - lam[i]*x[k1nscins]; // susceptible
	g[k1nscins+1] = dlam[k*nc+i]*x[ins] + lam[i]*x[k1nscins] - gamma_E*x[k1nscins+1];// Exposed 1
	// fixed flow boxes...
        g[k1nscins+2] = gamma_E*(x[k1nscins+1]-x[k1nscins+2]); // Exposed 2
        g[k1nscins+3] = (1-p_c)*gamma_E*x[k1nscins+2] - gamma_A*x[k1nscins+3]; // I_A infectious asymptomatic
        g[k1nscins+4] = p_c*gamma_E*x[k1nscins+2] - gamma_C*x[k1nscins+4]; // I_C infectious symptomatic
        g[k1nscins+36] = gamma_C*x[k1nscins+4] - gamma_ph*x[k1nscins+36]; // pre hospital
      }
      if (i>2&&i<13) { // add in direct effect of t0 (theta[nb+3])
        g[(nb+4)*nsc+ins+3] += dinoc*2*a2*(t-t0);
      }	
    }  
  } // class loop

  /* care home deaths */
  ins = (nc-1)*ns;
  c_1 = psi_h[nc-1]*gamma_C;
  g[ins+5] = p_h_max*c_1*p_g*x[ins+36] - gamma_G*x[ins+5]; // G_D^{i,1}
  g[ins+6] = gamma_G*x[ins+5] - gamma_G*x[ins+6];    // G_D^{i,2}
  if (deriv) { /* care home death sensitivities - up to nb+5 only matter */
    for (k=0;k<np;k++) { // parameter loop -- standard propagation
       k1=k+1;k1nsc=k1*nsc;k1nscins=k1nsc+ins;
       g[k1nscins+5] = p_h_max*c_1*p_g*x[k1nscins+36] - gamma_G*x[k1nscins+5]; 
       g[k1nscins+6] = gamma_G*x[k1nscins+5] - gamma_G*x[k1nscins+6];
    }
    k = nb+4;k1nsc=(k+1)*nsc; // p_h_max theta[nb+4]
    g[k1nsc+ins+5] += c_1*p_g*x[ins+36];
    k = nb+5;k1nsc=(k+1)*nsc; // p_g theta[nb+5]
    g[k1nsc+ins+5] += p_h_max*c_1*x[ins+36];
    k = nb+12;k1nsc=(k+1)*nsc; // gamma_G theta[nb+12]
    g[k1nsc+ins+5] += -x[ins+5];
    g[k1nsc+ins+6] += x[ins+5]-x[ins+6];
  }
  /* from here on only the first nc-1 classes are involved, not care home residents */

  /* the ICU sequence, states 7-16 */
  f_ic = egduf(&dmu_ic,mu_ic,t,f_start,f_end);
  f_icd = egduf(&dmu_d,mu_d,t,f_start,f_end);
  for (i=0;i<nc-1;i++) { // class loop
    ins = i * ns;
    c_0 =  p_h_max*psi_h[i]*psi_ic[i]*f_ic*p_ic_max*gamma_C;
    c_1 = (1-p_star)*c_0;
    g[ins+7] = c_1 * x[ins+36] -(gamma_IC_pre+gamma_U)*x[ins+7]; // ICU_pre^i
    c_2 = p_star*c_0;
    g[ins+8] = c_2 * x[ins+36] - gamma_IC_pre * x[ins+8] + gamma_U*x[ins+7]; // ICU_pre*^i
    c_0 = (1-p_icd_max*psi_icd[i]*f_icd)*gamma_IC_pre;      
    c_3 = c_0*(1-p_wd_max*psi_wd[i]*f_icd);
    g[ins+9] = c_3 * x[ins+7] - (gamma_IC_Wr + gamma_U)*x[ins+9];// ICU_Wr
    g[ins+10] = c_3 * x[ins+8] - gamma_IC_Wr * x[ins+10] + gamma_U*x[ins+9];// ICU_Wr*
    c_4 = c_0*p_wd_max*psi_wd[i]*f_icd;
    g[ins+11] = c_4*x[ins+7] - (gamma_IC_Wd+gamma_U)*x[ins+11]; // ICU_Wd
    g[ins+12] = c_4*x[ins+8] - gamma_IC_Wd * x[ins+12] + gamma_U* x[ins+11]; // ICU_Wd*
    c_5 = p_icd_max*psi_icd[i]*f_icd;
    g[ins+13] = c_5*x[ins+7] - (gamma_IC_d+gamma_U)*x[ins+13]; // ICU_D^{i,1}
    g[ins+14] = gamma_IC_d * x[ins+13] - (gamma_IC_d+gamma_U)*x[ins+14]; // ICU_D^{i,2}
    g[ins+15] = c_5*x[ins+8] - gamma_IC_d *x[ins+15] +gamma_U*x[ins+13]; // ICU_D*^{i,1}
    g[ins+16] = gamma_IC_d * x[ins+15] - gamma_IC_d*x[ins+16] +gamma_U*x[ins+14]; // ICU_D*^{i,2}								    
    if (deriv) {
      for (k=0;k<nb+14;k++) { // parameter loop -- standard propagation (nb+4/5 not needed actually)
         k1=k+1;k1nsc=k1*nsc;k1nscins=k1nsc+ins;
	 g[k1nscins+7] = c_1 * x[k1nscins+36] -(gamma_IC_pre+gamma_U)*x[k1nscins+7];
	 g[k1nscins+8] = c_2 * x[k1nscins+36] - gamma_IC_pre * x[k1nscins+8] + gamma_U*x[k1nscins+7];
	 g[k1nscins+9] = c_3 * x[k1nscins+7] - (gamma_IC_Wr + gamma_U)*x[k1nsc+ins+9];
	 g[k1nscins+10] = c_3 * x[k1nscins+8] - gamma_IC_Wr * x[k1nscins+10] + gamma_U*x[k1nscins+9];
	 g[k1nscins+11] = c_4*x[k1nscins+7] - (gamma_IC_Wd+gamma_U)*x[k1nscins+11];
	 g[k1nscins+12] = c_4*x[k1nscins+8] - gamma_IC_Wd * x[k1nscins+12] + gamma_U* x[k1nscins+11];
	 g[k1nscins+13] = c_5*x[k1nscins+7] - (gamma_IC_d+gamma_U)*x[k1nscins+13]; // ICU_D^{i,1}
         g[k1nscins+14] = gamma_IC_d * x[k1nscins+13] - (gamma_IC_d+gamma_U)*x[k1nscins+14]; // ICU_D^{i,2}
         g[k1nscins+15] = c_5*x[k1nscins+8] - gamma_IC_d *x[k1nscins+15] +gamma_U*x[k1nscins+13]; // ICU_D*^{i,1}
         g[k1nsc+ins+16] = gamma_IC_d * x[k1nscins+15] - gamma_IC_d*x[k1nsc+ins+16] +gamma_U*x[k1nsc+ins+14]; // ICU_D*^{i,2}
       }
       // direct dependencies on p_h_max (nb+4), mu_icu (nb+6) and p_icu_max (nb+7)
       k = nb+4;k1nsc=(k+1)*nsc; // p_h_max (nb+4): 7,8
       g[k1nsc+ins+7] += psi_h[i]*(1-p_star)*psi_ic[i]*f_ic*p_ic_max*gamma_C*x[ins+36];
       g[k1nsc+ins+8] += psi_h[i]*p_star*psi_ic[i]*f_ic*p_ic_max*gamma_C*x[ins+36];
       k = nb+6;k1nsc=(k+1)*nsc; //mu_ic (nb+6): 7,8
       g[k1nsc+ins+7] += p_h_max*psi_h[i]*(1-p_star)*psi_ic[i]*dmu_ic*p_ic_max*gamma_C*x[ins+36];
       g[k1nsc+ins+8] += p_h_max*psi_h[i]*p_star*psi_ic[i]*dmu_ic*p_ic_max*gamma_C * x[ins+36];
       k = nb+7;k1nsc=(k+1)*nsc; //p_ic_max (nb+7): 7,8
       g[k1nsc+ins+7] += p_h_max*psi_h[i]*(1-p_star)*psi_ic[i]*f_ic*gamma_C* x[ins+36];
       g[k1nsc+ins+8] += p_h_max*psi_h[i]*p_star*psi_ic[i]*f_ic*gamma_C* x[ins+36];
       // direct dependencies on p_icd_max=theta[nb+8] and mu_d=theta[nb+9] and p_wd_max=theta[nb+10]
       k = nb+8;k1nsc=(k+1)*nsc;k1nscins=k1nsc+ins;// p_icd_max (nb+8): 9,10,11,12,13,15
       g[k1nscins+9] += -psi_icd[i]*f_icd*(1-p_wd_max*psi_wd[i]*f_icd)*gamma_IC_pre*x[ins+7];
       g[k1nscins+10] += -psi_icd[i]*f_icd*(1-p_wd_max*psi_wd[i]*f_icd)*gamma_IC_pre*x[ins+8];
       g[k1nscins+11] += -psi_icd[i]*f_icd*p_wd_max*psi_wd[i]*f_icd*gamma_IC_pre*x[ins+7];
       g[k1nscins+12] += -psi_icd[i]*f_icd*p_wd_max*psi_wd[i]*f_icd*gamma_IC_pre*x[ins+8];
       g[k1nscins+13] += psi_icd[i]*f_icd*x[ins+7];
       g[k1nscins+15] += psi_icd[i]*f_icd*x[ins+8];
	 
       k = nb+9;k1nsc=(k+1)*nsc;k1nscins=k1nsc+ins;// mu_d (nb+9): 9,10,11,12,13,15
       g[k1nscins+9] += -p_icd_max*psi_icd[i]*dmu_d*(1-p_wd_max*psi_wd[i]*f_icd)*gamma_IC_pre* x[ins+7]
	                 -(1-p_icd_max*psi_icd[i]*f_icd)*p_wd_max*psi_wd[i]*dmu_d*gamma_IC_pre* x[ins+7];
       g[k1nscins+10] += -p_icd_max*psi_icd[i]*dmu_d*(1-p_wd_max*psi_wd[i]*f_icd)*gamma_IC_pre* x[ins+8]
	                  -(1-p_icd_max*psi_icd[i]*f_icd)*p_wd_max*psi_wd[i]*dmu_d*gamma_IC_pre* x[ins+8];
       g[k1nscins+11] += -p_icd_max*psi_icd[i]*dmu_d*p_wd_max*psi_wd[i]*f_icd*gamma_IC_pre*x[ins+7] +
	                   (1-p_icd_max*psi_icd[i]*f_icd)*p_wd_max*psi_wd[i]*dmu_d*gamma_IC_pre*x[ins+7];
       g[k1nscins+12] += -p_icd_max*psi_icd[i]*dmu_d*p_wd_max*psi_wd[i]*f_icd*gamma_IC_pre*x[ins+8] +
	                   (1-p_icd_max*psi_icd[i]*f_icd)*p_wd_max*psi_wd[i]*dmu_d*gamma_IC_pre*x[ins+8];
       g[k1nscins+13] += p_icd_max*psi_icd[i]*dmu_d*x[ins+7];
       g[k1nscins+15] += p_icd_max*psi_icd[i]*dmu_d*x[ins+8];
       
       k = nb+10;k1nsc=(k+1)*nsc;// p_wd_max (nb+10): 9,10,11,12 
       g[k1nsc+ins+9] += -(1-p_icd_max*psi_icd[i]*f_icd)*psi_wd[i]*f_icd*gamma_IC_pre* x[ins+7];
       g[k1nsc+ins+10] += -(1-p_icd_max*psi_icd[i]*f_icd)*psi_wd[i]*f_icd*gamma_IC_pre* x[ins+8];
       g[k1nsc+ins+11] += (1-p_icd_max*psi_icd[i]*f_icd)*psi_wd[i]*f_icd*gamma_IC_pre*x[ins+7];
       g[k1nsc+ins+12] += (1-p_icd_max*psi_icd[i]*f_icd)*psi_wd[i]*f_icd*gamma_IC_pre*x[ins+8];

       k = nb+13;k1nsc=(k+1)*nsc;//gamma_IC_Wr nb+13: 9, 10
       g[k1nsc+ins+9] +=  - x[ins+9];// ICU_Wr
       g[k1nsc+ins+10] +=  - x[ins+10] ;// ICU_Wr*
    }		
  }  // class loop ICU

  /* The step down from ICU to general ward sequence */
  for (i=0;i<nc-1;i++) { // class loop
    ins = i*ns;
    g[ins+17] = gamma_IC_Wr * x[ins+9] - (gamma_Wr+gamma_U)*x[ins+17]; // W_r^{i,1}
    g[ins+18] = gamma_Wr*x[ins+17] - (gamma_Wr+gamma_U)*x[ins+18]; // W_r^{i,2}
    g[ins+19] = gamma_IC_Wr * x[ins+10] - gamma_Wr * x[ins+19] + gamma_U*x[ins+17]; // W_r*^{i,1}
    g[ins+20] = gamma_Wr * x[ins+19] - gamma_Wr * x[ins+20] + gamma_U*x[ins+18]; // W_r*^{i,2}
    g[ins+21] = gamma_IC_Wd * x[ins+11] - (gamma_Wd + gamma_U)*x[ins+21]; // W_D^i
    g[ins+22] = gamma_IC_Wd * x[ins+12] - gamma_Wd * x[ins+22] + gamma_U*x[ins+21]; // W_D*^i
    if (deriv) {
      for (k=0;k<nb+15;k++) { // parameter loop -- standard propagation (nb+4/5 not needed actually)
         k1=k+1;k1nsc=k1*nsc;
	 g[k1nsc+ins+17] = gamma_IC_Wr * x[k1nsc+ins+9] - (gamma_Wr+gamma_U)*x[k1nsc+ins+17]; // W_r^{i,1}
         g[k1nsc+ins+18] = gamma_Wr*x[k1nsc+ins+17] - (gamma_Wr+gamma_U)*x[k1nsc+ins+18]; // W_r^{i,2}
         g[k1nsc+ins+19] = gamma_IC_Wr * x[k1nsc+ins+10] - gamma_Wr * x[k1nsc+ins+19] + gamma_U*x[k1nsc+ins+17]; // W_r*^{i,1}
         g[k1nsc+ins+20] = gamma_Wr * x[k1nsc+ins+19] - gamma_Wr * x[k1nsc+ins+20] + gamma_U*x[k1nsc+ins+18]; // W_r*^{i,2}
         g[k1nsc+ins+21] = gamma_IC_Wd * x[k1nsc+ins+11] - (gamma_Wd + gamma_U)*x[k1nsc+ins+21]; // W_D^i
         g[k1nsc+ins+22] = gamma_IC_Wd * x[k1nsc+ins+12] - gamma_Wd * x[k1nsc+ins+22] + gamma_U*x[k1nsc+ins+21]; // W_D*^i
      }
      // direct dependency...
      k = nb+13;k1nsc=(k+1)*nsc;//gamma_IC_Wr nb+13: 17,19
      g[k1nsc+ins+17] += x[ins+9]; // W_r^{i,1}	 
      g[k1nsc+ins+19] += x[ins+10]; // W_r*^{i,1}

      k = nb+14;k1nsc=(k+1)*nsc;//gamma_Wr nb+14: 17,18,19,20
      g[k1nsc+ins+17] +=  - x[ins+17]; // W_r^{i,1}
      g[k1nsc+ins+18] +=  x[ins+17] - x[ins+18]; // W_r^{i,2}
      g[k1nsc+ins+19] +=  - x[ins+19]; // W_r*^{i,1}
      g[k1nsc+ins+20] += x[ins+19] - x[ins+20]; // W_r*^{i,2}
    }
    
  } // class loop step down
  
  /* The general ward sequence */
  for (i=0;i<nc-1;i++) { // class loop
    ins = i*ns;
    c_a = p_h_max*psi_h[i]*(1-psi_ic[i]*f_ic*p_ic_max)*gamma_C;
    c_b = p_hd_max*psi_hd[i]*f_icd;
    c_0 = c_a*(1-c_b);  
    // c_0 = p_h_max*psi_h[i]*(1-psi_ic[i]*f_ic*p_ic_max)*(1-p_hd_max*psi_hd[i]*f_icd)*gamma_C;
     
    c_1 = (1-p_star)*c_0;

    g[ins+23] = c_1 * x[ins+36] -(gamma_Hr + gamma_U)*x[ins+23]; // H_r^i

    c_2 = p_star*c_0;

    g[ins+24] = c_2 * x[ins+36] + gamma_U *x[ins+23] - gamma_Hr * x[ins+24]; // H_{r*}^i
    c_0 = c_a*c_b;
    c_3 = c_0*(1-p_star);
    //c_3 = p_h_max*psi_h[i]*(1-p_star)*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*f_icd*gamma_C;

    g[ins+25] = c_3 * x[ins+36] - (gamma_Hd+gamma_U)*x[ins+25]; // H_D^{i,1} 
    g[ins+26] = gamma_Hd*x[ins+25] - (gamma_Hd+gamma_U)*x[ins+26];// H_D^{i,2}
    c_4 = c_0*p_star;
    //c_4 = p_h_max*psi_h[i]*p_star*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*f_icd*gamma_C;

    g[ins+27] = c_4 * x[ins+36] + gamma_U*x[ins+25] - gamma_Hd * x[ins+27]; //  H_{D*}^{i,1}
    g[ins+28] = gamma_Hd * x[ins+27] - gamma_Hd * x[ins+28] + gamma_U *x[ins+26];  //  H_{D*}^{i,2}
    if (deriv) {
      for (k=0;k<nb+16;k++) { // parameter loop -- standard propagation (nb+4/5 not needed actually)
         k1=k+1;k1nsc=k1*nsc;
	 g[k1nsc+ins+23] = c_1 * x[k1nsc+ins+36] -(gamma_Hr + gamma_U)*x[k1nsc+ins+23]; // H_r^i
	 g[k1nsc+ins+24] = c_2 * x[k1nsc+ins+36] + gamma_U *x[k1nsc+ins+23] - gamma_Hr * x[k1nsc+ins+24]; // H_{r*}^i
         g[k1nsc+ins+25] = c_3 * x[k1nsc+ins+36] - (gamma_Hd+gamma_U)*x[k1nsc+ins+25]; // H_D^{i,1}
         g[k1nsc+ins+26] = gamma_Hd*x[k1nsc+ins+25] - (gamma_Hd+gamma_U)*x[k1nsc+ins+26];// H_D^{i,2}
         g[k1nsc+ins+27] = c_4 * x[k1nsc+ins+36] + gamma_U*x[k1nsc+ins+25] - gamma_Hd * x[k1nsc+ins+27]; //  H_{D*}^{i,1}
         g[k1nsc+ins+28] = gamma_Hd * x[k1nsc+ins+27] - gamma_Hd * x[k1nsc+ins+28] + gamma_U *x[k1nsc+ins+26];  //  H_{D*}^{i,2}
       }
       // Direct depndencies on p_h_max (nb+4), mu_ic (nb+6), mu_d (nb+9), p_ic_max (nb+7), p_hd_max (nb+11)
       k=nb+4;k1nsc = (k+1)*nsc; // p_h_max=theta[nb+4]: 23,24,25,27
       g[k1nsc+ins+23] += psi_h[i]*(1-p_star)*(1-psi_ic[i]*f_ic*p_ic_max)*(1-p_hd_max*psi_hd[i]*f_icd)*gamma_C* x[ins+36];
       g[k1nsc+ins+24] += psi_h[i]*p_star*(1-psi_ic[i]*f_ic*p_ic_max)*(1-p_hd_max*psi_hd[i]*f_icd)*gamma_C * x[ins+36];
       g[k1nsc+ins+25] += psi_h[i]*(1-p_star)*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       g[k1nsc+ins+27] += psi_h[i]*p_star*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       k=nb+6;k1nsc = (k+1)*nsc; // mu_ic=theta[nb+6]: 23,24,25,27
       g[k1nsc+ins+23] += -p_h_max*psi_h[i]*(1-p_star)*psi_ic[i]*dmu_ic*p_ic_max*(1-p_hd_max*psi_hd[i]*f_icd)*gamma_C* x[ins+36];
       g[k1nsc+ins+24] += -p_h_max*psi_h[i]*p_star*psi_ic[i]*dmu_ic*p_ic_max*(1-p_hd_max*psi_hd[i]*f_icd)*gamma_C * x[ins+36];
       g[k1nsc+ins+25] += -p_h_max*psi_h[i]*(1-p_star)*psi_ic[i]*dmu_ic*p_ic_max*p_hd_max*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       g[k1nsc+ins+27] += -p_h_max*psi_h[i]*p_star*psi_ic[i]*dmu_ic*p_ic_max*p_hd_max*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       k=nb+9;k1nsc = (k+1)*nsc; // mu_d=theta[nb+9]: 23,24,25,27
       g[k1nsc+ins+23] += -p_h_max*psi_h[i]*(1-p_star)*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*dmu_d*gamma_C* x[ins+36];
       g[k1nsc+ins+24] += -p_h_max*psi_h[i]*p_star*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*dmu_d*gamma_C * x[ins+36];
       g[k1nsc+ins+25] += p_h_max*psi_h[i]*(1-p_star)*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*dmu_d*gamma_C * x[ins+36];
       g[k1nsc+ins+27] += p_h_max*psi_h[i]*p_star*(1-psi_ic[i]*f_ic*p_ic_max)*p_hd_max*psi_hd[i]*dmu_d*gamma_C * x[ins+36];
       k =nb+7;k1nsc = (k+1)*nsc; // p_ic_max=theta[nb+7]: 23,24,25,27
       g[k1nsc+ins+23] += -p_h_max*psi_h[i]*(1-p_star)*psi_ic[i]*f_ic*(1-p_hd_max*psi_hd[i]*f_icd)*gamma_C* x[ins+36];
       g[k1nsc+ins+24] += -p_h_max*psi_h[i]*p_star*psi_ic[i]*f_ic*(1-p_hd_max*psi_hd[i]*f_icd)*gamma_C * x[ins+36];
       g[k1nsc+ins+25] += -p_h_max*psi_h[i]*(1-p_star)*psi_ic[i]*f_ic*p_hd_max*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       g[k1nsc+ins+27] += -p_h_max*psi_h[i]*p_star*psi_ic[i]*f_ic*p_hd_max*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       k =nb+11;k1nsc = (k+1)*nsc; // p_hd_max=theta[nb+11]: 23,24,25,27
       g[k1nsc+ins+23] += -p_h_max*psi_h[i]*(1-p_star)*(1-psi_ic[i]*f_ic*p_ic_max)*psi_hd[i]*f_icd*gamma_C* x[ins+36];
       g[k1nsc+ins+24] += -p_h_max*psi_h[i]*p_star*(1-psi_ic[i]*f_ic*p_ic_max)*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       g[k1nsc+ins+25] += p_h_max*psi_h[i]*(1-p_star)*(1-psi_ic[i]*f_ic*p_ic_max)*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       g[k1nsc+ins+27] += p_h_max*psi_h[i]*p_star*(1-psi_ic[i]*f_ic*p_ic_max)*psi_hd[i]*f_icd*gamma_C * x[ins+36];
       k =nb+15;k1nsc = (k+1)*nsc; // gamma_Hr=theta[nb+15]: 23,24
       g[k1nsc+ins+23] +=  - x[ins+23]; // H_r^i
       g[k1nsc+ins+24] +=  - x[ins+24]; // H_{r*}^i
    }  
  } // class loop general ward
  /* Recovered */
  for (i=0;i<nc;i++) { // class loop
    ins = i*ns;
    g[ins+29] = gamma_A*x[ins+3] + (1 - p_h_max*psi_h[i])*gamma_C*x[ins+36] + gamma_Hr*(x[ins+23]+x[ins+24]) +
                gamma_Wr * (x[ins+18]+x[ins+20]); // R^i
    if (deriv) {
      for (k=0;k<np;k++) { // parameter loop -- standard propagation (nb+4/5 not needed actually)
         k1=k+1;k1nsc=k1*nsc;
	 g[k1nsc+ins+29] = gamma_A*x[k1nsc+ins+3] + (1 - p_h_max*psi_h[i])*gamma_C*x[k1nsc+ins+36] +
	        gamma_Hr*(x[k1nsc+ins+23]+x[k1nsc+ins+24]) + gamma_Wr * (x[k1nsc+ins+18]+x[k1nsc+ins+20]);
      }
      // direct dependence on p_h_max = theta[nb+4]
      k = nb+4;k1nsc = (k+1)*nsc;
      g[k1nsc+ins+29] += -psi_h[i]*gamma_C*x[ins+36];

      k = nb+14;k1nsc = (k+1)*nsc; // gamma_Wr
      g[k1nsc+ins+ins+29] += x[ins+18]+x[ins+20]; // R^i

      k = nb+15;k1nsc = (k+1)*nsc; // gamma_Hr
      g[k1nsc+ins+29] += x[ins+23]+x[ins+24];
    }    
  } // recovered class loop
  /* finally the 6 test compartments */
  for (i=0;i<nc;i++) { // class loop
    ins = i*ns;
    g[ins+30] = gamma_E*x[ins+2] - gamma_S*x[ins+30]; // T_S^i
    g[ins+31] = p_S*gamma_S*x[ins+30]; // T_{S+}^i
    g[ins+32] = (1-p_S)*gamma_S*x[ins+30]; // T_{S-}^i
    g[ins+33] = lam[i]*x[ins] - gamma_P*x[ins+33]; // T_P^i
    g[ins+34] = gamma_P*x[ins+33] - gamma_Po*x[ins+34]; // T_{P+}^i
    g[ins+35] = gamma_Po*x[ins+34]; // T_{P-}^i
    if (deriv) {
      for (k=0;k<nb+4;k++) { // parameter loop
	k1=k+1;k1nsc=k1*nsc;
	g[k1nsc+ins+30] = gamma_E*x[k1nsc+ins+2] - gamma_S*x[k1nsc+ins+30]; // T_S^i
        g[k1nsc+ins+31] = p_S*gamma_S*x[k1nsc+ins+30]; // T_{S+}^i
        g[k1nsc+ins+32] = (1-p_S)*gamma_S*x[k1nsc+ins+30]; // T_{S-}^i
        g[k1nsc+ins+33] = dlam[k*nc+i]*x[ins] + lam[i]*x[k1nsc+ins] - gamma_P*x[k1nsc+ins+33]; // T_P^i
        g[k1nsc+ins+34] = gamma_P*x[k1nsc+ins+33] - gamma_Po*x[k1nsc+ins+34]; // T_{P+}^i
        g[k1nsc+ins+35] = gamma_Po*x[k1nsc+ins+34]; // T_{P-}^i
      }
    }  
  }  
} // f  

void trans(double *y,double *x,double *theta,int *iaux,double *daux,int ny,double t) {
/* transformations y of state x required by output 
   y is an ny*(np+1) vector, ith transform in y[i]
   derivatives of y[i] w.r.t. theta[k] in y[i+(k+1)*ny] (0 based indices)
   On input ny contains ny*(np+1) 
*/
  double *C,*gamma,*baux,*lam,*dlam,*work,gamma_C,p_star,gamma_U,gamma_IC_d, gamma_IC_pre,gamma_IC_Wd,
    gamma_A,gamma_Wd,gamma_Hd,gamma_G,gamma_IC_Wr,gamma_Hr,gamma_Wr,*psi_h,p_h_max,Ia,Is,Itot,beta;
  int nb,np,nc,ns,ng,nk,deriv,nsc,i,k,ins,k1,k1nsc,k1ny,t_freeze;

  nb=iaux[0]; // number of beta params controlling b(t)
  np=iaux[1]; // length of theta - free param vector
  ny = ny/(np+1); /* length of y on input - reduce to number of output states */
  nc=iaux[2]; // number of classes
  ns=iaux[3]; // number of states per class
  ng=iaux[4]; // number of fixed parameters
  nk= iaux[5]; // number of knots
  deriv=iaux[6]; // compute derivatives?
  t_freeze = iaux[7]; // when to freeze b(t)
  nsc = ns*nc;
  C = daux; daux += (nc-1)*(nc-1);
  gamma = daux; daux += ng;
  baux = daux; daux += nk+1;
  lam = daux; daux += nc;
  dlam = daux; daux += nc*np;
  work = daux;
  // unpack fixed coeffs
  gamma_A = gamma[1];gamma_C=gamma[2];
  p_star = p_s(t);//gamma[6];
  gamma_IC_pre= gamma[9];gamma_U= gamma[10];gamma_IC_Wd = gamma[12];
  gamma_IC_d = gamma[13];gamma_Wd = gamma[15];gamma_Hd = gamma[17];
  //gamma_G=gamma[5];
 
  k = 22;
  psi_h = gamma + k;
  // unpack free parameters needed here...
  p_h_max = theta[nb+4]; gamma_G=theta[nb+12];
  gamma_IC_Wr=theta[nb+13];gamma_Wr=theta[nb+14];gamma_Hr=theta[nb+15]; // recovery rates in hospital - most clinical latitude
  for (i=0;i<(np+1)*ny;i++) y[i]=0.0; // clear the output 

  for (i=0;i<nc-1;i++) {
    ins = i*ns;
    y[0] += p_h_max*psi_h[i]*p_star*gamma_C*x[ins+36] + // X_adm
            gamma_U*(x[ins+7]+x[ins+9]+x[ins+11]+x[ins+13]+x[ins+14]+x[ins+17]+
	             x[ins+18]+x[ins+21]+x[ins+23]+x[ins+25]+x[ins+26]);
    y[1] += x[ins+8]+ x[ins+19]+ x[ins+20]+ x[ins+22]+ x[ins+24]+ x[ins+27]+ x[ins+28]; // X_hos
    y[2] += x[ins+10]+ x[ins+12]+ x[ins+15]+ x[ins+16]; // X_icu
    y[3] += gamma_IC_d*(x[ins+14]+x[ins+16])+gamma_Wd*(x[ins+21]+x[ins+22])+gamma_Hd*(x[ins+26]+x[ins+28]); // X_Hd
    y[10] += p_h_max*psi_h[i]*p_star*gamma_C*x[ins+36] + gamma_IC_Wr*(x[ins+10]+x[ins+12]) + // ward arrival
             gamma_U*(x[ins+7]+x[ins+17]+x[ins+18]+x[ins+21]+x[ins+23]+x[ins+25]+x[ins+26]);
    y[11] += gamma_IC_pre*x[ins+8] + gamma_Wr * x[ins+20] +  gamma_Wd * x[ins+22] + // ward departures
             gamma_Hr * x[ins+24] + gamma_Hd * x[ins+28];
    y[12] += gamma_IC_pre*x[ins+8] + gamma_U*(x[ins+9]+x[ins+11]+x[ins+13]+x[ins+14]); // ICU arrival
    y[13] +=  gamma_IC_Wr * x[ins+10] +  gamma_IC_Wd * x[ins+12] + gamma_IC_d*x[ins+16]; //ICU departure
    if (deriv) {
      for (k=0;k<np;k++) {
	k1 = k+1;k1nsc=k1*nsc;
	k1ny = k1*ny;
	y[k1ny] += p_h_max*psi_h[i]*p_star*gamma_C*x[k1nsc+ins+36] + 
            gamma_U*(x[k1nsc+ins+7]+x[k1nsc+ins+9]+x[k1nsc+ins+11]+x[k1nsc+ins+13]+x[k1nsc+ins+14]+x[k1nsc+ins+17]+
	             x[k1nsc+ins+18]+x[k1nsc+ins+21]+x[k1nsc+ins+23]+x[k1nsc+ins+25]+x[k1nsc+ins+26]);
	y[k1ny+1] += x[k1nsc+ins+8]+ x[k1nsc+ins+19]+ x[k1nsc+ins+20]+ x[k1nsc+ins+22]+ x[k1nsc+ins+24]+
	             x[k1nsc+ins+27]+ x[k1nsc+ins+28];
	y[k1ny+2] += x[k1nsc+ins+10]+ x[k1nsc+ins+12]+ x[k1nsc+ins+15]+ x[k1nsc+ins+16]; // X_icu
	y[k1ny+3] += gamma_IC_d*(x[k1nsc+ins+14]+x[k1nsc+ins+16])+gamma_Wd*(x[k1nsc+ins+21]+x[k1nsc+ins+22])
	             +gamma_Hd*(x[k1nsc+ins+26]+x[k1nsc+ins+28]);
	y[k1ny+10] += p_h_max*psi_h[i]*p_star*gamma_C*x[k1nsc+ins+36] + gamma_IC_Wr*(x[k1nsc+ins+10]+x[k1nsc+ins+12]) + // ward arrival
             gamma_U*(x[k1nsc+ins+7]+x[k1nsc+ins+17]+x[k1nsc+ins+18]+x[k1nsc+ins+21]+x[k1nsc+ins+23]+x[k1nsc+ins+25]+x[k1nsc+ins+26]);
	y[k1ny+11] += gamma_IC_pre*x[k1nsc+ins+8] + gamma_Wr * x[k1nsc+ins+20] +  gamma_Wd * x[k1nsc+ins+22] + // ward departures
             gamma_Hr * x[k1nsc+ins+24] + gamma_Hd * x[k1nsc+ins+28];
	y[k1ny+12] += gamma_IC_pre*x[k1nsc+ins+8] + gamma_U*(x[k1nsc+ins+9]+x[k1nsc+ins+11]+x[k1nsc+ins+13]+x[k1nsc+ins+14]); // ICU arrival
        y[k1ny+13] +=  gamma_IC_Wr * x[k1nsc+ins+10] +  gamma_IC_Wd * x[k1nsc+ins+12] + gamma_IC_d*x[k1nsc+ins+16]; //ICU departure
      }
      // direct effect of p_h_max
      k = nb+4;k1 = k+1; k1ny = k1*ny;
      y[k1ny] += psi_h[i]*p_star*gamma_C*x[ins+36];
      y[k1ny+10] += psi_h[i]*p_star*gamma_C*x[ins+36];
      // direct effect of gamma_IC_Wr (nb+13): 10,13
      k = nb+13;k1 = k+1; k1ny = k1*ny;
      y[k1ny+10] += x[ins+10]+x[ins+12];
      y[k1ny+13] += x[ins+10];
      // gamma_Wr (nb+14)
      k = nb+14;k1 = k+1; k1ny = k1*ny;
      y[k1ny+11] += x[ins+20];
      //gamma_Hr (nb+15)
      k = nb+15;k1 = k+1; k1ny = k1*ny;
      y[k1ny+11] += x[ins+24];
    }
  }
  ins = (nc-1)*ns; 
  y[4] = gamma_G*x[ins+6]; // X_Gd
  if (deriv) {
    for (k=0;k<np;k++) {
       k1 = k+1;k1nsc=k1*nsc;
       k1ny = k1*ny;
       y[k1ny+4] += gamma_G*x[k1nsc+ins+6];
    }
    k = nb+12;k1 = k+1;k1ny = k1*ny;
    y[k1ny+4] += x[ins+6]; // direct effect of gamma_G
  }
  for (i=3;i<12;i++) {
    ins = i*ns;
    y[5] += x[ins+31]; // X_S+ seropositives in testing
    if (deriv) {
      for (k=0;k<np;k++) {
	k1 = k+1;k1nsc=k1*nsc;k1ny = k1*ny;
	y[k1ny+5] += x[k1nsc+ins+31];
      }	
    }  
  }
  for (i=1;i<nc-1;i++) {
    ins = i*ns;
    y[6] += x[ins+34]; // X_R1+ - react 1 positives
    if (deriv) {
      for (k=0;k<np;k++) {
	k1 = k+1;k1nsc=k1*nsc;k1ny = k1*ny;
	y[k1ny+6] += x[k1nsc+ins+34];
      }	
    }  
  }
  // get incidence and R - derivatives not needed
  lambda(&beta,lam,dlam,x,theta,C,t,nb,np,nc,ns,nk,0,baux,work,1,t_freeze);
  Ia=Is=0; // infections caused by asymp and symp
  Itot=0; // total infections outside care homes
  for (i=0;i<nc-2;i++) {
    ins = ns*i;
    y[7] += lam[i]+dlam[i]; // daily infection rate (excluding from care homes)
    Ia += lam[i];Is +=dlam[i];
    Itot += x[ins+3] + x[ins+4];
  }
  if (Itot<1) y[8] = -1.0; else 
  y[8] = (Ia/gamma_A+Is/gamma_C)/Itot;
  y[9] = beta;
} // trans  

void RK4(double *yout,double *xout,double *x,double *theta,double *h,double *t0,int *n,
	 int *nx,int *ny,int *ostep,double *work,int *iaux,double *daux) {
/* basic RK4 solver for dx/dt = f(x,t,theta,iaux,daux)
   nx is dim of x; work should be nx*5 long; h is step length; t0 start time.
   ostep is steps between the n outputs, which are stored in columns of nx by n 
   matrix xout. yout is an ny*n matrix for transformed output.
   
*/
  double *k1,*k2,*k3,*k4,t,*xp,h2,h6;
  int i,j,k,ii;
  for (i=0;i < *nx * 5;i++) work[i] = 0.0;
  xp = work; work += *nx;
  k1 = work; work += *nx;k2 = work; work += *nx;
  k3 = work; work += *nx;k4 = work; work += *nx;
  t = *t0; h2 = 0.5 * *h; h6 = *h/6;
  for (k=i=0;i<*n;i++) {
    for (ii=0;ii<*ostep;ii++,k++) {
      f(k1,x,t,theta,iaux,daux);
      for (j=0;j<*nx;j++) xp[j] = x[j] + h2 * k1[j]; 
      f(k2,xp,t+h2,theta,iaux,daux);
      for (j=0;j<*nx;j++) xp[j] = x[j] + h2 * k2[j]; 
      f(k3,xp,t+h2,theta,iaux,daux);
      for (j=0;j<*nx;j++) xp[j] = x[j] + *h * k3[j];
      f(k4,xp,t + *h,theta,iaux,daux);
      for (j=0;j<*nx;j++) x[j] += h6*(k1[j]+2*(k2[j]+k3[j])+k4[j]);
      t += *h;
    }
    for (j=0;j<*nx;j++,xout++) *xout = x[j]; /* states output */
    trans(yout,x,theta,iaux,daux,*ny,t); /* transformed output and derivatives */
    yout += *ny;
  }
} /* RK4 */ 
