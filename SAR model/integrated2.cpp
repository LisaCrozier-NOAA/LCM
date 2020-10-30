// Separable covariance on lattice with AR1 structure in each direction.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(yShift); //the offset for the environmental and pitTotal start years
  DATA_IVECTOR(yr); //Vector of years
  DATA_IVECTOR(j); //vector of days
  DATA_VECTOR(s_n); //total number of fish observed on a day and year
  DATA_VECTOR(s_k); //total number of fish that survived for a given day and year
  DATA_IVECTOR(k); //Fish grouping based on hatchery or wild, in-river or transport
  DATA_INTEGER(nvar); //total number of env. vars.
  DATA_ARRAY(env); //environmental variables by year
  DATA_ARRAY(timing); //Compass timing data
  DATA_MATRIX(xm); //environmental variables for each fish group. This could be deprecated
  DATA_IVECTOR(marVars); //vector of which marine covariates to use
  DATA_SCALAR(sd); //observation error for the environmental state space model
  DATA_INTEGER(re_j); //day effect flag
  DATA_INTEGER(re_t); //year effect flag
  DATA_INTEGER(re_jt); //day X year effect flag
  DATA_INTEGER(cov_pars); //env. covariance flag
  DATA_INTEGER(retro); //flag for determining whether to get SDs for retrospective
  DATA_VECTOR(env_mu); //mean of the environmental data for the entire time series
  DATA_VECTOR(env_sc); //standard deviation for the environmental data for the entire time series
  DATA_VECTOR(env_mu_2000_2015); //mean of the env data from 2000-2015
  DATA_VECTOR(env_sc_2000_2015); //mean of the env data from 2000-2015
  DATA_INTEGER(calibration_flag); //flag for determining whether to get SDs for calibration
  DATA_MATRIX(calibration_timing); //compass timing for the calibration
  DATA_IVECTOR(calibration_years); //years the you want the calibration for 
  
  PARAMETER_VECTOR(mu_s); //mean survival fixed effect
  PARAMETER_VECTOR(frho_j); //temporal correlation for day effect
  PARAMETER_VECTOR(frho_t); //temporal correlation for year effect
  PARAMETER_VECTOR(frho1_jt); //temporal correlation for day in dayXyear interaction
  PARAMETER_VECTOR(frho2_jt); //temporal correlation for year in dayXyear interaction
  PARAMETER_VECTOR(fpsi_j); //SD of day of effect
  PARAMETER_VECTOR(fpsi_t); //SD of year effect
  PARAMETER_VECTOR(fpsi_jt); //SD of dayXyear interaction
  
  PARAMETER_MATRIX(eps_j);    //random effects for day, nj * nk; nk is the number of fish groups (H, W, in-river, trans)
  PARAMETER_MATRIX(eps_t);    //random effects for year, nj * nk
  PARAMETER_ARRAY(eps_jt);    //random effects for dayXyear interaction, nj * nt * nk
  
  PARAMETER_MATRIX(beta_mar); //coefficients for the marine variables


  PARAMETER(frho_x); //temporal correlation for env. variables
  PARAMETER_VECTOR(frho_Rx); //unstructured correlations
  PARAMETER_VECTOR(fpsi_x); //process error
  PARAMETER_ARRAY(eps_x);    //random effects for each env. variable, nt * ni, ni being the number of env covariates
  
  // 
  using namespace density;
  
  Type ff = 0.;  //joint likelihood

  //Get all of the requisite dimension data	
  int nk = mu_s.size();
  int ni = s_k.size();
  int nt = eps_t.rows();
  int nte = env.dim[0];
  int nj = eps_j.rows();
  int rb = beta_mar.rows();
  
  //Correlation coefficient and process variances for environmental data
  vector<Type> psi_x = exp(fpsi_x);

  Type rho_x = 1/(1+exp(-frho_x));
  UNSTRUCTURED_CORR_t<Type> Sigma_x(frho_Rx);
  if(cov_pars==1){
	  ff += AR1(rho_x,Sigma_x)(eps_x);
//     for(int t=1;t<nte;t++){
// 		  ff += Sigma_x(eps_x.col(t)-rho_x*eps_x.col(t)); //Don't do this. Doesn't use a sparse matrix
// 	  }	
  }

  //Estimated env. covariates
  matrix<Type> env_hat(nte,nvar);
  for(int j=0;j<nvar;j++){
    for(int t=0;t<nte;t++){
      env_hat(t,j) = eps_x(j,t)*psi_x(j);
      if(env(t,j)>-100 & cov_pars==1){
        ff -= dnorm(env_hat(t,j), env(t,j), sd, true);
      }
    }
  }

  //Correlation and SD parameters for random effects. Atan seems to work better than logit.
  vector<Type> rho_j = atan(frho_j)*2./3.154;//1/(1+exp(-frho_j));//
  vector<Type> rho_t = atan(frho_t)*2./3.154;//1/(1+exp(-frho_t));//
  vector<Type> rho1_jt = atan(frho1_jt)*2./3.154;//1/(1+exp(-frho1_jt));//
  vector<Type> rho2_jt = atan(frho2_jt)*2./3.154;//1/(1+exp(-frho2_jt));//
  vector<Type> psi_j = exp(fpsi_j);
  vector<Type> psi_t = exp(fpsi_t);
  vector<Type> psi_jt = exp(fpsi_jt);

  
  // Random effect for day*year
  for(int kk=0;kk<nk;kk++){
    array<Type> tmp_eps(nj,nt);
    for(int tt=0;tt<nt;tt++){
      for(int jj=0;jj<nj;jj++){
        tmp_eps(jj,tt) = eps_jt(kk,jj,tt);
      }
    }
    if(re_jt==1){
      ff += AR1(rho1_jt(kk),AR1(rho2_jt(kk)))(tmp_eps);		
    }
  }

  
  //Get the subset of the marine variable that you want to model
  matrix<Type> xm_tmp(ni,rb);
  matrix<Type> eMar(ni,nk);
  for(int rr=0;rr<rb;rr++){
    xm_tmp.col(rr) = xm.col(marVars(rr)); //subset the marine effect you need
  }
  //Model the temporal processes and get the marine effects
  for(int kk=0;kk<nk;kk++){
    if(re_j==1){
      ff += AR1(rho_j(kk))(eps_j.col(kk)); //day effect		
    }
    if(re_t==1){
      ff += AR1(rho_t(kk))(eps_t.col(kk)); //year effect
    }
    eMar.col(kk) = xm_tmp*beta_mar.col(kk); //Marine effects
  }
  
  vector<Type> s_hat(ni);
  for(int i=0;i<s_k.size();i++){
    Type nu = mu_s(k(i)) +        //Mean
      eMar(i) +         //Annual environmental effects
      eps_t(yr(i) - yShift,k(i)) * psi_t(k(i)) +        //1D AR1 year
      eps_j(j(i),k(i)) * psi_j(k(i)) +       //1D AR1 day
      eps_jt(k(i),j(i),yr(i)-yShift) * psi_jt(k(i)); //2D AR1XAR1 yearXday

    s_hat(i) = exp(nu)/(1+exp(nu)); //Logit
    ff -= dbinom(s_k(i),s_n(i),s_hat(i),true); //likelihood 
  }


  vector<Type> s_t(nt); //Retrospective survival
  matrix<Type> s_tj(nt,nj); //Retrospective survival
  matrix<Type> xm_tmp2(nte,rb); //Retrospective environmental variables - they have a different mean and SD
  for(int rr=0;rr<rb;rr++){
    for(int t=0;t<nte;t++){
      //Un-zscore
      xm_tmp2(t,rr) = env_hat(t,marVars(rr)) * env_sc(marVars(rr)) + env_mu(marVars(rr));
      //Zscore for the env_mu_2000_2015 to env_sc_2000_2015 mean and sd
      xm_tmp2(t,rr) = (xm_tmp2(t,rr) - env_mu_2000_2015(marVars(rr)))/env_sc_2000_2015(marVars(rr));
    }
  }
  
  vector<Type> eMar2(nt); //Marine effects for retrospective
  eMar2 = xm_tmp2*beta_mar.col(0);
  s_t = 0.;
  if(retro){
    for(int tt=0;tt<nte;tt++){
      for(int jj=0;jj<nj;jj++){
        Type nu = mu_s(0) +        //Mean
          eMar2(tt) +         //Annual environmental effects
          eps_t(tt,0) * psi_t(0) +        //1D AR1 year
          eps_j(jj,0) * psi_j(0) +       //1D AR1 day
          eps_jt(0,jj,tt) * psi_jt(0); //2D AR1XAR1 yearXday
        
        s_t(tt) += nu * timing(jj,tt); //annual estimate with weighted effect
        s_tj(tt,jj) = nu; //survival by day and year
      }
    }
    ADREPORT(s_t);   //yearly survival
  }

  vector<Type> s_t_calibration(calibration_years.size());
  s_t_calibration.fill(0.); //Calibration survival
  if(calibration_flag){
    for(int tt=0;tt<calibration_years.size();tt++){
      int t_tmp = calibration_years(tt);
      for(int jj=0;jj<nj;jj++){
        Type nu = mu_s(0) +        //Mean
          eMar2(t_tmp) +         //Annual environmental effects
          eps_t(t_tmp,0) * psi_t(0) +        //1D AR1 year
          eps_j(jj,0) * psi_j(0) +       //1D AR1 day
          eps_jt(0,jj,t_tmp) * psi_jt(0); //2D AR1XAR1 yearXday
        
        s_t_calibration(tt) += nu * calibration_timing(tt,jj);
      }
    }
    ADREPORT(s_t_calibration);   //yearly survival
  }
  
  //REport everything
  REPORT(mu_s);
  REPORT(s_t);
  REPORT(frho_j); 
  REPORT(frho_t); 
  REPORT(frho1_jt); 
  REPORT(frho2_jt); 
  REPORT(fpsi_j); 
  REPORT(fpsi_t); 
  REPORT(fpsi_jt); 

  REPORT(eps_j);    
  REPORT(eps_t);    
  REPORT(eps_jt);    
  REPORT(beta_mar);
  REPORT(s_hat);
  REPORT(nj);
  REPORT(nt);
  REPORT(nte);
  REPORT(nk);
  REPORT(xm);

  REPORT(calibration_timing);
  REPORT(psi_x);
  REPORT(rho_x);
  REPORT(Sigma_x.cov());
  REPORT(ff);
  REPORT(eps_x);
  REPORT(env);
  REPORT(env_hat);
  REPORT(s_t_calibration);
  REPORT(env_mu_2000_2015);
  REPORT(env_sc_2000_2015);
  
//  if(!retro){
    ADREPORT(mu_s);   //mean survival
    ADREPORT(beta_mar); //fixed marine effects
    ADREPORT(eps_j);    //random effects for day
    ADREPORT(eps_t);    //random effects for year
    ADREPORT(eps_jt);    //random effects for dayXyear
//  }

  
  return ff;
}
