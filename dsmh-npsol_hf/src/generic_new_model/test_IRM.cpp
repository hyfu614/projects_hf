
#include <fstream>
#include <iomanip>
#include <sstream>

#include <math.h>

//#include "TData.hpp"
#include "IRM.hpp"
#include "specification_io.hpp"
#include "sbvar.hpp"
#include "test_IRM_SBVAR.hpp"

using namespace std;
using namespace DM;
using namespace TS;

void GeneralVAR(void)
{
  // file id 
  unsigned int nfile=3; 
  string specification_file="IRM_SBVAR_3var_4lag_restricted.txt";
  // sizes
  unsigned int nvar=3, nlag_original=4, nobs=100, truncate=10, nlag;

  cout << "Number variables: " << nvar << "\nNumber lags: " << nlag_original << "\nNumber observations: " << nobs << endl;
  if (truncate > 0)
    cout << "Number decimal places: " << truncate/10 << endl;
  else
    cout << "No truncation\n";

  // generate random parameters for stationary VAR
  TDM A=RandomNormal(nvar,nvar), F=RandomNormal(nvar,nvar*nlag_original+1);
  if (nlag_original == 0)
    {
      nlag=1;
      F=HCat(Zeros(nvar,nvar),F);
    }
  else
    {
      nlag=nlag_original;
      TDM B=VCat(A % F(I(0,End),I(0,nvar*nlag-1)),HCat(Identity(nvar*(nlag-1)),Zeros(nvar*(nlag-1),nvar)));
      TDV re, im;
      Eig(re,im,B);
      double mae=0.0;
      for (unsigned int i=re.Dim(); i-- > 0; )
	if (re[i]*re[i]+im[i]*im[i] > mae) mae=re[i]*re[i]+im[i]*im[i];
      if (mae > 1.0)
	for (unsigned int i=0; i < nlag; i++)
	  F(I(0,End),I(i*nvar,End),pow(0.98/sqrt(mae),i+1)*F(I(0,End),I(i*nvar,(i+1)*nvar-1)));

      // test if roots are all less than or equal to one
      B=VCat(A % F(I(0,End),I(0,nvar*nlag-1)),HCat(Identity(nvar*(nlag-1)),Zeros(nvar*(nlag-1),nvar)));
      Eig(re,im,B);
      mae=0.0;
      for (unsigned int i=re.Dim(); i-- > 0; )
	if (re[i]*re[i]+im[i]*im[i] > mae) mae=re[i]*re[i]+im[i]*im[i];
      cout << "Maximum root of system: " << sqrt(mae) << endl;
    }

  // reduced form companion matrix
  TDM B=VCat(VCat(A % F,HCat(Identity(nvar*(nlag-1)),Zeros(nvar*(nlag-1),nvar+1))),One(nvar*nlag,nvar*nlag+1));
  //cout << "B: \n" << B << endl;
  //===============================================================================
  //=== generate data
  //===============================================================================
  // initial conditions
  TDV e=Zeros(nvar*nlag+1);
  TDV x0=RandomNormal(nvar*nlag+1);
  x0(nvar*nlag)=1.0;
  for (unsigned int t=0; t < nobs; t++)
    x0=B*x0+e(I(0,End),A % RandomNormal(nvar));

  // data and predetermined data
  TDM Y(nobs,nvar);
  TDM X(nobs,nvar*nlag+1);
  TDV x=x0;
  for (unsigned int t=0; t < nobs; t++)
    {
      X(t,I(0,End),x);
      x=B*x+e(I(0,End),A % RandomNormal(nvar));
      Y(t,I(0,End),x,I(0,nvar-1));
    }

  // trucate data and predetermined data
  if (truncate > 0)
    for (unsigned int t=0; t < nobs; t++)
      {
	for (unsigned int i=0; i < nvar; i++) Y(t,i)=round(Y(t,i)*truncate)/truncate;
	for (unsigned int i=0; i < nvar*nlag+1; i++) X(t,i)=round(X(t,i)*truncate)/truncate;
      }
  cout << "data generated\n";

  // write data
  stringstream filename;
  ofstream output;
  filename.str("");
  filename << "SBVAR_data_" << nvar << "var_" << nlag << "lag" << nfile << ".txt";
  output.open(filename.str().c_str());
  if (!output.is_open())
     throw dw_exception("GeneralVAR() - unable to open " + filename.str());
  TDV y(nvar);
  for (int t=nlag-1; t >= 0; t--)
    {
      y(I(0,End),X,0,I(t*nvar,(t+1)*nvar-1));
      output << y;
    }
  output << Y;
  output.close();

  unsigned int npre=nvar*nlag+1;
  // Setup SBVAR model
  TSpecification spec(specification_file);
  // get SBVAR contemporanous restriction matrix
  TDM AC;
  vector< vector<bool> > AR;
  spec.RestrictionMatrix(AR,AC,"//== SBVAR: Simple Restrictions A",' ');
  if ((AC.Rows() != nvar) || (AC.Cols() != nvar))
     throw dw_exception("TTimeSeries_SBVAR_Specification() - invalid size of contemporaneous restriction matrix");
  if (AC != 0.0)
     throw dw_exception("TTimeSeries_SBVAR_Specification() - restrictions must be linear");

  // get SBVAR predetermined restriction matrix
  TDM FC;
  vector< vector<bool> > FR;
  spec.RestrictionMatrix(FR,FC,"//== SBVAR: Simple Restrictions F",' ');
  if ((FC.Rows() != nvar) || (FC.Cols() != npre))
     throw dw_exception("TTimeSeries_SBVAR_Specification() - invalid size of predetermined restriction matrix");
  if (FC != 0.0)
     throw dw_exception("TTimeSeries_SBVAR_Specification() - restrictions must be linear");

  // create SBVAR equation mappings
  TPArray<TLinear> Fnct;
  TDM parameters=HCat(AC,FC);
  for (unsigned int i=0; i < nvar; i++)
    {
	I idx;
	for (unsigned int j=0; j < nvar; j++)
	  if (AR[i][j]) idx(j);
	for (unsigned int j=0; j < npre; j++)
	  if (FR[i][j]) idx(nvar+j);
	Fnct.push(new TLinear_insert(nvar+npre,idx));
    }
  TLinearProduct Fn(Fnct);

  TDV Hyperparameters=spec.Vector("//== SBVAR: Prior Hyperparameters");
  TLikelihood_SBVAR sbvar_likelihood(Y,X);
  //TTimeSeries_SBVAR sbvar(sbvar_likelihood,TFlatPrior(sbvar_likelihood.NumberParameters()),TIdentityFunction(sbvar_likelihood.NumberParameters()));
  TTimeSeries_ConjugateSBVAR sbvar(Y,X,TDensity_SimsZha(Y,X,Fn,1,Hyperparameters));

  // Setup IRM
  TPArray<TRealFunction> f0(nvar);
  for (unsigned int i=0; i < nvar; i++)
    f0.push(state_space_response(nvar*nlag+1,i));
  TPArray<TRealFunction> s(nvar*nvar);

  for (unsigned int j=0; j < nvar; j++)
    for (unsigned int i=0; i < nvar; i++)
      s.push(state_space_response(nvar*nlag+1,i));

  TLikelihood_IRM irm_likelihood(Y,f0,s,0);
  TTimeSeries_IRM irm(irm_likelihood,TFlatPrior(irm_likelihood.NumberParameters()),TIdentityFunction(irm_likelihood.NumberParameters()));

  // Set SBVAR parameters
  TDV sbvar_parameters=sbvar.GetLikelihoodParameters(A,F);
  //cout << "sbvar_parameters: " << sbvar_parameters << endl;
  cout << "log likelihood (sbvar): " << sbvar.LogLikelihood(sbvar_parameters) << endl;
  filename.str("");
  filename << "SBVAR_parameters_" << nvar << "var_" << nlag << "lag" << nfile << ".txt";
  output.open(filename.str().c_str());
  if (!output.is_open())
     throw dw_exception("GeneralVAR() - unable to open " + filename.str());
  output << "SBVAR_parameters:\n" << sbvar_parameters << endl;
  output << "SBVAR_logMDD: " << sbvar.LogMDD();
  output.close();
  
  // Set IRM parameters
  TDV irm_parameters(irm.NumberParameters());
  TDM invA=Inverse(A);
  TDV v=Cat(Vec(B),X(0,I(0,End)));
  int k=0;
  for (unsigned int i=0; i < nvar; k+=v.Dim(), i++) irm_parameters(I(k,End),v);
  for (unsigned int j=0; j < nvar; j++)
    {
      v=Cat(Vec(B),e(I(0,nvar-1),invA(I(0,End),j)));
      for (unsigned int i=0; i < nvar; k+=v.Dim(), i++) irm_parameters(I(k,End),v);
    }
  //cout << "Number of irm parameters: " << irm_parameters.Dim() << endl;
  //cout << "irm_parameters: " << irm_parameters << endl;
  cout << "log likelihood (irm): " << irm.LogLikelihood(irm_parameters) << endl;

  //test SBVAR to IRM parameters mapping function
  TTimeSeries_IRM irm2=TTimeSeries_IRM_SBVAR_Specification(specification_file);
  TDV irm2_parameters=irm2.LikelihoodParameters(sbvar_parameters);
  //cout << "irm2_parameters: " << irm2_parameters << endl;
  cout << "Diff between parameters: " << Sum(irm_parameters-irm2_parameters) << endl;

  // no dynamics?
  if (nlag_original == 0)
    {
      // A*y(t) = c + epsilon(t)
      //
      // initial forecast
      //
      //   c
      //
      // impulse response
      //
      //   inverse(A) if t = 0
      //       0      it t > 0
      //

      unsigned int sigma_type=1;

      TPArray<TRealFunction> f0_nd(nvar);
      for (unsigned int i=0; i < nvar; i++)
	f0_nd.push(new polynomial_exponential_linear_trend(0));
      TPArray<TRealFunction> s_nd(nvar*nvar);
      for (unsigned int i=0; i < nvar*nvar; i++)
	s_nd.push(new polynomial_exponential_linear_trend(0));
      TLikelihood_IRM irm_likelihood_nd(Y,f0_nd,s_nd,sigma_type);
      TTimeSeries_IRM irm_nd=TTimeSeries_IRM(irm_likelihood_nd,TFlatPrior(irm_likelihood_nd.NumberParameters()),TIdentityFunction(irm_likelihood_nd.NumberParameters()));

      TDV sigma;
      TDM Sigma;
      if (sigma_type == 0)
	Sigma=Identity(nvar);
      else if (sigma_type == 1)
	Sigma=Diag(sigma=RandomNormal(nvar));
      else
	Sigma=Reshape(sigma=RandomNormal(nvar*nvar),nvar,nvar);

      TDV irm_nd_parameters(irm_nd.NumberParameters());
      unsigned int k=0;
      TDV v(4,0.0);
      TDV c=A % F(I(0,End),nvar);
      for (unsigned int i=0; i < nvar; k+=4, i++)
	{
	  v(3)=c(i);
	  irm_nd_parameters(I(k,End),v);
	}
      v(2)=10000.0;
      TDM inv_sigmaA=Inverse(Sigma*A);
      for (unsigned int j=0; j < nvar; j++)
	for (unsigned int i=0; i < nvar; k+=4, i++)
	  {
	    v(3)=inv_sigmaA(i,j);
	    irm_nd_parameters(I(k,End),v);
	  }
      irm_nd_parameters(I(k,End),sigma);

      cout << "log likelihood (irm-nd): " << irm_nd.LogLikelihood(irm_nd_parameters) << endl;

      vector<TDV> irm_initial_forecast=irm.InitialForecast(irm_parameters,nobs+1);
      vector<TDV> irm_nd_initial_forecast=irm_nd.InitialForecast(irm_nd_parameters,nobs+1);
      vector<TDM> irm_impulse_response=irm.ImpulseResponse(irm_parameters,nobs);
      vector<TDM> irm_nd_impulse_response=irm_nd.ImpulseResponse(irm_nd_parameters,nobs);
      double diff=0.0;
      for (unsigned int t=0; t < nobs; t++)
	{
	  diff+=Norm(irm_nd_impulse_response[t]*Sigma - irm_impulse_response[t]);
	  if (t > 0) diff+=Norm(irm_nd_initial_forecast[t] - irm_initial_forecast[t]);
	  //cout << "t = " << t << endl << irm_nd_initial_forecast[t] << endl << irm_initial_forecast[t] << endl << irm_nd_impulse_response[t]*Sigma << endl << irm_impulse_response[t] << endl;
	}
      cout << "difference initial forecasts and impulse responses (irm/irm_nd): " << diff << endl;
    }
  
  // AR(1) process?
  if ((nvar == 1) && (nlag_original == 1)) 
    {    
      //
      // a*y(t) = f*y(t-1) + c + epsilon(t)
      //
      // f/a must be positive
      //
      // initial forecast:
      //
      //  (f/a)^t*y(0) + (1 + f/a + ... + (f/a)^(t-1))*(c/a)
      //    = (f/a)^t*y(0) + ((1-(f/a)^t)/(1-f/a))*(c/a)
      //    = (f/a)^t*(y(0) - (c/a)/(1-f/a)) + (c/a)/(1-f/a)
      //
      // impulse response:
      //
      //  (f/a)^t*(1/a)
      //

      if (A(0,0)*F(0,0) <= 0.0)
	cout << "A(0,0)*F(0,0) must be positive\n";
      else
	{
	  TPArray<TRealFunction> f0_ar(nvar);
	  for (unsigned int i=0; i < nvar; i++)
	    f0_ar.push(new polynomial_exponential_linear_trend(0));
	  TPArray<TRealFunction> s_ar(nvar*nvar);
	  for (unsigned int i=0; i < nvar*nvar; i++)
	    s_ar.push(new polynomial_exponential_linear_trend(0));
	  TLikelihood_IRM irm_likelihood_ar(Y,f0_ar,s_ar,true);
	  TTimeSeries_IRM irm_ar=TTimeSeries_IRM(irm_likelihood_ar,TFlatPrior(irm_likelihood_ar.NumberParameters()),TIdentityFunction(irm_likelihood_ar.NumberParameters()));

	  if (irm_ar.NumberParameters() != 9)
	    {
	      cout << "number of irm_ar parameters: " << irm_ar.NumberParameters() << endl;
	      throw dw_exception("Invalid number of irm_ar parameters");
	    }
	  TDV irm_ar_parameters(irm_ar.NumberParameters());
	  double d=(F(0,1)/A(0,0))/(1.0 - F(0,0)/A(0,0)), beta=-log(F(0,0)/A(0,0)), sigma=3.3345; //1.0;
	  irm_ar_parameters(0)=d;
	  irm_ar_parameters(1)=0.0;
	  irm_ar_parameters(2)=beta;
	  irm_ar_parameters(3)=X(0,0);
      
	  irm_ar_parameters(4)=0.0;
	  irm_ar_parameters(5)=0.0;
	  irm_ar_parameters(6)=beta;
	  irm_ar_parameters(7)=1.0/(A(0,0)*sigma);
      
	  irm_ar_parameters(8)=sigma;

	  cout << "log likelihood (irm-ar): " << irm_ar.LogLikelihood(irm_ar_parameters) << endl;

	  vector<TDV> irm_initial_forecast=irm.InitialForecast(irm_parameters,nobs+1);
	  vector<TDV> irm_ar_initial_forecast=irm_ar.InitialForecast(irm_ar_parameters,nobs+1);
	  vector<TDM> irm_impulse_response=irm.ImpulseResponse(irm_parameters,nobs);
	  vector<TDM> irm_ar_impulse_response=irm_ar.ImpulseResponse(irm_ar_parameters,nobs);
	  double diff=0.0;
	  for (unsigned int t=0; t < nobs; t++)
	    {
	      diff+=Norm(irm_ar_impulse_response[t]*sigma - irm_impulse_response[t]);
	      if (t > 0) diff+=Norm(irm_ar_initial_forecast[t] - irm_initial_forecast[t]);
	      //cout << "t = " << t << endl << irm_ar_initial_forecast[t] << endl << irm_initial_forecast[t] << endl << irm_ar_impulse_response[t]*sigma << endl << irm_impulse_response[t] << endl;
	    }
	  cout << "difference initial forecasts and impulse responses (irm/irm_ar): " << diff << endl;
	}
    }

  // cout << "sbvar log conditional likelihoods:\n" << sbvar.LogConditionalLikelihoodVector(sbvar_parameters) << endl;
  // cout << "irm log conditional likelihoods:\n" << irm.LogConditionalLikelihoodVector(irm_parameters) << endl;
  
  // This block of code should at some point go in sbvar.*, but the problem is
  // with the deterministic component.
  // Get SBVAR initial forecast and impulse responses
  vector<TDV> sbvar_initial_forecast(nobs+1);
  sbvar_initial_forecast[0]=Zeros(nvar);
  x=X(0,I(0,End));
  for (unsigned int t=0; t < nobs; t++)
    {
      sbvar_initial_forecast[t]=x(I(0,nvar-1));
      x=B*x;
    }
  sbvar_initial_forecast[nobs]=x(I(0,nvar-1));

  vector<TDM> sbvar_impulse_response(nobs);
  TDM ir=VCat(Inverse(A),Zeros(nvar*(nlag-1)+1,nvar));
  for (unsigned int t=0; t < nobs; t++)
    {
      sbvar_impulse_response[t]=ir(I(0,nvar-1),I(0,nvar-1));
      ir=B*ir;
    }

  vector<TDV> irm_initial_forecast=irm.InitialForecast(irm_parameters,nobs+1);
  vector<TDM> irm_impulse_response=irm.ImpulseResponse(irm_parameters,nobs);
  double diff=0.0;
  for (unsigned int t=0; t < nobs; t++)
    {
      diff+=Norm(sbvar_impulse_response[t] - irm_impulse_response[t]);
      if (t > 0) diff+=Norm(sbvar_initial_forecast[t] - irm_initial_forecast[t]);
      // cout << "t = " << t << endl << sbvar_initial_forecast[t] - irm_initial_forecast[t] << endl << sbvar_impulse_response[t] - irm_impulse_response[t] << endl;
    }
  cout << "difference initial forecasts and impulse responses (irm/sbvar): " << diff << endl;
}

/* ==========================================
void TestGamma(void)
{
  unsigned int N=10000;
  char ch;
  do
    {
      cout << "gamma distribution\n";
      double alpha=3.0*dw_uniform_rnd(), beta=5.0*dw_uniform_rnd();
      //TDensity_gamma_transform gamma_dst(1,alpha,beta);
      TDensity_gamma_transform gamma_dst(TDV(1,alpha),TDV(1,beta));
      TFunction_halfline f(0.0,1,1);

      // draws
      TDV draws(N), g_draws(N);
      for (unsigned int i=N; i-- > 0; )
	{
	  TDV x=gamma_dst.Draw();
	  draws[i]=x(0);
	  g_draws[i]=f(x)(0);
	}

      // cumulative distributions
      double x0, x1, c, g_x0, g_x1, g_c, g_mid;
      double constant=1.0/(pow(beta,alpha)*exp(dw_log_gamma(alpha)));
      TDV sorted_draws=Sort(draws), g_sorted_draws=Sort(g_draws);
      TDV x(1);
      x1=sorted_draws[0];
      g_x1=0.5*g_sorted_draws[0];
      g_mid=0.5*g_x1;
      c=g_c=g_x1*constant*pow(g_mid,alpha-1)*exp(-g_mid/beta);
      for (unsigned int i=1; i < N; i++)
	{
	  x0=x1;
	  x1=0.5*(sorted_draws[i-1]+sorted_draws[i]);
	  x[0]=0.5*(x0+x1);
	  c+=(x1-x0)*exp(gamma_dst.LogDensity(x));

	  g_x0=g_x1;
	  g_x1=0.5*(g_sorted_draws[i-1]+g_sorted_draws[i]);
	  g_mid=0.5*(g_x0+g_x1);
	  g_c+=(g_x1-g_x0)*constant*pow(g_mid,alpha-1)*exp(-g_mid/beta);

	  if ((i % 1000) == 0)
	    cout << g_mid << " - " << x[0] << " : " << (double)i/(double)N << "  " << g_c << "  " << c << endl;
	}


      cout << "computed mean: " << Mean(g_draws) << "  analytic mean: " << alpha*beta << endl;
      cout << "computed standard deviation: " << StdDev(g_draws) << "   analytic standard deviation: " << sqrt(alpha)*beta << endl;
      cout << "alpha: " << alpha << "  beta: " << beta << endl;

      cout << "enter q to quit.\n";
      cin >> ch;
    }
  while (ch != 'q');

  do
    {
      // setup
      cout << "inverse gamma distribution\n";
      double alpha=5.0*dw_uniform_rnd(), beta=5.0*dw_uniform_rnd();
      //TDensity_inverse_gamma_transform inverse_gamma_dst(1,alpha,beta);
      TDensity_inverse_gamma_transform inverse_gamma_dst(TDV(1,alpha),TDV(1,beta));
      TFunction_halfline f(0.0,1,1);

      // draws
      TDV draws(N), ig_draws(N);
      for (unsigned int i=N; i-- > 0; )
	{
	  TDV x=inverse_gamma_dst.Draw();
	  draws[i]=x(0);
	  ig_draws[i]=f(x)(0);
	}

      // cumulative distributions
      double x0, x1, c, ig_x0, ig_x1, ig_c, ig_mid;
      double constant=1.0/(pow(beta,alpha)*exp(dw_log_gamma(alpha)));
      TDV sorted_draws=Sort(draws), ig_sorted_draws=Sort(ig_draws);
      TDV x(1);
      x1=sorted_draws[0];
      ig_x1=0.5*ig_sorted_draws[0];
      ig_mid=0.5*ig_x1;
      c=ig_c=ig_x1*constant*pow(ig_mid,-alpha-1)*exp(-1.0/(ig_mid*beta));
      for (unsigned int i=1; i < N; i++)
	{
	  x0=x1;
	  x1=0.5*(sorted_draws[i-1]+sorted_draws[i]);
	  x[0]=0.5*(x0+x1);
	  c+=(x1-x0)*exp(inverse_gamma_dst.LogDensity(x));

	  ig_x0=ig_x1;
	  ig_x1=0.5*(ig_sorted_draws[i-1]+ig_sorted_draws[i]);
	  ig_mid=0.5*(ig_x0+ig_x1);
	  ig_c+=(ig_x1-ig_x0)*constant*pow(ig_mid,-alpha-1)*exp(-1.0/(ig_mid*beta));

	  if ((i % 1000) == 0)
	    cout << ig_mid << " - " << x[0] << " : " << (double)i/(double)N << "  " << ig_c << "  " << c << endl;
	}

      // moments
      if (alpha > 1)
	cout << "computed mean: " << Mean(ig_draws) << "  analytic mean: " << 1.0/(beta*(alpha-1.0)) << endl;
      if (alpha > 2)
	cout << "computed standard deviation: " << StdDev(ig_draws) << "   analytic standard deviation: " << 1.0/(beta*(alpha-1.0)*sqrt(alpha-2.0)) << endl;
      cout << "alpha: " << alpha << "  beta: " << beta << endl;

      cout << "enter q to quit.\n";
      cin >> ch;
    }
  while (ch != 'q');
}
*/
//void GenerateData_alt(void)
//{
  // unsigned int nvar=2;
  // unsigned int nobs=100;
  // unsigned int truncate=100;

  // //=== setup irm model =======================================================
  // unsigned int degree=1;
  // bool is_sigma=false; TDM sigma=Identity(nvar);
  // TPArray<TRealFunction> f0(nvar);
  // for (unsigned int i=0; i < nvar; i++)
  //   f0.push(new polynomial_exponential_linear_trend(0));
  // TPArray<TRealFunction> s(nvar*nvar);
  // for (unsigned int i=0; i < nvar*nvar; i++)
  //   s.push(new polynomial_exponential_linear_trend(degree));
  // TLikelihood_IRM likelihood(RandomNormal(nobs,nvar),f0,s,is_sigma);
  // TTimeSeries_IRM irm(likelihood,TFlatPrior(likelihood.NumberParameters()),TIdentityFunction(likelihood.NumberParameters()));

  // //=== set parameters ========================================================
  // vector<TDV> initial_forecast_parameters(nvar);
  // vector< vector<TDV> > impulse_response_parameters(nvar);
  // for (unsigned int i=0; i < nvar; i++)
  //   {
  //     initial_forecast_parameters[i]=Zeros(4);
  //     for (unsigned int j=0; j < nvar; j++)
  // 	impulse_response_parameters[i].push_back(Zeros(4+degree));
  //   }

  // //---------------------------------------------------------------------------
  // // nvar=2   degree=1
  // initial_forecast_parameters[0](0)=2.5;  // constant trend
  // initial_forecast_parameters[0](3)=initial_forecast_parameters[0](0);
  
  // initial_forecast_parameters[1](1)=2;    // linear trend
  // initial_forecast_parameters[1](3)=initial_forecast_parameters[1](0);

  // impulse_response_parameters[0][0](2)=3.5;   // not so persistent
  // impulse_response_parameters[0][0](3)=1.0;
  // impulse_response_parameters[1][1](2)=0.05;  // relatively persistent
  // impulse_response_parameters[1][1](3)=-1.5;

  // // impulse_response_parameters[0][1](2)=0.9;  
  // // impulse_response_parameters[0][1](3)=-1.3;
  // // impulse_response_parameters[1][0](2)=1.1;  
  // // impulse_response_parameters[1][0](3)=-0.2;
  
  // impulse_response_parameters[0][1](2)=3*dw_uniform_rnd();  
  // impulse_response_parameters[0][1](3)=dw_gaussian_rnd();
  // impulse_response_parameters[1][0](2)=3*dw_uniform_rnd();   
  // impulse_response_parameters[1][0](3)=dw_gaussian_rnd();
  
  // //---------------------------------------------------------------------------
  
  // TDV parameters(irm.NumberParameters());
  // unsigned int k=0;
  // for (unsigned int i=0; i < nvar; k+=4, i++)
  //   parameters(I(k,End),initial_forecast_parameters[i]);
  // for (unsigned int j=0; j < nvar; j++)
  //   for (unsigned int i=0; i < nvar; k+=4+degree, i++)
  //     parameters(I(k,End),impulse_response_parameters[i][j]);
  // if (irm.IsSigma()) parameters(I(k,End),Vec(sigma));
							     
  // //=== generate data =========================================================
  // TDM Y=irm.GenerateData(parameters,1,nobs);

  // //=== truncate data =========================================================
  // if (truncate > 0)
  //   for (unsigned int t=0; t < nobs; t++)
  //     for (unsigned int i=0; i < nvar; i++)
  // 	Y(t,i)=round(Y(t,i)*truncate)/truncate;

  // //=== create new model with generated data ==================================
  // TLikelihood_IRM likelihood_data(Y,f0,s,is_sigma);
  // TTimeSeries_IRM irm_data(likelihood_data,TFlatPrior(likelihood_data.NumberParameters()),TIdentityFunction(likelihood_data.NumberParameters()));

  // //== print data =============================================================
  // cout << Y << endl;
  
  // cout << "IRM log likelihood: " << irm_data.LogLikelihood(parameters) << endl;
  // cout << "IRM log conditional likelihood:\n" << irm_data.LogConditionalLikelihoodVector(parameters) << endl;

  // cout << "parameters:\n" << parameters << endl;
  // vector<TDM> ir=irm.ImpulseResponse(parameters,1);
  // cout << "condition number ir[0]*sigma: " << ConditionNumber(ir[0]*sigma) << endl;

  // //===========================================================================
//}

void GenerateData(void)
{
   unsigned int nvar=4;
   unsigned int nobs=100;
   stringstream filename;
   ofstream output;

   //=== generate data set =====================================================
   filename << "IRM_data_" << nvar << "var.txt";
   output.open(filename.str().c_str());
   if (!output.is_open())
     throw dw_exception("GenerateData() - unable to open " + filename.str());
   output << RandomNormal(nobs,nvar);
   output.close();

   filename.str("");
   filename << "IRM_" << nvar << "var_restricted.txt";
   TTimeSeries_IRM irm=TTimeSeries_IRM_Specification(filename.str());

   // parameters
   TSpecification spec(filename.str());
   vector<string> series=spec.Strings("//== IRM: Series");
   unsigned int n_var=series.size();
   // degree
   unsigned int degree=spec.NonNegativeInteger("//== IRM: Degree");
   // get initial forecast restriction
   TDM f0C;
   vector< vector<bool> > f0R;
   spec.RestrictionMatrix(f0R,f0C,"//== IRM: Initial Forecast Restrictions",' ');
   // get impulse response restrictions
   TDM sC;
   vector< vector<bool> > sR;
   spec.RestrictionMatrix(sR,sC,"//== IRM: Impulse Response Restrictions",' ');
   // create mapping
   I idx, idx_iG;
   unsigned int k=0;
   // initial forecasts
   for (unsigned int i=0; i < n_var; i++)
     for (unsigned int j=0; j < 4; k++, j++)
       if (f0R[i][j])
         {
	   if (j == 2)
	     idx_iG(k);
	   else 
	     idx(k);
	 }
   // impulse responses
   for (unsigned int j=0; j < n_var; j++)
     for (unsigned int i=0; i < n_var; i++)
       for (unsigned int m=0; m < degree+4; k++, m++)
         if (sR[i][j*(degree+4)+m])
    	   {
	      if (m == 2)
	        idx_iG(k);
	      else
		idx(k);
	   }
    
   TDV raw_parameters=spec.Vector("//== IRM: parameters",' ');
   TDV parameters(idx.Size()+idx_iG.Size());
   parameters(I(0,idx.Size()-1),raw_parameters,idx);
   parameters(I(idx.Size(),End),raw_parameters,idx_iG);
  
   // TDV parameters(irm.NumberParameters());
  
   // parameters(I(0,End),raw_parameters(I(0,nvar*4-1)));
   // unsigned int k=nvar*4;
   // for (unsigned int i=0; i < nvar; i++)
   //   for (unsigned int j=0; j < nvar; j++)
   //     parameters(I(k+(i+j*nvar)*(degree+4),End),raw_parameters(I(k+(j+i*nvar)*(degree+4),k+(j+i*nvar+1)*(degree+4)-1)));

   if (parameters.Dim() != irm.NumberParameters())
     {
       cout << "dimension: " << parameters.Dim() << "   expected dimension: " << irm.NumberParameters() << endl;
       throw dw_exception("GenerateData() - invalid number of parameters");
     }

   cout << parameters << endl;
   cout << irm.LikelihoodParameters(parameters) << endl;
   vector<TDV> f0=irm.InitialForecast(parameters,100);
   vector<TDM> s=irm.ImpulseResponse(parameters,100);
   for (unsigned int i=0; i < 2; i++)
     cout << f0[i] << endl << s[i] << endl;
   // write data
   filename.str("");
   filename << "IRM_data_" << nvar << "var.txt";
   output.open(filename.str().c_str());
   if (!output.is_open())
     throw dw_exception("GenerateData() - unable to open " + filename.str());
   output << setprecision(14);
   output << irm.GenerateData(parameters,1,nobs);
   output.close();
   filename.str("");
   filename << "IRM_" << nvar << "var_restricted.txt";
   TTimeSeries_IRM irm2=TTimeSeries_IRM_Specification(filename.str());
   f0=irm.InitialForecast(parameters,irm2.NumberObservations());
   s=irm.ImpulseResponse(parameters,irm2.NumberObservations());
   cout << "IRM log likelihood: " << irm2.LogLikelihood(parameters) << endl;
   cout << "IRM log conditional likelihood:\n" << irm2.LogConditionalLikelihoodVector(parameters) << endl;

   for (unsigned int i=0; i < 100; i++)
     {
	parameters = irm2.DrawPrior();
	cout << "LogPrior: " << irm2.LogPrior(parameters) << endl;
	cout << "LogLikelihood: " << irm2.LogLikelihood(parameters) << endl;
     }
   //===========================================================================
}

int main(void)
{
  try
    {
      dw_initialize_generator(0);

      //TestGamma();

      GeneralVAR();

      //GenerateData_alt();

      //GenerateData();

      // vector<TDV> f0(irm.NumberObservations()+1);
      // vector<TDM> s(irm.NumberObservations());
      // TDV m=MeanCols(irm.Data());
      // for (unsigned int i=0; i < f0.size(); i++) f0[i]=m;
      // s[0]=Identity(irm.NumberVariables());
      // for (unsigned int i=1; i < irm.NumberObservations(); i++)
      // 	s[i]=Zeros(irm.NumberVariables(),irm.NumberVariables());
      // irm.FixedPoint(f0,s,1000);
      // for (unsigned int t=0; t < s.size(); t++) cout << Vec(s[t]);
      

    }
  catch (dw_exception &e)
    {
      std::cout << "Exception: " << e.what() << std::endl;
    }
  return 0;
}
