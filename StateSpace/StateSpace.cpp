#include "specification_io.hpp"
#include "StateSpace.hpp"
#include "TFunction.hpp"
#include "TData.hpp"
#include "DM.hpp"
#include "dw_rand.h"
#include "dw_exception.hpp"

using namespace std;
using namespace DM;


namespace TS
{
  //===============================================================================
  //== TLikelihood_StateSpace
  //===============================================================================
  TLikelihood_StateSpace::TLikelihood_StateSpace(const TLikelihood_StateSpace &likelihood)
    : ny(likelihood.ny),
      nx(likelihood.nx),
      nepsilon(likelihood.nepsilon),
      n_obs(likelihood.n_obs),
      n_pre(likelihood.n_pre),
      Ex(likelihood.Ex),
      VARx(likelihood.VARx),
      lcl(likelihood.lcl),
      y(likelihood.y),
      y_init(likelihood.y_init),
      ln_constant(likelihood.ln_constant)
  { }

  TLikelihood_StateSpace::TLikelihood_StateSpace(const TDM &Y, const TDM &initial_Y, unsigned int n_var, unsigned n_state, unsigned int n_shock, unsigned int n_init)
    : ny(n_var),
      nx(n_state),
      nepsilon(n_shock),
      n_obs(Y.Rows()),
      n_pre(n_init),
      Ex(nx,0.0),
      VARx(Identity(nx)),
      lcl(n_obs),
      y(n_obs),
      y_init(n_pre),
      ln_constant(0.918938533204673*ny)  // 0.918938533204673 = 0.5*log(2*pi)
  {
    for (unsigned int t=0; t < n_obs; t++)
      {
	y[t]=Y(t,I(0,End));
      }
    if (n_pre > 0)
      {
        for (unsigned int t=0; t < n_pre; t++)
          {
	    y_init[t]=initial_Y(t,I(0,End));
          }
      }  
  }

  void TLikelihood_StateSpace::KalmanFilter(const TDV &parameters) const
  {
    TDV a=Geta(parameters), b=Getb(parameters);
    TDM H=GetH(parameters), F=GetF(parameters), Phi=GetPhi(parameters);
   
    TDV v, Ey, Ext, Extc=Ex;
    TDM K, N, Temp, VARxt, VARxtc=VARx;

    // state initialization
    if (n_pre > 0)
      {
        for (unsigned int t=0; t<n_pre; t++)
    	{
      	  // Kalman Filter updating
      	  K = VARxtc * T(H);
      	  N = H * K;
      	  Ey = a + H*Extc;
      	  v = y_init[t] - Ey;
      	  Temp = MultiplyInverse(K, N);
      	  Ext = Extc + Temp*v;
      	  VARxt = VARxtc - Temp*T(K); 
      	  Extc = b + F*Ext;
      	  VARxtc = (F*VARxt)*T(F) + Phi*T(Phi);
	}
      } 
    
    for (unsigned int t=0; t<n_obs; t++)
    {
      // Kalman Filter updating
      K = VARxtc * T(H);
      N = H * K;
      Ey = a + H*Extc;
      v = y[t] - Ey;
      Temp = MultiplyInverse(K, N);
      Ext = Extc + Temp*v;
      VARxt = VARxtc - Temp*T(K); 
      Extc = b + F*Ext;
      VARxtc = (F*VARxt)*T(F) + Phi*T(Phi); 
	
      // Log conditional likelihood 
      lcl[t] = -ln_constant - 0.5*(LogAbsDet(N) + InnerProduct(v, InverseMultiply(N,v)));
    }
  }


  double TLikelihood_StateSpace::LogConditionalLikelihood(const DM::TDV &parameters, unsigned int t) const
  {
    if (t >= n_obs) throw dw_exception("TLikelihood_StateSpace::LogConditionalLikelihood() - time index out of range");
    
    KalmanFilter(parameters);  
    
    return lcl[t];
  }

  double TLikelihood_StateSpace::LogLikelihood(const DM::TDV &parameters) const
  {
    double ll=0.0;

    KalmanFilter(parameters);

    for (unsigned int t=0; t<n_obs; t++) ll += lcl[t];
    return ll;
  }

  TDV TLikelihood_StateSpace::LogConditionalLikelihoodVector(const DM::TDV &parameters) const
  {
    KalmanFilter(parameters);
    return lcl;
  }



  TDM TLikelihood_StateSpace::Data(void) const
  {
    TDM Y(n_obs,ny,false);
    for (unsigned int t=0; t < n_obs; t++)
      Y(t,I(0,End),y[t]);
    return Y;
  }

  // extracts measurement equation parameters from parameter vector
  TDV TLikelihood_StateSpace::Geta(const TDV &parameters) const
  {
    TDV a(ny); 
    
    a(I(0,ny-1),parameters,I(0,ny-1));
    
    return a;
  }

  TDM TLikelihood_StateSpace::GetH(const TDV &parameters) const
  { 
    TDM H(ny,nx,false);

    for (unsigned int i=0; i < ny; i++)
    	H(i,I(0,nx-1),parameters,I(ny+i*nx,ny+(i+1)*nx-1));
    
    return H;
  }

  // extracts state transition equation parameters from parameter vector
  TDV TLikelihood_StateSpace::Getb(const TDV &parameters) const
  { 
    TDV b(nx);
    
    b(I(0,nx-1),parameters,I(ny*(1+nx),ny*(1+nx)+nx-1));
    
    return b;
  }

  TDM TLikelihood_StateSpace::GetF(const TDV &parameters) const
  {
    TDM F(nx,nx,false);

    for (unsigned int i=0; i < nx; i++)
    	F(i,I(0,nx-1),parameters,I(ny*(1+nx)+(1+i)*nx,ny*(1+nx)+(2+i)*nx-1));
    
    return F;   
  }

  TDM TLikelihood_StateSpace::GetPhi(const TDV &parameters) const
  {
    TDM Phi(nx,nepsilon,false);
  
    for (unsigned int i=0; i < nx; i++)
    	Phi(i,I(0,nepsilon-1),parameters,I((ny+nx)*(1+nx)+i*nepsilon,(ny+nx)*(1+nx)+(i+1)*nepsilon-1));
    
    return Phi;  
  }
  
  // forms parameter vector from all parameter matrices
  TDV TLikelihood_StateSpace::GetParameters(const TDV &a, const TDM &H, const TDV &b, const TDM &F, const TDM &Phi) const
  {
    if ((unsigned int)a.Dim() != ny)
      throw dw_exception("TLikelihood_StateSpace::GetParameters() - invalid dimensions of a");	
    if ((unsigned int)b.Dim() != nx)
      throw dw_exception("TLikelihood_StateSpace::GetParameters() - invalid dimensions of b");	
    if (((unsigned int)H.Rows() != ny) || ((unsigned int)H.Cols() != nx))
      throw dw_exception("TLikelihood_StateSpace::GetParameters() - invalid dimensions of H");
    if (((unsigned int)F.Rows() != nx) || ((unsigned int)F.Cols() != nx))
      throw dw_exception("TLikelihood_StateSpace::GetParameters() - invalid dimensions of F");
    if (((unsigned int)Phi.Rows() != nx) || ((unsigned int)Phi.Cols() != nepsilon))
      throw dw_exception("TLikelihood_StateSpace::GetParameters() - invalid dimensions of Phi");
    
    TDV parameters(NumberParameters());  
	
    parameters(I(0,ny-1),a);
    parameters(I(ny,ny*(1+nx)-1),Vec(T(H)));
    parameters(I(ny*(1+nx),ny*(1+nx)+nx-1),b);
    parameters(I(ny*(1+nx)+nx,(ny+nx)*(1+nx)-1),Vec(T(F)));
    parameters(I((ny+nx)*(1+nx),End),Vec(T(Phi)));

    return parameters;
  }

  // extracts information from specification file
  TLikelihood_StateSpace TLikelihood_StateSpace_Specification(const std::string &filename)
  {
    vector<string> specs;
    vector< vector<string> > values;
    
    //get specification file
    ParseSpecificationFile(filename,specs,values);

    // number of variables
    int n_var=ParsePositiveInteger(specs,values,"//== StateSpace: Number Variables");

    // number of states
    int n_state=ParsePositiveInteger(specs,values,"//== StateSpace: Number States");

    // number of shocks
    int n_shock=ParseNonNegativeInteger(specs,values,"//== StateSpace: Number Shocks");

    // number of state initialization
    unsigned int n_init=ParseNonNegativeInteger(specs,values,"//== StateSpace: Number State Initialization");
	
    // series
    vector<string> series=ParseStrings(specs,values,"//== StateSpace: Series");

    // begin and end dates
    string begin_date=ParseString(specs,values,"//== StateSpace: Begin Date");
    string end_date=ParseString(specs,values,"//== StateSpace: End Date");

    // get raw data
    TTimeSeriesDataRaw raw_data=TTimeSeriesDataRaw_Specification(filename);

    // Setup Y
    TDM Y=raw_data.Data(series,begin_date,end_date);
    
    // setup initial data for state initialization
    unsigned int idx=raw_data.DateIndex(begin_date);
    if (idx < n_init)
      n_init=idx;
    TDM initial_Y;
    if (n_init > 0) 
      initial_Y=raw_data.Data(series,idx-n_init,idx-1);
     
    return TLikelihood_StateSpace(Y,initial_Y,n_var,n_state,n_shock,n_init);    
  }

  //===============================================================================
  //== TNormalPrior 
  //===============================================================================
  double TNormalPrior::LogDensity(const DM::TDV &x) const
  {
    double log_prior=-0.918938533204673*dim;      //0.918938533204673 = 0.5*log(2*pi)
    
    log_prior -= Sum(Log(deviation))+0.5*InnerProduct(eMultiplyInverse(x-mean, deviation));
    
    return log_prior;
    
  }
  
  TDV TNormalPrior::Draw(void) const
  {
    TDV parameters(dim);
      
    return parameters = eMultiply(deviation, RandomNormal(dim))+mean;
  }
      
  //===============================================================================
  //== TTimeSeries_StateSpace 
  //===============================================================================
  TTimeSeries_StateSpace TTimeSeries_StateSpace_Specification(const string &filename)
  {
    TLikelihood_StateSpace likelihood=TLikelihood_StateSpace_Specification(filename);
    unsigned int ny=likelihood.NumberMeasurements(), nx=likelihood.NumberStates(), nepsilon=likelihood.NumberShocks();

    vector<string> specs;
    vector< vector<string> > values;
 
    //get specification file
    ParseSpecificationFile(filename,specs,values);

    // get measurement equation constant restriction vector
    TDM AC;
    vector< vector<bool> > ar;
    ParseRestrictionMatrix(specs,values,"//== StateSpace: Restrictions a",' ',ar,AC);
    TDV ac=Vec(AC);

    // get measurement model restriction matrix
    TDM HC;
    vector< vector<bool> > HR;
    ParseRestrictionMatrix(specs,values,"//== StateSpace: Restrictions H",' ',HR,HC);

    // get state transition equation constant restriction vector
    TDM BC;
    vector< vector<bool> > br;
    ParseRestrictionMatrix(specs,values,"//== StateSpace: Restrictions b",' ',br,BC);
    TDV bc=Vec(BC);

    // get state transition restriction matrix
    TDM FC;
    vector< vector<bool> > FR;
    ParseRestrictionMatrix(specs,values,"//== StateSpace: Restrictions F",' ',FR,FC);

    // get state transition shock variance restriction matrix
    TDM PhiC;
    vector< vector<bool> > PhiR;
    ParseRestrictionMatrix(specs,values,"//== StateSpace: Restrictions Phi",' ',PhiR,PhiC);

    // create parameter mapping
    TDV c(likelihood.GetParameters(ac,HC,bc,FC,PhiC));
    I idx;
    
    for (unsigned int i=0; i < ny; i++)
      if (ar[i][0]) idx(i);
    for (unsigned int i=0; i < ny; i++)
      {
	for (unsigned int j=0; j < nx; j++)
	  if (HR[i][j]) idx(ny+i*nx+j);
      }
    for (unsigned int i=0; i < nx; i++)
      if (br[i][0]) idx(ny*(1+nx)+i);
    for (unsigned int i=0; i < nx; i++)
      {
        for (unsigned int j=0; j < nx; j++)
	  if (FR[i][j]) idx(ny*(1+nx)+nx*(1+i)+j);
      }
    for (unsigned int i=0; i < nx; i++)
      {
	for (unsigned int j=0; j < nepsilon; j++)  
	  if (PhiR[i][j]) idx((ny+nx)*(1+nx)+i*nepsilon+j);  
      }
    TLinear_inclusion f(c,idx);


    // prior
    string prior_type=ParseString(specs,values,"//== StateSpace: Prior Type");
    if (prior_type.compare("Flat") == 0)
      {
	return TTimeSeries_StateSpace(likelihood,TFlatPrior(f.Domain()),f);
      }
    else if (prior_type.compare("Normal") == 0)
      {
        return TTimeSeries_StateSpace(likelihood,TNormalPrior(f.Domain()),f);
      }
    else
      throw dw_exception("TTimeSeries_StateSpace_Specification() - unknown prior type (" + prior_type + ")");
  }

}





 
