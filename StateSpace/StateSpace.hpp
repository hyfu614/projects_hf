
#ifndef __HF_STATESPACE__
#define __HF_STATESPACE__

#include "DM.hpp"
#include "TTimeSeries.hpp"

namespace TS 
{
  //
  // The STATESPACE equations are 
  //
  //     y(t) = a + H*x(t) 
  //	 x(t) = b + F*x(t-1) + Phi*epsilon(t)
  //
  // y(t):	 ny vector, measurement
  // x(t):	 nx vector, states  
  // epsilon(t): nepsilon vector, iid and standard normal 
  // Phi:	 nx * nepsilon matrix, covariance matrix of epsilon(t)
  // a:		 ny vector
  // H:		 ny * nx matrix
  // b:		 nx vector
  // F:		 nx * nx matrix 
  // 
	
  class TLikelihood_StateSpace : public TLikelihood
  {
  protected:
    unsigned int ny, nx, nepsilon, n_obs, n_pre;
    DM::TDV Ex;				// States mean
    DM::TDM VARx;			// States variance
    mutable DM::TDV lcl;		// Log conditional likelihood vector 
    std::vector<DM::TDV> y, y_init;	// Data 
    double ln_constant;

  public:
    // constructors/destructors
    TLikelihood_StateSpace(const TLikelihood_StateSpace &likelihood);
    TLikelihood_StateSpace(const DM::TDM &Y, const DM::TDM &initial_Y, unsigned int n_var, unsigned int n_state, unsigned int n_shock, unsigned int n_init);
    ~TLikelihood_StateSpace() {};
    TLikelihood_StateSpace* Clone(void) const { return new TLikelihood_StateSpace(*this); };

    // likelihood
    unsigned int NumberParameters(void) const { return ny*(1+nx)+nx*(1+nx+nepsilon); };
    double LogConditionalLikelihood(const DM::TDV &parameters, unsigned int t) const; 
    double LogLikelihood(const DM::TDV &parameters) const;
    DM::TDV LogConditionalLikelihoodVector(const DM::TDV &parameters) const;
        
    // Kalman filter
    void KalmanFilter(const DM::TDV &parameters) const;

    // data
    unsigned int NumberObservations(void) const { return n_obs; };
    unsigned int NumberVariables(void) const {  return ny; };
    unsigned int NumberMeasurements(void) const { return ny; };
    unsigned int NumberStates(void) const { return nx; };
    unsigned int NumberShocks(void) const { return nepsilon; };
    DM::TDM Data(void) const;
    
    // info
    DM::TDV Geta(const DM::TDV &parameters) const;
    DM::TDM GetH(const DM::TDV &parameters) const;
    DM::TDV Getb(const DM::TDV &parameters) const;
    DM::TDM GetF(const DM::TDV &parameters) const;
    DM::TDM GetPhi(const DM::TDV &parameters) const;
    DM::TDV GetParameters(const DM::TDV &a, const DM::TDM &H, const DM::TDV &b, const DM::TDM &F, const DM::TDM &Phi) const;
  };
  TLikelihood_StateSpace TLikelihood_StateSpace_Specification(const std::string &filename);

  
  class TNormalPrior : public TDensity
  {
  protected:
    unsigned int dim;
    DM::TDV mean;
    DM::TDV deviation;

  public:
    TNormalPrior(const TNormalPrior &prior) : dim(prior.dim), mean(prior.mean), deviation(prior.deviation) {};
    TNormalPrior(const unsigned int Dim, const DM::TDV &Mean, const DM::TDV &Deviation) : dim(Dim), mean(Mean), deviation(Deviation) {};
    TNormalPrior(const unsigned int Dim) : dim(Dim), mean(DM::Zeros(Dim)), deviation(DM::Ones(Dim)) {};
    virtual ~TNormalPrior() {};

    virtual TNormalPrior* Clone(void) const { return new TNormalPrior(*this); };

    virtual unsigned int Dim(void) const { return dim; };
    virtual double LogDensity(const DM::TDV &x) const;
    virtual DM::TDV Draw(void) const;
  };

  
  class TTimeSeries_StateSpace : public TTimeSeries
  {
  protected:

  public:
    TTimeSeries_StateSpace(const TTimeSeries_StateSpace &Model) : TTimeSeries(Model) {};
    TTimeSeries_StateSpace(const TLikelihood_StateSpace &Likelihood, const TDensity &Prior, const TFunction &Function) : TTimeSeries(Likelihood,Prior,Function) {};
    ~TTimeSeries_StateSpace() {};

    // data
    unsigned int NumberObservations(void) const { return static_cast<TLikelihood_StateSpace&>(likelihood).NumberObservations(); };
    unsigned int NumberMeasurements(void) const { return static_cast<TLikelihood_StateSpace&>(likelihood).NumberMeasurements(); };
    unsigned int NumberStates(void) const { return static_cast<TLikelihood_StateSpace&>(likelihood).NumberStates(); };
    unsigned int NumberShocks(void) const { return static_cast<TLikelihood_StateSpace&>(likelihood).NumberShocks(); };

    // info
    DM::TDV Geta(const DM::TDV &parameters) const { return static_cast<TLikelihood_StateSpace&>(likelihood).Geta(f(parameters)); };
    DM::TDV Getb(const DM::TDV &parameters) const { return static_cast<TLikelihood_StateSpace&>(likelihood).Getb(f(parameters)); };
    DM::TDM GetH(const DM::TDV &parameters) const { return static_cast<TLikelihood_StateSpace&>(likelihood).GetH(f(parameters)); };
    DM::TDM GetF(const DM::TDV &parameters) const { return static_cast<TLikelihood_StateSpace&>(likelihood).GetF(f(parameters)); };
    DM::TDM GetPhi(const DM::TDV &parameters) const { return static_cast<TLikelihood_StateSpace&>(likelihood).GetPhi(f(parameters)); };
    DM::TDV GetLikelihoodParameters(const DM::TDV &a, const DM::TDM &H, const DM::TDV &b, const DM::TDM &F, const DM::TDM &Phi) const { return static_cast<TLikelihood_StateSpace&>(likelihood).GetParameters(a,H,b,F,Phi); };
    
    //// StateSpace specific posterior functions
    //DM::TDV MaximizePosterior(double tolerance=1.0e-6, bool verbose=false);
    //DM::TDV SimulatePosterior(unsigned int N, int type=0, int thin=1, int burn_in=10);
  };
  TTimeSeries_StateSpace TTimeSeries_StateSpace_Specification(const std::string &filename);
}

#endif
