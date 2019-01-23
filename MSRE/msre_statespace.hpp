#ifndef __MSRE_STATESPACE__
#define __MSRE_STATESPACE__

#include "TTimeSeries.hpp"
#include "TRegimeProcesses.hpp"
#include "tv_state_space.hpp"
#include "msv_msre.hpp"

namespace TS
{

/********************************************************************************
Markov-switching rational expectations (MSRE) model with minimal state variable (MSV)
solutions. If r(t) denotes the regime at time t, for 0 <= t < T, then the MSRE model
has the form:

    A(r(t)) x(t) = B(r(t)) x(t-1) + Psi(r(t)) epsilon(t) + Pi eta(t)

where
    
    r(t)	: takes values between 0 and nr-1 inclusive
    x(t)	: n vector of endogenous and predetermined variables
    epsilon(t)  : nepsilon vector of iid stationary exogenous shocks
    eta(t)	: s vector of expectational errors
    A(i)	: n x n matrix
    B(i)	: n x n matrix
    Psi(i)	: n x nepsilon matrix
    Pi		: n x s matrix  

The parameters of the model are A(i), B(i), Psi(i), and P, where P are the parameters 
controling the regime process (transition probabilities). The parameters are arranged as 

    vec(A(0))
    vec(B(0))
    vec(Psi(0))
    .
    .
    .
    vec(A(nr-1))
    vec(B(nr-1))
    vec(Psi(nr-1))
    vec(P)
 

Assumes that Pi' = [Zeros(s,n-s)  Identity(s)] and that A(i) is invertible.  
P is the regime transition matrix and P(i,j) is the probability that r(t)=j 
given that r(t-1)=i. Note that the rows of P must sum to one. 

The MSV solutions is of the form

    x(t) = V{r(t)}*F1{r(t)}*x(t-1) + V{r(t)}*G1{r(t)}*epsilon(t)
    eta(t) = F2{r(t)}*x(t-1) + G2{r(t)}*epsilont(t)
 
The measurement equation of the linear state space model is given by

    y(t) = a(r(t)) + H(r(t))*x(t) + Phiy(r(t))*u(t)

The parameters of the measurement equation are arranged as

    a(0)
    vec(H(0))
    vec(Phiy(0))
    .
    .
    .
    a(nr-1)
    vec(H(nr-1))
    vec(Phiy(nr-1))

**************************************************************************************/
  class TLikelihood_MSRE : public TLikelihood
  {
  protected: 
    unsigned int nr;			// number of regimes
    unsigned int n;			// number of states
    unsigned int s;			// number of expectational errors
    unsigned int nepsilon;		// number of shocks
    unsigned int ny;			// number of variables in measurement equation
    unsigned int nerrors;		// number of errors in measurement equation 	
    unsigned int nparameters;  		// number of MSRE parameters

    // regime process
    unsigned int rho;			// number of lags for kim filter
    mutable DM::TDM Ptr;		// transition matrix for regime switching
    TRegimeProcess_markov regime_process;	

    // data
    DM::TDM data;

    // log-likelihood vector 
    mutable DM::TDV lcli;	

    // tv state space
    TLikelihood_TimeVaryingLinearStateSpace& tv_statespace;

  public:
    virtual ~TLikelihood_MSRE() { delete &tv_statespace; };
    TLikelihood_MSRE(const TLikelihood_MSRE &Msre_SS);
    TLikelihood_MSRE(const DM::TDM Data, unsigned int nRegimes, unsigned int nErrors, unsigned int nStates, unsigned int nShocks, unsigned int nMeaErrors, unsigned int nLags);
    virtual TLikelihood_MSRE* Clone(void) const { return new TLikelihood_MSRE(*this); };
  
    // likelihood
    virtual unsigned int NumberParameters(void) const { return nparameters; };
    virtual double LogConditionalLikelihood(const DM::TDV &parameters, unsigned int t) const;
    virtual double LogLikelihood(const DM::TDV &parameters) const;
    virtual DM::TDV LogConditionalLikelihoodVector(const DM::TDV &parameters) const;

    // data: NumberObservations() x NumberVariables() matrix of data used to compute log-likelihood
    virtual unsigned int NumberObservations(void) const { return data.Rows(); };
    virtual unsigned int NumberVariables(void) const { return data.Cols(); };
    virtual DM::TDM Data(void) const { return data; };   
 
    // state space form parameters
    DM::TDV Get_StateSpace_parameters(const DM::TDV &parameters) const;		// parameters=[msre_parameters' measurement_parameters']'

  };

  class TTimeSeries_MSRE : public TTimeSeries
  {
  public:
    TTimeSeries_MSRE(const TTimeSeries_MSRE &TimeSeries) : TTimeSeries(TimeSeries) {};
    TTimeSeries_MSRE(const TLikelihood_MSRE &Likelihood, const TDensity &Prior, const TFunction &Function) : TTimeSeries(Likelihood,Prior,Function) {};
    virtual ~TTimeSeries_MSRE() { };
    virtual TTimeSeries_MSRE* Clone(void) const { return new TTimeSeries_MSRE(*this); };

  };
  
}

/************************************************************************************/

#endif
