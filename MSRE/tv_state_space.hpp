#ifndef __DW_TIME_VARYING_STATE_SPACE__
#define __DW_TIME_VARYING_STATE_SPACE__

#include "TTimeSeries.hpp"
#include "TRegimeProcesses.hpp"

namespace TS
{

/********************************************************************************
Linear state space model with time varying parameters. If r(t) denotes the regime
at time t, for 0 <= t < T, then the model is given by a measurment equation
   
    y(t) = a(r(t)) + H(r(t))*z(t) + Psiy(r(t))*u(t)

and a state equation

    z(t) = b(r(t)) + F(r(t))*z(t-1) + Phiz(r(t))*epsilon(t)

where

    r(t)       : takes values between 0 and nr-1 inclusive
    y(t)       : ny vector
    u(t)       : nu vector (nu could be zero)
    z(t)       : nz vector
    epsilon(t) : nepsilon vector (nepsilon >= ny)
    a(i)       : ny vector
    H(i)       : ny x nz matrix
    Phiy(i)    : ny x nu matrix
    b(i)       : nz vector
    F(i)       : nz x nz matrix
    Phiz(i)    : nz x nepsilon matrix

the parameter of the model are a(i), H(i), Phiy(i), b(i), F(i), Phiz(i) and q, 
where q are the parameters controling the regime process. We denote the state
space parameters by theta. The parameters are arranged as 

    a(0)
    vec(H(0))
    vec(Phiy(0))
    b(0)
    vec(F(0))
    vec(Phiz(0))
    .
    .
    .
    a(nr-1)
    vec(H(nr-1))
    vec(Phiy(nr-1))
    b(nr-1)
    vec(F(nr-1))
    vec(Phiz(nr-1))
    q
 
To implement the Kim filter, define

    zeta(t) = (s(t), ... ,s(t-rho+1))  
    xi(t)   = (s(t), ... ,s(t-rho+1),s(t-rho)) = (s(t),zeta(t-1)) = (zeta(t),s(t-rho))

The index rho >= 0 is fixed.  If rho == 0, then zeta(t) is empty.  Both xi and zeta 
are encoded as integers.  Under this encodeding, 

    xi(t) = zeta(t)*h + s(t-rho) = s(t)*(h^rho) + zeta(t-1)

********************************************************************************/
class TLikelihood_TimeVaryingLinearStateSpace : public TLikelihood
{
protected:
  unsigned int nr;
  unsigned int ny;
  unsigned int nu;
  unsigned int nz;
  unsigned int nepsilon;
  unsigned int nparameters;

  // regime process
  TRegimeProcess_markov regime_process;
  unsigned int rho;			// number of lags
  mutable std::vector<DM::TDM> Q;	// transition matrix, Q[t](i,j)=P[xi(t)=i | zeta(t-1)=j]

  // data
  DM::TDM data;

  // Kim filter - not all of these may be used or others may be needed
  mutable std::vector<DM::TDV> Pxi1;  			// Pxi1[t](k)  = P[xi(t) = k | Y(t-1), theta, q]
  mutable std::vector<DM::TDV> Pxi;   			// Pxi[t](k)   = P[xi(t) = k | Y(t), theta, q]
  mutable std::vector< std::vector<DM::TDV> > Pr;	// Pr[t][j](r) = P[s(t-rho) = r | Y(t), zeta(t)=j, theta, q]
  mutable std::vector< std::vector<DM::TDV> > Ez;    	// Ez[t][k]    = E[z(t) | Y(t), xi(t)=k, theta, q]
  mutable std::vector< std::vector<DM::TDM> > Ezz;   	// Ezz[t][k]   = E[(z(t) - Ez(t,k))(z(t) - Ez(t,k))' | Y(t), xi(t)=k, theta, q] 
  mutable std::vector< std::vector<DM::TDV> > IEz;   	// IEz[t][j]   = E[z(t) | Y(t), zeta(t)=j, theta, q]
  mutable std::vector< std::vector<DM::TDM> > IEzz;  	// IEzz[t][j]  = E[(z(t) - Ez(t,zeta(t)))(z(t) - Ez(t,zeta(t)))' | Y(t), zeta(t)=j, theta, q] 
  mutable std::vector< std::vector<DM::TDV> > Ez1;   	// Ez1[t][k]   = E[z(t) | Y(t-1), xi(t)=k, theta, q]
  mutable std::vector< std::vector<DM::TDM> > Ezz1;  	// Ezz1[t][k]  = E[(z(t) - Ez1(t,xi(t)))(z(t) - Ez1(t,xi(t)))' | Y(t-1), xi(t)=k, theta, q]  
  mutable std::vector< std::vector<DM::TDV> > Ey1;   	// Ey1[t][k]   = E[y(t) | Y(t-1), xi(t)=k, theta, q]
  mutable std::vector< std::vector<DM::TDM> > Eyy1;  	// Eyy1[t][k]  = E[(y(t) - Ey1(t,k))(y(t) - Ey1(t,k))' | Y(t-1), xi(t)=k, theta, q] 
  mutable std::vector< std::vector<DM::TDV> > ey1;   	// ey1[t][k]   = y(t) - E[y(t) | Y(t-1), xi(t)=k, theta, q]
  mutable std::vector<DM::TDV> L;     			// L[t](k)     = P[y(t) | Y(t-1), xi(t)=k, theta, q]
  mutable DM::TDV lcli;			     		// lcli(t)     = Sum(L[t][k])


public:
  virtual ~TLikelihood_TimeVaryingLinearStateSpace() {};
  TLikelihood_TimeVaryingLinearStateSpace(const TLikelihood_TimeVaryingLinearStateSpace &TV_StateSpace);
  TLikelihood_TimeVaryingLinearStateSpace(const DM::TDM Data, unsigned int nRegimes, unsigned int nErrors, unsigned int nStates, unsigned int nShocks, unsigned int nLags, const TRegimeProcess_markov &Regime);
  virtual TLikelihood_TimeVaryingLinearStateSpace* Clone(void) const { return new TLikelihood_TimeVaryingLinearStateSpace(*this); };
  
  // Kim filter
  void Filter(const DM::TDV &parameters, unsigned int tau) const;

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
  std::vector<DM::TDV> Get_a(const DM::TDV &parameters) const;
  std::vector<DM::TDM> Get_H(const DM::TDV &parameters) const;
  std::vector<DM::TDM> Get_Phiy(const DM::TDV &parameters) const;
  std::vector<DM::TDV> Get_b(const DM::TDV &parameters) const;
  std::vector<DM::TDM> Get_F(const DM::TDV &parameters) const;
  std::vector<DM::TDM> Get_Phiz(const DM::TDV &parameters) const;
  DM::TDM Get_Ptr(const DM::TDV &parameters) const;
  
  // smoother
  std::vector<DM::TDV> SmoothedProbabilities(unsigned int t0, unsigned int t1) const;
  std::vector< std::vector<DM::TDV> > SmoothedMean(unsigned int t0, unsigned int t1, const DM::TDV &parameters, const std::vector<DM::TDV> &SPxi) const;     

};
  
}

#endif
