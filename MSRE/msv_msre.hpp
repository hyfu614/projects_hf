#ifndef __MSV_MSRE__
#define __MSV_MSRE__

#include "TTimeSeries.hpp"
#include "TRegimeProcesses.hpp"

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

In most applications, x(t) is partitioned as 
    x(t)' = [y(t)' z(t)' Et[y(t+1)']
where the first pair [y(t)' z(t)'] is of dimension n-s and
    y(t) = E(t-1)[y(t)] + eta(t)

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

    x(t) = V{s(t)}*F1{s(t)}*x(t-1) + V{s(t)}*G1{s(t)}*epsilon(t)
    eta(t) = F2{s(t)}*x(t-1) + G2{s(t)}*epsilont(t)
 

********************************************************************************/
class TMSV_MSRE
{
protected:
  unsigned int nr;		// number of regimes
  unsigned int n;		// number of states
  unsigned int s;		// number of expectational errors
  unsigned int nepsilon;	// number of shocks
  unsigned int n_parameters;   	// number of MSRE parameters
  mutable std::vector<DM::TDM> V;
  mutable std::vector<DM::TDM> F1;
  mutable std::vector<DM::TDM> F2;
  mutable std::vector<DM::TDM> G1; 
  mutable std::vector<DM::TDM> G2;

  // regime process
  TRegimeProcess_markov regime_process;
  mutable DM::TDM Ptr;			// transition matrix

public:
  ~TMSV_MSRE() {};
  TMSV_MSRE(const TMSV_MSRE &Msre);
  TMSV_MSRE(unsigned int nRegimes, unsigned int nExpErrors, unsigned int nStates, unsigned int nShocks);
  TMSV_MSRE* Clone(void) const { return new TMSV_MSRE(*this); };
  
  // msv solutions
  // max_count: the maximum number of iterations on Newton's method before failing 
  // tol: the convergence criterion
  // A positive return value is the number of iterations needed to obtain convergence 
  // and indicates success. A negitive return value is the number of iterations 
  // before the method terminated without convergence and indicates failure.
  int MSVsolution(const DM::TDV &parameters, const DM::TDV &initial_values, unsigned int max_count=5000, double tol=1.0e-6) const;

  // MSRE form parameters
  unsigned int NumberParameters(void) const { return n_parameters; };
  std::vector<DM::TDM> Get_A(const DM::TDV &parameters) const;
  std::vector<DM::TDM> Get_B(const DM::TDV &parameters) const;
  std::vector<DM::TDM> Get_Psi(const DM::TDV &parameters) const;
  DM::TDM Get_Ptr(const DM::TDV &parameters) const;

  // get solutions
  std::vector<DM::TDM> Get_V(void) const { return V; };
  std::vector<DM::TDM> Get_F1(void) const { return F1; };
  std::vector<DM::TDM> Get_F2(void) const { return F2; };
  std::vector<DM::TDM> Get_G1(void) const { return G1; };
  std::vector<DM::TDM> Get_G2(void) const { return G2; };

  // starting values
  DM::TDM StartingValules(const std::vector<DM::TDM> &A, const std::vector<DM::TDM> &B, bool all=true) const;
  DM::TDM AllInvariantSubspaces(const DM::TDM &A0, const DM::TDM &B0, unsigned int &n_explosive) const;
  
};

/************************************************************************************/
DM::TDV SortIdx(const DM::TDV &v, DM::TDV &sorted_id, bool ascending=true);
unsigned int nchoosek(unsigned int n, unsigned int k);
unsigned int factorial(unsigned int n);
  
}

#endif
