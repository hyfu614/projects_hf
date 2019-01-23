#include <fstream>
#include <sstream>
#include <math.h>

#include "dw_rand.h"
#include "dw_math.h"
#include "TRegimeProcesses.hpp"
#include "specification_io.hpp"
#include "tv_state_space.hpp"

using namespace std;
using namespace DM;

namespace TS
{
  TLikelihood_TimeVaryingLinearStateSpace::TLikelihood_TimeVaryingLinearStateSpace(const TLikelihood_TimeVaryingLinearStateSpace &TV_StateSpace)
    : nr(TV_StateSpace.nr),
      ny(TV_StateSpace.ny),
      nu(TV_StateSpace.nu),
      nz(TV_StateSpace.nz),
      nepsilon(TV_StateSpace.nepsilon),
      nparameters(TV_StateSpace.nparameters),
      regime_process(TV_StateSpace.regime_process),
      rho(TV_StateSpace.rho),
      Q(rho+1),
      data(TV_StateSpace.data),
      Pxi1(data.Rows()),
      Pxi(data.Rows()),
      Pr(data.Rows()),
      Ez(data.Rows()),
      Ezz(data.Rows()),
      IEz(data.Rows()),
      IEzz(data.Rows()),
      Ez1(data.Rows()),
      Ezz1(data.Rows()),
      Ey1(data.Rows()),
      Eyy1(data.Rows()),
      ey1(data.Rows()),
      L(data.Rows()),
      lcli(data.Rows())
  { }

  TLikelihood_TimeVaryingLinearStateSpace::TLikelihood_TimeVaryingLinearStateSpace(const TDM Data, unsigned int nRegimes, unsigned int nErrors, unsigned int nStates, unsigned int nShocks, unsigned int nLags, const TRegimeProcess_markov &Regime)
    : data(Data),
      nr(nRegimes),
      ny(Data.Cols()),
      nu(nErrors),
      nz(nStates),
      nepsilon(nShocks),
      nparameters(nr*(ny*(1+nz+nu)+nz*(1+nz+nepsilon)+nr)),
      regime_process(Regime),
      rho(nLags), 
      Q(rho+1),
      Pxi1(Data.Rows()),
      Pxi(Data.Rows()),
      Pr(Data.Rows()),
      Ez(Data.Rows()),
      Ezz(Data.Rows()),
      IEz(Data.Rows()),
      IEzz(Data.Rows()),
      Ez1(Data.Rows()),
      Ezz1(Data.Rows()),
      Ey1(Data.Rows()),
      Eyy1(Data.Rows()),
      ey1(Data.Rows()),
      L(Data.Rows()),
      lcli(Data.Rows())
  { }

  void TLikelihood_TimeVaryingLinearStateSpace::Filter(const TDV &parameters, unsigned int tau) const
  {
    if (tau < rho || tau >= NumberObservations()) 
      throw dw_exception("TLikelihood_TimeVaryingLinearStateSpace::Filter() - invalid time tau");

    if (parameters.Dim() != nparameters)
      throw dw_exception("TLikelihood_TimeVaryingLinearStateSpace::Filter() - invalid number of likelihood parameters");

    double ln_constant=0.918938533204673*ny; 	 // 0.918938533204673 = 0.5*log(2*pi)

    vector<TDV> a=Get_a(parameters), b=Get_b(parameters); 
    vector<TDM> H=Get_H(parameters), Phiy=Get_Phiy(parameters), F=Get_F(parameters), Phiz=Get_Phiz(parameters);
    TDM Ptr=Get_Ptr(parameters);

    lcli.UniqueMemory(tau+1);

    // state initialization 
    TDV Ez0=Zeros(nz);
    TDM Ezz0=Identity(nz);
    TDV Pxi0=regime_process.InitialProbabilities(Vec(Ptr));
    unsigned int nxi=nr, k;
    TDM K, Temp;

    for (unsigned int t=0; t<rho; t++)
      {
        Pxi1[t].UniqueMemory(nxi*nr); 	
	Pxi[t].UniqueMemory(nxi*nr); 
	L[t].UniqueMemory(nxi*nr);
	Q[t]=Zeros(nxi*nr,nxi);
	Ez1[t].resize(nxi*nr);
	Ezz1[t].resize(nxi*nr);
	Ey1[t].resize(nxi*nr);
	Eyy1[t].resize(nxi*nr);
	ey1[t].resize(nxi*nr);
	Ez[t].resize(nxi*nr);
	Ezz[t].resize(nxi*nr);

	for (unsigned int i=0; i<nr; i++)	// current time t state index
    	  {
	    for (unsigned int j=0; j<nxi; j++)	// time t-1 state index
 	      {
		k = i*nxi+j;			// state index
		
		// Kalman filter
		// Ez1[t][k] = E[z(t) | Y(t-1), xi(t)=k, theta, q]
	        Ez1[t][k] = (t == 0) ? b[i]+F[i]*Ez0 : b[i]+F[i]*Ez[t-1][j];
		// Ezz1[t][k] = E[(z(t) - Ez1(t,xi(t)))(z(t) - Ez1(t,xi(t)))' | Y(t-1), xi(t)=k, theta, q] 
		Ezz1[t][k] = (t == 0) ? (F[i]*Ezz0)*T(F[i])+Phiz[i]*T(Phiz[i]) : (F[i]*Ezz[t-1][j])*T(F[i])+Phiz[i]*T(Phiz[i]);
	 	Ezz1[t][k] = 0.5*(Ezz1[t][k]+T(Ezz1[t][k]));	// force symmetric
		// Ey1[t][k] = E[y(t) | Y(t-1), xi(t)=k, theta, q]
		Ey1[t][k] = a[i]+H[i]*Ez1[t][k];
		K = Ezz1[t][k]*T(H[i]);
		// Eyy1[t][k] = E[(y(t) - Ey1(t,k))(y(t) - Ey1(t,k))' | Y(t-1), xi(t)=k, theta, q]   	  	
		Eyy1[t][k] = H[i]*K+Phiy[i]*T(Phiy[i]);
		Eyy1[t][k] = 0.5*(Eyy1[t][k]+T(Eyy1[t][k]));
		// ey1[t][k] = y(t) - E[y(t) | Y(t-1), xi(t)=k, theta, q]
      	  	ey1[t][k] = data(t,I(0,End)) - Ey1[t][k];
      	  	Temp = MultiplyInverse(K, Eyy1[t][k]);
		// Ez[t][k] = E[z(t) | Y(t), xi(t)=k, theta, q]
      	  	Ez[t][k] = Ez1[t][k] + Temp*ey1[t][k];
		// Ezz[t][k] = E[(z(t) - Ez(t,k))(z(t) - Ez(t,k))' | Y(t), xi(t)=k, theta, q] 
      	  	Ezz[t][k] = Ezz1[t][k] - Temp*T(K); 
  		Ezz[t][k] = 0.5*(Ezz[t][k]+T(Ezz[t][k]));
		// Log conditional likelihood L[t][k] = P[y(t) | Y(t-1), xi(t)=k, theta, q]
      		L[t](k) = -ln_constant - 0.5*(LogAbsDet(Eyy1[t][k]) + InnerProduct(ey1[t][k], InverseMultiply(Eyy1[t][k],ey1[t][k])));

	      }
	    
	    Q[t](I(i*nxi,(i+1)*nxi-1),I(0,End),Kron(Diag(Ptr(I(0,End),i)),Identity(nxi/nr)));
	    
	  }

	// number of states at time t
	nxi *= nr;

	// Pxi1[t](k) = P[xi(t) = k | Y(t-1), theta, q]
	Pxi1[t] = (t == 0) ? Q[t]*Pxi0 : Q[t]*Pxi[t-1];
        
	lcli(t)=MINUS_INFINITY;
	// Pxi[t](k) = P[xi(t) = k | Y(t), theta, q]        
 	for (k=0; k<nxi; k++)
	  {
	    if (Pxi1[t](k) > 0.0)
	      {
	        lcli(t)=AddScaledLogs(1.0,lcli(t),Pxi1[t](k),L[t](k));
	        Pxi[t](k)=log(Pxi1[t](k)) + L[t](k);
	      }
	    else
	      Pxi[t](k)=MINUS_INFINITY;
	  }
	
        for (k=0; k<nxi; k++)
	  if (Pxi[t](k) != MINUS_INFINITY)
	    Pxi[t](k)=exp(Pxi[t](k) - lcli(t));
	  else
	    Pxi[t](k)=0.0;
	
      } 

  vector<TDV> Pzeta(tau+1);
  unsigned int ns=pow(nr,rho+1);	// total number of states
  unsigned int nzeta=ns/nr;		// number of states after collapsion
  
  // transition base matrix
  TDM Qbase=Zeros(nzeta/nr,nzeta);
  for (unsigned i=0; i<nzeta/nr; i++)
    Qbase(i,I(i*nr,(i+1)*nr-1),Ones(nr));
  if (regime_process.IsInvariant())
    {
      Q[rho]=Zeros(ns,ns);
      for (unsigned int i=0; i<nr; i++)
        Q[rho](I(i*nzeta,(i+1)*nzeta-1),I(0,End),Kron(Diag(Ptr(I(0,End),i)),Qbase));
    }

  // Kim filter
  for (unsigned int t=rho; t <= tau; t++)
    {
      // Sets IEz[t-1] and IEzz[t-1] 
      // Uses Ez[t-1], Ezz[t-1], and Pxi[t-1] = i | Y[t-1], Theta)
      // IEz[t][j] = E[z(t) | Y(t), zeta(t)=j, theta, q]
      // IEzz[t][j] = E[(z(t) - Ez(t,zeta(t)))(z(t) - Ez(t,zeta(t)))' | Y(t), zeta(t)=j, theta, q] 
      Pzeta[t-1].UniqueMemory(nzeta);
      IEz[t-1].resize(nzeta);
      IEzz[t-1].resize(nzeta);
      Pr[t-1].resize(nzeta);
      L[t].UniqueMemory(ns);
      Pxi1[t].UniqueMemory(ns);
      Pxi[t].UniqueMemory(ns);
      if (!regime_process.IsInvariant()) 
        Q.push_back(Zeros(ns,ns));
      Ez1[t].resize(ns);
      Ezz1[t].resize(ns);
      Ey1[t].resize(ns);
      Eyy1[t].resize(ns);
      ey1[t].resize(ns);
      Ez[t].resize(ns);
      Ezz[t].resize(ns);
      
      for (unsigned int j=0; j<nzeta; j++)
	{
	  Pzeta[t-1](j) = Sum(Pxi[t-1](I(j*nr,(j+1)*nr-1)));	// Pzeta[t](j) = P[zeta(t) = j | Y(t), theta, q] 
	  IEz[t-1][j] = Pxi[t-1](j*nr)*Ez[t-1][j*nr];
	  for (unsigned int i=1; i<nr; i++)
	    IEz[t-1][j] += Pxi[t-1](j*nr+i)*Ez[t-1][j*nr+i];
	  IEz[t-1][j] /= Pzeta[t-1](j);  
          
	  IEzz[t-1][j] = Pxi[t-1](j*nr)*(Ezz[t-1][j*nr]+OuterProduct(IEz[t-1][j]-Ez[t-1][j*nr]));
 	  for (unsigned int i=1; i<nr; i++)
	    IEzz[t-1][j] += Pxi[t-1](j*nr+i)*(Ezz[t-1][j*nr+i]+OuterProduct(IEz[t-1][j]-Ez[t-1][j*nr+i]));
	  IEzz[t-1][j] /= Pzeta[t-1](j);

    	  TDV prb(nr);
	  for (unsigned int i=0; i<nr; i++)
	    prb(i)=Pxi[t-1](j*nr+i) / Pzeta[t-1](j);
	  Pr[t-1][j]=prb;
      	}

      // Kalman filter
      for (unsigned int i=0; i<nr; i++)		// current time t state index
    	{
	  for (unsigned int j=0; j<nzeta; j++)	// time t-1 state index
 	    {
	      k = i*nzeta+j;			// state index
	      // Ez1[t][k] = E[z(t) | Y(t-1), xi(t)=k, theta, q]
	      Ez1[t][k] = b[i]+F[i]*IEz[t-1][j];
	      // Ezz1[t][k] = E[(z(t) - Ez1(t,xi(t)))(z(t) - Ez1(t,xi(t)))' | Y(t-1), xi(t)=k, theta, q] 
	      Ezz1[t][k] = (F[i]*IEzz[t-1][j])*T(F[i])+Phiz[i]*T(Phiz[i]);
	      Ezz1[t][k] = 0.5*(Ezz1[t][k]+T(Ezz1[t][k]));
	      // Ey1[t][k] = E[y(t) | Y(t-1), xi(t)=k, theta, q]
	      Ey1[t][k] = a[i]+H[i]*Ez1[t][k];
	      K = Ezz1[t][k]*T(H[i]);
	      // Eyy1[t][k] = E[(y(t) - Ey1(t,k))(y(t) - Ey1(t,k))' | Y(t-1), xi(t)=k, theta, q]   	  	
	      Eyy1[t][k] = H[i]*K+Phiy[i]*T(Phiy[i]);	
	      Eyy1[t][k] = 0.5*(Eyy1[t][k]+T(Eyy1[t][k]));
	      // ey1[t][k] = y(t) - E[y(t) | Y(t-1), xi(t)=k, theta, q]
      	      ey1[t][k] = data(t,I(0,End)) - Ey1[t][k];
      	      Temp = MultiplyInverse(K, Eyy1[t][k]);
	      // Ez[t][k] = E[z(t) | Y(t), xi(t)=k, theta, q]
      	      Ez[t][k] = Ez1[t][k] + Temp*ey1[t][k];
	      // Ezz[t][k] = E[(z(t) - Ez(t,k))(z(t) - Ez(t,k))' | Y(t), xi(t)=k, theta, q] 
      	      Ezz[t][k] = Ezz1[t][k] - Temp*T(K); 
	      Ezz[t][k] = 0.5*(Ezz[t][k]+T(Ezz[t][k]));
	      // Log conditional likelihood L[t](k) = P[y(t) | Y(t-1), xi(t)=k, theta, q]
      	      L[t](k) = -ln_constant - 0.5*(LogAbsDet(Eyy1[t][k]) + InnerProduct(ey1[t][k], InverseMultiply(Eyy1[t][k],ey1[t][k])));
		
	    }
    
	  if (!regime_process.IsInvariant()) 
	    Q[t](I(i*nzeta,(i+1)*nzeta-1),I(0,End),Kron(Diag(Ptr(I(0,End),i)),Qbase));
	}

      // Pxi1[t](k) = P[xi(t) = k | Y(t-1), theta, q]
      if (regime_process.IsInvariant()) 
	Pxi1[t]=Q[rho] * Pxi[t-1];
      else
        Pxi1[t]=Q[t] * Pxi[t-1];
 	 
      lcli(t)=MINUS_INFINITY;
      // Pxi[t](k) = P[xi(t) = k | Y(t), theta, q]
      for (k=0; k<ns; k++)
	{
	  if (Pxi1[t](k) > 0.0)
	    {
	      lcli(t)=AddScaledLogs(1.0,lcli(t),Pxi1[t](k),L[t](k));
	      Pxi[t](k)=log(Pxi1[t](k)) + L[t](k);
	    }
	  else
	    Pxi[t](k)=MINUS_INFINITY;
	}
  
      for (k=0; k<ns; k++)
	if (Pxi[t](k) != MINUS_INFINITY)
	  Pxi[t](k)=exp(Pxi[t](k) - lcli(t));
	else
	  Pxi[t](k)=0.0;	

    }
  	
  // IEz[tau] and IEzz[tau]
  Pzeta[tau].UniqueMemory(nzeta);
  IEz[tau].resize(nzeta);
  IEzz[tau].resize(nzeta);
  Pr[tau].resize(nzeta);
  for (unsigned int j=0; j<nzeta; j++)
    {
      Pzeta[tau](j) = Sum(Pxi[tau](I(j*nr,(j+1)*nr-1)));	// Pzeta[t](j) = P[zeta(t) = j | Y(t), theta, q]

      IEz[tau][j] = Pxi[tau](j*nr)*Ez[tau][j*nr];
      for (unsigned int i=1; i<nr; i++)
	IEz[tau][j] += Pxi[tau](j*nr+i)*Ez[tau][j*nr+i];
      IEz[tau][j] /= Pzeta[tau](j);  
          
      IEzz[tau][j] = Pxi[tau](j*nr)*(Ezz[tau][j*nr]+OuterProduct(IEz[tau][j]-Ez[tau][j*nr]));
      for (unsigned int i=1; i<nr; i++)
	IEzz[tau][j] += Pxi[tau](j*nr+i)*(Ezz[tau][j*nr+i]+OuterProduct(IEz[tau][j]-Ez[tau][j*nr+i]));
      IEzz[tau][j] /= Pzeta[tau](j);
      
      TDV prb(nr);
      for (unsigned int i=0; i<nr; i++)
	prb(i)=Pxi[tau](j*nr+i) / Pzeta[tau](j);
      Pr[tau][j]=prb;
    }
	
  }

  vector<TDV> TLikelihood_TimeVaryingLinearStateSpace::Get_a(const TDV &parameters) const
  {
    vector<TDV> a(nr);
    unsigned int stride=ny*(1+nz+nu)+nz*(1+nz+nepsilon);
    for(unsigned int i=0, k=0; i<nr; i++, k+=stride)
      a[i]=parameters(I(k,k+ny-1));

    return a;
  }

  vector<TDM> TLikelihood_TimeVaryingLinearStateSpace::Get_H(const TDV &parameters) const
  {
    vector<TDM> H(nr);
    unsigned int stride=ny*(1+nz+nu)+nz*(1+nz+nepsilon);
    for(unsigned int i=0, k=ny; i<nr; i++, k+=stride)
      H[i]=Reshape(parameters(I(k,k+ny*nz-1)),ny,nz);

    return H;
  }

  vector<TDM> TLikelihood_TimeVaryingLinearStateSpace::Get_Phiy(const TDV &parameters) const 
  {
    vector<TDM> Phiy(nr);
    unsigned int stride=ny*(1+nz+nu)+nz*(1+nz+nepsilon);
    for(unsigned int i=0, k=ny*(1+nz); i<nr; i++, k+=stride)
      Phiy[i]=Reshape(parameters(I(k,k+ny*nu-1)),ny,nu);

    return Phiy;
  }

  vector<TDV> TLikelihood_TimeVaryingLinearStateSpace::Get_b(const TDV &parameters) const
  {
    vector<TDV> b(nr);
    unsigned int stride=ny*(1+nz+nu)+nz*(1+nz+nepsilon);
    for(unsigned int i=0, k=ny*(1+nz+nu); i<nr; i++, k+=stride)
      b[i]=parameters(I(k,k+nz-1));

    return b;
  }

  vector<TDM> TLikelihood_TimeVaryingLinearStateSpace::Get_F(const TDV &parameters) const
  {
    vector<TDM> F(nr);
    unsigned int stride=ny*(1+nz+nu)+nz*(1+nz+nepsilon);
    for(unsigned int i=0, k=ny*(1+nz+nu)+nz; i<nr; i++, k+=stride)
      F[i]=Reshape(parameters(I(k,k+nz*nz-1)),nz,nz);

    return F;
  }

  vector<TDM> TLikelihood_TimeVaryingLinearStateSpace::Get_Phiz(const TDV &parameters) const
  {
    vector<TDM> Phiz(nr);
    unsigned int stride=ny*(1+nz+nu)+nz*(1+nz+nepsilon);
    for(unsigned int i=0, k=ny*(1+nz+nu)+nz*(1+nz); i<nr; i++, k+=stride)
      Phiz[i]=Reshape(parameters(I(k,k+nz*nepsilon-1)),nz,nepsilon);

    return Phiz;
  }

  TDM TLikelihood_TimeVaryingLinearStateSpace::Get_Ptr(const TDV &parameters) const
  {
    TDM Ptr(nr,nr);
    unsigned int offset=nr*(ny*(1+nz+nu)+nz*(1+nz+nepsilon));
    Ptr=Reshape(parameters(I(offset,offset+nr*nr-1)),nr,nr);

    return Ptr;
  }

  double TLikelihood_TimeVaryingLinearStateSpace::LogConditionalLikelihood(const TDV &parameters, unsigned int t) const
  {
    if (t < 0 || t > NumberObservations()) 
      throw dw_exception("TLikelihood_TimeVaryingLinearStateSpace::LogConditionalLikelihood() - invalid time t");

    if (parameters.Dim() != nparameters)
      throw dw_exception("TLikelihood_TimeVaryingLinearStateSpace::LogConditionalLikelihood() - invalid number of likelihood parameters");

    Filter(parameters,t);

    return lcli(t);
  }

  TDV TLikelihood_TimeVaryingLinearStateSpace::LogConditionalLikelihoodVector(const TDV &parameters) const
  {
    if (parameters.Dim() != nparameters)
      throw dw_exception("TLikelihood_TimeVaryingLinearStateSpace::LogConditionalLikelihoodVector() - invalid number of likelihood parameters");

    Filter(parameters,NumberObservations()-1);
 
    return lcli;
  }

  double TLikelihood_TimeVaryingLinearStateSpace::LogLikelihood(const TDV &parameters) const
  {
    return Sum(LogConditionalLikelihoodVector(parameters));
  }

  // Smoothed Probabilities via backward recursion 
  vector<TDV> TLikelihood_TimeVaryingLinearStateSpace::SmoothedProbabilities(unsigned int t0, unsigned int t1) const
  {
    unsigned int ns=Pxi[t1].Dim();
    vector<TDV> SPxi(t1+1);
    SPxi[t1]=Pxi[t1];  
    for (unsigned int t=t1-1; t >= t0; t--) 
      {
	SPxi[t].Resize(ns);
      	for (unsigned int k=0; k<ns; k++)
	  SPxi[t](k) = (Pxi1[t+1](k) > 0.0) ? SPxi[t+1](k)/Pxi1[t+1](k) : 0.0;

	if (regime_process.IsInvariant())
	  SPxi[t] = eMultiply(SPxi[t], Q[rho]*Pxi[t]);
	else
	  SPxi[t] = eMultiply(SPxi[t], Q[t+1]*Pxi[t]);
      }

    return SPxi;
  } 

  // Smoothed state mean via backward recursion 
  vector< vector<TDV> > TLikelihood_TimeVaryingLinearStateSpace::SmoothedMean(unsigned int t0, unsigned int t1, const TDV &parameters, const vector<TDV> &SPxi) const
  {
    TDV z;
    unsigned int r, zeta, zeta_t;
    unsigned int nzeta=IEz[t1].size(), ns=nzeta*nr;

    vector<TDM> F=Get_F(parameters);
    vector<TDV> SPzeta(t1);
    vector< vector<TDV> > SEz(t1);
    vector< vector<TDV> > ISEz(t1+1);
    
    // Computes ISEz[t1][j]
    ISEz[t1].resize(nzeta);
    for (unsigned int j=0; j<nzeta; j++)
      ISEz[t1][j]=IEz[t1][j];

    for (unsigned int t=t1-1; t >= t0; t--)
      {
 	SEz[t].resize(ns);
	ISEz[t].resize(nzeta);
        // k = xi(t+1)
        for (unsigned int k=0; k<ns; k++)
	  {
	    // r = r(t+1)
	    r=(unsigned int)k/nzeta;

	    // zeta = zeta(t+1)
	    zeta=(unsigned int)k/nr;

	    // zeta_t = zeta(t)
	    zeta_t = k - r*nzeta; 

	    z=ISEz[t+1][zeta]-Ez1[t+1][k];
	
            // Use GeneralizedInverse(Ezz1[t+1][k]) because it can be singular
	    z = GeneralizedInverse(Ezz1[t+1][k])*z;
	    z = T(F[r]) * z;
	    z = IEzz[t][zeta_t] * z;
	    SEz[t][k]=IEz[t][zeta_t]+z;	      
	  }
     
        SPzeta[t].Resize(nzeta);
        for (unsigned int j=0; j<nzeta; j++)
	  {
	    SPzeta[t](j) = Sum(SPxi[t+1](I(j*nr,(j+1)*nr-1)));	
	    ISEz[t][j] = SPxi[t+1](j*nr)*SEz[t][j*nr];
	    for (unsigned int i=1; i<nr; i++)
              ISEz[t][j] += SPxi[t+1](j*nr+i)*SEz[t][j*nr+i]; 
	    ISEz[t][j] /= SPzeta[t](j);
          }

      }

    return ISEz;
  }

}
