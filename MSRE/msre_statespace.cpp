#include <fstream>
#include <sstream>
#include <math.h>

#include "dw_rand.h"
#include "dw_math.h"
#include "TRegimeProcesses.hpp"
#include "specification_io.hpp"
#include "msre_statespace.hpp"

using namespace std;
using namespace DM;

namespace TS
{
  TLikelihood_MSRE::TLikelihood_MSRE(const TLikelihood_MSRE &Msre_SS)
    : nr(Msre_SS.nr),
      n(Msre_SS.n),
      s(Msre_SS.s),
      nepsilon(Msre_SS.nepsilon),
      ny(Msre_SS.ny),
      nerrors(Msre_SS.nerrors),
      rho(Msre_SS.rho),
      nparameters(Msre_SS.nparameters),
      Ptr(Msre_SS.Ptr),
      regime_process(Msre_SS.regime_process),
      data(Msre_SS.data),
      lcli(data.Rows()),
      tv_statespace(Msre_SS.tv_statespace)
  { }

  TLikelihood_MSRE::TLikelihood_MSRE(const TDM Data, unsigned int nRegimes, unsigned int nExpErrors, unsigned int nStates, unsigned int nShocks, unsigned int nMeaErrors, unsigned int nLags) 
    : nr(nRegimes),
      n(nStates),
      s(nExpErrors),
      nepsilon(nShocks),
      ny(Data.Cols()),
      nerrors(nMeaErrors),
      rho(nLags),
      nparameters(nr*(n*(2*n+nepsilon)+ny*(1+n+nerrors)+nr)),
      Ptr(nr,nr),
      regime_process(nRegimes),
      data(Data),
      lcli(Data.Rows()),
      tv_statespace(*(TLikelihood_TimeVaryingLinearStateSpace(Data,nr,nerrors,n,nepsilon,rho,regime_process).Clone()))
  { }
 
  double TLikelihood_MSRE::LogConditionalLikelihood(const TDV &parameters, unsigned int t) const
  {
    if (t < 0 || t > NumberObservations()) 
      throw dw_exception("TLikelihood_MSRE::LogConditionalLikelihood() - invalid time t");

    TDV statespace_parameters=Get_StateSpace_parameters(parameters);
    lcli=tv_statespace.LogConditionalLikelihoodVector(statespace_parameters);

    return lcli(t);
  }

  TDV TLikelihood_MSRE::LogConditionalLikelihoodVector(const TDV &parameters) const
  {
    TDV statespace_parameters=Get_StateSpace_parameters(parameters);
    lcli=tv_statespace.LogConditionalLikelihoodVector(statespace_parameters);
    
    return lcli;
  }

  double TLikelihood_MSRE::LogLikelihood(const TDV &parameters) const
  {
    TDV statespace_parameters=Get_StateSpace_parameters(parameters);
    lcli=tv_statespace.LogConditionalLikelihoodVector(statespace_parameters);

    return Sum(lcli);
  }

  TDV TLikelihood_MSRE::Get_StateSpace_parameters(const TDV &parameters) const
  {
    if (parameters.Dim() != nparameters)
      throw dw_exception("TLikelihood_MSRE::Get_StateSpace_parameters() - invalid number of parameters");
    
    TMSV_MSRE msv(nr,s,n,nepsilon);
    TDV msre_parameters=parameters(I(0,msv.NumberParameters()-1));
    TDV me_parameters=parameters(I(msv.NumberParameters(),End));		// measurement equation parameters
    TDV initial_values=Zeros(nr*s*(n-s));
  
    int err=msv.MSVsolution(msre_parameters,initial_values);
    if (err < 0)
      throw dw_exception("TLikelihood_MSRE::Get_StateSpace_parameters() - cannot find the MSV solutions of MSRE model");

    vector<TDM> V,F1,F2,G1,G2;
    V=msv.Get_V();
    F1=msv.Get_F1();
    F2=msv.Get_F2();
    G1=msv.Get_G1();
    G2=msv.Get_G2();

    unsigned int n_meparameters=ny*(1+n+nerrors);
    unsigned int stride=ny*(1+n+nerrors)+n*(1+n+nepsilon);
    TDV statespace_parameters(nr*(stride+nr));

    TDM F, Phiz;
    for (unsigned int i=0, offset=0; i<nr; i++)
      {
  	statespace_parameters(I(offset,offset+n_meparameters-1),me_parameters(I(i*n_meparameters,(i+1)*n_meparameters-1)));
	offset += n_meparameters;
	F=V[i]*F1[i];
	Phiz=V[i]*G1[i];
        statespace_parameters(I(offset,offset+n-1),Zeros(n));
	offset += n;
 	statespace_parameters(I(offset,offset+n*n-1),Vec(F));
	offset += n*n;
	statespace_parameters(I(offset,offset+n*nepsilon-1),Vec(Phiz));
	offset += n*nepsilon;
      } 
    statespace_parameters(I(nr*stride,nr*(stride+nr)-1),msre_parameters(I(nr*n*(2*n+nepsilon),End)));

    return statespace_parameters;
  }

}  
