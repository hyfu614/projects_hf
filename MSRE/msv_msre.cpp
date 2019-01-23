#include <fstream>
#include <sstream>
#include <math.h>

#include "dw_rand.h"
#include "dw_math.h"
#include "TRegimeProcesses.hpp"
#include "msv_msre.hpp"

using namespace std;
using namespace DM;

namespace TS
{
  TMSV_MSRE::TMSV_MSRE(const TMSV_MSRE &Msre)
    : nr(Msre.nr),
      n(Msre.n),
      s(Msre.s),
      nepsilon(Msre.nepsilon),
      n_parameters(Msre.n_parameters),
      Ptr(Msre.Ptr),
      regime_process(Msre.regime_process),
      V(Msre.V),
      F1(Msre.F1),
      F2(Msre.F2),
      G1(Msre.G1),
      G2(Msre.G2)
  { }

  TMSV_MSRE::TMSV_MSRE(unsigned int nRegimes, unsigned int nExpErrors, unsigned int nStates, unsigned int nShocks) 
    : nr(nRegimes),
      n(nStates),
      s(nExpErrors),
      nepsilon(nShocks),
      n_parameters(nr*(n*(2*n+nepsilon)+nr)),
      Ptr(nr,nr),
      regime_process(nRegimes),
      V(nr),
      F1(nr),
      F2(nr),
      G1(nr),
      G2(nr)
  { }
   
  int TMSV_MSRE::MSVsolution(const TDV &parameters, const TDV &initial_values, unsigned int max_count, double tol) const
  {
    if (parameters.Dim() != n_parameters)
      throw dw_exception("TMSV_MSRE::MSVsolution() - invalid number of MSRE parameters");

    if (initial_values.Dim() != nr*s*(n-s))
      throw dw_exception("TMSV_MSRE::MSVsolution() - invalid number of initial values");

    bool cont=true;

    Ptr=Get_Ptr(parameters);
    vector<TDM> A=Get_A(parameters), B=Get_B(parameters), Psi=Get_Psi(parameters);
    vector< vector<TDM> > C(nr);
    for (unsigned int i=0; i<nr; i++)
      {
  	C[i].resize(nr);
        for (unsigned int j=0; j<nr; j++)
	  {
            C[i][j].Resize(n,n);
	    C[i][j]=Ptr(i,j)*MultiplyInverse(B[j],A[i]);
    	  }
      }
    
    TDV x(initial_values);    
    vector<TDM> X(nr);    
    TDM D=Zeros(nr*s*(n-s),nr*s*(n-s));
    TDV f(nr*s*(n-s));

    unsigned int count=1; 
    while (cont)
      {
 	for (unsigned int i=0; i<nr; i++)
	  X[i]=Reshape(x(I(i*s*(n-s),(i+1)*s*(n-s)-1)),s,n-s);
	
	for (unsigned int i=0; i<nr; i++)
	  for (unsigned int j=0; j<nr; j++)
	    {
	      TDM W1,W2;
	      W1=C[i][j]*VCat(Identity(n-s), -X[i]);
	      W2=W1(I(0,n-s-1),I(0,n-s-1));
	      D(I(i*s*(n-s),(i+1)*s*(n-s)-1),I(j*s*(n-s),(j+1)*s*(n-s)-1),Kron(T(W2),Identity(s)));
	      if (i == j)
		{
		  TDM W3,W4,W5;
		  W3=Zeros(s,n);
		  for (unsigned int k=0; k<nr; k++)
		    W3 += HCat(X[k], Identity(s)) * C[i][k];
		  W4=-W3(I(0,s-1),I(n-s,n-1));
		  W5=D(I(i*s*(n-s),(i+1)*s*(n-s)-1),I(j*s*(n-s),(j+1)*s*(n-s)-1))+Kron(Identity(n-s),W4);
		  D(I(i*s*(n-s),(i+1)*s*(n-s)-1),I(j*s*(n-s),(j+1)*s*(n-s)-1),W5);
		}
	    }

	TDV y;
	for (unsigned int i=0; i<nr; i++)
	  {
	    TDM Mf=Zeros(s,n-s);
	    for (unsigned int j=0; j<nr; j++)
	      Mf += (HCat(X[j], Identity(s))*C[i][j]) * VCat(Identity(n-s), -X[i]);
	    f(I(i*s*(n-s),(i+1)*s*(n-s)-1),Vec(Mf));
	  }

  	y=InverseMultiply(D,f);
	x -= y;
	
	if ((count > max_count) || (Norm(f) < tol))
	  cont=false;
	else
	  count++;

      }

    int err;
    if (Norm(f) < tol)
      err = count;
    else
      err = -count;

    TDM W,F,G;
    for (unsigned int i=0; i<nr; i++)
      {
	X[i]=Reshape(x(I(i*s*(n-s),(i+1)*s*(n-s)-1)),s,n-s);
 	V[i]=InverseMultiply(A[i], VCat(Identity(n-s), -X[i])); 
	W=HCat(VCat(Identity(n-s), X[i]), VCat(Zeros(n-s,s), Identity(s)));
	F=W*B[i];
	F1[i]=F(I(0,n-s-1),I(0,n-1));
     	F2[i]=F(I(n-s,n-1),I(0,n-1));
	G=W*Psi[i];
	G1[i]=G(I(0,n-s-1),I(0,nepsilon-1));
	G2[i]=G(I(n-s,n-1),I(0,nepsilon-1));
      }

    return err;
  }
 
  vector<TDM> TMSV_MSRE::Get_A(const TDV &parameters) const
  {
    vector<TDM> A(nr);
    unsigned int stride=n*(2*n+nepsilon);
    for(unsigned int i=0, k=0; i<nr; i++, k+=stride)
      A[i]=Reshape(parameters(I(k,k+n*n-1)),n,n);

    return A;
  }

  vector<TDM> TMSV_MSRE::Get_B(const TDV &parameters) const
  {
    vector<TDM> B(nr);
    unsigned int stride=n*(2*n+nepsilon);
    for(unsigned int i=0, k=n*n; i<nr; i++, k+=stride)
      B[i]=Reshape(parameters(I(k,k+n*n-1)),n,n);

    return B;
  }

  vector<TDM> TMSV_MSRE::Get_Psi(const TDV &parameters) const
  {
    vector<TDM> Psi(nr);
    unsigned int stride=n*(2*n+nepsilon);
    for(unsigned int i=0, k=2*n*n; i<nr; i++, k+=stride)
      Psi[i]=Reshape(parameters(I(k,k+n*nepsilon-1)),n,nepsilon);

    return Psi;
  }

  // transition matrix
  TDM TMSV_MSRE::Get_Ptr(const TDV &parameters) const
  {
    TDM P;
    unsigned int offset=nr*(n*(2*n+nepsilon));
    P=Reshape(parameters(I(offset,offset+nr*nr-1)),nr,nr);
  
    return P;
  }

  
  TDM TMSV_MSRE::StartingValules(const vector<TDM> &A, const vector<TDM> &B, bool all) const
  {
    vector<TDM> Xi(nr);
    TDV n_total=Zeros(nr);
    unsigned int n_starting_values=1, n_explosive;
   
    for (unsigned int i=0; i<nr; i++)
      {
	Xi[i]=AllInvariantSubspaces(A[i],B[i],n_explosive);
     
	if (all)
	  n_total(i)=Xi[i].Rows();
	else
	  n_total(i)=Xi[i].Rows()-n_explosive;
 	n_starting_values *= n_total(i);	
      }
   
    unsigned int n_parameters=nr*s*(n-s);
    TDM X=Zeros(n_starting_values,n_parameters);
    TDV y=Zeros(n_parameters);
    TDV sel=Zeros(nr);
   
    for (unsigned int i=0; i<n_starting_values; i++)
      {
	for (unsigned int j=0, k=0; j<nr; j++, k += s*(n-s))
	  y(I(k,k+s*(n-s)-1),Xi[j],sel(j),I(0,Xi[j].Cols()-1)); 
	  
	X(i,I(0,End),y);
	unsigned int j=nr;
	while (j > 0)
	  {
	    j--;
	    if (sel(j) == n_total(j)-1)
	      sel(j)=0;
	    else
	      {
		sel(j)++;
 		break;
	      }
	  }
      }
     
      return X;
    
  }
  
  TDM TMSV_MSRE::AllInvariantSubspaces(const TDM &A0, const TDM &B0, unsigned int &n_explosive) const
  {
    TDM A1, A2, T, U, TT, UU;
    A1=A0(I(0,n-s-1),I(0,A0.Cols()-1));
    A2=A0(I(n-s,n-1),I(0,A0.Cols()-1));

    TDV evRe, evIm, evAbs, sortedIdx; 
    I idxReal, idxComp;
    Schur(TT,UU,evRe,evIm,InverseMultiply(A0,B0));
    // sort real and complex indices
    evAbs=eMultiply(evRe,evRe)+eMultiply(evIm,evIm);
    SortIdx(evAbs, sortedIdx);
 
    unsigned int i=0, nReal=0, nComp=0;
    while (i < n)
      {
	if (fabs(evIm((unsigned int)sortedIdx(i))) < 1e-10)
	  {
	    idxReal(sortedIdx(i));
	    nReal++;
	    i++;
	  }
	else
	  {
	    idxComp(sortedIdx(i));
	    nComp++;
	    i += 2;
	  }
      }
 
    // minimum and maximum number of complex conjugate pairs in invariant subspace    
    unsigned int nMinComp=fmax(0, floor(((n-s)-nReal+1)/2));
    unsigned int nMaxComp=fmin(nComp, floor((n-s)/2));
    if ((n-s) < 2*nMinComp)
      throw dw_exception("TMSV_MSRE::AllInvariantSubspaces() - No real solution");

    // number of invariant subspaces
    unsigned int nInvariant=0;
    for (unsigned int ii=nMinComp; ii<nMaxComp+1; ii++)
       nInvariant += nchoosek(nComp,ii)*nchoosek(nReal,n-s-2*ii);

    unsigned int nParameters, nExplosive=0, nNonExplosive=0;
    nParameters=s*(n-s);
    TDM ExplosivePara, NonExplosivePara;
    ExplosivePara=Zeros(nInvariant,nParameters);
    NonExplosivePara=Zeros(nInvariant,nParameters);
 
    for (unsigned int ii=nMinComp; ii<nMaxComp+1; ii++)
      {
	bool cont=true;
    	unsigned int jj=n-s-2*ii;
    	TDV complexChosen, realChosen;
	complexChosen=Cat(Ones(ii), Zeros(nComp-ii));
	realChosen=Cat(Ones(jj), Zeros(nReal-jj));

	while (cont)
	  {
	    bool explosive=false;
	    vector<int> select(n,0);
	    for (unsigned int i=0; i<nComp; i++)
	      {
		if (complexChosen(i))
		  {
		    select[idxComp[i]]=1;
		    select[idxComp[i]+1]=1;
		    if (evAbs(idxComp[i]) >= 1)
		      explosive=true;
		  }
	      } 

	    for (unsigned int i=0; i<nReal; i++)
	      {
		if (realChosen(i))
		  {
		    select[idxReal[i]]=1;
		    if (evAbs(idxReal[i]) >= 1)
		      explosive=true;
		  }
	      } 

	    // reorder schur and save
	    ReorderSchur(T,U,evRe,evIm,TT,UU,select);	
	    TDM U1=U(I(0,End),I(0,n-s-1));
	    TDM Z=-A2 * MultiplyInverse(U1,A1*U1);
	    if (explosive)
	      {
		ExplosivePara(nExplosive,I(0,End),Vec(Z));
		nExplosive++;
	      }
	    else
	      {
		NonExplosivePara(nNonExplosive,I(0,End),Vec(Z));
		nNonExplosive++;
	      }
   
	    // increment selection
	    unsigned int k=0;
	    while ((nReal > k) && (realChosen(nReal-k-1) == 1))
	      {
		k++;
		realChosen(nReal-k)=0;
	      }
	    unsigned int i=nReal-k;
	    while ((i > 0) && (realChosen(i-1) == 0))
	      i--;
	    if (i > 0)
	      {
		realChosen(i-1)=0;
		realChosen(I(i,i+k),Ones(k+1));
	      }
	    else
	      {
		realChosen(I(i,i+k-1),Ones(k));
	 	k=0;
		while ((nComp > k) && (complexChosen(nComp-k-1) == 1))
		  {
		    k++;
		    complexChosen(nComp-k)=0;
		  }
		i=nComp-k;
		while ((i > 0) && (complexChosen(i-1) == 0))
	      	  i--;

	    	if (i > 0)
	          {
		    complexChosen(i-1)=0;
		    complexChosen(I(i,i+k),Ones(k+1));
	          }
	        else
		  cont=false;
	      }
 
	  }
      }
  	    
    TDM X=VCat(NonExplosivePara(I(0,nNonExplosive-1),I(0,End)), ExplosivePara(I(0,nExplosive-1),I(0,End)));
    
    return X;
  }

//===============================================================================

  TDV SortIdx(const TDV &v, TDV &sorted_id, bool ascending)
  {
    TDV w(v);
    unsigned int n=v.Dim();
    sorted_id=Trend(n,1);

    for (unsigned int i=0; i<n-1; i++)
    {
      for (unsigned int j=0; j<n-i-1; j++)
      {
        if (ascending)
        {
          if(w(j) > w(j+1))
	  {
	    swap(w(j), w(j+1));
	    swap(sorted_id(j), sorted_id(j+1));
	  }
        }
        else
        {
	  if(w(j) < w(j+1))
	  {
	    swap(w(j), w(j+1));
	    swap(sorted_id(j), sorted_id(j+1));
	  }
        }
      }
    }

    return w;
  }

  unsigned int nchoosek(unsigned int n, unsigned int k)
  {
    if (n < k)
      dw_exception("nchoosek() - n can not be less than k");
 
    unsigned int b;
    b = factorial(n) / (factorial(n-k)*factorial(k));

    return b;
  }

  unsigned int factorial(unsigned int n)
  {
    unsigned int f=1;
    for (unsigned int i=0; i<n; i++)
      f *= i+1;

    return f; 
  }

  
}
