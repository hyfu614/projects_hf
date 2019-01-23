#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "dw_rand.h"
#include "DM.hpp"
#include "msv_msre.hpp"
#include "msre_statespace.hpp"

using namespace std;
using namespace DM;
using namespace TS;

double verify_solution(vector<TDM> &A, vector<TDM> &B, vector<TDM> &Psi, vector<TDM> &V, vector<TDM> &F1, vector<TDM> &F2, vector<TDM> &G1, vector<TDM> &G2, TDM &Ptr, unsigned int s)
{
  unsigned int n=A[0].Rows();  
  unsigned int nr=Ptr.Rows();
  double diff=0.0, tmp;
  TDM X=Zeros(n,n), Y, F, G, W;
  
  X(I(n-s,n-1),I(n-s,n-1),Identity(s));
  for (unsigned int i=0; i<nr; i++)
    {
      Y=A[i]*V[i];
      X(I(0,n-1),I(0,n-s-1),Y);
      F=VCat(F1[i],F2[i])-Inverse(X)*B[i];
      tmp=Norm(F);
      cout << "Norm(F[" << i << "]): " << tmp << endl;
      if (tmp > diff)
	diff=tmp;

      G=VCat(G1[i],G2[i])-Inverse(X)*Psi[i];
      tmp=Norm(G);
      cout << "Norm(G[" << i << "]): " << tmp << endl;
      if (tmp > diff)
	diff=tmp;
    }

  for (unsigned int i=0; i<nr; i++)
    {
      TDM Z=Zeros(s,n);
      for (unsigned int j=0; j<nr; j++)
        Z += Ptr(i,j)*F2[j];
      W=Z*V[i]; 
      tmp=Norm(W);
      cout << "Norm(W[" << i << "]): " << tmp << endl;
      if (tmp > diff)
	diff=tmp;
    }

  return diff;
}

void testMsvMsre(void)
{
  unsigned int n,s,nr,nepsilon,nu,rho,nobs,ny;
  n=9;
  s=2;
  nr=2;
  nepsilon=4;
  nu=5;
  rho=3;
  nobs=100;
  ny=6;

  vector<TDV> a(nr);
  vector<TDM> A(nr),B(nr),Psi(nr),H(nr),Phiy(nr);
  TDM P;
  
  for (unsigned int i=0; i<nr; i++)
    {
      a[i]=4*RandomUniform(ny)-2;
      A[i]=RandomNormal(n,n);
      B[i]=RandomNormal(n,n);
      Psi[i]=VCat(RandomNormal(n-s,nepsilon), Zeros(s,nepsilon));
      H[i]=RandomNormal(ny,n);
      Phiy[i]=RandomNormal(ny,nu);
    }

  P=RandomUniform(nr,nr);
  for (unsigned int i=0; i<nr; i++)
    P(i,I(0,End),P(i,I(0,End))/Sum(P(i,I(0,End))));
  
  TDM Data=RandomNormal(nobs,ny);

  TMSV_MSRE msre(nr, s, n, nepsilon);
  
  TDV msre_parameters(msre.NumberParameters());
  TDV me_parameters(nr*ny*(1+n+nu));
  unsigned int offset=0, offset_me=0;
  for (unsigned int i=0; i<nr; i++)
    {
      msre_parameters(I(offset,offset+n*n-1),Vec(A[i]));
      offset += n*n;
      msre_parameters(I(offset,offset+n*n-1),Vec(B[i]));
      offset += n*n;
      msre_parameters(I(offset,offset+n*nepsilon-1),Vec(Psi[i]));
      offset += n*nepsilon;
      me_parameters(I(offset_me,offset_me+ny-1),a[i]);
      offset_me += ny;
      me_parameters(I(offset_me,offset_me+ny*n-1),Vec(H[i]));
      offset_me += ny*n;
      me_parameters(I(offset_me,offset_me+ny*nu-1),Vec(Phiy[i]));
      offset_me += ny*nu;
    }
  msre_parameters(I(offset,offset+nr*nr-1),Vec(P));

  
  TDV initial_values=Zeros(nr*s*(n-s));
  int nItrn=msre.MSVsolution(msre_parameters,initial_values);
  cout << "Number of iterations is: " << nItrn << endl; 
  vector<TDM> V,F1,F2,G1,G2;
  V=msre.Get_V();
  F1=msre.Get_F1();
  F2=msre.Get_F2();
  G1=msre.Get_G1();
  G2=msre.Get_G2();
  double diff=verify_solution(A, B, Psi, V, F1, F2, G1, G2, P, s);
  cout << "The largest difference is: " << diff << endl;
 
  /*
  TDM X=msre.StartingValules(A, B);
  cout << "Number of starting values: " << X.Rows() << endl;
  for (unsigned int i=0; i<X.Rows(); i++)
    {
      int nItrn=msre.MSVsolution(msre_parameters,X(i,I(0,End)));
      cout << "Number of iterations is: " << nItrn << endl; 
     
      V=msre.Get_V();
      F1=msre.Get_F1();
      F2=msre.Get_F2();
      G1=msre.Get_G1();
      G2=msre.Get_G2();
      double diff=verify_solution(A, B, Psi, V, F1, F2, G1, G2, P, s);
      cout << "The largest difference is: " << diff << endl;
    }
  */

  TDV parameters=Cat(msre_parameters,me_parameters);
  TLikelihood_MSRE msre_ss(Data,nr,s,n,nepsilon,nu,rho);
  cout << "LogLikilihood: " << msre_ss.LogLikelihood(parameters) << endl;

}

int main(void)
{
  try
    {
      dw_initialize_generator(0);
	
      testMsvMsre();
      
    }
  catch (dw_exception &e)
    {
      cout << "terminated with exception: " << e.what() << endl;
    }

  return 0;
}

