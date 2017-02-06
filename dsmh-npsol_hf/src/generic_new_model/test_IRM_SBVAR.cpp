#include <fstream>
#include <iomanip>
#include <sstream>

#include <math.h>

#include "IRM.hpp"
#include "specification_io.hpp"
#include "sbvar.hpp"
#include "TFunction.hpp"
#include "TData.hpp"
#include "test_IRM_SBVAR.hpp"

using namespace std;
using namespace DM;

namespace TS
{
  TTimeSeries_IRM TTimeSeries_IRM_SBVAR_Specification(const std::string &filename)
  {
	TSpecification spec(filename);

    	string type=spec.String("//== IRM: Type");
    	if (!type.compare("state_space_response"))
      	{
		// get raw data
		TTimeSeriesDataRaw raw_data(spec);
	
        	// series
		vector<string> series=spec.Strings("//== IRM: Series");
		unsigned int n_var=series.size();
		if (n_var == 0)
	  		throw dw_exception("TTimeSeries_IRM_SBVAR_Specification() - at least one data series must be included");

		// get SBVAR number of lags
		unsigned int n_lag=spec.NonNegativeInteger("//== SBVAR: Number Lags");

		// begin and end index
		unsigned int begin_index=spec.NonNegativeInteger("//== IRM: Begin Index");
		unsigned int end_index=spec.NonNegativeInteger("//== IRM: End Index");

		// get data
		TDM Y=raw_data.Data(series,begin_index,end_index);

    		// Setup X - lagged variables
    		TDM X=raw_data.LaggedData(n_lag,series,begin_index,end_index);
 		X=HCat(X,Ones(X.Rows()));

		unsigned int n_pre=X.Cols();

    		// get SBVAR contemporanous restriction matrix
    		TDM AC;
    		vector< vector<bool> > AR;
    		spec.RestrictionMatrix(AR,AC,"//== SBVAR: Simple Restrictions A",' ');
    		if ((AC.Rows() != n_var) || (AC.Cols() != n_var))
      			throw dw_exception("TTimeSeries_SBVAR_Specification() - invalid size of contemporaneous restriction matrix");
    		if (AC != 0.0)
      			throw dw_exception("TTimeSeries_SBVAR_Specification() - restrictions must be linear");

    		// get SBVAR predetermined restriction matrix
    		TDM FC;
    		vector< vector<bool> > FR;
    		spec.RestrictionMatrix(FR,FC,"//== SBVAR: Simple Restrictions F",' ');
    		if ((FC.Rows() != n_var) || (FC.Cols() != n_pre))
      			throw dw_exception("TTimeSeries_SBVAR_Specification() - invalid size of predetermined restriction matrix");
    		if (FC != 0.0)
      			throw dw_exception("TTimeSeries_SBVAR_Specification() - restrictions must be linear");

		// create SBVAR equation mappings
    		TPArray<TLinear> Fnct;
    		TDM parameters=HCat(AC,FC);
    		for (unsigned int i=0; i < n_var; i++)
      		{
			I idx;
			for (unsigned int j=0; j < n_var; j++)
	  			if (AR[i][j]) idx(j);
			for (unsigned int j=0; j < n_pre; j++)
	  			if (FR[i][j]) idx(n_var+j);
			Fnct.push(new TLinear_insert(n_var+n_pre,idx));
      		}
    		TLinearProduct F(Fnct);

		// get SBVAR prior type
		string prior=spec.String("//== SBVAR: Prior Type");
		if (!prior.compare("Sims-Zha"))
		{
			TDV Hyperparameters=spec.Vector("//== SBVAR: Prior Hyperparameters");
    		
			// Setup IRM
  			TPArray<TRealFunction> f0(n_var);
  			for (unsigned int i=0; i < n_var; i++)
    				f0.push(state_space_response(n_var*n_lag+1,i));
  			TPArray<TRealFunction> s(n_var*n_var);
			for (unsigned int j=0; j < n_var; j++)
    				for (unsigned int i=0; i < n_var; i++)
      					s.push(state_space_response(n_var*n_lag+1,i));

			
			// variance matrix type
			string SigmaType=spec.String("//== IRM: Sigma");
			unsigned int sigma_type;
			if (!SigmaType.compare("identity"))
	  			sigma_type=0;
			else if (!SigmaType.compare("diagonal"))
	  			sigma_type=1;
			else if (!SigmaType.compare("full"))
	  			sigma_type=2;
			else
	  			throw dw_exception("//== IRM: Sigma -- unknown type (" + SigmaType + ")");
		
			TLikelihood_IRM irm_likelihood(Y,f0,s,0);
		
			// set IRM restriction matrix
			TDV c(irm_likelihood.NumberParameters(),0.0);
			I idx_irm;
			TDM B=VCat(VCat(Zeros(n_var,n_var*n_lag+1),HCat(Identity(n_var*(n_lag-1)),Zeros(n_var*(n_lag-1),n_var+1))),One(n_var*n_lag,n_var*n_lag+1));
  			TDV v=Cat(Vec(B),X(0,I(0,End)));
  			int k=0;
  			for (unsigned int i=0; i < n_var; k+=v.Dim(), i++) 
			{
				c(I(k,End),v);
				for (unsigned int j=0; j<n_var*n_lag+1; j++) 
				{
					for (unsigned m=0; m<n_var; m++) idx_irm(k+j*(n_var*n_lag+1)+m);
				}
			}
			TDV e=Zeros(n_var*n_lag+1);
			v=Cat(Vec(B),e);
  			for (unsigned int j=0; j < n_var; j++)
    			{
      				for (unsigned int i=0; i < n_var; k+=v.Dim(), i++) 
				{	
					c(I(k,End),v);
					for (unsigned int l=0; l<n_var*n_lag+1; l++)
						for (unsigned m=0; m<n_var; m++) idx_irm(k+l*(n_var*n_lag+1)+m);
				}
    			}
			k=n_var*(n_var*n_lag+1)*(n_var*n_lag+2)+(n_var*n_lag+1)*(n_var*n_lag+1);
			for (unsigned int j=0; j < n_var; j++)
    			{
      				for (unsigned int i=0; i < n_var; k+=v.Dim(), i++) 
				{	
					for (unsigned m=0; m<n_var; m++) idx_irm(k+m);
				}
    			}
	
			// set mapping functions from SBVAR prior parameters to IRM likelihood parameters
			TBvar2IrmFunction g(n_var,n_pre); 
			TAffine_insert f(c,idx_irm);
			TCompositeFunction ParaFunction(f,g);
	
  			return TTimeSeries_IRM(irm_likelihood,TDensity_SimsZha(Y,X,F,1,Hyperparameters),ParaFunction);
		}
		else
      			throw dw_exception("TTimeSeries_IRM_SBVAR_Specification() - unknown prior type (" + prior + ")");
	
	}
	else
      		throw dw_exception("TTimeSeries_IRM_SBVAR_Specification() - unknown type (" + type + ")");
  }

//=========================================================================================================

  TDV TBvar2IrmFunction::operator()(const TDV &x) const
  {
    if (x.Dim() != Domain())
      throw dw_exception("TBvar2IrmFunction::operator() - invalid dimension of argument"); 

    TDM A(n_var,n_var,false), F(n_var,n_pre,false);
    for (unsigned int i=0, k=0; i < n_var; A(i++,I(0,End),x(I(k,k+n_var-1))), k+=n_var+n_pre);
    for (unsigned int i=0, k=n_var; i < n_var; F(i++,I(0,End),x(I(k,k+n_pre-1))), k+=n_var+n_pre);

    TDV v1=Vec(A%F);
    TDM invA=Inverse(A);
    TDV parameters(Range());
    unsigned int k=0;
    for (unsigned int i=0; i<n_var*(n_var+1); i++, k+=v1.Dim())
        parameters(I(k,End),v1);
    for (unsigned int j=0; j<n_var; j++)
    {
      TDV v2=invA(I(0,End),j);
      for (unsigned int i=0; i<n_var; k+=v2.Dim(), i++) parameters(I(k,End),v2);
    }
    
    return parameters;
  }
}


