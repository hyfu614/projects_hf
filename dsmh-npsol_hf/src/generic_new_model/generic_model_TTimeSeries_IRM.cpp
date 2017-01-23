#include  <string>
#include "DM.hpp"
#include "generic_model.hpp"
#include "IRM.hpp"
#include "generic_model_TTimeSeries_IRM.hpp"

using namespace std;
using namespace DM;
using namespace TS;

double Generic_Model_TTimeSeries_IRM::log_posterior_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return irm.LogPosterior(parameters);
}

double Generic_Model_TTimeSeries_IRM::log_likelihood_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return irm.LogLikelihood(parameters);
}

double Generic_Model_TTimeSeries_IRM::log_prior_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return irm.LogPrior(parameters);
}

bool Generic_Model_TTimeSeries_IRM::DrawParametersFromPrior(double *x, int n)
{
        TDV parameters;
        
        parameters=irm.DrawPrior();
        for (int i=0; i<n; i++)
            *(x+i) = parameters[i];
                
        return true;
}
