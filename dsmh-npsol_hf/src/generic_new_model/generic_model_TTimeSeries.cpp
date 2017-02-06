#include  <string>
#include "DM.hpp"
#include "generic_model.hpp"
#include "TTimeSeries.hpp"
#include "generic_model_TTimeSeries.hpp"

using namespace std;
using namespace DM;
using namespace TS;

double Generic_Model_TTimeSeries::log_posterior_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return timeseries.LogPosterior(parameters);
}

double Generic_Model_TTimeSeries::log_likelihood_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return timeseries.LogLikelihood(parameters);
}

double Generic_Model_TTimeSeries::log_prior_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return timeseries.LogPrior(parameters);
}

bool Generic_Model_TTimeSeries::DrawParametersFromPrior(double *x, int n)
{
        TDV parameters;
        
        parameters=timeseries.DrawPrior();
        for (int i=0; i<n; i++)
            *(x+i) = parameters[i];
                
        return true;
}
