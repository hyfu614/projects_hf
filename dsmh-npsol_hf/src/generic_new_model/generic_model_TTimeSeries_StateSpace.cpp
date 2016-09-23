#include  <string>
#include "DM.hpp"
#include "generic_model.hpp"
#include "StateSpace.hpp"
#include "TTimeSeries.hpp"
#include "generic_model_TTimeSeries_StateSpace.hpp"

using namespace std;
using namespace DM;
using namespace TS;

double Generic_Model_TTimeSeries_StateSpace::log_posterior_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return statespace.LogPosterior(parameters);
}

double Generic_Model_TTimeSeries_StateSpace::log_likelihood_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return statespace.LogLikelihood(parameters);
}

double Generic_Model_TTimeSeries_StateSpace::log_prior_function(const double *x, int n)
{
        DM::TDV parameters(n);
        
        for (int i=0; i<n; i++)
            parameters[i] = *(x+i);
        
        return statespace.LogPrior(parameters);
}

bool Generic_Model_TTimeSeries_StateSpace::DrawParametersFromPrior(double *x, int n)
{
        TDV parameters(x,n);
        
        parameters=statespace.DrawPrior();
                
        return true;
}