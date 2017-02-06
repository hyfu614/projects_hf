#ifndef GENERIC_MODEL_TIMESERIES_HEADER
#define GENERIC_MODEL_TIMESERIES_HEADER

#include  <string>
#include "generic_model.hpp"
#include "DM.hpp"
#include "TTimeSeries.hpp"


class Generic_Model_TTimeSeries : public Generic_Model
{
protected: 
	TS::TTimeSeries &timeseries; 
        
public:
        // constructors/destructors
        Generic_Model_TTimeSeries(const Generic_Model_TTimeSeries &Model) : timeseries(Model.timeseries) {};
        Generic_Model_TTimeSeries(TS::TTimeSeries &TimeSeries) : timeseries(TimeSeries) {};
	~Generic_Model_TTimeSeries() {}; 
	
        virtual double log_posterior_function(const double *x, int n);
	virtual double log_likelihood_function(const double *x, int n);
	virtual double log_prior_function(const double *x, int n); 
	virtual bool DrawParametersFromPrior(double *x, int n);
	virtual int GetNumberParameters() const { return timeseries.NumberParameters(); }; 
                
};


#endif
