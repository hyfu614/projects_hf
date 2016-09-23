#ifndef GENERIC_MODEL_STATESPACE_HEADER
#define GENERIC_MODEL_STATESPACE_HEADER

#include  <string>
#include "generic_model.hpp"
#include "DM.hpp"
#include "StateSpace.hpp"
#include "TTimeSeries.hpp"


class Generic_Model_TTimeSeries_StateSpace : public Generic_Model
{
protected: 
	TS::TTimeSeries_StateSpace &statespace; 
        
public:
        // constructors/destructors
        Generic_Model_TTimeSeries_StateSpace(const Generic_Model_TTimeSeries_StateSpace &Model) : statespace(Model.statespace) {};
        Generic_Model_TTimeSeries_StateSpace(TS::TTimeSeries_StateSpace &Statespace) : statespace(Statespace) {};
	~Generic_Model_TTimeSeries_StateSpace() {}; 
	
        virtual double log_posterior_function(const double *x, int n);
	virtual double log_likelihood_function(const double *x, int n);
	virtual double log_prior_function(const double *x, int n); 
	virtual bool DrawParametersFromPrior(double *x, int n);
	virtual int GetNumberParameters() const { return statespace.NumberParameters(); }; 
                
};


#endif
