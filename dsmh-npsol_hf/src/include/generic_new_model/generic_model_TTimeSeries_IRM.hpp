#ifndef GENERIC_MODEL_IRM_HEADER
#define GENERIC_MODEL_IRM_HEADER

#include  <string>
#include "generic_model.hpp"
#include "DM.hpp"
#include "IRM.hpp"


class Generic_Model_TTimeSeries_IRM : public Generic_Model
{
protected: 
	TS::TTimeSeries_IRM &irm; 
        
public:
        // constructors/destructors
        Generic_Model_TTimeSeries_IRM(const Generic_Model_TTimeSeries_IRM &Model) : irm(Model.irm) {};
        Generic_Model_TTimeSeries_IRM(TS::TTimeSeries_IRM &Irm) : irm(Irm) {};
	~Generic_Model_TTimeSeries_IRM() {}; 
	
        virtual double log_posterior_function(const double *x, int n);
	virtual double log_likelihood_function(const double *x, int n);
	virtual double log_prior_function(const double *x, int n); 
	virtual bool DrawParametersFromPrior(double *x, int n);
	virtual int GetNumberParameters() const { return irm.NumberParameters(); };   
};


#endif
