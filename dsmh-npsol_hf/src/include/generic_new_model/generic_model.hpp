#ifndef GENERIC_NEW_MODEL_HEADER
#define GENERIC_NEW_MODEL_HEADER

class Generic_Model
{
public:
	virtual double log_posterior_function(const double *x, int n) = 0; 
	virtual double log_likelihood_function(const double *x, int n) = 0; 
	virtual double log_prior_function(const double *x, int n) = 0; 
	virtual bool DrawParametersFromPrior(double *x, int n)=0; 
	virtual int GetNumberParameters() const = 0; 
	//virtual std::vector<DM::TDM> ImpulseResponse(const DM::TDV &parameters, unsigned int horizon) =0;
	//virtual DM::TDV LogConditionalLikelihoodVector(const DM::TDV &parameters)=0; 
	Generic_Model(){}
	virtual ~Generic_Model(){}
}; 

#endif
