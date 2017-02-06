#ifndef __IRM_SBVAR_MAPPING__
#define __IRM_SBVAR_MAPPING__

#include "sbvar.hpp"
#include "IRM.hpp"
#include "TFunction.hpp"

namespace TS 
{
  TTimeSeries_IRM TTimeSeries_IRM_SBVAR_Specification(const std::string &filename);

  class TBvar2IrmFunction : public TFunction
  {
    protected:
      unsigned int n_var;
      unsigned int n_pre;

    public:
      TBvar2IrmFunction(const TBvar2IrmFunction &F) : n_var(F.n_var), n_pre(F.n_pre) {};
      TBvar2IrmFunction(unsigned int nvar, unsigned int npre) : n_var(nvar), n_pre(npre) {};
      virtual ~TBvar2IrmFunction() {};
      virtual TBvar2IrmFunction* Clone(void) const { return new TBvar2IrmFunction(*this); };

      virtual unsigned int Domain(void) const { return n_var*(n_var+n_pre); };
      virtual unsigned int Range(void) const { return n_var*n_var*((n_var+1)*n_pre+n_var); };

      // function evaluation
      virtual DM::TDV operator()(const DM::TDV &x) const;
 
  };
}

#endif
