#include <fstream>

#include "DM.hpp"
#include "TData.hpp"
#include "StateSpace.hpp"
#include "specification_io.hpp"
//#include "dw_rand.h"

using namespace std;
using namespace DM;
using namespace TS;

int main(void)
{
  try
    {
      string specification_file="statespace_specification_test.txt";

      TTimeSeries_StateSpace statespace=TTimeSeries_StateSpace_Specification(specification_file);
      cout << "n_obs=" << statespace.NumberObservations() << "\n";
     
      vector<double> b={0.1, 0.2, 0.3, 0.4, 0.5};
      vector<double> F={1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1};
      
      TLinear_inclusion &f=dynamic_cast<TLinear_inclusion&>(statespace.Function());	
      cout << "domain=" << f.Domain() << endl;

      TDV parameters(f.Domain());
      unsigned int nx=statespace.NumberStates();
      
      for (unsigned int i=0; i<nx; i++)
	parameters(i)=b[i];
      for (unsigned int i=0; i<nx*nx; i++)
	parameters(nx+i)=F[i];

      cout << "a=" << statespace.Geta(parameters) << endl;
      cout << "H=\n" << statespace.GetH(parameters) << endl;
      cout << "b=" << statespace.Getb(parameters) << endl;
      cout << "F=\n" << statespace.GetF(parameters) << endl;
      cout << "Phi=\n" << statespace.GetPhi(parameters) << endl;

      double likelihood;
      likelihood = statespace.LogLikelihood(parameters);
      cout << "likelihood: " << likelihood << endl;

      cout << "likelihood vector:\n" << statespace.LogConditionalLikelihoodVector(parameters) << endl;

      // cout << "F =\n" << F << endl;

      // TDV likelihood_parameters=Vec(T(HCat(A,F)));
      // TLinear &f=dynamic_cast<TLinear&>(sbvar.Function());
      // TLinear_general f_inv=f.GeneralizedInverse();
      // TDV parameters=f_inv(likelihood_parameters);
      // cout << "parameters =\n" << parameters;

      // cout << "\nlog likelihood = " << sbvar.LogLikelihood(parameters) << endl;
      // cout << "\nlog prior = " << sbvar.LogPrior(parameters) << endl;
      // cout << "\nlog posterior = " << sbvar.LogPosterior(parameters) << endl;




      // TDV parameters=RandomNormal(sbvar.NumberParameters());
      // TDM A=sbvar.GetA(sbvar.LikelihoodParameters(parameters)), F=sbvar.GetF(sbvar.LikelihoodParameters(parameters));
      // cout << "\nparameters =\n" << parameters;
      // cout << "\nlikelihood parameters =\n" << sbvar.LikelihoodParameters(parameters) << "\nA =\n" << A << "\nF=\n" << F;
      // cout << "\ndiff parameters = " << Norm(sbvar.LikelihoodParameters(parameters)-sbvar.GetParameters(A,F)) << endl;

      // parameters=sbvar.DrawPrior();
      // A=sbvar.GetA(sbvar.LikelihoodParameters(parameters));
      // F=sbvar.GetF(sbvar.LikelihoodParameters(parameters));
      // cout << "\nrandom prior draw of parameters =\n" << parameters;
      // cout << "\nlikelihood parameters =\n" << sbvar.LikelihoodParameters(parameters) << "\nA =\n" << A << "\nF=\n" << F;
      // cout << "\ndiff parameters = " << Norm(sbvar.LikelihoodParameters(parameters)-sbvar.GetParameters(A,F)) << endl;

      // cout << "\nlog likelihood = " << sbvar.LogLikelihood(parameters) << endl;
      // cout << "\nlog prior = " << sbvar.LogPrior(parameters) << endl;
      // cout << "\nlog posterior = " << sbvar.LogPosterior(parameters) << endl;

      //string specification_file="tvsbvar_model_file.txt";
      //dw_initialize_generator(0);

      // // Test raw data 
      // TTimeSeriesDataRaw raw_data=TTimeSeriesDataRaw_Specification(specification_file);
      // vector< vector<string> > data_strings=raw_data.Data();
      // for (unsigned int t=0; t < data_strings.size(); cout << endl, t++) 
      // 	for (unsigned int i=0; i < data_strings[t].size(); i++) cout << data_strings[t][i] << ',';
      // cout << "\n\n";
      // vector<string> dates=raw_data.Dates();
      // for (unsigned int t=0; t < raw_data.NumberObservations(); t++) cout << dates[t] << endl; 
      // cout << "\n\n";
      // vector<string> series=raw_data.SeriesNames();
      // for (unsigned int j=0; j < raw_data.NumberSeries(); j++) cout << series[j] << endl;
      // cout << "\n\n";

      //TLikelihood_SBVAR likelihood=TLikelihood_SBVAR_Specification(specification_file);

      //TTimeSeries_SBVAR sbvar=TTimeSeries_SBVAR_Specification(specification_file);

      // // Test Dirichlet density
      // int dirichlet_dim=10;
      // TDV Hyperparameters=10*RandomUniformVector(dirichlet_dim);
      // TDirichletDensity density(Hyperparameters);
      // TDV sum(dirichlet_dim,0.0);
      // int N=1000000;
      // for (int i=0; i < N; i++)
      // 	sum+=density.Probabilities(density.Draw());
      // cout << density.Mean() - (1.0/(double)N)*sum;
      // exit(0);

      // TRegimeProcess& RP=*TRegimeProcess_Specification(specification_file);
      // TFunction& f=*TRegimeProcess_LikelihoodFunction_Specification(specification_file);
      // TDensity& prior=*TRegimeProcess_Prior_Specification(specification_file);

      // TIndependentRegimeProcesses& IRP=dynamic_cast<TIndependentRegimeProcesses&>(RP);

      // TDV parameters=prior.Draw();
      // TDV likelihood_parameters=f(parameters);
      // cout << "parameters: " << parameters << endl;
      // cout << "likelihood parameters: " << likelihood_parameters << endl;
      // cout << RP.InitialProbabilities(likelihood_parameters) << endl;
      // cout << RP.TransitionMatrix(likelihood_parameters,0) << endl;
 
      // for (unsigned int i=0; i < IRP.NumberProcesses(); i++)
      // 	{
      // 	  cout << IRP[i].InitialProbabilities(IRP.Parameters(i,likelihood_parameters)) << endl;
      // 	  cout << IRP[i].TransitionMatrix(IRP.Parameters(i,likelihood_parameters),0) << endl;
      // 	}

      // vector< vector<unsigned int> > idx=IRP.RegimeIndexes();
      // for (unsigned int i=0; i < IRP.NumberRegimes(); i++)
      // 	{
      // 	  for (unsigned int j=0; j < IRP.NumberProcesses(); j++)
      // 	    cout << idx[i][j] << " ";
      // 	  cout << endl;
      // 	}



      // vector<string> specs;
      // vector< vector<string> > values;
      // ParseSpecificationFile(specification_file,specs,values);
      // TDV parameters=ParseVector(specs,values,"//== SBVAR: Default Parameters");
      // // TDV parameters=RandomNormalVector(sbvar.NumberParameters());
      // cout << "log likelihood: " << sbvar.LogLikelihood(parameters) << endl;
      // cout << "parameters:\n" << parameters << endl;
      // cout << "A:\n" << sbvar.GetA(parameters) << endl;
      // cout << "F:\n" << sbvar.GetF(parameters) << endl;

      //TRegimeProcess_endogenous_linear process;



      //TRegimeProcessArray a=SetupArrayRegimeProcess("restriction_rp_endo_test.txt");

    }
  catch (dw_exception &e)
    {
      std::cout << "Exception: " << e.what() << std::endl;
    }
  return 0;
}
