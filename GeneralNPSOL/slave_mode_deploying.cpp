#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <mpi.h>

#include "specification_io.hpp"
#include "mpi_constant.hpp"
#include "General_Maximization.hpp"

using namespace std;
using namespace DM;
using namespace TS;


TTimeSeries* MinusLogPosterior_NPSOL::model;

void MinusLogPosterior_NPSOL::Maximize_Posterior(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{ 
  TDV xx(*n);
  memcpy(xx.Vector(),x,(*n)*sizeof(double));
  *f = -model->LogPosterior(xx);
}

vector<TDV> HillClimb_NPSOL(TTimeSeries *TSModel, const TDV &start_point, unsigned int niteration, bool cold_start)
{
    const string COLD_START = string("Cold Start");
    const string WARM_START = string("Warm Start");
    const string NO_PRINT_OUT = string("Major print level = 0");
    const string MAJOR_PRINT_LEVEL = string("Major print level = 10");
    const string DERIVATIVE_LEVEL = string("Derivative level = 0");
    const string HESSIAN = string("Hessian = Yes"); 

    // npsol 
    MinusLogPosterior_NPSOL::model=TSModel; 
	
    int n = MinusLogPosterior_NPSOL::model->NumberParameters();

    // linear constraints
    int nclin = 0; 	// number of linear constraints
    int ldA = nclin > 1 ? nclin : 1; 	// number of rows of A 
    double A[ldA*n]; 
	
    // nonlinear constraints	
    int ncnln = 0; 	// number of nonlinear constraints	
    int ldJ = ncnln > 1 ? ncnln : 1; 	// number of rows of cJac; 
    double cJac[ldJ*n]; 
	
    // R provides the upper-triangular Cholesky factor of the initial approximation
    // of the Hessian of the Lagrangian
    int ldR = n; 	// number of rows of R;
    double R[ldR*n];   
	
    int nctotal = n + nclin + ncnln; 
    int istate[nctotal]; 	// array indicating whether the constaits are effective
    double bl[nctotal]; 	// lower bounds for samples, A and cJac
    double bu[nctotal]; 	// upper bounds for samples, A and cJac
    for (int i=0; i<n; i++)
      {
	bl[i] = MINUS_INFINITY; 
        bu[i] = PLUS_INFINITY; 
      }
    double clambda[nctotal];	// lagragin parameters of the constraints	
  
    for (int i=0; i<nctotal; i++)
      clambda[i] = 0; 

    // work space
    int leniw = 3*n + nclin + 2*ncnln; 
    int iw[leniw]; 		// integer work space
    int lenw; 
    if (nclin == 0 && ncnln == 0)
      lenw = 20*n; 
    else if (ncnln == 0)
      lenw = 2*n*n + 20*n + 11*nclin; 
    else 
      lenw = 2*n*n + n*nclin + 2*n*ncnln + 20*n + 11*nclin + 21*ncnln;
    double w[lenw]; 	// floating work space

    // returning value of npsol
    int inform, iter;
    double c[1];		// because ncnln = 0
    double f; 			// value of objective f(x) at the final iterate
    double g[n]; 		// objective gradient	   

    npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length());
    npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());
    npoptn_((char*)COLD_START.c_str(), COLD_START.length());
    npoptn_((char*)HESSIAN.c_str(), HESSIAN.length());

    vector<TDV> solution; 
    TDV maxima_point(n+2);
    TDV hessian_cholesky(n*n,0.0);     
    double LogPosterior_maxima=start_point(0);	
    f = -1.0e300;		
    TDV parameters(start_point(I(2,End)));
    unsigned int inpsol=0;

    while ((f+LogPosterior_maxima)<-1e-2 || inpsol<niteration) 
      {  
        npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, MinusLogPosterior_NPSOL::Maximize_Posterior, &inform, &iter, istate, c, cJac, clambda, &f, g, R, parameters.Vector(), iw, &leniw, w, &lenw);

        if (f < -LogPosterior_maxima)
	{
  	  LogPosterior_maxima = -f;
	  maxima_point(0)=TSModel->LogPosterior(parameters);
	  maxima_point(1)=TSModel->LogLikelihood(parameters);
	  maxima_point(I(2,End),parameters);
	  memcpy(hessian_cholesky.Vector(), R, sizeof(double)*n*n); 
	}
	
	inpsol++;
	//cout << "Iter " << inpsol <<  " Done! LogPosterior: " << -f << endl;
	if (cold_start)
	{
	  for (int i=0; i<nctotal; i++)
      	    clambda[i] = 0;   
	  for (int i=0; i<ldR*n; i++)
            R[i] = 0; 
	}
	else
	  npoptn_((char*)WARM_START.c_str(), WARM_START.length());
      }
	
    solution.push_back(maxima_point);	
    solution.push_back(hessian_cholesky);
    solution.push_back(TDV(1,(double)inpsol));

    return solution; 
}

void slave_mode_deploying(const int N_MESSAGE, TTimeSeries *model, const string &output_dir, bool cold_start)
{	
  int my_rank, nNode;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  double *sPackage= new double [N_MESSAGE];
  double *rPackage= new double [N_MESSAGE];   
  MPI_Status status;

  while (1)
    {
      MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
      if (status.MPI_TAG == HILL_CLIMB_TAG)
	{
	  unsigned int niteration = rPackage[NITERATION_INDEX]; 
	  unsigned int start_index = rPackage[GROUP_INDEX];
	  unsigned int length = rPackage[LENGTH_INDEX];
	  vector<TDV> all_samples, maxima_points;
	  TDV start_point;
 
 	  all_samples = slave_node_initial_points(output_dir, my_rank, start_index, length);
	  
	  for (unsigned int i=0; i<all_samples.size(); i++)
	    {
	      start_point=all_samples[i](I(1,all_samples[i].Dim()-1));
	      vector<TDV> solution = HillClimb_NPSOL(model, start_point, niteration, cold_start);
	      cout << "Node " << my_rank << " point " << i << " Done! LogPosterior: " << solution[0](0) << ", LogLikelihood: " << solution[0](1) << endl;
	      TDV one_point(all_samples[i].Dim());
	      one_point(I(0),solution[2]);
	      one_point(I(1,End),solution[0]);
	      maxima_points.push_back(one_point);
	    }
		
	  stringstream output_file;
      	  output_file.str("");
      	  output_file << output_dir << "/all_node_results/maxima_points_node" << my_rank << ".txt";
	  ofstream output;
  	  output.open(output_file.str().c_str());
  	  if (!output.is_open())
     	    throw dw_exception("slave_mode_deploying() - unable to open " + output_file.str());  
	  
	  for (unsigned int i=0; i<maxima_points.size(); i++)	  
	      output << maxima_points[i];

  	  output.close(); 
   	  
	}
      
      else
    	{
	  delete []rPackage; 
	  delete []sPackage; 
	  exit(0); 
	}  

      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 	  
   	    
    }
}

vector<TDV> HillClimb_NPSOL_iter(TTimeSeries *TSModel, const vector<TDV> &start_point)
{
	const string COLD_START = string("Cold Start");
        const string NO_PRINT_OUT = string("Major print level = 0");
	const string MAJOR_PRINT_OUT = string("Major print level = 10");
        const string DERIVATIVE_LEVEL = string("Derivative level = 0");
	const string HESSIAN = string("Hessian = Yes"); 

	// npsol 
	MinusLogPosterior_NPSOL::model=TSModel; 
	
	int n = MinusLogPosterior_NPSOL::model->NumberParameters();

	// linear constraints
	int nclin = 0; 	// number of linear constraints
	int ldA = nclin > 1 ? nclin : 1; 	// number of rows of A 
	double A[ldA*n]; 
	
	// nonlinear constraints	
	int ncnln = 0; 	// number of nonlinear constraints	
	int ldJ = ncnln > 1 ? ncnln : 1; 	// number of rows of cJac; 
	double cJac[ldJ*n]; 
	
	// R provides the upper-triangular Cholesky factor of the initial approximation
	// of the Hessian of the Lagrangian
	int ldR = n; 	// number of rows of R;
	double R[ldR*n];   
	
	int nctotal = n + nclin + ncnln; 
	int istate[nctotal]; 	// array indicating whether the constaits are effective
	double bl[nctotal]; 	// lower bounds for samples, A and cJac
	double bu[nctotal]; 	// upper bounds for samples, A and cJac
	for (int i=0; i<n; i++)
	{
	  bl[i] = MINUS_INFINITY; 
	  bu[i] = PLUS_INFINITY; 
	}
	double clambda[nctotal];	// lagragin parameters of the constraints	
	  
	for (int i=0; i<nctotal; i++)
	  clambda[i] = 0; 

	// work space
	int leniw = 3*n + nclin + 2*ncnln; 
	int iw[leniw]; 		// integer work space
	int lenw; 
	if (nclin == 0 && ncnln == 0)
	  lenw = 20*n; 
	else if (ncnln == 0)
	  lenw = 2*n*n + 20*n + 11*nclin; 
	else 
	  lenw = 2*n*n + n*nclin + 2*n*ncnln + 20*n + 11*nclin + 21*ncnln;
	double w[lenw]; 	// floating work space

	// returning value of npsol
	int inform, iter;
	double c[1];		// because ncnln = 0
	double f; 		// value of objective f(x) at the final iterate
	double g[n]; 		// objective gradient	   

	npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length());
        npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());
        npoptn_((char*)COLD_START.c_str(), COLD_START.length());
        npoptn_((char*)HESSIAN.c_str(), HESSIAN.length());

	int nSolution = start_point.size(); 

	vector<TDV> solutions; 
    
	for (int i=0; i<nSolution; i++)
	{
	  TDV parameters(start_point[i](I(2,End)));
	  npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, MinusLogPosterior_NPSOL::Maximize_Posterior, &inform, &iter, istate, c, cJac, clambda, &f, g, R, parameters.Vector(), iw, &leniw, w, &lenw);
	  TDV one_solution(n+2);
	  one_solution(0)=TSModel->LogPosterior(parameters);
	  one_solution(1)=TSModel->LogLikelihood(parameters);
	  one_solution(I(2,End),parameters);
	  solutions.push_back(one_solution);
    
	}
		
	return solutions; 

}

vector<TDV> slave_node_initial_points(const string &output_dir, int my_rank, unsigned int start_index, unsigned int length) 
{
  vector<TDV> all_samples;
  ifstream input;
  string line;
  vector<string> fields;
      
  stringstream samples_file;
  samples_file.str("");
  samples_file << output_dir << "/maxima_points_iter0.txt";
  input.open(samples_file.str().c_str());
  if (!input.is_open())
    throw dw_exception("slave_node_initial_points() - unable to open " + samples_file.str());

  input.seekg(0);
  if (my_rank > 1)
  {
    for (int i=0; i<start_index-1; i++)
    {
      getline(input,line);
      if (!input.good())
        throw dw_exception("slave_node_initial_points() - there is error in file " + samples_file.str());  
    }
  }

  for (int i=0; i<length; i++)
  {
    getline(input,line);
    if (input.good())
    {  
      ParseLine(line,' ',fields);
      TDV one_sample(fields.size()+1);
      one_sample(0)=1.0;
      for (unsigned int j=0; j<fields.size(); j++)
	one_sample(j+1) = atof(fields[j].c_str());
      all_samples.push_back(one_sample);
    }
    else if (!input.eof())
      throw dw_exception("slave_node_initial_points() - there is error in file " + samples_file.str());  
  }
  input.close();

  return all_samples;
}

vector<TDV> slave_node_start_points(const string &output_dir, int my_rank) 
{
  vector<TDV> all_samples;
  ifstream input;
  string line;
  vector<string> fields;
      
  stringstream samples_file;
  samples_file.str("");
  samples_file << output_dir << "/all_node_results/maxima_points_node" << my_rank << ".txt";
  input.open(samples_file.str().c_str());
  if (!input.is_open())
    throw dw_exception("slave_node_start_points() - unable to open " + samples_file.str());
  
  getline(input,line);
  while (input.good())
  {  
    ParseLine(line,' ',fields);
    TDV one_sample(fields.size());
    for (unsigned int j=0; j<fields.size(); j++)
	one_sample(j) = atof(fields[j].c_str());
    all_samples.push_back(one_sample);
    getline(input,line);
  } 
  if (!input.eof())
      throw dw_exception("slave_node_start_points() - there is error in file " + samples_file.str());  
  
  input.close();

  return all_samples;
}


void slave_mode_deploying_iter(const int N_MESSAGE, TTimeSeries *model, const string &output_dir)
{	
  int my_rank, nNode;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  double *sPackage= new double [N_MESSAGE];
  double *rPackage= new double [N_MESSAGE];   
  MPI_Status status;

  while (1)
    {
      MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
      if (status.MPI_TAG == HILL_CLIMB_TAG)
	{
	  unsigned int iteration = rPackage[ITERATION_INDEX]; 
	  unsigned int start_index = rPackage[GROUP_INDEX];
	  unsigned int length = rPackage[LENGTH_INDEX];
	  vector<TDV> all_samples, start_points;
	  TDV one_point;
 
 	  if (iteration == 0)   
	    all_samples = slave_node_initial_points(output_dir, my_rank, start_index, length);
	  else
	    all_samples = slave_node_start_points(output_dir, my_rank);
	  for (unsigned int i=0; i<all_samples.size(); i++)
	    {
	      if ((int)all_samples[i](0))
 	      {
	        one_point=all_samples[i](I(1,all_samples[i].Dim()-1));
		start_points.push_back(one_point);
	      }
	    }
  	    
  
	  vector<TDV> solutions = HillClimb_NPSOL_iter(model, start_points);

	  stringstream output_file;
      	  output_file.str("");
      	  output_file << output_dir << "/all_node_results/maxima_points_node" << my_rank << ".txt";
	  ofstream output;
  	  output.open(output_file.str().c_str());
  	  if (!output.is_open())
     	    throw dw_exception("slave_mode_deploying() - unable to open " + output_file.str());  
	  unsigned int j=0;	// npsol points index
	  unsigned int num_climbing=0;	// number of climbing points
	  for (unsigned int i=0; i<all_samples.size(); i++)
	  {
	    if ((int)all_samples[i](0))
	    {
	      if ((solutions[j](0)-start_points[j](0)) > 0.01)
	      {
		output << 1 << " " << solutions[j];
		num_climbing++;
	      }
	      else
		output << 0 << " " << start_points[j];	
	      j++;
	    }
	    else
	      output << all_samples[i];
	  }

  	  output.close(); 
   	  sPackage[NUMCLIMBING_INDEX]=num_climbing;
	}
      
      else
    	{
	  delete []rPackage; 
	  delete []sPackage; 
	  exit(0); 
	}  

      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 	  
   	    
    }
}
