#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>

#include "specification_io.hpp"
#include "mpi_constant.hpp"
#include "General_Maximization.hpp"
#include "dw_rand.h"

using namespace std;
using namespace DM;
using namespace TS;

void master_mode_deploying(const int N_MESSAGE, int nNode, unsigned int nMaxima, int nIterations, TTimeSeries *model, const string &draws_file, int interval, const string &output_dir)
{	
    double *sPackage= new double [N_MESSAGE];
    double *rPackage= new double [N_MESSAGE];   
    MPI_Status status;

    vector<TDV> samples;

    if (draws_file.empty())	// no draws file, drawing prior
    {
      for (unsigned int i=0; i < 10*nMaxima; i++)
        {
	  TDV parameters = model->DrawPrior();
	  TDV one_sample(parameters.Dim()+2);
	  one_sample(0) = model->LogPosterior(parameters);
	  one_sample(1) = model->LogLikelihood(parameters);
	  one_sample(I(2,End), parameters);
	  samples.push_back(one_sample);
	}
    }
    else			// read random samples from draws file 
    {
      ifstream input;
      string line;
      vector<string> fields;
      
      input.open(draws_file.c_str());
      if (!input.is_open())
	throw dw_exception("master_mode_deploying() - unable to open " + draws_file);
      input.seekg(0);
      for (unsigned int i=0; i < nMaxima; i++)
	{
	  getline(input,line);
	  if (input.good())
     	    {
	      ParseLine(line,'\t',fields);
              TDV one_sample(fields.size());
              for (unsigned int j=0; j<fields.size(); j++)
		one_sample(j) = atof(fields[j].c_str());
      	      samples.push_back(one_sample);
	    }
	  unsigned int inc = (i == 0) ? interval+dw_uniform_int(interval) : interval;
	  for (unsigned int j=0; j<inc; j++)
  	    {
	      getline(input,line);
              if (input.eof())
		input.seekg(0);
	    }
	}
      input.close();
    }
    // sort the samples in descending LogPosterior order
    TDV_Sort comparator(0);
    sort(samples.begin(), samples.end(), comparator);
      
    stringstream start_point_filename;
    start_point_filename.str("");
    start_point_filename << output_dir << "/maxima_points_iter0.txt";
    ofstream output;
    output.open(start_point_filename.str().c_str());
    if (!output.is_open())
      throw dw_exception("master_mode_deploying() - unable to open " + start_point_filename.str());
    for (unsigned int i=0; i<nMaxima; i++)
      output << samples[i];
    output.close();	
	
    unsigned int nsamples_per_node = floor((double)nMaxima/(double)(nNode-1));

    // send out hill_climb jobs
    for (int j=1; j<nNode; j++)
      {
          unsigned int start_index = (j-1)*nsamples_per_node;
          unsigned int end_index = j*nsamples_per_node;
          end_index = end_index < nMaxima ? end_index : nMaxima; 
  	  sPackage[NITERATION_INDEX] = nIterations; 
          sPackage[GROUP_INDEX] = start_index; 
          sPackage[LENGTH_INDEX] = end_index-start_index; 
          MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, j, HILL_CLIMB_TAG, MPI_COMM_WORLD);
      }

    for (int j=1; j<nNode; j++)
      MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);

    output_to_file_final(output_dir, nNode, nMaxima);
      
    // tell all the slaves to exit by sending an empty messag with 0 simulation length 
    for (int i=1; i<nNode; i++)
      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
    delete [] sPackage;
    delete [] rPackage;

}

void output_to_file_final(const string &output_dir, unsigned int nNode, unsigned int nMaxima)
{
    vector<TDV> allpoints(0);
    for (unsigned int i=1; i<nNode; i++)
    {
      stringstream input_filename;
      input_filename.str("");
      input_filename << output_dir << "/all_node_results/maxima_points_node" << i << ".txt";
      ifstream input;
      input.open(input_filename.str().c_str());
      if (!input.is_open())
        throw dw_exception("output_to_file() - unable to open " + input_filename.str());   

      vector<TDV> points;
      string line;
      vector<string> fields;
      getline(input,line);
      while (input.good() && !input.eof())
	{
	  ParseLine(line,' ',fields);
          TDV one_point(fields.size()-1);
          for (unsigned int j=1; j<fields.size(); j++)
	    one_point(j-1) = atof(fields[j].c_str());
      	  points.push_back(one_point);	
	  getline(input,line);
        }
      if (!input.good() && !input.eof())
        throw dw_exception("output_to_file() - there is error in file " + input_filename.str());  
      if (input.eof()) 
	input.close();
      allpoints.insert(allpoints.end(), points.begin(), points.end());
    }

    stringstream output_filename;
    output_filename.str("");
    output_filename << output_dir << "/maxima_points_final.txt";
    ofstream output;
    output.open(output_filename.str().c_str());
    if (!output.is_open())
      throw dw_exception("output_to_file() - unable to open " + output_filename.str());  
    for (unsigned int i=0; i<allpoints.size(); i++)
      output << allpoints[i];
    output.close();

    TDV_Sort comparator(0);
    sort(allpoints.begin(), allpoints.end(), comparator);
    output_filename.str("");
    output_filename << output_dir << "/sorted_maxima_points_final.txt";
    output.open(output_filename.str().c_str());
    if (!output.is_open())
      throw dw_exception("output_to_file() - unable to open " + output_filename.str());  
    for (unsigned int i=0; i<allpoints.size(); i++)
      output << allpoints[i];
    output.close();
}


void master_mode_deploying_iter(const int N_MESSAGE, int nNode, unsigned int nMaxima, int nIterations, TTimeSeries *model, const string &draws_file, int interval, const string &output_dir)
{	
    double *sPackage= new double [N_MESSAGE];
    double *rPackage= new double [N_MESSAGE];   
    MPI_Status status;

    vector<TDV> samples;

    if (draws_file.empty())	// no draws file, drawing prior
    {
      for (unsigned int i=0; i < 100*nMaxima; i++)
        {
	  TDV parameters = model->DrawPrior();
	  TDV one_sample(parameters.Dim()+2);
	  one_sample(0) = model->LogPosterior(parameters);
	  one_sample(1) = model->LogLikelihood(parameters);
	  one_sample(I(2,End), parameters);
	  samples.push_back(one_sample);
	}
    }
    else			// read random samples from draws file 
    {
      ifstream input;
      string line;
      vector<string> fields;
      
      input.open(draws_file.c_str());
      if (!input.is_open())
	throw dw_exception("master_mode_deploying() - unable to open " + draws_file);
      input.seekg(0);
      for (unsigned int i=0; i < nMaxima; i++)
	{
	  getline(input,line);
	  if (input.good())
     	    {
	      ParseLine(line,'\t',fields);
              TDV one_sample(fields.size());
              for (unsigned int j=0; j<fields.size(); j++)
		one_sample(j) = atof(fields[j].c_str());
      	      samples.push_back(one_sample);
	    }
	  unsigned int inc = (i == 0) ? interval+dw_uniform_int(interval) : interval;
	  for (unsigned int j=0; j<inc; j++)
  	    {
	      getline(input,line);
              if (input.eof())
		input.seekg(0);
	    }
	}
      input.close();
    }
    // sort the samples in descending LogPosterior order
    TDV_Sort comparator(0);
    sort(samples.begin(), samples.end(), comparator);
      
    stringstream start_point_filename;
    start_point_filename.str("");
    start_point_filename << output_dir << "/maxima_points_iter0.txt";
    ofstream output;
    output.open(start_point_filename.str().c_str());
    if (!output.is_open())
      throw dw_exception("master_mode_deploying() - unable to open " + start_point_filename.str());
    for (unsigned int i=0; i<nMaxima; i++)
      output << samples[i];
    output.close();	
	
    unsigned int nsamples_per_node = floor((double)nMaxima/(double)(nNode-1));
    unsigned int total_climbing_points = nMaxima;

    // send out hill_climb jobs
    //for (int i=0; i<nIterations; i++)
    int i = 0;		// iteration index
    while (total_climbing_points > 0)
    {
      for (int j=1; j<nNode; j++)
        {
          unsigned int start_index = (j-1)*nsamples_per_node;
          unsigned int end_index = j*nsamples_per_node;
          end_index = end_index < nMaxima ? end_index : nMaxima; 
  	  sPackage[ITERATION_INDEX] = i; 
          sPackage[GROUP_INDEX] = start_index; 
          sPackage[LENGTH_INDEX] = end_index-start_index; 
          MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, j, HILL_CLIMB_TAG, MPI_COMM_WORLD);
        }

      total_climbing_points = 0;
      for (int j=1; j<nNode; j++)
      {
        MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);
	total_climbing_points += rPackage[NUMCLIMBING_INDEX];
      }

      output_to_file(output_dir, i+1, nNode, nMaxima);
      cout << "Iter " << i+1 << ": Number of climbing points is " << total_climbing_points << " (out of " << nMaxima << ").\n";
      i++;
    }

    // tell all the slaves to exit by sending an empty messag with 0 simulation length 
    for (int i=1; i<nNode; i++)
      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
    delete [] sPackage;
    delete [] rPackage;

}

void output_to_file(const string &output_dir, unsigned int iIteration, unsigned int nNode, unsigned int nMaxima)
{
    vector<TDV> allpoints(0);
    for (unsigned int i=1; i<nNode; i++)
    {
      stringstream input_filename;
      input_filename.str("");
      input_filename << output_dir << "/all_node_results/maxima_points_node" << i << ".txt";
      ifstream input;
      input.open(input_filename.str().c_str());
      if (!input.is_open())
        throw dw_exception("output_to_file() - unable to open " + input_filename.str());   

      vector<TDV> points;
      string line;
      vector<string> fields;
      getline(input,line);
      while (input.good() && !input.eof())
	{
	  ParseLine(line,' ',fields);
          TDV one_point(fields.size()-1);
          for (unsigned int j=1; j<fields.size(); j++)
	    one_point(j-1) = atof(fields[j].c_str());
      	  points.push_back(one_point);	
	  getline(input,line);
        }
      if (!input.good() && !input.eof())
        throw dw_exception("output_to_file() - there is error in file " + input_filename.str());  
      if (input.eof()) 
	input.close();
      allpoints.insert(allpoints.end(), points.begin(), points.end());
    }

    stringstream output_filename;
    output_filename.str("");
    output_filename << output_dir << "/maxima_points_iter" << iIteration << ".txt";
    ofstream output;
    output.open(output_filename.str().c_str());
    if (!output.is_open())
      throw dw_exception("output_to_file() - unable to open " + output_filename.str());  
    for (unsigned int i=0; i<allpoints.size(); i++)
      output << allpoints[i];
    output.close();

    TDV_Sort comparator(0);
    sort(allpoints.begin(), allpoints.end(), comparator);
    output_filename.str("");
    output_filename << output_dir << "/sorted_maxima_points_iter" << iIteration << ".txt";
    output.open(output_filename.str().c_str());
    if (!output.is_open())
      throw dw_exception("output_to_file() - unable to open " + output_filename.str());  
    for (unsigned int i=0; i<allpoints.size(); i++)
      output << allpoints[i];
    output.close();
}



