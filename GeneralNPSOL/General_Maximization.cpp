#include <fstream>
#include <sstream>
#include <mpi.h>
#include <getopt.h>
#include <sys/stat.h>

#include "IRM.hpp"
#include "General_Maximization.hpp"
#include "mpi_constant.hpp"

using namespace std;
using namespace DM;
using namespace TS;

bool makedir(const string &results_dir, const string &sub_dir)
{
  errno = 0;   // clear error number
  int status = mkdir(results_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (status !=0 && errno != EEXIST) 
    { 
      cerr << results_dir << endl;
      return false;
    }
 
  stringstream dir_name; 
  dir_name.str(""); 
  dir_name << results_dir << "/" << sub_dir;
  errno = 0;   // clear error number
  status = mkdir(dir_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (status !=0 && errno != EEXIST) 
    { 
      cerr << dir_name << endl;
      return false;
    } 

  dir_name << "/all_node_results";
  errno = 0;   // clear error number
  status = mkdir(dir_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (status !=0 && errno != EEXIST) 
    { 
      cerr << dir_name << endl;
      return false;
    } 
  return true;
}


int main(int argc, char **argv)
{
  try
    {
      static struct option long_options[] =
        {
		// Simulation options
		{"output sub_directory", required_argument, 0, 'R'},
		{"prior draws", no_argument, 0, 'p'},
		{"file for starting point", required_argument, 0, 'F'},		 
		{"number of points in draws file", required_argument, 0, 'M'},
		{"number of optimization points", required_argument, 0, 'N'},
                {"number of NPSOL calls for optimization", required_argument, 0,  'O'},
		{"warm start for NPSOL call", no_argument, 0, 'w'},
		{0, 0, 0, 0}
	}; 
	
      string specification_file = "IRM_6var_4shock_restricted_specification_generalf0_deg21_np.txt";
      //string specification_file = "IRM_6var_4shock_unrestricted_specification_generalf0_deg21.txt";
      string draws_file = "./results/general_4shock_deg21_unr/general_4shock_deg21_unr.draw.text";
      string results_dir = "./maxresults"; 	// default directory for saving results
      string output_subdir = "posterior_draws";	// default sub_directory for saving results
      bool cold_start = true;		// default: cold start for each NPSOL call

      unsigned int nDraws = 200000;
      unsigned int nMaxima = 20;
      unsigned int nIterations = 100;

      // Command line option processing
      int option_index = 0;	
      while (1)
        {
          int c = getopt_long(argc, argv, "R:pF:M:N:O:w", long_options, &option_index);
          if (c == -1)
            break;
	  switch(c)
            {
	      case 'R':
		output_subdir = string(optarg); break;
	      case 'p':
		draws_file = string(""); break; 
	      case 'F':
                draws_file = string(optarg); break;
	      case 'M':
		nDraws = atoi(optarg); break;
	      case 'N':
                nMaxima = atoi(optarg); break;
              case 'O':
                nIterations = atoi(optarg); break;
	      case 'w':
		cold_start = false; break;
              default: 
		break; 
	    }
	}

      TTimeSeries_IRM irm=TTimeSeries_IRM_Specification(specification_file);    	
      TTimeSeries *target_model=irm.Clone();   
      unsigned int interval=(unsigned int)floor(nDraws/nMaxima);

      MPI_Init(&argc, &argv);
      int my_rank, nNode;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &nNode);
	
      dw_initialize_generator(time(NULL));
      
      stringstream output_dir; 
      output_dir.str(""); 
      output_dir << results_dir << "/" << output_subdir;
  
      const int N_MESSAGE = 100; 
      if (my_rank == 0)
        {
          if (!makedir(results_dir, output_subdir))
            {
              cerr << "Error in making directory for " << output_dir << endl;
              double *sMessage= new double [N_MESSAGE];
              for (int i=1; i<nNode; i++)
                MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
              delete [] sMessage;
              exit(1);
            }
	  master_mode_deploying(N_MESSAGE, nNode, nMaxima, nIterations, target_model, draws_file, interval, output_dir.str());
        }
      else
	slave_mode_deploying(N_MESSAGE, target_model, output_dir.str(), cold_start);

      MPI_Finalize();
    }
  catch (dw_exception &e)
    {
      cout << "exception: " << e.what() << endl;
    }
  return 0;
}

      
