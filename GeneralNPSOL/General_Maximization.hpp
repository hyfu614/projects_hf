#ifndef __GENERAL_MAXIMIZATION__
#define __GENERAL_MAXIMIZATION__

#include <string>
#include <vector>
#include "DM.hpp"
#include "TTimeSeries.hpp"

class MinusLogPosterior_NPSOL : public TS::TTimeSeries
{
public:
  static TS::TTimeSeries *model;
  MinusLogPosterior_NPSOL(const MinusLogPosterior_NPSOL &Model) : TTimeSeries(Model) {};
  static void Maximize_Posterior(int *mode, int *n, double *x, double *f, double *g, int *nstate);  
};

class TDV_Sort
{
private:
  int index;
public:
  TDV_Sort(int Index) : index(Index) {};
  bool operator() (const DM::TDV &previous, const DM::TDV &next) { return previous(index) > next(index); };
};

extern "C"
{
  // npsol defines
  void npoptn_(char *string, unsigned int length);
  void npsol_(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A, double *bl, double *bu,
	      void funcon(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate),
	      void funobj(int *mode, int *n, double *x, double *f, double *g, int *nstate), int *inform, int *iter, int *istate,
	      double *c, double *cJac, double *clamda, double *f, double *g, double *R, double *x, int *iw, int *leniw, double *w, int *lenw);

}

void master_mode_deploying(const int N_MESSAGE, int nNode, unsigned int nMaxima, int nIterations, TS::TTimeSeries *model, const std::string &draws_file, int max_interval, const std::string &output_dir);
void output_to_file_final(const std::string &output_dir, unsigned int nNode, unsigned int nMaxima);
void master_mode_deploying_iter(const int N_MESSAGE, int nNode, unsigned int nMaxima, int nIterations, TS::TTimeSeries *model, const std::string &draws_file, int interval, const std::string &output_dir);
void output_to_file(const std::string &output_dir, unsigned int iIteration, unsigned int nNode, unsigned int nMaxima);
void slave_mode_deploying(const int N_MESSAGE, TS::TTimeSeries *model, const std::string &output_dir, bool cold_start);
std::vector<DM::TDV> HillClimb_NPSOL(TS::TTimeSeries *TSModel, const DM::TDV &start_point, unsigned int niteration, bool cold_start);
void slave_mode_deploying_iter(const int N_MESSAGE, TS::TTimeSeries *model, const std::string &output_dir);
std::vector<DM::TDV> HillClimb_NPSOL_iter(TS::TTimeSeries *TSModel, const std::vector<DM::TDV> &start_point);
std::vector<DM::TDV> slave_node_initial_points(const std::string &output_dir, int my_rank, unsigned int start_index, unsigned int length);
std::vector<DM::TDV> slave_node_start_points(const std::string &output_dir, int my_rank); 
#endif
