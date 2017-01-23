#include <cmath>
#include <fstream>
#include "dw_rand.h"
#include "CSampleIDWeight.hpp"
#include "CEquiEnergyModel.hpp"
#include "CMetropolis.hpp"
#include "DM.hpp"

using namespace std; 
using namespace DM;

bool CMetropolis::SetScale(const vector<double> &_scale)
{
	if (_scale.size() != blocks.size() )
		return false; 
	for (int i=0; i<(int)blocks.size(); i++)
		blocks[i] = blocks[i] * _scale[i];  
	return true; 
}

bool CMetropolis::SetScale(double _scale)
{
	for (int i=0; i<(int)blocks.size(); i++)
		blocks[i] = blocks[i] * _scale; 
	return true; 
}

void CMetropolis::SetBlocks(const vector<DM::TDM> &_blocks)
{
	blocks = _blocks; 
	for (int i=0; i<(int)blocks.size(); i++)
		//blocks[i].CopyContent(_blocks[i]);
                blocks[i]=_blocks[i];		// changed by HF 
		
}

void CMetropolis::SetBlockScheme(const vector<DM::I> &_index)
{
	block_scheme = _index;
	//for (int i=0; i<(int)block_scheme.size(); i++)
		//block_scheme[i](_index[i]);	// changed by HF
}

bool CMetropolis:: BlockRandomWalkMetropolis(double &log_posterior_y, CSampleIDWeight &y, const CSampleIDWeight &initial_v, int thin)
// Given
//	initial_v:	initial value, vector
// Results:
// 	y:	new vector according to random-walk metropolis
// 	log_posterior_y:	log_posterior(y)
// Returns:
// 	true:	if y is a new sample (different than x)
// 	false:	if y is the same as x
// Involving:
// 	blocks:	k matrices, each of which is n*b(i) consisting of the direction vectors of the i-th block
// 	model:	contains target distribution	
{
	int k=blocks.size(); 	// number of blocks	
	vector<int> b(k); 	// size of each block
        
	for (int i=0; i<k; i++)
		//b[i] = blocks[i].cols;
                b[i] = blocks[i].Cols();	// changed by HF
                 
	CSampleIDWeight x = initial_v;  
	double log_previous = model->log_posterior_function(x), log_current; 
	bool if_new_sample = false;  
	for (int i_thin=0; i_thin<thin; i_thin++)
	{
		for (int i=0; i<k; i++)
		{
			//y.data = x.data+ blocks[i]*RandomNormalVector(b[i]); 
			y.data = x.data + blocks[i]*RandomNormal(b[i]);		// changed by HF
			y.DataChanged(); 
			log_current = model->log_posterior_function(y); 
                        //cout << "Old Posterior=" << log_previous << " New Posterior=" << log_current << endl;
			if (log_current - log_previous >= log(dw_uniform_rnd()) )
			{
                                //cout << "new sample here!\n"; 
				x = y; 
				log_previous = log_current; 
				if (!if_new_sample)
					if_new_sample = true; 
			}
		}
	}
	
	log_posterior_y = log_previous; 
	y = x; 
	return if_new_sample; 
}

bool CMetropolis::ReadBlocks(const string &file_name)
{
	fstream iFile(file_name.c_str(), ios::in|ios::binary); 
	if (!iFile)
		return false; 
	int n_blocks; 
	iFile.read((char*)(&n_blocks),sizeof(int) ); 
	blocks = vector<DM::TDM>(n_blocks); 
	for (int i=0; i<n_blocks; i++)
	{
		int m, n; 
		iFile.read((char*)(&m),sizeof(int) ); 
		iFile.read((char*)(&n),sizeof(int) ); 
		blocks[i]=TDM(m,n); 
		for (int j=0; j<n; j++)
		{
			for (int k=0; k<m; k++)
			{
				double data; 
				iFile.read((char*)(&data),sizeof(double)); 
				blocks[i](k,j)=data; 
			}	
		}
	}

	iFile.close(); 
	return true; 	
}

bool CMetropolis::WriteBlocks(const string &file_name)
{
	fstream oFile(file_name.c_str(), ios::out|ios::binary); 
	if (!oFile)
		return false;
	int n_blocks = blocks.size(); 
	oFile.write((char*)(&n_blocks), sizeof(int));	// number of blocks
	for (int i=0; i<n_blocks; i++)
	{
		int m = blocks[i].Rows();	// blocks[i].rows; changed by HF 
		int n = blocks[i].Cols();	// blocks[i].cols; changed by HF
		oFile.write((char*)(&m), sizeof(int)); 
		oFile.write((char*)(&n), sizeof(int)); 
		for (int j=0; j<n; j++)
		{
			for (int k=0; k<m; k++)
			{
				double data = blocks[i](k,j); 
				oFile.write((char*)(&data), sizeof(double)); 
			}
		}
	}
	oFile.close(); 
	return true;
}

void CMetropolis::GetBlockMatrix_WeightedSampling(const std::vector<CSampleIDWeight> &Y, const std::vector<double> &weight)
{
	blocks.resize(block_scheme.size()); 
        DM::TDV mean, d_vector;
        DM::TDM variance, U_matrix, V_matrix, D_matrix;
        //cout << "Size of block scheme:" << block_scheme.size() << endl;
        for (int i_block=0; i_block<(int)block_scheme.size(); i_block++)
        {
                mean = Zeros(block_scheme[i_block].Size());	// Zeros(block_scheme[i_block].size); changed by HF
                //variance = Zeros(block_scheme[i_block].size, block_scheme[i_block].size);
		variance = Zeros(block_scheme[i_block].Size(), block_scheme[i_block].Size());	// changed by HF
                for (int i=0; i<(int)Y.size(); i++)
                {
                        //mean = mean + weight[i] * Y[i].data.SubVector(block_scheme[i_block]);
			//variance = variance + weight[i] * Multiply(Y[i].data.SubVector(block_scheme[i_block]),Y[i].data.SubVector(block_scheme[i_block]));
			//cout << "weight[i]: " << weight[i] << endl;
			mean = mean + weight[i] * Y[i].data(block_scheme[i_block]);	//changed by HF
			variance = variance + weight[i] * OuterProduct(Y[i].data(block_scheme[i_block]));	// changed by HF
                 }
                //cout << "mean: " << mean << endl;
                //variance = variance - Multiply(mean, mean);
                //variance = 0.5 * (variance+Transpose(variance));
		variance = variance - OuterProduct(mean);	// changed by HF
		variance = 0.5 * (variance+T(variance));	// changed by HF
		//cout << "Variance: " << variance << endl;
                SVD(U_matrix, d_vector, V_matrix, variance);
                //for (int i=0; i<d_vector.dim; i++)
		//	d_vector[i] = sqrt(d_vector[i]);
		//D_matrix = DiagonalMatrix(d_vector);
		for (unsigned int i=0; i<d_vector.Dim(); i++)
                        d_vector(i) = sqrt(d_vector(i));	// changed by HF
                D_matrix = Diag(d_vector);
                U_matrix = U_matrix * D_matrix;
                //cout << "U_matrix:" << U_matrix << endl;
                //blocks[i_block]=Zeros(Y[0].data.dim, block_scheme[i_block].size);
		//blocks[i_block].Insert(block_scheme[i_block], I(0,block_scheme[i_block].size-1), U_matrix);
		blocks[i_block]=Zeros(Y[0].data.Dim(), block_scheme[i_block].Size());	// changed by HF
                blocks[i_block](block_scheme[i_block], I(0,End), U_matrix);
       }
}

bool CMetropolis::ReadBlockScheme(const string &file_name)
{
	ifstream input; 
	input.open(file_name.c_str(), ios::in|ios::binary); 
	if (!input.is_open())
		return false; 
	
	int n_blocks; 
	input.read((char*)&(n_blocks), sizeof(int)); 
	block_scheme.resize(n_blocks);
	
	for (int i_block=0; i_block<n_blocks; i_block++)
	{
		int n_indices; 
		input.read((char*)&(n_indices), sizeof(int)); 
		int index; 
		for (int i_index=0; i_index<n_indices; i_index++)
		{
			input.read((char*)&(index), sizeof(int)); 
			//block_scheme[i_block] += index;
			block_scheme[i_block](index);	// changed by HF 
		}
	} 
	input.close(); 
	return true; 
}
