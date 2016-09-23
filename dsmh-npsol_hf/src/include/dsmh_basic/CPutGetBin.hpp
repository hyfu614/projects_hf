#ifndef CLASS_PUT_GET_BIN_HEADER
#define CLASS_PUT_GET_BIN_HEADER

#include <string>
#include <vector>

class CSampleIDWeight; 

using namespace std; 

class CPutGetBin					 
{
private: 
	bool check_file_fetch_consolidate; 
	std::vector<string> filename_fetch; 
	std::vector<string> filename_consolidate; 
protected:
	int suffix;  
	string id; 
	int nDumpFile;			// total number of samples generated by far 
	int capacity; 			// capacity of this bin for put and get;
	int nPutUsed; 			// available space for put (will be reset as capacity after each dump, and will decrease as a sample is depoisited)
	int nGetUsed; 			// available data for get (will be reset as capacity after each fetch, and will decrease as a sample is drawn)	
	string filename_prefix; 		// run_id.id.indxx (index = nSampleGeneratedByFar/capacity); 
	vector <CSampleIDWeight> dataPut; 	// space for put
	vector <CSampleIDWeight> dataGet; 	// data for get 

	string GetFileNameForDump() const; 
	bool Dump(const string & =string()); 			// dump the current materials to file;
	int GetNumberFileForDump() const; 

	bool Fetch(); 
 	bool ReadFromOneFile(const string &, int &, const vector<int> &index);
	int NumberRecord(const string &) const; 
	void GetFileNameForFetchConsolidate();

public: 
	CPutGetBin(const string & _id=string(), int _nDumpFile=0, int _capacity=0, string _grandPrefix=string(), int _suffix=0); 
	CPutGetBin(const CPutGetBin &); 
	CPutGetBin & operator=(const CPutGetBin &); 
	~CPutGetBin() ; 

	void SetBinID(string _id, int _suffix=0) { id = _id; suffix=_suffix; }
	const string & GetBinID() const { return id; }

	void SetCapacity (int _capacity) { capacity = _capacity; }
	int GetCapacity() const { return capacity; }
	
	void SetFileNamePrefix(const string &_grandPrefix) { filename_prefix = _grandPrefix; } 
	string GetFileNamePrefix() const { return filename_prefix;}	

	int DepositSample(const CSampleIDWeight &); 
	bool DrawSample(CSampleIDWeight &); 

	void finalize(); 	// save unsaved data
	void consolidate(); 	// conslidate partial sample files into complete sample files
	void restore();	// load data from a partial file
	void RestoreForFetch(); // load data from a partial file but will not update it later. This is used for single-thread mpi version so that for each stage it will load partial files for its higher stage for ee draw later
	bool empty() const; 

	/* for reassigning samples into different bins */ 
	void DisregardHistorySamples(); 
	void ClearDepositDrawHistory(); 

	/* Reset check_file_fetch_consolidate etc */
	void ClearStatus() { 
		check_file_fetch_consolidate = false; 
		filename_fetch.clear(); 
		filename_consolidate.clear(); 
	}; 
}; 

#endif
