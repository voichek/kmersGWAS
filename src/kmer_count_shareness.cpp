#include <sys/time.h>
#include "kmer_general.h"
#include <bitset>
#include <fstream>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>
#include <cstdlib>

#include <sparsehash/dense_hash_map>
//#include <sparsepp/sparsepp/spp.h>


#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#define HASH_TABLE_SIZE 1300000000
using namespace std;

//using spp::sparse_hash_map;
using google::dense_hash_map;
class CKmerAPI_YV: public CKmerAPI {
	public:
		CKmerAPI_YV (uint32 length = 0): CKmerAPI(length) {}
		uint64 to_uint() {
			return (uint64)kmer_data[0];
		}
		void infoYV() {
			cerr << "kmer_length = " <<  kmer_length << endl;				// Kmer's length, in symbols
			cerr << "byte_alignment = " <<  (int)byte_alignment << endl;			// A number of "empty" symbols placed before prefix to let sufix's symbols to start with a border of a byte
			cerr << "no_of_rows = " <<  no_of_rows << endl;				// A number of 64-bits words allocated for kmer_data 	
		}
};

struct Hash64 {
	size_t operator()(uint64_t key) const { 
		key ^= key >> 33;
		key *= 0xff51afd7ed558ccd;
		key ^= key >> 33;
		key *= 0xc4ceb9fe1a85ec53;
		key ^= key >> 33;
		return key;
	}
};

double get_time(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + (tv.tv_usec / 1000000.0);
}

typedef dense_hash_map<uint64, uint64, Hash64> my_hash; 
//typedef sparse_hash_map<uint64, uint64, Hash64> my_hash; 

vector<string> read_accession_db_list(char *filename) {
	ifstream fin(filename);
	vector<string> res;
	string name;
	while(fin >> name) {res.push_back(name);}
	return res;
}


//typedef sparse_hash_map<uint64, uint64 > my_hash; 
int main(int argc, char *argv[]) {
	cerr << "DB list file: " << argv[1] << endl;
	cerr << "going to run on the first: " << argv[2] << endl;
	uint64 DB_ro_run = atoi(argv[2]);

	vector<string> db_names = read_accession_db_list(argv[1]);
	CKmerAPI_YV	kmer_obj(KMER_LEN);
	uint c1;

	my_hash main_db(HASH_TABLE_SIZE);
	main_db.set_empty_key(-1);

	vector<uint64> k_mers(0);
	double t0,t1;

	// Going over K-mers DBs
	for(uint64 i=0; i<DB_ro_run; i++) {
		cerr << i << "\tloading: " << db_names[i] << endl;

		CKMCFile kmer_db;
		kmer_db.OpenForListing(db_names[i]);

		t0 = get_time();
		while (kmer_db.ReadNextKmer(kmer_obj, c1))  { 
			if((k_mers.size()%5000000) == 0) {cerr << "."; cerr.flush();}
			k_mers.push_back(kmer_obj.to_uint());
		}
		t1 = get_time();
		cerr << "\t\tLoaded " << k_mers.size() << " k-mers timing = " << (t1-t0) << endl;
		
		t0 = get_time();
		my_hash::iterator it_hash;
		for(uint k=0; k<k_mers.size(); k++) {
			if((k%5000000) == 0) {cerr << "|"; cerr.flush();}
			it_hash = main_db.find(k_mers[k]);
			if(it_hash == main_db.end()) 
				{ main_db.insert(my_hash::value_type(k_mers[k], 1));}
			else {
				it_hash->second++;			 
			}
		}
		t1 = get_time();
		cerr << "\t\tfor 100,000,000 time: " << (t1 - t0) << " size of hash table(M) " <<
			(main_db.size()/1000000) <<	endl;
		cout << main_db.size() << endl;
		k_mers.resize(0);	
	}
	vector<uint64> shareness(DB_ro_run	+ 1, 0);
	for(auto it : main_db) {
		shareness.at(it.second)++;
	}
	for(unsigned int i=0; i<shareness.size(); i++) {
		cout << i << "\t" << shareness[i] << endl;
	}
	return 0;
}
