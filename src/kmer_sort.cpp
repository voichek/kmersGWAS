#include <sys/time.h>
#include "kmer_general.h"
#include <bitset>
#include <fstream>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>
#include <cstdlib>

using namespace std;

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

double get_time(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (tv.tv_usec / 1000000.0);
}


vector<string> read_accession_db_list(char *filename) {
	ifstream fin(filename);
	vector<string> res;
	string name;
	while(fin >> name) {res.push_back(name);}
	return res;
}

int main(int argc, char *argv[]) {
	cerr << "reading accession DB list" << endl;
	vector<string> db_names = read_accession_db_list(argv[1]);

	CKmerAPI_YV	kmer_obj(KMER_LEN);
	uint c1;

	vector<uint64> temp;
	double t0,t1;
	for(uint64 i=0; i<(uint64)atoi(argv[2]); i++) {
		cerr << i << "\tloading: " << db_names[i] << endl;

		CKMCFile kmer_db;
		kmer_db.OpenForListing(db_names[i]);

		uint64 j =0;
		t0 = get_time();
		while (kmer_db.ReadNextKmer(kmer_obj, c1)) { 
			if((j%5000000) == 0) {cerr << "."; cerr.flush();}
			temp.push_back(kmer_obj.to_uint());
			j++;
		}
		t1 = get_time();
		cerr << "\t\tvector size = " << temp.size() << " timing = " << (t1-t0) << endl;
		t0 = get_time();
		std::sort(temp.begin(), temp.end());
		t1 = get_time();
		cerr << "sorting " << temp.size() << " elements, took: " << (t1-t0) << endl;
		//uint64 cnt_uniq = 0;
		uint64 last = temp[0];
		uint64 i1,i2;
		i1 = 1;
		i2 = 1;
		while(i2<temp.size()){
			if(temp[i2] != last) {
				last = temp[i2];
				temp[i1] = last;
				i1++;
			}
			i2++;
		}
		temp.resize(i1);
		//for(uint64 i = 1; i<temp.size(); i++) {
		//	if(temp[i] !=temp[i-1]) {cnt_uniq++;}
		//}	
		cout << temp.size() << "\t" << i2 << endl;
	}
	return 0;
}
