#include "kmer_DB.h"
using namespace std;

int main(int argc, char *argv[]) {
	cerr << "DB list file: " << argv[1] << endl;
	cerr << "going to run on the first: " << argv[2] << endl;

	uint64 DB_to_run = atoi(argv[2]); //the first * of accession to run on

	/* Read accessions to use list */
	vector<string> db_names = read_accession_db_list(argv[1]);
	string dir_dbs(argv[3]);	

	kmer_set kmer_list_to_use(800000000);
	kmer_list_to_use.set_empty_key(-1); // need to define empty value for google dense hash table

	double t0,t1;

	ifstream kmer_file("kmers_more_than_once", std::ifstream::binary);
	if(kmer_file) { // if file could be open
		kmer_file.seekg(0, kmer_file.end); // 
		uint64 kmers_in_file = (kmer_file.tellg()) >> 3;
		kmer_file.seekg(0, kmer_file.beg);
		uint64 kmer_uint;
		t0 = get_time();
		for(uint64 i=0; i<(kmers_in_file); i++) {
			kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
			kmer_list_to_use.insert(kmer_uint);
		}
		kmer_file.close();
		t1 = get_time();
		std::cerr << "time to load: " << t1-t0 << std::endl;
		cerr << "Hash set size = " << kmer_list_to_use.size() << endl;
	} else { cerr << "could not open file?" << endl;}
	for(uint64 i=0; i<DB_to_run; i++) {
		Kmer_DB cur_acc(dir_dbs + "/" + db_names[i] , db_names[i]);
		cur_acc.intersect_kmers(kmer_list_to_use, "order_kmers_appear_more_than_once");
	}
	return 0;
}
