#include "kmer_general.h"

using namespace std;

int main(int argc, char *argv[]) {
	cerr << "DB list file: " << argv[1] << endl;
	cerr << "going to run on the first: " << argv[2] << endl;

	uint64 DB_ro_run = atoi(argv[2]); //the first * of accession to run on

	/* Read accessions to use list */
	vector<string> db_names = read_accession_db_list(argv[1]);
	CKmerAPI_YV	kmer_obj(KMER_LEN);


	/* Build the hash table that will contain all the kmers counts */
	my_hash main_db(HASH_TABLE_SIZE);
	main_db.set_empty_key(-1); // need to define empty value for google dense hash table

	vector<uint64> k_mers(0);
	uint kmer_counter;

	// Going over K-mers DBs
	for(uint64 i=0; i<DB_ro_run; i++) {
		cerr << i << "\tloading: " << db_names[i] << endl;

		CKMCFile kmer_db;
		kmer_db.OpenForListing(db_names[i]);

		while (kmer_db.ReadNextKmer(kmer_obj, kmer_counter))  { // Reading k-mers in file
			if((k_mers.size()%5000000) == 0) {cerr << "."; cerr.flush();}
			k_mers.push_back(kmer_obj.to_uint());
		}
		cerr << endl;

		my_hash::iterator it_hash;
		for(uint k=0; k<k_mers.size(); k++) {
			if((k%5000000) == 0) {cerr << "."; cerr.flush();}
			it_hash = main_db.find(k_mers[k]);
			if(it_hash == main_db.end()) 
				{ main_db.insert(my_hash::value_type(k_mers[k], 1));}
			else {
				it_hash->second++;			 
			}
		}
		cerr << endl;
		k_mers.resize(0);	
	}

	/* output all the k-mers apearing more than once to a file */
	ofstream file ("kmers_more_than_once.txt", ios::binary);

	vector<uint64> shareness(DB_ro_run+1, 0);
	for(auto it : main_db) {
		shareness[it.second]++;
		if(it.second>1) {
			file.write(reinterpret_cast<const char *>(&it.first), sizeof(it.first));
		}
	}
	
	file.close();
	// Output shareness
	for(unsigned int i=0; i<shareness.size(); i++) {
		cout << i << "\t" << shareness[i] << endl;
	}
	return 0;
}
