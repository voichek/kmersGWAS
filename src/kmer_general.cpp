#include "kmer_general.h"

using namespace std;

vector<KMC_db_handle> read_accession_db_list(string filename) {
	ifstream fin(filename);
	vector<KMC_db_handle> res;
	KMC_db_handle db_info;

	while(fin >> db_info.dir_path) {
		fin >> db_info.name;
		res.push_back(db_info);
	}
	fin.close();
	return res;
}

string bits2kmer31(uint64 w) {
	const static char dict_bp[] = {'A','C','G','T'};
	const static uint64 mask2bits = 0x0000000000000003;

	string res(31,'X');
	for(std::size_t i=0; i<31; i++) {
		res[30-i] = dict_bp[w & mask2bits];
		w = (w>>2);
	}
	return res;
}

double get_time(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + (tv.tv_usec / 1000000.0);
}

kmer_set load_kmer_raw_file(std::string filename, size_t set_initial_size) {
	kmer_set kmer_list_to_use(set_initial_size);
	kmer_list_to_use.set_empty_key(NULL_KEY); // need to define empty value for google dense hash table

	ifstream kmer_file(filename, std::ifstream::binary);
	if(kmer_file) { // if file could be open
		kmer_file.seekg(0, kmer_file.end); // 
		uint64 kmers_in_file = (kmer_file.tellg()) >> 3;
		kmer_file.seekg(0, kmer_file.beg);
		uint64 kmer_uint;
		for(uint64 i=0; i<(kmers_in_file); i++) {
			kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
			kmer_list_to_use.insert(kmer_uint);
		}
		kmer_file.close();
	}
	return kmer_list_to_use;
}

kmer_set load_kmer_and_score_raw_file(std::string filename, size_t set_initial_size) {
	kmer_set kmer_list_to_use(set_initial_size);
	kmer_list_to_use.set_empty_key(NULL_KEY); // need to define empty value for google dense hash table

	ifstream kmer_file(filename, std::ifstream::binary);
	if(kmer_file) { // if file could be open
		kmer_file.seekg(0, kmer_file.end); // 
		uint64 kmers_in_file = (kmer_file.tellg()) >> (3+1);
		kmer_file.seekg(0, kmer_file.beg);
		uint64 kmer_uint;
		double score;
		for(uint64 i=0; i<(kmers_in_file); i++) {
			kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
			kmer_file.read(reinterpret_cast<char *>(&score), sizeof(score));
			kmer_list_to_use.insert(kmer_uint);
		}
		kmer_file.close();
	}
	std::cerr << "loaded set of k-mers #" << kmer_list_to_use.size() << std::endl;
	return kmer_list_to_use;
}
