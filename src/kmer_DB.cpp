#include "kmer_DB.h"
#include <fstream>

using google::dense_hash_set;
using namespace std;


Kmer_DB::Kmer_DB(std::string dir_path, std::string db_name):
	m_db_name(db_name),
	m_dir_path(dir_path) {}


bool lookup_x(const kmer_set& Set, const uint64& kmer)
{
	kmer_set::const_iterator it  = Set.find(kmer);
	return (it != Set.end()); 
}


// go over the KMC DB and output to file only the kmers in the kmers set inputed
void Kmer_DB::intersect_kmers(const kmer_set& kmers_to_use, std::string file_name) {
	string output_file = m_dir_path + "/" +  file_name;
	string output_log_file = output_file + ".log";

	ofstream of(output_file, ios::binary);
	ofstream of_log(output_log_file);

	cerr << output_file << endl;
	CKMCFile kmer_database = get_KMC_handle(); //open DB
	CKmerAPI_YV kmer_obj(KMER_LEN);

	double t0,t1;
	t0 = get_time();
	unsigned int counter;
	uint64 cnt_f = 0, cnt_nf = 0, kmer;

	vector<uint64> kmers;
	while(kmer_database.ReadNextKmer(kmer_obj, counter)) {
		kmer = kmer_obj.to_uint();
		if(lookup_x(kmers_to_use, kmer)) {
			cnt_f++;
			kmers.push_back(kmer);
		} else { cnt_nf++;}
	}
	// Sort kmers 
	sort(kmers.begin(), kmers.end());
	for(vector<uint64>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
		of.write(reinterpret_cast<const char *>(&(*it)), sizeof(*it));
	}
	t1 = get_time();
	of_log << m_db_name << "\ttime\t" << t1-t0 << "\tfound\t" << cnt_f << "\tnot_found\t" << cnt_nf << endl;
	of.close();
	of_log.close();
}


// Counts how many times each k-mer appeared and plot to std::cout
vector<std::size_t> Kmer_DB::calculate_kmers_counts_histogram() {
	CKMCFile kmer_database = get_KMC_handle();
	CKmerAPI kmer_obj(KMER_LEN);

	vector<std::size_t> counters(0, 0);
	unsigned int counter;
	while(kmer_database.ReadNextKmer(kmer_obj, counter)) {
		if(counter >= counters.size()) { counters.resize(counter+1,0);}
		counters[counter] += 1;		
	}
	return counters;
}

// get a hanle to the KMC DB to run over its k-mers 
CKMCFile Kmer_DB::get_KMC_handle() {
	CKMCFile kmer_database;
	kmer_database.OpenForListing(m_dir_path + "/" +  m_db_name);

	return kmer_database;
}


//int main(int argc, char *argv[]) {
//	cerr << "DB list file: " << argv[1] << endl;
//	cerr << "going to run on the first: " << argv[2] << endl;
//
//	uint64 DB_to_run = atoi(argv[2]); //the first * of accession to run on
//
//	/* Read accessions to use list */
//	vector<string> db_names = read_accession_db_list(argv[1]);
//	string dir_dbs(argv[3]);	
//
//	kmer_set kmer_list_to_use(800000000);
//	kmer_list_to_use.set_empty_key(-1); // need to define empty value for google dense hash table
//
//	double t0,t1;
//
//	ifstream kmer_file("kmers_more_than_once", std::ifstream::binary);
//	if(kmer_file) { // if file could be open
//		kmer_file.seekg(0, kmer_file.end); // 
//		uint64 kmers_in_file = (kmer_file.tellg()) >> 3;
//		kmer_file.seekg(0, kmer_file.beg);
//		uint64 kmer_uint;
//		t0 = get_time();
//		for(uint64 i=0; i<(kmers_in_file); i++) {
//			kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
//			kmer_list_to_use.insert(kmer_uint);
//		}
//		kmer_file.close();
//		t1 = get_time();
//		std::cerr << "time to load: " << t1-t0 << std::endl;
//		cerr << "Hash set size = " << kmer_list_to_use.size() << endl;
//	} else { cerr << "could not open file?" << endl;}
//	for(uint64 i=0; i<DB_to_run; i++) {
//		Kmer_DB cur_acc(dir_dbs + "/" + db_names[i] , db_names[i]);
//		cur_acc.intersect_kmers(kmer_list_to_use, "order_kmers_appear_more_than_once");
//	}
//	return 0;
//}
