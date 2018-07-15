#include "kmer_general.h"

using namespace std;

int main(int argc, char *argv[]) {
	CKMCFile kmer_database;
	kmer_database.OpenForListing(argv[1]);
	CKmerAPI kmer_obj(KMER_LEN);
	
	vector<unsigned int> counters(0, 0);
	unsigned int counter, cnt;
	string str;
	cnt = 0;
	while(kmer_database.ReadNextKmer(kmer_obj, counter)) {
		if(counter >= counters.size()) { counters.resize(counter+1,0);}
		counters[counter] += 1;		
		cnt++;
	}
	
	for(unsigned int i=0; i<counters.size(); i++) {
		cout << i << "\t" << counters[i] << endl;											
	}

	return 0;
}
