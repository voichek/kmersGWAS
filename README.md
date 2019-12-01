# A library for genome-wide association using k-mers
A library for running k-mers based GWAS

Binaries for beta version can be downloaded [from here](https://github.com/voichek/kmersGWAS/releases/download/v0.1-beta/v0_1_beta.zip).

Details on how to create a k-mers table and how to run GWA can be found in [the manual](/manual.pdf).

## key functionalities
#### Run pipeline:
+ **python2.7 src/py/pipeline.py** - Running the k-mers-GWAS pipeline
#### Build k-mers table:
+ **kmers_add_strand_information** - Create k-mers list from KMC output
+ **list_kmers_found_in_multiple_samples** - Create a list of k-mers to use in the k-mers table
+ **build_kmers_table** - Build the k-mers table 
#### other:
+ **kmers_table_to_bed** - Convert the k-mers table to PLINK binary format
+ **filter_kmers** - Output the presence/absence patterns of a set of k-mers
+ **emma_kinship_kmers** - Calculate the kinship matrix from PLINK binary format
+ **emma_kinship** - Calcualte the kinship matrix of the k-mers table


## Citing
If you used our library for published work, please cite us:

[Finding genetic variants in plants without complete genomes](https://www.biorxiv.org/content/10.1101/818096v2)

Yoav Voichek and Detlef Weigel (2019).
