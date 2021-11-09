# A library for genome-wide associations using *k*-mers
A library for running *k*-mers-based GWAS

Binaries for beta version can be downloaded [from here](https://github.com/voichek/kmersGWAS/releases/download/v0.2-beta/v0_2_beta.zip).

Details on how to create a *k*-mers table and how to run GWA can be found in [the manual](/manual.pdf).

## key functionalities
#### Run pipeline:
+ **python2.7 src/py/pipeline.py** - Running the *k*-mers-GWAS pipeline
#### Build the *k*-mers table:
+ **kmers_add_strand_information** - Create *k*-mers list from KMC output
+ **list_kmers_found_in_multiple_samples** - Create a list of *k*-mers to use in the *k*-mers table
+ **build_kmers_table** - Build the *k*-mers table 
#### Other functionalities:
+ **kmers_table_to_bed** - Convert the *k*-mers table to PLINK binary format
+ **filter_kmers** - Textual output the presence/absence patterns of a set of *k*-mers
+ **emma_kinship_kmers** - Calculate the kinship matrix from PLINK binary format
+ **emma_kinship** - Calcualte the kinship matrix of the *k*-mers table

## Prerequisites:
+ Linux system with a 64 bit CPU
+ R
+ python**2.7**
+ KMC (part of the release under external_programs/ directory)
+ GEMMA (part of the release under external_programs/ directory)
+ R packages (if not present the pipeline will try to automatically install):
  + MASS
  + Mvnpermute
  + matrixcalc


## Examples:
In the examples directory there are two examples of how to use the library:
1. **Pre-existing _k_-mers table, only run the _k_-mers-base GWA** - Using the _k_-mers table for _A. thaliana_ 1001G to run _k_-mers-based GWA on flowering time (same as Fig. 1 in our paper). This table can be used also to run GWA on other phenotypes measured on accessions part of the 1001G.
2. **Building the k-mers table and using it to run GWAS** - Building the _k_-mers table on 241 _E. coli_ accessions (from Earle et al. 2016) and then running the GWA using this table.

## Citing
If you used our library for published work, please cite us:

[Identifying genetic variants underlying phenotypic variation in plants without complete genomes](https://www.nature.com/articles/s41588-020-0612-7)

Yoav Voichek and Detlef Weigel (2019).
