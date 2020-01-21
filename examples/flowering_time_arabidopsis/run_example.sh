#!/bin/bash

# Run k-mers GWAS on flowering time in A. thaliana using the pre-built k-mers table

#1. Download the kmers_table matrix for the 1001G
wget https://zenodo.org/record/3596255/files/A_thaliana_1001G_kmers_table.tar.gz?download=1 -O 1001G_kmers_table.tar.gz

#2. Extract the kmers_table from the .tar.gz
tar -xf 1001G_kmers_table.tar.gz

#3. Run the kmers GWAS (use 8 threads)
python2.7 ../../kmers_gwas.py --pheno FT10.pheno --kmers_table A_thaliana_1001G_kmers_table/kmers_table -l 31 -p 8 --outdir run_GWAS_FT10
