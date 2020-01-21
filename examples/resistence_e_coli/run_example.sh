#! /bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This is an example of how to use the k-mers GWAS library for building the initial k-mers table and then use it to conduct GWAS
# Please change the paths to be correct in your system. you can also change the parameters used to to run the different parts of the procedure.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library_dir="../../" # The path to the directory with the library
base_dir="./example_dir" # This will be our working directory
THREADS=8 #Max number of threads to use 
SRATOOLS="/ebio/abt6/yvoichek/smallproj/programs/sratoolkit.2.9.2-ubuntu64" # path to sratools (for fasterq-dump)
KMC="$library_dir/external_programs/kmc_v3"
# Parameters for building the k-mers table
K=31 #k-mer size
THRESHOLD_COUNT=2 #threshold to count k-mers per sample

# Sub-folders to use
samples_dir="$base_dir/samples"

# Create directories
mkdir -p $base_dir
mkdir -p $samples_dir

# Go over all strain names in the phenotype file and create separate k-mers list
line_index=0
while read name value 
do
	if [ $line_index -gt 0 ]
	then
		echo "Working on sample #$line_index : $name"
		sample_dir="$samples_dir/$name"
		mkdir -p $sample_dir # Create each samples directory
		# 1. Download sequence files
		CMD="$SRATOOLS/bin/fasterq-dump $name --outdir $sample_dir --temp $sample_dir --threads $THREADS --split-files >>$sample_dir/fasterq-dump.log.out 2>&1"
		fastq_r1="$sample_dir/${name}_1.fastq"
		fastq_r2="$sample_dir/${name}_2.fastq"
		download_try=1 # In case download doesn't succeed we will try multiple times
		while [[ ! -f $fastq_r1 ]] || [[ ! -f $fastq_r2 ]]
		do
			if [ $download_try -gt 1 ]
			then
				echo "Couldn't download sequencing file, wait two minutes and try again(#${download_try})"
				sleep 2m
			fi
			echo "$CMD"; eval $CMD
			let download_try++
		done
		# 2. Count k-mers twice using KMC
		input_file="$sample_dir/input_files.txt" # File with list of fastq files
		CMD="ls -l $sample_dir/*.fastq | awk '{print \$NF}' > $input_file"
		echo "$CMD"; eval $CMD
		kmc_canon="$sample_dir/kmc_canon" # KMC output of canonical counting
		kmc_non_canon="$sample_dir/kmc_non_canon" #KMC output of non-canonical counting
		# 2.1. KMC with canonization
		CMD="$KMC -t$THREADS -k$K -ci$THRESHOLD_COUNT @$input_file $kmc_canon $sample_dir > $sample_dir/kmc_canon.log.out"
		echo "$CMD"; eval $CMD
		# 2.2. KMC without canonization
		CMD="$KMC -t$THREADS -k$K -ci0 -b @$input_file $kmc_non_canon $sample_dir > $sample_dir/kmc_non_canon.log.out"
		echo "$CMD"; eval $CMD
		# 3. create a list of k-mers from the two KMC DBs
		CMD="$library_dir/bin/kmers_add_strand_information -c $kmc_canon -n $kmc_non_canon -k $K -o $sample_dir/kmers_with_strand > $sample_dir/add_strand.log.out"
		echo "$CMD"; eval $CMD
		# 4. delete KMC DBs and fastq files
		CMD="rm $sample_dir/*.kmc*"
		echo "$CMD"; eval $CMD
		CMD="rm $sample_dir/*.fastq"
		echo "$CMD"; eval $CMD
	fi
	((line_index++))
done < "./resistence.pheno"

# Create a file with the list of k-mers lists per sample
list_samples_file="$base_dir/kmers_list_paths.txt"
CMD="cat ./resistence.pheno | tail -n +2 | awk '{printf \"$samples_dir/%s/kmers_with_strand\\t%s\\n\",\$1,\$1}' > $list_samples_file"
echo "$CMD"; eval $CMD

# Filter k-mers from separate lists to one list with all k-mers to use
kmers_to_use_file="$base_dir/kmers_to_use"
CMD="$library_dir/bin/list_kmers_found_in_multiple_samples -l $list_samples_file -k $K --mac 5 -p 0.2 -o $kmers_to_use_file"
echo "$CMD"; eval $CMD

# Create the k-mers table
kmers_table_prefix="$base_dir/kmers_table"
CMD="$library_dir/bin/build_kmers_table -l $list_samples_file -k $K -a $kmers_to_use_file -o $kmers_table_prefix"
echo "$CMD"; eval $CMD

# At this stage all the k-mers information needed for GWAS is in the k-mers table. Therefore to save space we can remove all the individual sample
# k-mers list
CMD="rm $samples_dir/*/kmers_with_strand -f"
echo "$CMD"; eval $CMD


# Calculate the kinship matrix of the k-mers table (can also be calculated from a SNP matrix, making sure the order is the same)
CMD="$library_dir/bin/emma_kinship_kmers -t $kmers_table_prefix -k $K --maf 0.05 > $kmers_table_prefix.kinship"
echo "$CMD"; eval $CMD

# run GWAS
gwas_outdir="$base_dir/gwas_results"
CMD="python2.7 $library_dir/kmers_gwas.py --pheno resistence.pheno --kmers_table $kmers_table_prefix -l $K -p $THREADS --outdir $gwas_outdir >>$base_dir/gwas_run.log.txt 2>&1"
echo "$CMD"; eval $CMD

echo "k-mers passing 5% threshold are found in $gwas_outdir/kmers/pass_threshold_5per"
