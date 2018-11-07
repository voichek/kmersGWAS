import argparse # Needs python2.7+

def check_positive(value):
    "Check if a variable entered to argparse is positive"
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

parser = argparse.ArgumentParser(description='Run k-mer associations and then GEMMA')

# Files:
# phenotype file
parser.add_argument("--pheno", dest = "fn_phenotype", type=str, required=True,
        help='phenotype file name Format:sample name[TAB]phenotype val[NEW-LINE]...)')

# 
# output directory
parser.add_argument("--outdir", dest = "outdir", type=str, required=True,
        help='Directory to output pipeline results (shouldnt exist)')

# out names
parser.add_argument("-o", "--out", dest = "name", type=str, default="results",
        help='base name for all output files')

# DBs list
parser.add_argument("-d", "--db_list", dest = "db_list", type=str, required=True,
        help='List of DBs')

# k-mers presence/absence table
parser.add_argument("--kmers_table", dest = "kmers_table", type=str, required=True,
        help='presence/absence table')

# path for GEMMA
parser.add_argument("--gemma_path", dest = "gemma_path", type=str,
        default = "/ebio/abt6/yvoichek/smallproj/prefix/bin/gemma",
        help='path to GEMMA')

# Parallel
parser.add_argument("-p", "--parallel", dest = "parallel", default=1, type=check_positive, 
        help='Maximal number of threads to use')

# k-mer length
parser.add_argument("-l", "--kmer_len", dest = "kmers_len", type=int,choices=range(15,32), metavar="[15-31]", 
        help='Length of K-mers in the database table')

# number of k-mers to take
parser.add_argument("-k", "--kmers_number", dest = "n_kmers", type=int, default=1000001,
        help='Numbers of k-mers to filter from first step (due to bug in GEMMA 0.98 number shouldnt be a multiplication of 20K)')

# MAF (for k-mers and also for SNPs if used)
parser.add_argument("--maf", dest = "maf", type=float, default=0.05,
        help='Minor allele frequency')

# Number of permutation
parser.add_argument("--permutations", dest = "n_permutations", type=int, default=0,
        help='number of permutation for permutation test')

# file with kinship matrix (regular IBS)
parser.add_argument("--kinship_IBS", dest = "kinship_IBS", type=str,
        default = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/tests/build_kinship_matrix_of_all_1001G/kinship_emma.maf0_05",
        help='IBS kinship matrix, should be positive positive semidefinite (used for phenotype transformation & permutations)')


# Run GEMMA on SNPs
parser.add_argument("--run_on_snps", dest = "run_snps", help="run pipeline with the same parameters on SNPs (for comparisons)",
        action="store_true")

# Percent of missing values of SNPs to tolerate
parser.add_argument("--miss_gemma", dest = "miss_gemma", type=float, default=0.5,
        help='Tolerance for missing values in SNPs table')

# file with kinship matrix (created by GEMMA)
parser.add_argument("--kinship_gemma", dest = "kinship_GEMMA", type=str,
        default = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/tests/build_kinship_matrix_of_all_1001G/kinship_emma.maf0_05",
        help='kinship produced by GEMMA for running GEMMA in the second step')
        #        default = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/tests/build_kinship_matrix_of_all_1001G/kinshipname_all_miss_0_75_gk2.sXX.txt",

# SNP files (bed/bim/fam)
parser.add_argument("--snp_matrix", dest = "snps_matrix", type=str,
        default = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/tests/build_kinship_matrix_of_all_1001G/1001genomes_snp-short-indel_only_ACGTN.vcf.plink",
        help='base name for snps bed/bim/fam files')


parser.add_argument("-v", "--verbose", dest = "verbose",  help="increase output verbosity",
        action="store_true")


args = parser.parse_args()
