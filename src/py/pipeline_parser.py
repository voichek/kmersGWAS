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

# output directory
parser.add_argument("--outdir", dest = "outdir", type=str, required=True,
        help='Directory to output pipeline results (shouldnt exist)')

# out names
parser.add_argument("-o", "--out", dest = "name", type=str, default="results",
        help='base name for all output files')

# k-mers presence/absence table
parser.add_argument("--kmers_table", dest = "kmers_table", type=str, required=True,
        help='Base for presence/absence table and accessions list')

# Parallel
parser.add_argument("-p", "--parallel", dest = "parallel", default=1, type=check_positive, 
        help='Maximal number of threads to use')

# k-mer length
parser.add_argument("-l", "--kmer_len", dest = "kmers_len", type=int,choices=range(15,32), metavar="[15-31]", 
        help='Length of K-mers in the database table')

# number of k-mers to take
parser.add_argument("-k", "--kmers_number", dest = "n_kmers", type=int, default=100001,
        help='Numbers of k-mers to filter from first step (due to bug in GEMMA 0.98 number shouldnt be a multiplication of 20K)')

# number of k-mers to take
parser.add_argument("--kmers_for_no_perm_phenotype", dest = "n_extra_phenotype_kmers", type=int, 
        help='Save a different number of k-mers in the heap for the non-permuted phenotype')

# number of snps to take (if a two step method is used)
parser.add_argument("--snps_number", dest = "n_snps", type=int, default=10001,
        help='Numbers of snps to filter from first step (used only if there using a two step snps approximation)')


# Number of permutation
parser.add_argument("--permutations", dest = "n_permutations", type=int, default=0,
        help='number of permutation for permutation test')

# Use kinship matrix from the kmers table
parser.add_argument("--kinship_kmers", dest = "use_kinship_from_kmers",
        help="Use the kinship matrix from kmers_table", action="store_true")

# Run SNPs associations in ONE step - only run GEMMA
parser.add_argument("--run_on_kmers", dest = "run_kmers", help="run pipeline on k-mers",
        action="store_true")

# Run SNPs associations in ONE step - only run GEMMA
parser.add_argument("--run_on_snps_one_step", dest = "run_one_step_snps", help="run pipeline with the same parameters on SNPs",
        action="store_true")

# RUN SNPs association in TWO steps - for permutations, first filter likley snps and then run GEMMA on them
parser.add_argument("--run_on_snps_two_steps", dest = "run_two_steps_snps", 
        help="run pipeline with the same parameters on SNPs - first filtering using GRAMMAR-Gamma and then using GEMMA for the top ones",
        action="store_true")

### Percent of missing values of SNPs to tolerate
##parser.add_argument("--miss_gemma", dest = "miss_gemma", type=float, default=0.5,
##        help='Tolerance for missing values in SNPs table')

## MAF (for k-mers and also for SNPs if used)
parser.add_argument("--maf", dest = "maf", type=float, default=0.05,
        help='Minor allele frequency')

## MAC (for k-mers and also for SNPs if used)
parser.add_argument("--mac", dest = "mac", type=float, default=5,
        help='Minor allele count')

## Min data poinrt
parser.add_argument("--min_data_points", dest = "min_data_points", type=float, default=30,
        help='Stop running if there is less data points than this threshold')

# SNP files (bed/bim/fam)
parser.add_argument("--snp_matrix", dest = "snps_matrix", type=str,
        help='base name for snps bed/bim/fam files')

# Control the verbosity of the program
parser.add_argument("-v", "--verbose", dest = "verbose",  help="increase output verbosity",
        action="store_true")

## count patterns of presence absence
parser.add_argument("--pattern_counter", dest = "kmers_pattern_counter",  
        help="Count the number of presence absence patterns k-mers pattern (has a large RAM usage)",
        action="store_true")

## an option specifying if to calculate or not a a qq plot for the k-mers (default is yes)
parser.add_argument("--no_qq_plot", dest = "qq_plot",  
        help="Don't calculate a qq plot (less computations)", action="store_false")

## Keep the intermediate files
parser.add_argument("--dont_remove_intermediates", dest = "remove_intermediate",  
        help="Mostly for debugging, to keep the intermediate files", action="store_false")

# path for GEMMA
parser.add_argument("--gemma_path", dest = "gemma_path", type=str,
        default = "/ebio/abt6/yvoichek/smallproj/prefix/bin/gemma",
        help='path to GEMMA')

args = parser.parse_args()
