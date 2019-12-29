import argparse # Needs python2.7+
import sys

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
        help='path to phenotype file (format:sample name[TAB]phenotype value[NEW-LINE]...)')

# output directory
parser.add_argument("--outdir", dest = "outdir", type=str, required=True,
        help='output directory')

# k-mers presence/absence table
parser.add_argument("--kmers_table", dest = "kmers_table", type=str, required=True,
        help='path to the k-mers table')

# k-mer length
parser.add_argument("-l", "--kmer_len", dest = "kmers_len", type=int,choices=range(15,32), metavar="[15-31]",required=True, 
        help='length of the k-mers')

# Parallel
parser.add_argument("-p", "--parallel", dest = "parallel", default=1, type=check_positive, 
        help='maximal number of threads')

# number of k-mers to take (notice: GEMMA 0.98 had a bug which didn't work with variant number which are multiplies of 20K)
parser.add_argument("-k", "--kmers_number", dest = "n_kmers", type=int, default=10001,
        help='numbers of k-mers to filter from first step') 

# number of k-mers to take
parser.add_argument("--kmers_for_no_perm_phenotype", dest = "n_extra_phenotype_kmers", type=int, 
        help='number of k-mers to filter from first step for the non-permuted phenotype (if different)')

# Number of permutation
parser.add_argument("--permutations", dest = "n_permutations", type=int, default=100,
        help='number of permutation for defining the threshold (should be at least 20)')

## MAF (for k-mers and also for SNPs if used)
parser.add_argument("--maf", dest = "maf", type=float, default=0.05,
        help='minor allele frequency')

## MAC (for k-mers and also for SNPs if used)
parser.add_argument("--mac", dest = "mac", type=float, default=5,
        help='minor allele count')

## Min data poinrt
parser.add_argument("--min_data_points", dest = "min_data_points", type=float, default=30,
        help='stop running with less phenotypic individuals than this threshold')

## count patterns of presence absence
parser.add_argument("--pattern_counter", dest = "kmers_pattern_counter",  
        help="count number of unique presence/absence k-mers patterns",
        action="store_true")

# Run SNPs associations in ONE step - only run GEMMA
parser.add_argument("--run_on_snps_one_step", dest = "run_one_step_snps", 
        help="run pipeline with the same parameters on SNPs (for comparison)",
        action="store_true")

# RUN SNPs association in TWO steps - for permutations, first filter likley snps and then run GEMMA on them
parser.add_argument("--run_on_snps_two_steps", dest = "run_two_steps_snps", 
        help="run pipeline with the same parameters on SNPs - first approximate for initial filtering and then using GEMMA top hits",
        action="store_true")

# SNP files (bed/bim/fam)
parser.add_argument("--snp_matrix", dest = "snps_matrix", type=str,
        help='base name for snps bed/bim/fam files')

# number of snps to take (if a two step method is used)
parser.add_argument("--snps_number", dest = "n_snps", type=int, default=10001,
        help='numbers of snps to filter from first step (used when running a two step snps approximation)')

# If you run it only on SNPs 
parser.add_argument("--dont_run_on_kmers", dest = "run_kmers", help="don't run pipeline on k-mers",
        action="store_false")

# Use kinship matrix from the kmers table
parser.add_argument("--kinship_snps", dest = "use_kinship_from_kmers",
        help="use the kinship matrix from the SNPs table", action="store_false")

## Keep the intermediate files
parser.add_argument("--dont_remove_intermediates", dest = "remove_intermediate",  
        help="for debugging, keep the intermediate files", action="store_false")

# path for GEMMA
parser.add_argument("--gemma_path", dest = "gemma_path", type=str,
        default = sys.path[0] + "/external_programs/gemma_0_96",
        help='path to GEMMA (tested with version 0.96)')

args = parser.parse_args()
